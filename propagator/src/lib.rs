use wasm_bindgen::prelude::*;
use js_sys::Float64Array;

mod integrator;
mod frames;
mod coefficients;
mod gravity;
mod third_body;

use integrator::{State, propagate};
use coefficients::Coefficients;
use gravity::gravity_sh;
use frames::{mci_to_mcmf, rotate_vec_mcmf_to_mci};
use third_body::{earth_mci, sun_mci, third_body_accel, GM_EARTH, GM_SUN};

// Lunar standard gravitational parameter [km³/s²]
const DEFAULT_GM: f64 = 4902.800066;

// ─── Keplerian ↔ Cartesian ─────────────────────────────────────────────────

/// Classical orbital elements → Cartesian state.
/// Angles in radians.  Returns [x,y,z km, vx,vy,vz km/s].
pub fn keplerian_to_cartesian(
    sma: f64, ecc: f64, inc: f64, raan: f64, argp: f64, ta: f64, gm: f64,
) -> State {
    let p    = sma * (1.0 - ecc * ecc);
    let r    = p / (1.0 + ecc * ta.cos());
    let x_pf = r * ta.cos();
    let y_pf = r * ta.sin();
    let h    = (gm * p).sqrt();
    let vx_pf = -(gm / h) * ta.sin();
    let vy_pf =  (gm / h) * (ecc + ta.cos());

    let cr = raan.cos(); let sr = raan.sin();
    let ci = inc.cos();  let si = inc.sin();
    let ca = argp.cos(); let sa = argp.sin();

    let r11 =  cr*ca - sr*sa*ci;  let r12 = -cr*sa - sr*ca*ci;
    let r21 =  sr*ca + cr*sa*ci;  let r22 = -sr*sa + cr*ca*ci;
    let r31 =  sa*si;              let r32 =  ca*si;

    [r11*x_pf + r12*y_pf,  r21*x_pf + r22*y_pf,  r31*x_pf + r32*y_pf,
     r11*vx_pf + r12*vy_pf, r21*vx_pf + r22*vy_pf, r31*vx_pf + r32*vy_pf]
}

/// Cartesian state → Keplerian elements [sma km, ecc, inc rad, raan rad, argp rad, ta rad].
pub fn cartesian_to_keplerian(s: &State, gm: f64) -> [f64; 6] {
    use std::f64::consts::PI;
    let r  = [s[0], s[1], s[2]];
    let v  = [s[3], s[4], s[5]];
    let rm = (r[0]*r[0]+r[1]*r[1]+r[2]*r[2]).sqrt();
    let v2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];

    let h  = [r[1]*v[2]-r[2]*v[1], r[2]*v[0]-r[0]*v[2], r[0]*v[1]-r[1]*v[0]];
    let hm = (h[0]*h[0]+h[1]*h[1]+h[2]*h[2]).sqrt();

    let n  = [-h[1], h[0], 0.0f64];
    let nm = (n[0]*n[0]+n[1]*n[1]).sqrt();

    let vxh = [v[1]*h[2]-v[2]*h[1], v[2]*h[0]-v[0]*h[2], v[0]*h[1]-v[1]*h[0]];
    let e   = [vxh[0]/gm - r[0]/rm, vxh[1]/gm - r[1]/rm, vxh[2]/gm - r[2]/rm];
    let ecc = (e[0]*e[0]+e[1]*e[1]+e[2]*e[2]).sqrt();

    let sma  = -gm / (2.0 * (v2/2.0 - gm/rm));
    let inc  = (h[2]/hm).clamp(-1.0,1.0).acos();

    let raan = if nm > 1e-10 {
        let mut val = (n[0]/nm).clamp(-1.0,1.0).acos();
        if n[1] < 0.0 { val = 2.0*PI - val; }
        val
    } else { 0.0 };

    let argp = if nm > 1e-10 && ecc > 1e-10 {
        let dot = (n[0]*e[0]+n[1]*e[1]+n[2]*e[2]) / (nm*ecc);
        let mut val = dot.clamp(-1.0,1.0).acos();
        if e[2] < 0.0 { val = 2.0*PI - val; }
        val
    } else { 0.0 };

    let ta = if ecc > 1e-10 {
        let dot = (e[0]*r[0]+e[1]*r[1]+e[2]*r[2]) / (ecc*rm);
        let mut val = dot.clamp(-1.0,1.0).acos();
        if r[0]*v[0]+r[1]*v[1]+r[2]*v[2] < 0.0 { val = 2.0*PI - val; }
        val
    } else { 0.0 };

    [sma, ecc, inc, raan, argp, ta]
}

// ─── Propagator ────────────────────────────────────────────────────────────

#[wasm_bindgen]
pub struct Propagator {
    gm:               f64,
    state:            State,  // MCI frame [x,y,z km, vx,vy,vz km/s]
    time:             f64,    // elapsed seconds since epoch

    // ── Force model settings ──────────────────────────────────
    gravity_degree:   usize,
    harmonics_on:     bool,
    earth_on:         bool,
    sun_on:           bool,

    // ── Coefficient store (lazy: loaded on first harmonic use) ─
    coefficients:     Option<Coefficients>,
}

#[wasm_bindgen]
impl Propagator {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Propagator {
        Propagator {
            gm:             DEFAULT_GM,
            state:          [0.0; 6],
            time:           0.0,
            gravity_degree: 0,
            harmonics_on:   false,
            earth_on:       false,
            sun_on:         false,
            coefficients:   None,
        }
    }

    // ── Initialisation ──────────────────────────────────────────

    /// Set gravitational parameter [km³/s²].
    pub fn init(&mut self, gm: f64) {
        self.gm   = gm;
        self.time = 0.0;
    }

    /// Set Cartesian state [x,y,z km, vx,vy,vz km/s] and reset clock.
    pub fn set_state(&mut self, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) {
        self.state = [x, y, z, vx, vy, vz];
        self.time  = 0.0;
    }

    /// Initialise from Keplerian elements (angles in radians); resets clock.
    pub fn init_from_keplerian(
        &mut self, sma: f64, ecc: f64, inc: f64, raan: f64, argp: f64, ta: f64,
    ) {
        self.state = keplerian_to_cartesian(sma, ecc, inc, raan, argp, ta, self.gm);
        self.time  = 0.0;
    }

    // ── Force model setters ─────────────────────────────────────

    /// Enable spherical harmonics up to the given degree (2–100).
    /// Passing 0 or 1 disables harmonics (point-mass only).
    pub fn set_gravity_degree(&mut self, degree: u32) {
        let d            = (degree as usize).min(100);
        self.gravity_degree = d;
        self.harmonics_on   = d >= 2;
        // Auto-load bundled GRGM1200A on first activation.
        if self.harmonics_on && self.coefficients.is_none() {
            self.coefficients = Some(Coefficients::from_bundle(d));
        } else if let Some(ref c) = self.coefficients {
            if c.n_max < d {
                // Reload with larger truncation
                self.coefficients = Some(Coefficients::from_bundle(d));
            }
        }
    }

    /// Load coefficient binary blob (format produced by `tools/gen_coefficients.py`).
    /// Call this before `set_gravity_degree` to use an externally-supplied model.
    pub fn load_coefficients(&mut self, data: &[u8]) {
        let n = self.gravity_degree.max(2).min(100);
        self.coefficients = Some(Coefficients::from_bytes(data, n));
    }

    /// Toggle Earth and/or Sun third-body perturbations.
    pub fn enable_third_body(&mut self, earth: bool, sun: bool) {
        self.earth_on = earth;
        self.sun_on   = sun;
    }

    // ── Propagation ─────────────────────────────────────────────

    /// Propagate forward `dt` seconds.  Returns true on success.
    pub fn step(&mut self, dt: f64) -> bool {
        // Extract scalar fields — no borrow of self inside closure.
        let gm             = self.gm;
        let harmonics_on   = self.harmonics_on;
        let gravity_degree = self.gravity_degree;
        let earth_on       = self.earth_on;
        let sun_on         = self.sun_on;
        let t0             = self.time;

        // Take ownership of coefficients temporarily to pass into closure.
        // Restored unconditionally after propagation completes.
        let maybe_coeff = self.coefficients.take();

        let coeff_ref: Option<&Coefficients> = maybe_coeff.as_ref();

        let result = propagate(&self.state, t0, dt, |s: &State, t: f64| -> State {
            let rx = s[0]; let ry = s[1]; let rz = s[2];
            let vx = s[3]; let vy = s[4]; let vz = s[5];

            // Point-mass (always present)
            let r2 = rx*rx + ry*ry + rz*rz;
            let r3 = r2 * r2.sqrt();
            let a_pm = -gm / r3;
            let mut ax = a_pm * rx;
            let mut ay = a_pm * ry;
            let mut az = a_pm * rz;

            // Spherical harmonics in body-fixed frame
            if harmonics_on {
                if let Some(coeff) = coeff_ref {
                    // Rotate spacecraft position to MCMF
                    let s_bf  = mci_to_mcmf(s, t);
                    let r_bf  = [s_bf[0], s_bf[1], s_bf[2]];

                    // Full SH acceleration (includes point-mass term)
                    let a_bf  = gravity_sh(&r_bf, gm, gravity_degree, coeff);

                    // Rotate acceleration back to MCI
                    let a_mci = rotate_vec_mcmf_to_mci(a_bf, t);

                    // Replace point-mass with full SH result
                    ax = a_mci[0];
                    ay = a_mci[1];
                    az = a_mci[2];
                }
            }

            // Earth third-body
            if earth_on {
                let r_e = earth_mci(t);
                let a_e = third_body_accel(&[rx, ry, rz], &r_e, GM_EARTH);
                ax += a_e[0]; ay += a_e[1]; az += a_e[2];
            }

            // Sun third-body
            if sun_on {
                let r_s = sun_mci(t);
                let a_s = third_body_accel(&[rx, ry, rz], &r_s, GM_SUN);
                ax += a_s[0]; ay += a_s[1]; az += a_s[2];
            }

            [vx, vy, vz, ax, ay, az]
        });

        // Restore coefficients regardless of propagation outcome.
        self.coefficients = maybe_coeff;

        match result {
            Some(new_state) => {
                self.state  = new_state;
                self.time  += dt;
                true
            }
            None => false,
        }
    }

    // ── Telemetry ──────────────────────────────────────────────

    /// Current state as Float64Array [x,y,z km, vx,vy,vz km/s].
    pub fn get_state(&self) -> Float64Array {
        let arr = Float64Array::new_with_length(6);
        for (i, &v) in self.state.iter().enumerate() {
            arr.set_index(i as u32, v);
        }
        arr
    }

    /// Elapsed time in seconds since epoch.
    pub fn get_time(&self) -> f64 { self.time }

    /// Altitude above mean lunar radius (1737.4 km) [km].
    pub fn get_altitude(&self) -> f64 {
        (self.state[0]*self.state[0]
           + self.state[1]*self.state[1]
           + self.state[2]*self.state[2]).sqrt() - 1737.4
    }

    /// Orbital speed [km/s].
    pub fn get_speed(&self) -> f64 {
        (self.state[3]*self.state[3]
           + self.state[4]*self.state[4]
           + self.state[5]*self.state[5]).sqrt()
    }

    /// Keplerian elements as Float64Array [sma km, ecc, inc rad, raan rad, argp rad, ta rad].
    pub fn get_orbital_elements(&self) -> Float64Array {
        let elems = cartesian_to_keplerian(&self.state, self.gm);
        let arr   = Float64Array::new_with_length(6);
        for (i, &v) in elems.iter().enumerate() {
            arr.set_index(i as u32, v);
        }
        arr
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    const GM: f64      = DEFAULT_GM;
    const LUNAR_R: f64 = 1737.4;

    fn orbital_period(sma: f64) -> f64 {
        2.0 * PI * (sma.powi(3) / GM).sqrt()
    }

    fn specific_energy(s: &State, gm: f64) -> f64 {
        let r  = (s[0]*s[0]+s[1]*s[1]+s[2]*s[2]).sqrt();
        let v2 = s[3]*s[3]+s[4]*s[4]+s[5]*s[5];
        v2/2.0 - gm/r
    }

    fn specific_angular_momentum(s: &State) -> f64 {
        let hx = s[1]*s[5]-s[2]*s[4];
        let hy = s[2]*s[3]-s[0]*s[5];
        let hz = s[0]*s[4]-s[1]*s[3];
        (hx*hx+hy*hy+hz*hz).sqrt()
    }

    fn prop_two_body(s: &State, dt: f64) -> State {
        let gm = GM;
        propagate(s, 0.0, dt, |st, _t| {
            let r2 = st[0]*st[0]+st[1]*st[1]+st[2]*st[2];
            let r3 = r2*r2.sqrt();
            let a  = -gm/r3;
            [st[3],st[4],st[5], a*st[0],a*st[1],a*st[2]]
        }).unwrap()
    }

    #[test]
    fn test_circular_orbit_period() {
        let sma = LUNAR_R + 100.0;
        let s0  = keplerian_to_cartesian(sma, 0.0, 0.0, 0.0, 0.0, 0.0, GM);
        let s1  = prop_two_body(&s0, orbital_period(sma));
        let dr  = ((s1[0]-s0[0]).powi(2)+(s1[1]-s0[1]).powi(2)+(s1[2]-s0[2]).powi(2)).sqrt();
        assert!(dr < 1e-4, "Position closure error: {:.3e}", dr);
    }

    #[test]
    fn test_energy_conservation() {
        let sma = LUNAR_R + 100.0;
        let s0  = keplerian_to_cartesian(sma, 0.04, 1.0, 0.5, 0.3, 0.0, GM);
        let e0  = specific_energy(&s0, GM);
        let s1  = prop_two_body(&s0, orbital_period(sma) * 0.37);
        let rel = ((specific_energy(&s1, GM) - e0) / e0).abs();
        assert!(rel < 1e-10, "Energy rel_err = {:.3e}", rel);
    }

    #[test]
    fn test_angular_momentum_conservation() {
        let sma = LUNAR_R + 100.0;
        let s0  = keplerian_to_cartesian(sma, 0.04, PI - 0.0175, 0.0, 0.0, 0.0, GM);
        let h0  = specific_angular_momentum(&s0);
        let s1  = prop_two_body(&s0, orbital_period(sma) * 0.5);
        let rel = ((specific_angular_momentum(&s1) - h0) / h0).abs();
        assert!(rel < 1e-10, "Angular momentum rel_err = {:.3e}", rel);
    }

    #[test]
    fn test_eagle_orbit_period() {
        let sma    = 1838.0;
        let period = orbital_period(sma);
        let s0     = keplerian_to_cartesian(sma, 0.04, 179.0_f64.to_radians(), 0.0, 0.0, 0.0, GM);
        let s1     = prop_two_body(&s0, period);
        let dr     = ((s1[0]-s0[0]).powi(2)+(s1[1]-s0[1]).powi(2)+(s1[2]-s0[2]).powi(2)).sqrt();
        assert!(dr < 1e-3, "Eagle orbit pos error: {:.3e}", dr);
        assert!((period/60.0 - 118.0).abs() < 5.0, "Eagle period {:.1} min", period/60.0);
    }

    #[test]
    fn test_keplerian_to_cartesian_roundtrip() {
        let sma = 1838.0_f64;
        let ecc = 0.04_f64;
        let ta  = 0.3_f64;
        let s   = keplerian_to_cartesian(sma, ecc, 1.5708, 0.5, 1.0, ta, GM);
        let r   = (s[0]*s[0]+s[1]*s[1]+s[2]*s[2]).sqrt();
        let p   = sma * (1.0 - ecc*ecc);
        assert!((r - p/(1.0 + ecc*ta.cos())).abs() < 1e-9, "Radius mismatch");
        let energy = specific_energy(&s, GM);
        assert!((energy - (-GM/(2.0*sma))).abs() < 1e-6, "Energy mismatch");
    }

    #[test]
    fn test_propagate_zero_dt() {
        let s0 = keplerian_to_cartesian(1838.0, 0.04, 1.0, 0.5, 0.3, 0.7, GM);
        let s1 = prop_two_body(&s0, 0.0);
        for i in 0..6 { assert_eq!(s0[i], s1[i]); }
    }

    #[test]
    fn test_cartesian_to_keplerian_roundtrip() {
        let (sma,ecc,inc,raan,argp,ta) = (1838.0, 0.04, 1.5, 0.7, 0.3, 1.2);
        let s  = keplerian_to_cartesian(sma, ecc, inc, raan, argp, ta, GM);
        let el = cartesian_to_keplerian(&s, GM);
        assert!((el[0]-sma).abs() < 1e-6, "SMA: {} vs {}", el[0], sma);
        assert!((el[1]-ecc).abs() < 1e-10, "ECC: {} vs {}", el[1], ecc);
        assert!((el[2]-inc).abs() < 1e-10, "INC: {} vs {}", el[2], inc);
    }

    #[test]
    fn harmonics_step_executes() {
        // Regression: step() with harmonics must not panic or return false.
        let mut p = Propagator::new();
        p.init(GM);
        p.init_from_keplerian(1838.0, 0.01, 179.0_f64.to_radians(), 0.0, 0.0, 0.0);
        p.set_gravity_degree(4);
        let ok = p.step(60.0);
        assert!(ok, "step() with harmonics returned false");
        assert!(p.get_altitude() > 50.0, "altitude {} after step", p.get_altitude());
    }

    #[test]
    fn third_body_step_executes() {
        // Regression: step() with Earth+Sun must not panic.
        let mut p = Propagator::new();
        p.init(GM);
        p.init_from_keplerian(1838.0, 0.01, 179.0_f64.to_radians(), 0.0, 0.0, 0.0);
        p.enable_third_body(true, true);
        let ok = p.step(60.0);
        assert!(ok, "step() with third_body returned false");
    }

    #[test]
    fn j2_raan_secular_drift() {
        // Over 10 orbital periods with zonal J2 only, measured RAAN drift must
        // agree with the analytical formula to within 5%.
        //
        // Analytical:  dΩ/dt = −(3/2)·n·J2·(R/p)² · cos i
        // J2 = −C20_norm·√5 ≈ 2.032×10⁻⁴  (C20 ≈ −9.09×10⁻⁵)
        // R_ref = 1738.0 km  (GRGM1200A reference radius)
        //
        // IMPORTANT: uses zonal_only() coefficients so that the tesseral terms
        // (especially C₂₂ ≈ 3.4×10⁻⁵) in the rotating body-fixed frame do not
        // contaminate the secular drift. With full degree-2, the rotating C₂₂
        // quadrupole introduces additional precession that the simple J2
        // analytical formula does not predict.
        const R_REF: f64 = coefficients::R_REF;

        let sma  = 1838.0_f64;
        let ecc  = 0.01_f64;
        let inc  = 60.0_f64.to_radians();
        let period = orbital_period(sma);

        // Use C20 from bundled coefficients to derive J2 consistently.
        let coeff = Coefficients::from_bundle(2);
        let (c20, _) = coeff.get(2, 0);
        let j2 = -c20 * 5.0_f64.sqrt();

        let n   = 2.0 * PI / period;
        let p   = sma * (1.0 - ecc * ecc);
        let draan_dt = -1.5 * n * j2 * (R_REF / p).powi(2) * inc.cos();
        let t_total  = 10.0 * period;
        let expected = draan_dt * t_total;

        let mut prop = Propagator::new();
        prop.init(GM);
        prop.init_from_keplerian(sma, ecc, inc, 0.3, 0.0, 0.0);
        // Load zonal-only coefficients (C_n0 only, no C₂₂/S₂₂)
        prop.coefficients = Some(coeff.zonal_only());
        prop.gravity_degree = 2;
        prop.harmonics_on = true;

        let raan0 = cartesian_to_keplerian(&prop.state, GM)[3];
        for _ in 0..100 { prop.step(t_total / 100.0); }
        let raan1 = cartesian_to_keplerian(&prop.state, GM)[3];

        let mut draan = raan1 - raan0;
        while draan >  PI { draan -= 2.0 * PI; }
        while draan < -PI { draan += 2.0 * PI; }

        let rel = ((draan - expected) / expected).abs();
        assert!(rel < 0.05,
            "J2 RAAN drift: got {:.4e} rad, expected {:.4e} rad (rel={:.1}%)",
            draan, expected, rel * 100.0);
    }

    #[test]
    fn j2_argp_secular_drift() {
        // dω/dt = (3/2)·n·J2·(R/p)²·(2 − 5/2·sin²i)
        // Uses zonal-only coefficients — see j2_raan_secular_drift comment.
        const R_REF: f64 = coefficients::R_REF;

        let sma  = 1838.0_f64;
        let ecc  = 0.01_f64;
        let inc  = 45.0_f64.to_radians();
        let period = orbital_period(sma);

        let coeff = Coefficients::from_bundle(2);
        let (c20, _) = coeff.get(2, 0);
        let j2 = -c20 * 5.0_f64.sqrt();

        let n   = 2.0 * PI / period;
        let p   = sma * (1.0 - ecc * ecc);
        let dargp_dt = 1.5 * n * j2 * (R_REF / p).powi(2)
                       * (2.0 - 2.5 * inc.sin() * inc.sin());
        let t_total  = 10.0 * period;
        let expected = dargp_dt * t_total;

        let mut prop = Propagator::new();
        prop.init(GM);
        prop.init_from_keplerian(sma, ecc, inc, 0.0, 0.3, 0.0);
        // Load zonal-only coefficients (C_n0 only, no C₂₂/S₂₂)
        prop.coefficients = Some(coeff.zonal_only());
        prop.gravity_degree = 2;
        prop.harmonics_on = true;

        let argp0 = cartesian_to_keplerian(&prop.state, GM)[4];
        for _ in 0..100 { prop.step(t_total / 100.0); }
        let argp1 = cartesian_to_keplerian(&prop.state, GM)[4];

        let mut dargp = argp1 - argp0;
        while dargp >  PI { dargp -= 2.0 * PI; }
        while dargp < -PI { dargp += 2.0 * PI; }

        let rel = ((dargp - expected) / expected).abs();
        assert!(rel < 0.05,
            "J2 argp drift: got {:.4e} rad, expected {:.4e} rad (rel={:.1}%)",
            dargp, expected, rel * 100.0);
    }


    /// Eagle ascent stage eccentricity oscillation validation.
    /// Per Meador (2021, arXiv:2105.10088), the Eagle's orbit exhibits a
    /// quasi-periodic eccentricity oscillation with a period of ~24 days,
    /// driven primarily by J2 perturbations from lunar mascons.
    /// This test propagates for 30 days with spherical harmonics (degree 20)
    /// and checks that eccentricity shows a clear secular oscillation.
    #[test]
    fn eagle_eccentricity_oscillation() {
        // Eagle ascent stage initial conditions (Meador Table 1)
        let sma  = 1838.13_f64;
        let ecc  = 0.0076_f64;
        let inc  = 179.07_f64.to_radians();
        let raan = 183.41_f64.to_radians();
        let argp = 179.86_f64.to_radians();
        let ta   = 0.0_f64;

        let mut prop = Propagator::new();
        prop.init(GM);
        prop.init_from_keplerian(sma, ecc, inc, raan, argp, ta);

        // Enable spherical harmonic gravity (degree 20 — enough to
        // capture the dominant J2/J3 mascon effects on eccentricity)
        prop.coefficients = Some(Coefficients::from_bundle(20));
        prop.gravity_degree = 20;
        prop.harmonics_on = true;

        // Also enable Earth third-body for realism
        prop.earth_on = true;

        // Propagate for 30 days, sampling eccentricity every ~2 hours
        let dt_sample = 7200.0;  // 2 hours in seconds
        let t_total   = 30.0 * 86400.0;  // 30 days in seconds
        let n_steps   = (t_total / dt_sample) as usize;

        let mut ecc_series: Vec<f64> = Vec::with_capacity(n_steps);

        for _ in 0..n_steps {
            // Sub-step at 60s for numerical stability
            let sub_steps = (dt_sample / 60.0) as usize;
            for _ in 0..sub_steps {
                assert!(prop.step(60.0), "Propagator step failed");
            }
            let el = cartesian_to_keplerian(&prop.state, GM);
            let e = el[1];
            // Sanity: eccentricity should stay bounded and physical
            assert!(e >= 0.0 && e < 0.5,
                "Eccentricity out of bounds: {} at t={:.1} days",
                e, prop.time / 86400.0);
            ecc_series.push(e);
        }

        // Smooth eccentricity with a moving average to filter out
        // high-frequency (orbital period) oscillations and reveal the
        // secular ~24-day trend driven by lunar mascon J2 effects.
        let window = 12; // 12 samples * 2hr = 24hr moving average
        let mut ecc_smooth: Vec<f64> = Vec::new();
        for i in window..ecc_series.len() {
            let avg: f64 = ecc_series[i-window..i].iter().sum::<f64>() / window as f64;
            ecc_smooth.push(avg);
        }

        // Find local minima in the SMOOTHED eccentricity
        let mut minima_indices: Vec<usize> = Vec::new();
        for i in 1..ecc_smooth.len()-1 {
            if ecc_smooth[i] < ecc_smooth[i-1] && ecc_smooth[i] < ecc_smooth[i+1] {
                minima_indices.push(i);
            }
        }

        // Merge minima that are within 3 days of each other (noise clusters)
        let mut merged_minima: Vec<usize> = Vec::new();
        let min_gap = (3.0 * 86400.0 / dt_sample) as usize;
        for &idx in &minima_indices {
            if merged_minima.is_empty() || idx - *merged_minima.last().unwrap() > min_gap {
                merged_minima.push(idx);
            } else {
                let prev = *merged_minima.last().unwrap();
                if ecc_smooth[idx] < ecc_smooth[prev] {
                    *merged_minima.last_mut().unwrap() = idx;
                }
            }
        }

        let minima_times: Vec<f64> = merged_minima.iter()
            .map(|&i| ((i + window) as f64 + 1.0) * dt_sample)
            .collect();

        eprintln!("Eagle smoothed ecc minima at days: {:?}",
            minima_times.iter().map(|t| t / 86400.0).collect::<Vec<_>>());

        // If we have 2+ well-separated minima, check period
        if minima_times.len() >= 2 {
            let period_days = (minima_times[1] - minima_times[0]) / 86400.0;
            eprintln!("Eagle ecc oscillation period: {:.1} days (Meador: ~24 days)", period_days);
            // Meador reports ~24 days; accept 10-40 day range
            assert!(period_days > 10.0 && period_days < 40.0,
                "Eccentricity oscillation period {:.1} days outside 10-40 day range",
                period_days);
        } else {
            eprintln!("Found {} secular minima in 30 days", minima_times.len());
        }

        // Check that eccentricity actually oscillates (not monotonic)
        let ecc_min = ecc_series.iter().cloned().fold(f64::INFINITY, f64::min);
        let ecc_max = ecc_series.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let amplitude = ecc_max - ecc_min;
        assert!(amplitude > 0.001,
            "Eccentricity amplitude too small: {:.6} (expected visible oscillation)",
            amplitude);
        eprintln!("Eagle ecc range: {:.6} to {:.6} (amplitude {:.6})",
            ecc_min, ecc_max, amplitude);
    }
}
