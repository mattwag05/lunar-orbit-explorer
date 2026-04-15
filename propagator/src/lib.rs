use wasm_bindgen::prelude::*;
use js_sys::Float64Array;

mod integrator;
use integrator::{State, propagate};

// Lunar standard gravitational parameter (km³/s²)
const DEFAULT_GM: f64 = 4902.800066;

// ─── Keplerian orbit elements → Cartesian state ───────────────────────────

/// Convert classical orbital elements to Cartesian state in the same frame.
/// sma: semi-major axis (km)
/// ecc: eccentricity
/// inc: inclination (radians)
/// raan: right ascension of ascending node (radians)
/// argp: argument of perigee (radians)
/// ta: true anomaly (radians)
/// Returns [x, y, z, vx, vy, vz] in km, km/s
pub fn keplerian_to_cartesian(sma: f64, ecc: f64, inc: f64, raan: f64, argp: f64, ta: f64, gm: f64) -> State {
    // Semi-latus rectum
    let p = sma * (1.0 - ecc * ecc);
    let r = p / (1.0 + ecc * ta.cos());

    // Position in perifocal frame
    let x_pf = r * ta.cos();
    let y_pf = r * ta.sin();

    // Velocity in perifocal frame
    let h = (gm * p).sqrt();
    let vx_pf = -(gm / h) * ta.sin();
    let vy_pf =  (gm / h) * (ecc + ta.cos());

    // Rotation matrix from perifocal to inertial
    let cos_raan = raan.cos();
    let sin_raan = raan.sin();
    let cos_inc  = inc.cos();
    let sin_inc  = inc.sin();
    let cos_argp = argp.cos();
    let sin_argp = argp.sin();

    // R = Rz(-raan) * Rx(-inc) * Rz(-argp)
    let r11 =  cos_raan * cos_argp - sin_raan * sin_argp * cos_inc;
    let r12 = -cos_raan * sin_argp - sin_raan * cos_argp * cos_inc;
    let r21 =  sin_raan * cos_argp + cos_raan * sin_argp * cos_inc;
    let r22 = -sin_raan * sin_argp + cos_raan * cos_argp * cos_inc;
    let r31 =  sin_argp * sin_inc;
    let r32 =  cos_argp * sin_inc;

    let x  = r11 * x_pf + r12 * y_pf;
    let y  = r21 * x_pf + r22 * y_pf;
    let z  = r31 * x_pf + r32 * y_pf;
    let vx = r11 * vx_pf + r12 * vy_pf;
    let vy = r21 * vx_pf + r22 * vy_pf;
    let vz = r31 * vx_pf + r32 * vy_pf;

    [x, y, z, vx, vy, vz]
}

/// Cartesian state [x,y,z km, vx,vy,vz km/s] → Keplerian elements.
/// Returns [sma km, ecc, inc rad, raan rad, argp rad, ta rad].
pub fn cartesian_to_keplerian(s: &State, gm: f64) -> [f64; 6] {
    use std::f64::consts::PI;

    let r = [s[0], s[1], s[2]];
    let v = [s[3], s[4], s[5]];

    let r_mag = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
    let v2    = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

    // Specific angular momentum h = r × v
    let h = [
        r[1]*v[2] - r[2]*v[1],
        r[2]*v[0] - r[0]*v[2],
        r[0]*v[1] - r[1]*v[0],
    ];
    let h_mag = (h[0]*h[0] + h[1]*h[1] + h[2]*h[2]).sqrt();

    // Ascending node vector n = k × h (k = [0,0,1])
    let n = [-h[1], h[0], 0.0f64];
    let n_mag = (n[0]*n[0] + n[1]*n[1]).sqrt();

    // Eccentricity vector e = (v×h)/GM - r/|r|
    let vxh = [
        v[1]*h[2] - v[2]*h[1],
        v[2]*h[0] - v[0]*h[2],
        v[0]*h[1] - v[1]*h[0],
    ];
    let e = [
        vxh[0]/gm - r[0]/r_mag,
        vxh[1]/gm - r[1]/r_mag,
        vxh[2]/gm - r[2]/r_mag,
    ];
    let ecc = (e[0]*e[0] + e[1]*e[1] + e[2]*e[2]).sqrt();

    // Specific energy → SMA
    let energy = v2 / 2.0 - gm / r_mag;
    let sma    = -gm / (2.0 * energy);

    // Inclination
    let inc = (h[2] / h_mag).clamp(-1.0, 1.0).acos();

    // RAAN
    let raan = if n_mag > 1e-10 {
        let mut val = (n[0] / n_mag).clamp(-1.0, 1.0).acos();
        if n[1] < 0.0 { val = 2.0 * PI - val; }
        val
    } else { 0.0 };

    // Argument of perigee
    let argp = if n_mag > 1e-10 && ecc > 1e-10 {
        let dot = (n[0]*e[0] + n[1]*e[1] + n[2]*e[2]) / (n_mag * ecc);
        let mut val = dot.clamp(-1.0, 1.0).acos();
        if e[2] < 0.0 { val = 2.0 * PI - val; }
        val
    } else { 0.0 };

    // True anomaly
    let ta = if ecc > 1e-10 {
        let dot = (e[0]*r[0] + e[1]*r[1] + e[2]*r[2]) / (ecc * r_mag);
        let mut val = dot.clamp(-1.0, 1.0).acos();
        let rdotv = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
        if rdotv < 0.0 { val = 2.0 * PI - val; }
        val
    } else { 0.0 };

    [sma, ecc, inc, raan, argp, ta]
}

// ─── WASM-bindgen API ─────────────────────────────────────────────────────

#[wasm_bindgen]
pub struct Propagator {
    gm:    f64,
    state: State,   // MCI frame: [x,y,z km, vx,vy,vz km/s]
    time:  f64,     // elapsed seconds since epoch
}

#[wasm_bindgen]
impl Propagator {
    /// Create a new propagator (point-mass, no force model active).
    #[wasm_bindgen(constructor)]
    pub fn new() -> Propagator {
        Propagator {
            gm:    DEFAULT_GM,
            state: [0.0; 6],
            time:  0.0,
        }
    }

    /// Set gravitational parameter [km³/s²].
    pub fn init(&mut self, gm: f64) {
        self.gm   = gm;
        self.time = 0.0;
    }

    /// Set Cartesian state [x, y, z (km), vx, vy, vz (km/s)].
    pub fn set_state(&mut self, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) {
        self.state = [x, y, z, vx, vy, vz];
        self.time  = 0.0;
    }

    /// Initialise from classical Keplerian elements (angles in radians).
    pub fn init_from_keplerian(
        &mut self,
        sma: f64, ecc: f64, inc: f64, raan: f64, argp: f64, ta: f64,
    ) {
        self.state = keplerian_to_cartesian(sma, ecc, inc, raan, argp, ta, self.gm);
        self.time  = 0.0;
    }

    /// Propagate forward by `dt` seconds. Returns true on success.
    pub fn step(&mut self, dt: f64) -> bool {
        let gm = self.gm;
        let result = propagate(&self.state, dt, |s: &State| -> State {
            let r2 = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
            let r3 = r2 * r2.sqrt();
            let a  = -gm / r3;
            [s[3], s[4], s[5], a*s[0], a*s[1], a*s[2]]
        });
        match result {
            Some(new_state) => {
                self.state = new_state;
                self.time += dt;
                true
            }
            None => false,
        }
    }

    /// Get current state as Float64Array [x, y, z, vx, vy, vz].
    pub fn get_state(&self) -> Float64Array {
        let arr = Float64Array::new_with_length(6);
        for (i, &v) in self.state.iter().enumerate() {
            arr.set_index(i as u32, v);
        }
        arr
    }

    /// Get elapsed time in seconds.
    pub fn get_time(&self) -> f64 { self.time }

    /// Get current altitude above mean lunar radius (1737.4 km).
    pub fn get_altitude(&self) -> f64 {
        let r = (self.state[0].powi(2) + self.state[1].powi(2) + self.state[2].powi(2)).sqrt();
        r - 1737.4
    }

    /// Get current orbital speed in km/s.
    pub fn get_speed(&self) -> f64 {
        (self.state[3].powi(2) + self.state[4].powi(2) + self.state[5].powi(2)).sqrt()
    }

    /// Get Keplerian elements [sma km, ecc, inc rad, raan rad, argp rad, ta rad].
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

    const GM: f64     = DEFAULT_GM;
    const LUNAR_R: f64 = 1737.4;

    fn orbital_period(sma: f64) -> f64 {
        2.0 * PI * (sma.powi(3) / GM).sqrt()
    }

    fn specific_energy(state: &State, gm: f64) -> f64 {
        let r  = (state[0].powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
        let v2 = state[3].powi(2) + state[4].powi(2) + state[5].powi(2);
        v2 / 2.0 - gm / r
    }

    fn specific_angular_momentum(state: &State) -> f64 {
        let hx = state[1]*state[5] - state[2]*state[4];
        let hy = state[2]*state[3] - state[0]*state[5];
        let hz = state[0]*state[4] - state[1]*state[3];
        (hx.powi(2) + hy.powi(2) + hz.powi(2)).sqrt()
    }

    fn propagate_two_body(s: &State, dt: f64) -> State {
        let gm = GM;
        propagate(s, dt, |st: &State| {
            let r2 = st[0]*st[0] + st[1]*st[1] + st[2]*st[2];
            let r3 = r2 * r2.sqrt();
            let a  = -gm / r3;
            [st[3], st[4], st[5], a*st[0], a*st[1], a*st[2]]
        }).unwrap()
    }

    #[test]
    fn test_circular_orbit_period() {
        let sma    = LUNAR_R + 100.0;
        let period = orbital_period(sma);
        let s0     = keplerian_to_cartesian(sma, 0.0, 0.0, 0.0, 0.0, 0.0, GM);
        let s1     = propagate_two_body(&s0, period);
        let dr     = ((s1[0]-s0[0]).powi(2) + (s1[1]-s0[1]).powi(2) + (s1[2]-s0[2]).powi(2)).sqrt();
        let dv     = ((s1[3]-s0[3]).powi(2) + (s1[4]-s0[4]).powi(2) + (s1[5]-s0[5]).powi(2)).sqrt();
        assert!(dr < 1e-4, "Position error after one period: {:.3e} km", dr);
        assert!(dv < 1e-7, "Velocity error after one period: {:.3e} km/s", dv);
    }

    #[test]
    fn test_energy_conservation() {
        let sma = LUNAR_R + 100.0;
        let ecc = 0.04;
        let s0  = keplerian_to_cartesian(sma, ecc, 1.0, 0.5, 0.3, 0.0, GM);
        let e0  = specific_energy(&s0, GM);
        let s1  = propagate_two_body(&s0, orbital_period(sma) * 0.37);
        let e1  = specific_energy(&s1, GM);
        let rel_err = ((e1 - e0) / e0).abs();
        // DOP853 at 1e-11 tolerance — much tighter than Phase 1's 1e-8
        assert!(rel_err < 1e-10, "Energy not conserved: rel_err = {:.3e}", rel_err);
    }

    #[test]
    fn test_angular_momentum_conservation() {
        let sma = LUNAR_R + 100.0;
        let ecc = 0.04;
        let inc = PI - 0.0175;
        let s0  = keplerian_to_cartesian(sma, ecc, inc, 0.0, 0.0, 0.0, GM);
        let h0  = specific_angular_momentum(&s0);
        let s1  = propagate_two_body(&s0, orbital_period(sma) * 0.5);
        let h1  = specific_angular_momentum(&s1);
        let rel_err = ((h1 - h0) / h0).abs();
        assert!(rel_err < 1e-10, "Angular momentum not conserved: rel_err = {:.3e}", rel_err);
    }

    #[test]
    fn test_eagle_orbit_period() {
        let sma    = 1838.0;
        let ecc    = 0.04;
        let inc    = 179.0_f64.to_radians();
        let period = orbital_period(sma);
        let s0     = keplerian_to_cartesian(sma, ecc, inc, 0.0, 0.0, 0.0, GM);
        let s1     = propagate_two_body(&s0, period);
        let dr     = ((s1[0]-s0[0]).powi(2) + (s1[1]-s0[1]).powi(2) + (s1[2]-s0[2]).powi(2)).sqrt();
        assert!(dr < 1e-3, "Eagle orbit position error: {:.3e} km", dr);
        let period_min = period / 60.0;
        assert!((period_min - 118.0).abs() < 5.0, "Eagle period {:.1} min", period_min);
    }

    #[test]
    fn test_keplerian_to_cartesian_roundtrip() {
        let sma = 1838.0_f64;
        let ecc = 0.04_f64;
        let inc = 1.5708_f64;
        let ta  = 0.3_f64;
        let state = keplerian_to_cartesian(sma, ecc, inc, 0.5, 1.0, ta, GM);
        let r = (state[0].powi(2) + state[1].powi(2) + state[2].powi(2)).sqrt();
        let p = sma * (1.0 - ecc * ecc);
        let r_expected = p / (1.0 + ecc * ta.cos());
        assert!((r - r_expected).abs() < 1e-9, "Radius mismatch: got {:.6} expected {:.6}", r, r_expected);
        let energy = specific_energy(&state, GM);
        let expected_energy = -GM / (2.0 * sma);
        assert!((energy - expected_energy).abs() < 1e-6, "Energy mismatch");
    }

    #[test]
    fn test_propagate_zero_dt() {
        let s0 = keplerian_to_cartesian(1838.0, 0.04, 1.0, 0.5, 0.3, 0.7, GM);
        let s1 = propagate_two_body(&s0, 0.0);
        for i in 0..6 {
            assert_eq!(s0[i], s1[i], "State changed on zero-dt at index {}", i);
        }
    }

    #[test]
    fn test_cartesian_to_keplerian_roundtrip() {
        let sma  = 1838.0_f64;
        let ecc  = 0.04_f64;
        let inc  = 1.5_f64;
        let raan = 0.7_f64;
        let argp = 0.3_f64;
        let ta   = 1.2_f64;
        let s = keplerian_to_cartesian(sma, ecc, inc, raan, argp, ta, GM);
        let el = cartesian_to_keplerian(&s, GM);
        assert!((el[0] - sma).abs() < 1e-6, "SMA roundtrip: {} vs {}", el[0], sma);
        assert!((el[1] - ecc).abs() < 1e-10, "ECC roundtrip: {} vs {}", el[1], ecc);
        assert!((el[2] - inc).abs() < 1e-10, "INC roundtrip: {} vs {}", el[2], inc);
    }
}
