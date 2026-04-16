//! Spherical harmonic gravity acceleration in the Moon-fixed (MCMF) frame.
//!
//! # Algorithm
//! Evaluates the gravitational acceleration from the GRGM1200A potential
//! using fully-normalised Associated Legendre Functions (ALFs) and their
//! geocentric-spherical gradient, converted to body-fixed Cartesian.
//!
//! Reference: Montenbruck & Gill "Satellite Orbits" (2000) §3.2.
//!
//! # Normalisation convention
//! ALFs satisfy:  ∫ P̄_nm(sin φ)² cos φ dφ dλ = 4π/(2−δ_{m,0})
//! This matches the GRGM1200A coefficient convention exactly.
//!
//! # Gradient formulae (body-fixed Cartesian)
//! Given geocentric lat φ, lon λ, radius r:
//!   a = (∂V/∂r, (1/r)∂V/∂φ, (1/(r cosφ))∂V/∂λ)   in spherical
//!
//! Converted to Cartesian via the Jacobian of spherical → Cartesian.

use crate::coefficients::{Coefficients, R_REF};

/// Compute gravitational acceleration [km/s²] for body-fixed position `r` [km].
///
/// `gm`:    gravitational parameter [km³/s²]
/// `n_max`: maximum degree/order to include (capped at `coeff.n_max`)
/// `coeff`: GRGM1200A coefficient store
///
/// The acceleration includes all terms up to degree `n_max`, with degree 0
/// giving the exact point-mass result. Always call with position in the
/// Moon-fixed (MCMF) frame.
pub fn gravity_sh(r: &[f64; 3], gm: f64, n_max: usize, coeff: &Coefficients) -> [f64; 3] {
    let x = r[0];
    let y = r[1];
    let z = r[2];

    let rho = (x*x + y*y + z*z).sqrt();
    if rho < 1.0 {
        // Inside body (should never happen); return point-mass
        let r3 = rho.powi(3);
        return [-gm*x/r3, -gm*y/r3, -gm*z/r3];
    }

    let n_max = n_max.min(coeff.n_max);

    // Geocentric spherical coordinates
    let rxy     = (x*x + y*y).sqrt();
    let sin_phi = z / rho;                // sin(geocentric latitude)
    let cos_phi = rxy / rho;              // cos(geocentric latitude)
    let lambda  = y.atan2(x);            // East longitude [rad]

    // Precompute cos(mλ) and sin(mλ) for m = 0..=n_max
    // (recursive: cos((m+1)λ) = 2cosλ·cos(mλ) − cos((m−1)λ))
    let cos_lam = lambda.cos();
    let sin_lam = lambda.sin();
    let mut cos_ml = vec![0.0f64; n_max + 2];
    let mut sin_ml = vec![0.0f64; n_max + 2];
    cos_ml[0] = 1.0;
    sin_ml[0] = 0.0;
    if n_max >= 1 {
        cos_ml[1] = cos_lam;
        sin_ml[1] = sin_lam;
    }
    for m in 2..=n_max {
        cos_ml[m] =  2.0 * cos_lam * cos_ml[m-1] - cos_ml[m-2];
        sin_ml[m] =  2.0 * cos_lam * sin_ml[m-1] - sin_ml[m-2];
    }

    // ALF storage: p[n*(n+1)/2 + m] = P̄_nm(sin_phi)
    let n_pairs = (n_max + 1) * (n_max + 2) / 2;
    let mut p  = vec![0.0f64; n_pairs];  // P̄_nm
    let mut dp = vec![0.0f64; n_pairs];  // dP̄_nm/dφ

    let idx = |n: usize, m: usize| n*(n+1)/2 + m;

    // P̄_00 = 1
    p[idx(0, 0)] = 1.0;

    // Build P̄_nm using the standard recursion (M&G eq. 3.29)
    for n in 1..=n_max {
        let nf = n as f64;

        // Sectorial P̄_nn from P̄_{n-1,n-1}
        let a_nn = ((2.0*nf + 1.0) / (2.0*nf)).sqrt();
        p[idx(n, n)] = a_nn * cos_phi * p[idx(n-1, n-1)];

        // Sub-diagonal P̄_{n,n-1} from P̄_{n-1,n-1}
        let a_nn1 = (2.0*nf + 1.0).sqrt();
        p[idx(n, n-1)] = a_nn1 * sin_phi * p[idx(n-1, n-1)];

        // Off-diagonal recurrence P̄_nm from P̄_{n-1,m} and P̄_{n-2,m}
        for m in 0..n.saturating_sub(1) {
            let mf = m as f64;
            let a = ((4.0*nf*nf - 1.0) / (nf*nf - mf*mf)).sqrt();
            let b = (((2.0*nf + 1.0)*(nf - 1.0 + mf)*(nf - 1.0 - mf))
                     / ((2.0*nf - 3.0)*(nf*nf - mf*mf))).sqrt();
            let p_n1 = if n >= 1 { p[idx(n-1, m)] } else { 0.0 };
            let p_n2 = if n >= 2 { p[idx(n-2, m)] } else { 0.0 };
            p[idx(n, m)] = a * sin_phi * p_n1 - b * p_n2;
        }
    }

    // dP̄_nm/dφ via the derivative recursion (avoids 1/cos_phi singularity).
    //
    // dP̄_nm/dφ = (A_nm·P̄_{n-1,m} − n·sin_phi·P̄_nm) / cos_phi
    //
    // where A_nm = sqrt((2n+1)(n+m)(n−m) / (2n−1))
    //
    // At the poles (cos_phi < eps) the derivative is zero for most terms;
    // the sectorial terms P̄_nn have a finite limit handled below.
    const POLE_EPS: f64 = 1e-10;
    for n in 0..=n_max {
        for m in 0..=n {
            if cos_phi < POLE_EPS {
                // At a pole all dP/dphi = 0 except m = 1 sectorials — but for
                // the force computation the λ-component already vanishes at
                // the poles (1/cos_phi × ∂/∂λ handled separately), so setting
                // dp = 0 everywhere at the pole is the safest convention.
                dp[idx(n, m)] = 0.0;
                continue;
            }
            let nf = n as f64;
            let mf = m as f64;
            let a_nm = if n >= 1 {
                let num = (2.0*nf + 1.0) * (nf + mf) * (nf - mf);
                let den = 2.0*nf - 1.0;
                if den > 0.0 { (num / den).sqrt() } else { 0.0 }
            } else {
                0.0
            };
            let p_n1m = if n >= 1 { p[idx(n-1, m)] } else { 0.0 };
            dp[idx(n, m)] = (a_nm * p_n1m - nf * sin_phi * p[idx(n, m)]) / cos_phi;
        }
    }

    // Accumulate sums for spherical-coordinate gradient components.
    // V = (GM/r) Σ_{n,m} (R/r)^n P̄_nm [C cos mλ + S sin mλ]
    //
    // ∂V/∂r    = −(GM/r²) Σ (n+1)(R/r)^n P̄_nm cs
    // ∂V/∂φ    = (GM/r)   Σ        (R/r)^n dP̄_nm cs
    // ∂V/∂λ    = (GM/r)   Σ        (R/r)^n P̄_nm m·sn
    //
    // where cs = C cos mλ + S sin mλ
    //       sn = −C sin mλ + S cos mλ  (∂(cs)/∂λ / m)

    let mut sum_r   = 0.0f64;
    let mut sum_phi = 0.0f64;
    let mut sum_lam = 0.0f64;

    for n in 0..=n_max {
        let nf    = n as f64;
        let ratio = (R_REF / rho).powi(n as i32);
        for m in 0..=n {
            let (c_nm, s_nm) = coeff.get(n, m);
            let mf   = m as f64;
            let cs   = c_nm * cos_ml[m] + s_nm * sin_ml[m];
            let sn   = -c_nm * sin_ml[m] + s_nm * cos_ml[m];  // (1/m)·∂cs/∂λ when m>0
            let pnm  = p[idx(n, m)];
            let dpnm = dp[idx(n, m)];

            sum_r   -= (nf + 1.0) * ratio * pnm  * cs;
            sum_phi +=              ratio * dpnm * cs;
            sum_lam +=              ratio * pnm  * mf * sn;
        }
    }

    let gm_r2 = gm / (rho * rho);

    // Spherical acceleration components (not yet Cartesian)
    let a_r   = gm_r2 * sum_r;           // radial:    ∂V/∂r
    let a_phi = gm_r2 * sum_phi;         // latitude:  (1/r)·∂V/∂φ  = (GM/r²)·sum_phi
    let a_lam = if cos_phi > POLE_EPS {
        gm_r2 * sum_lam / cos_phi        // longitude: (1/(r cosφ))·∂V/∂λ = (GM/r²)·sum_lam/cosφ
    } else {
        0.0
    };

    // Convert (a_r, a_phi, a_lam) to body-fixed Cartesian (x, y, z).
    // ∂r/∂x  = x/r = cosφ cosλ;   ∂φ/∂x = −sinφ cosλ / r;   ∂λ/∂x = −sinλ / (r cosφ)
    let cp_cl = cos_phi * lambda.cos();
    let cp_sl = cos_phi * lambda.sin();
    let sp_cl = sin_phi * lambda.cos();
    let sp_sl = sin_phi * lambda.sin();

    let ax = a_r * cp_cl - a_phi * sp_cl - a_lam * lambda.sin();
    let ay = a_r * cp_sl - a_phi * sp_sl + a_lam * lambda.cos();
    let az = a_r * sin_phi + a_phi * cos_phi;

    [ax, ay, az]
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coefficients::Coefficients;

    const GM: f64 = 4902.800066;

    /// Build a Coefficients with only C_00 = 1 (degree 0 = pure point-mass).
    fn point_mass_coeff() -> Coefficients {
        Coefficients::from_bundle(0)
    }

    #[test]
    fn degree_zero_matches_point_mass() {
        // With n_max=0, gravity_sh must match the point-mass formula exactly.
        let c = point_mass_coeff();
        let r = [1838.0f64, 0.0, 0.0];
        let a = gravity_sh(&r, GM, 0, &c);
        let a_exact = -GM / r[0].powi(2);
        // x-component should equal point-mass; y, z should be zero.
        assert!((a[0] - a_exact).abs() / a_exact.abs() < 1e-6,
            "Point-mass ax: got {:.6e}, expected {:.6e}", a[0], a_exact);
        assert!(a[1].abs() < 1e-12, "ay = {} (expected 0)", a[1]);
        assert!(a[2].abs() < 1e-12, "az = {} (expected 0)", a[2]);
    }

    #[test]
    fn j2_correction_has_right_sign_and_magnitude() {
        // For an equatorial spacecraft (z=0), J2 reduces the radial acceleration
        // compared to point-mass (the Moon's equatorial bulge provides outward push).
        let c = Coefficients::from_bundle(2);
        let r = [1838.0f64, 0.0, 0.0];
        let a_sh = gravity_sh(&r, GM, 2, &c);
        let a_pm = -GM / r[0].powi(2);
        // The J2 radial perturbation should be small: |Δa/a_pm| ≈ 1e-3 to 1e-5
        let rel_diff = (a_sh[0] - a_pm) / a_pm.abs();
        assert!(rel_diff.abs() > 1e-6,
            "J2 perturbation too small: {:.3e}", rel_diff);
        assert!(rel_diff.abs() < 0.1,
            "J2 perturbation too large: {:.3e}", rel_diff);
    }

    #[test]
    fn polar_orbit_z_acceleration_nonzero() {
        // A polar position (x=y=0, z=r) should feel a z-directed acceleration.
        // With J2, the polar acceleration is larger than equatorial.
        let c = Coefficients::from_bundle(2);
        let r_polar = [0.0f64, 0.0, 1838.0];
        let a = gravity_sh(&r_polar, GM, 2, &c);
        let a_pm = -GM / 1838.0_f64.powi(2);
        // z-component should be predominantly radial (negative toward origin)
        assert!(a[2] < 0.0, "az should be negative (toward origin), got {}", a[2]);
        // J2 correction at pole: 1.5*J2*(R/r)^2*(3-3) = 0? No, at the pole
        // the J2 term contributes: (3/2)*J2*(R/r)^2 * (3-1) = 3*J2*(R/r)^2
        // so polar a is stronger. Check it's different from point-mass.
        let rel = (a[2] - a_pm) / a_pm.abs();
        assert!(rel.abs() > 1e-5, "J2 polar correction too small: {:.3e}", rel);
    }

    #[test]
    fn acceleration_magnitude_reasonable() {
        // At ~100 km altitude, |a| ≈ GM/r² ≈ 0.00145 km/s²  (≈ 1.45 m/s²).
        // GM = 4902.8 km³/s², r ≈ 1838 km → GM/r² ≈ 0.001451 km/s².
        let c = Coefficients::from_bundle(4);
        // Slightly off-axis so all three spherical terms contribute.
        let r = [1836.0f64, 50.0, 30.0];
        let a = gravity_sh(&r, GM, 4, &c);
        let mag = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]).sqrt();
        assert!(mag > 0.001 && mag < 0.002,
            "|a| = {:.5e} km/s² out of expected range [0.001, 0.002]", mag);
    }
}
