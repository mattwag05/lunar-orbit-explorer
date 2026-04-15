//! Moon-Centred Inertial (MCI) ↔ Moon-Centred Moon-Fixed (MCMF) transforms.
//!
//! Model: uniform rotation about +Z at the IAU 2009 sidereal rate.
//!
//! Convention:
//!   r_MCMF = R_z(-ω·t) · r_MCI          (position: rotate backward by elapsed angle)
//!   v_MCMF = R_z(-ω·t) · v_MCI - ω⃗ × r_MCMF   (velocity: remove frame spin)
//!
//! Inverse (MCMF → MCI):
//!   r_MCI  = R_z(+ω·t) · r_MCMF
//!   v_MCI  = R_z(+ω·t) · (v_MCMF + ω⃗ × r_MCMF)

use crate::integrator::State;

/// Moon sidereal rotation rate [rad/s].
/// 13.1763582°/day = 2π/(27.3217×86400) rad/s  (IAU 2009).
pub const OMEGA_MOON: f64 = 2.661_709_79e-6;

/// Rotate 3-vector by R_z(theta): x' = c·x − s·y, y' = s·x + c·y, z' = z.
#[inline(always)]
fn rz(v: [f64; 3], theta: f64) -> [f64; 3] {
    let c = theta.cos();
    let s = theta.sin();
    [c*v[0] - s*v[1],
     s*v[0] + c*v[1],
     v[2]]
}

/// MCI state → MCMF state at elapsed time `t` [s] since epoch.
///
/// Steps:
///   1. θ = -ω·t  (body-fixed frame has rotated +ω·t, so we rotate the vector -ω·t)
///   2. r_mcmf = R_z(θ) · r_mci
///   3. v_mcmf = R_z(θ) · v_mci  −  ω⃗ × r_mcmf
///      where ω⃗ × r = ω·[−ry, rx, 0]  (ω⃗ = [0,0,ω])
pub fn mci_to_mcmf(s: &State, t: f64) -> State {
    let theta = -OMEGA_MOON * t;
    let r_mcmf = rz([s[0], s[1], s[2]], theta);
    let v_rot  = rz([s[3], s[4], s[5]], theta);
    // Subtract frame spin: ω × r_mcmf = [-ω·ry, ω·rx, 0]
    let v_mcmf = [
        v_rot[0] - (-OMEGA_MOON * r_mcmf[1]),
        v_rot[1] - ( OMEGA_MOON * r_mcmf[0]),
        v_rot[2],
    ];
    [r_mcmf[0], r_mcmf[1], r_mcmf[2],
     v_mcmf[0], v_mcmf[1], v_mcmf[2]]
}

/// MCMF state → MCI state at elapsed time `t` [s] since epoch.
pub fn mcmf_to_mci(s: &State, t: f64) -> State {
    let theta = OMEGA_MOON * t;
    // Add frame spin back before rotating to inertial
    let v_with_spin = [
        s[3] + (-OMEGA_MOON * s[1]),
        s[4] + ( OMEGA_MOON * s[0]),
        s[5],
    ];
    let r_mci = rz([s[0], s[1], s[2]], theta);
    let v_mci = rz(v_with_spin, theta);
    [r_mci[0], r_mci[1], r_mci[2],
     v_mci[0], v_mci[1], v_mci[2]]
}

/// Rotate only the position part of a state from MCMF to MCI at time `t`.
/// Used when rotating an acceleration vector back to inertial frame.
#[inline(always)]
pub fn rotate_vec_mcmf_to_mci(v: [f64; 3], t: f64) -> [f64; 3] {
    rz(v, OMEGA_MOON * t)
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_identity() {
        // Forward then inverse must return the original state exactly.
        let s0: State = [1838.0, 100.0, 50.0, 0.1, -1.63, 0.05];
        let t = 3600.0 * 24.0 * 7.0; // 7 days
        let s_mcmf = mci_to_mcmf(&s0, t);
        let s_back = mcmf_to_mci(&s_mcmf, t);
        for i in 0..6 {
            assert!(
                (s_back[i] - s0[i]).abs() < 1e-10,
                "Round-trip error at index {}: {} vs {}", i, s_back[i], s0[i]
            );
        }
    }

    #[test]
    fn position_magnitude_preserved() {
        // Rotation preserves vector magnitude.
        let s0: State = [1838.0, 200.0, -100.0, 0.0, 1.63, 0.0];
        let r0 = (s0[0].powi(2) + s0[1].powi(2) + s0[2].powi(2)).sqrt();
        for &t in &[0.0, 3600.0, 86400.0, 27.3217 * 86400.0] {
            let s1 = mci_to_mcmf(&s0, t);
            let r1 = (s1[0].powi(2) + s1[1].powi(2) + s1[2].powi(2)).sqrt();
            assert!((r1 - r0).abs() < 1e-9,
                "Position magnitude not preserved at t={}: {} vs {}", t, r1, r0);
        }
    }

    #[test]
    fn one_sidereal_period_returns_to_start() {
        // After exactly one sidereal period θ = 2π → R_z(-2π) = I → round-trip is identity.
        let t_sid = 2.0 * std::f64::consts::PI / OMEGA_MOON; // ≈ 27.32 days
        let s0: State = [1838.0, 0.0, 0.0, 0.0, 1.63, 0.0];
        let s1 = mci_to_mcmf(&s0, t_sid);
        for i in 0..3 {
            assert!((s1[i] - s0[i]).abs() < 1e-7,
                "Position not preserved after one period at index {}", i);
        }
    }
}
