//! Third-body point-mass perturbations: Earth and Sun.
//!
//! Both perturbing bodies use simplified Keplerian circular-orbit ephemerides
//! in the Moon-Centred Inertial (MCI) frame. The inclinations match the mean
//! orbital planes relative to the Moon's equatorial plane.
//!
//! Earth: circular, a = 384400 km, T = 27.3217 days, i ≈ 6.68° to Moon equator
//! Sun:   the Moon orbits the Earth-Moon barycentre ~1 AU from the Sun; for
//!        the indirect perturbation we use the Sun's direction as opposite to
//!        the Earth-Moon direction projected onto the ecliptic.
//!
//! Reference: Battin "An Introduction to the Mathematics and Methods of
//!            Astrodynamics" (1987), §9.2; and Meador (1994) for the 24-day
//!            eccentricity oscillation induced by Earth + J2.

/// Earth GM [km³/s²]  (IAU 2009)
pub const GM_EARTH: f64 = 398_600.4418;

/// Sun GM [km³/s²]
pub const GM_SUN: f64 = 1.327_124_400e11;

/// Earth–Moon mean semi-major axis [km]
pub const EARTH_MOON_SMA: f64 = 384_400.0;

/// Earth–Moon sidereal period [s]  (27.3217 days)
pub const EARTH_MOON_T: f64 = 27.3217 * 86400.0;

/// Earth orbit inclination to Moon equatorial plane [rad]  (≈ 6.68°)
pub const EARTH_MOON_INC: f64 = 0.116_587; // rad

/// Earth–Sun mean semi-major axis (Moon from Sun) [km]  (1 AU)
pub const EARTH_SUN_SMA: f64 = 1.495_978_707e8;

/// Earth–Sun sidereal period [s]  (365.25 days)
pub const EARTH_SUN_T: f64 = 365.25 * 86400.0;

/// Sun direction inclination to Moon equatorial plane [rad]  (≈ 1.54°)
pub const SUN_INC: f64 = 0.026_878; // rad

// ─── Simplified analytical ephemerides ───────────────────────────────────────

/// Position of Earth in MCI frame at elapsed time `t` [s] since epoch [km].
///
/// Circular Keplerian orbit. At t=0 Earth is on the +x axis.
#[inline]
pub fn earth_mci(t: f64) -> [f64; 3] {
    let n     = std::f64::consts::TAU / EARTH_MOON_T;
    let theta = n * t;
    let ci    = EARTH_MOON_INC.cos();
    let si    = EARTH_MOON_INC.sin();
    let a     = EARTH_MOON_SMA;
    [
         a * theta.cos(),
         a * theta.sin() * ci,
         a * theta.sin() * si,
    ]
}

/// Position of Sun in MCI frame at elapsed time `t` [s] since epoch [km].
///
/// Treats the Moon as orbiting the Earth–Moon barycentre (≈ Earth) which
/// itself orbits the Sun. The Sun direction is opposite to the Earth direction
/// in the heliocentric frame, projected for the Moon.
/// Approximation: Sun position ≈ −(1 AU) × (unit vector toward Earth in
/// ecliptic plane, rotated for Sun's annual motion).
#[inline]
pub fn sun_mci(t: f64) -> [f64; 3] {
    let n     = std::f64::consts::TAU / EARTH_SUN_T;
    let theta = n * t;
    let ci    = SUN_INC.cos();
    let si    = SUN_INC.sin();
    let a     = EARTH_SUN_SMA;
    // Sun is in the opposite direction of the Earth's annual motion as seen
    // from the Moon; we place it at an initial phase of π to separate it
    // from the Earth direction at t=0.
    [
         a * (theta + std::f64::consts::PI).cos(),
         a * (theta + std::f64::consts::PI).sin() * ci,
         a * (theta + std::f64::consts::PI).sin() * si,
    ]
}

// ─── Point-mass perturbation acceleration ────────────────────────────────────

/// Third-body indirect + direct point-mass acceleration on a spacecraft.
///
/// Battin formula:  a = GM₃ × ( (r₃ − r_sc)/|r₃ − r_sc|³ − r₃/|r₃|³ )
///
/// `r_sc`: spacecraft position in MCI [km]
/// `r_3`:  perturbing body position in MCI [km]
/// `gm_3`: GM of perturbing body [km³/s²]
#[inline]
pub fn third_body_accel(r_sc: &[f64; 3], r_3: &[f64; 3], gm_3: f64) -> [f64; 3] {
    let d = [
        r_3[0] - r_sc[0],
        r_3[1] - r_sc[1],
        r_3[2] - r_sc[2],
    ];
    let d2     = d[0]*d[0] + d[1]*d[1] + d[2]*d[2];
    let d_mag3 = d2 * d2.sqrt();

    let r3_sq  = r_3[0]*r_3[0] + r_3[1]*r_3[1] + r_3[2]*r_3[2];
    let r3_mag3 = r3_sq * r3_sq.sqrt();

    [
        gm_3 * (d[0] / d_mag3 - r_3[0] / r3_mag3),
        gm_3 * (d[1] / d_mag3 - r_3[1] / r3_mag3),
        gm_3 * (d[2] / d_mag3 - r_3[2] / r3_mag3),
    ]
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_at_correct_sma() {
        let r = earth_mci(0.0);
        let dist = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
        assert!((dist - EARTH_MOON_SMA).abs() < 1.0,
            "Earth distance {} km != expected {} km", dist, EARTH_MOON_SMA);
    }

    #[test]
    fn sun_at_correct_sma() {
        let r = sun_mci(0.0);
        let dist = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
        assert!((dist - EARTH_SUN_SMA).abs() / EARTH_SUN_SMA < 0.01,
            "Sun distance {:.3e} km off from 1 AU", dist);
    }

    #[test]
    fn earth_perturbation_order_of_magnitude() {
        // Earth's tidal acceleration on a 100 km lunar orbit:
        //   a ≈ 2·GM_E·r_sc / d_E³ ≈ 2 × 398600 × 1838 / 384400³
        //     ≈ 2.58×10⁻⁷ km/s²  (literature: ~3×10⁻⁷ km/s² typical)
        let r_sc = [1838.0f64, 0.0, 0.0];
        let r_e  = earth_mci(0.0);
        let a    = third_body_accel(&r_sc, &r_e, GM_EARTH);
        let mag  = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]).sqrt();
        assert!(mag > 1e-9 && mag < 1e-4,
            "Earth third-body |a| = {:.3e} km/s² out of range [1e-9, 1e-4]", mag);
    }

    #[test]
    fn sun_perturbation_order_of_magnitude() {
        // Sun is ≈ 27× weaker than Earth for lunar orbiters.
        let r_sc  = [1838.0f64, 0.0, 0.0];
        let r_sun = sun_mci(0.0);
        let a     = third_body_accel(&r_sc, &r_sun, GM_SUN);
        let mag   = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]).sqrt();
        assert!(mag > 1e-12 && mag < 1e-5,
            "Sun third-body |a| = {:.3e} km/s² out of range [1e-12, 1e-5]", mag);
    }

    #[test]
    fn earth_sun_not_collinear_at_t0() {
        // At t=0, Earth is at +x; Sun is at −x. They should be in opposite
        // hemispheres (dot product of unit vectors < −0.5).
        let re = earth_mci(0.0);
        let rs = sun_mci(0.0);
        let re_mag = (re[0]*re[0]+re[1]*re[1]+re[2]*re[2]).sqrt();
        let rs_mag = (rs[0]*rs[0]+rs[1]*rs[1]+rs[2]*rs[2]).sqrt();
        let dot = (re[0]*rs[0]+re[1]*rs[1]+re[2]*rs[2]) / (re_mag * rs_mag);
        assert!(dot < -0.5, "Earth and Sun should be in opposite hemispheres; dot = {}", dot);
    }
}
