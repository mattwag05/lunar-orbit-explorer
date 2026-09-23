/**
 * Class B defined constants (PRD 5.0).
 *
 * Each value mirrors the constant the Rust propagator actually uses, so copy
 * and physics cannot disagree. web/tests/constants.test.mjs parses the Rust
 * sources and fails if any of these drift. This module is the only place in
 * web/ where these figures appear as literals.
 */

/** Lunar GM [km³/s²]: `DEFAULT_GM` in propagator/src/lib.rs. */
export const LUNAR_GM_KM3_S2 = 4902.800066;

/** Mean lunar radius [km]: the constant in `get_altitude` in propagator/src/lib.rs. */
export const MOON_RADIUS_KM = 1737.4;

/** Earth–Moon mean distance [km]: `EARTH_MOON_SMA` in propagator/src/third_body.rs. */
export const EARTH_MOON_SMA_KM = 384_400.0;

/** Earth–Sun mean distance [km]: `EARTH_SUN_SMA` in propagator/src/third_body.rs. */
export const EARTH_SUN_SMA_KM = 1.495_978_707e8;

/** Lunar sidereal rotation rate [rad/s]: `OMEGA_MOON` in propagator/src/frames.rs. */
export const OMEGA_MOON_RAD_S = 2.661_709_79e-6;

/**
 * Mean surface gravity [km/s²] from the two constants above. Used only to
 * express a returned anomaly as a fraction of surface gravity when drawing
 * Act 4's exaggerated relief.
 */
export const SURFACE_GRAVITY_KM_S2 = LUNAR_GM_KM3_S2 / (MOON_RADIUS_KM * MOON_RADIUS_KM);

/** 1 mGal in km/s². */
export const MGAL_KM_S2 = 1e-8;
