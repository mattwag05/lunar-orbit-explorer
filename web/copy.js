/**
 * Copy block (PRD 5.0).
 *
 * Class C mission facts are the only physics-bearing literals in the
 * presentation layer, and they all live here with a source. Class B figures
 * are templated from physics-constants.js; Class A figures are filled in at
 * render time from the propagator.
 */

import { MOON_RADIUS_KM, EARTH_MOON_SMA_KM, EARTH_SUN_SMA_KM } from './physics-constants.js';
import { num, words } from './format.js';

const DEG = Math.PI / 180;

// ─── Class C: mission facts ──────────────────────────────────────────────

/**
 * Apollo 11 lunar module ascent stage "Eagle" after liftoff.
 * Source: this repository's Phase 1 default orbit, an approximation of the
 * post-ascent orbit (about 1.25° to the lunar equator, retrograde). It is not
 * a reconstructed ephemeris. Retrograde is written as 180° − 0.93° here
 * because the propagator measures inclination from +Z.
 */
export const EAGLE = {
  mission: 'Apollo 11',
  name: 'Eagle',
  elements: { sma: 1838.13, ecc: 0.0076, inc: 179.07 * DEG, raan: 183.41 * DEG, argp: 179.86 * DEG, ta: 0 },
};

/**
 * Apollo 16 Particles and Fields Subsatellite, PFS-2.
 * Source set, used consistently in copy and validation:
 *   NASA NSSDCA, "Apollo 16 Subsatellite" (1972-031D),
 *   https://nssdc.gsfc.nasa.gov/nmc/spacecraft/display.action?id=1972-031D
 *   periselene 90 km, aposelene 130 km, 10° to the lunar equator, "clockwise
 *   as viewed from north" (retrograde), released 24 April 1972, impacted
 *   29 May 1972 "after 34 days (425 revolutions)".
 *   NASA Science, "Apollo 16 Subsatellite",
 *   https://science.nasa.gov/mission/apollo-16-subsatellite/
 *   "34 days in orbit rather than the planned one year".
 * Node and periapsis direction are not in either source. The assumed values
 * below come from the validation sweep in docs/validation.md, and the UI
 * says they are assumed.
 */
export const PFS2 = {
  mission: 'Apollo 16',
  name: 'PFS-2',
  year: 1972,
  periseleneKm: 90,
  aposeleneKm: 130,
  incPublishedDeg: 10,
  lifetimeDays: 34,
  revolutions: 425,
  plannedLifetime: 'a year',
  // Median case of the 16-direction sweep (impact at 34.0 days there).
  assumedRaanDeg: 180,
  assumedArgpDeg: 180,
};

/** PFS-2 in propagator elements: retrograde 10° is 170° from +Z. */
export function pfs2Elements() {
  const rp = MOON_RADIUS_KM + PFS2.periseleneKm;
  const ra = MOON_RADIUS_KM + PFS2.aposeleneKm;
  return {
    sma: (rp + ra) / 2,
    ecc: (ra - rp) / (ra + rp),
    inc: (180 - PFS2.incPublishedDeg) * DEG,
    raan: PFS2.assumedRaanDeg * DEG,
    argp: PFS2.assumedArgpDeg * DEG,
    ta: 0,
  };
}

/**
 * Inclinations where low lunar orbits stay stable ("frozen orbits").
 * Source: Science@NASA, Trudy E. Bell, "Bizarre Lunar Orbits" (2006),
 * republished at https://phys.org/news/2006-11-bizarre-lunar-orbits.html
 */
export const FROZEN_INCLINATIONS_DEG = [27, 50, 76, 86];

// ─── Act copy (PRD section 5, verbatim after templating) ─────────────────

export const APP_NAME = 'Lunar Orbit Explorer';

export const ACT_TITLES = [
  'The hook', 'Where you are', 'The orbit', 'The assumption',
  'The lumps', 'The drift', 'The other pullers', 'Try it yourself',
];

export const COPY = {
  act0: {
    eyebrow: 'NO GPS AT THE MOON',
    headline: 'The Moon has no satellites to guide you. So how did Eagle find its way home?',
  },
  act1: {
    eyebrow: 'WHERE YOU ARE',
    headline: () => `The Moon is ${num(MOON_RADIUS_KM, 0)} km in radius and has no air. Nothing up here can hear you.`,
    chip: 'Distances shown to scale. Orbit shown at 1x.',
  },
  act2: {
    eyebrow: 'THE ORBIT',
    headline: (periodMin, altKm) =>
      `Eagle circled the Moon every ${num(periodMin, 0)} minutes, ${num(altKm, 0)} km above the ground.`,
    body: 'Drag to look around it.',
  },
  act3: {
    eyebrow: 'IF THE MOON WERE SMOOTH',
    headline: 'Treat the Moon as a single point of mass and this orbit never changes. Not in a day. Not in a year.',
  },
  act4: {
    eyebrow: 'THE MOON IS NOT SMOOTH',
    headline: 'The Moon is lumpy. There is extra mass buried under the maria, and it pulls.',
    stepLabels: {
      0: 'A point of mass.',
      2: 'Two: the Moon is not round.',
      20: 'Twenty: the buried mass shows up.',
    },
    fullLabel: (degree, count) => `${degree}: ${num(count, 0)} measured coefficients.`,
    chip: (factor) => `Gravity field exaggerated ${num(factor, 0)}x so you can see it.`,
  },
  act5: {
    eyebrow: 'LET A WEEK PASS',
    headline: 'Now let the real gravity act, and watch what a week does to the same orbit.',
    toggle: ['POINT MASS', 'REAL GRAVITY'],
    checkpoint: (days, km, direction) =>
      `${days} ${days === 1 ? 'day' : 'days'} on: ${num(km, 0)} km off, mostly ${direction}`,
    conclusion: () =>
      `${PFS2.mission} released ${PFS2.name} into a low lunar orbit in ${PFS2.year}. ` +
      `It was expected to last ${PFS2.plannedLifetime}. It fell after ${PFS2.lifetimeDays} days.`,
    chipDegree: (degree, days, chunkS) =>
      `Drift computed at degree ${degree} for speed, over ${days} days in ${chunkS} s chunks.`,
  },
  act6: {
    eyebrow: 'AND TWO MORE THINGS PULL',
    headline: () =>
      `Earth is ${num(EARTH_MOON_SMA_KM, 0)} km away and the Sun is ${words(EARTH_SUN_SMA_KM)} km. Both still change this orbit.`,
    body: 'The lumps are switched off here, so Earth and Sun act alone against a point-mass Moon.',
  },
  act7: {
    eyebrow: 'NOW TRY IT YOURSELF',
    headline: 'Move the orbit. Break it. Ride along.',
    actions: ['Move the orbit', 'Break it', 'Ride along'],
    footer: 'Esc to return.',
  },
  cockpit: {
    eyebrow: 'RIDING ALONG',
    headline: 'You are aboard. The numbers below are this spacecraft, right now.',
  },
  provenance: (degree) =>
    `True simulation. Gravity from GRGM1200A (${degree}x${degree}, NASA GSFC). ` +
    'Propagated with DOP853. Not a live tracking feed.',
};
