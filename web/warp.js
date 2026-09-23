/**
 * Time-warp ladders. Every step carries its own label and multiplier so a
 * label can never drift from its rate (PRD 5.9).
 */

/** Cockpit ladder, PRD 5.9. Replaces the old WARP_LEVELS array. */
export const COCKPIT_WARP = [
  { label: 'Pause',    mult: 0 },
  { label: 'Real',     mult: 1 },
  { label: '1 min/s',  mult: 60 },
  { label: '4 min/s',  mult: 240 },
  { label: '20 min/s', mult: 1200 },
  { label: '1 h/s',    mult: 3600 },
];

/** Rates the guided acts use. Same shape, disclosed the same way. */
export const ACT_WARP = {
  real:     { label: 'Real',     mult: 1 },
  fourMin:  { label: '4 min/s',  mult: 240 },
  twentyMin:{ label: '20 min/s', mult: 1200 },
  sixHour:  { label: '6 h/s',    mult: 21600 },
  day:      { label: '1 day/s',  mult: 86400 },
  halfDay:  { label: '12 h/s',   mult: 43200 },
};

const grouped = new Intl.NumberFormat('en-US', { maximumFractionDigits: 0 });

/**
 * Honesty chip for a warp step, or null when time is not distorted.
 * @param {{label:string,mult:number}} step
 * @param {number} [achieved] measured sim-seconds per wall-second, when the
 *   propagator cannot keep up with the requested rate
 */
export function warpChip(step, achieved) {
  if (!step || step.mult === 0 || step.mult === 1) return null;
  if (achieved !== undefined && achieved < step.mult * 0.9) {
    return `Time shown about ${grouped.format(achieved)}x faster than real. ` +
      `Asked for ${step.label}; the propagator is running as fast as it can.`;
  }
  return `Time shown ${grouped.format(step.mult)}x faster than real (${step.label}).`;
}
