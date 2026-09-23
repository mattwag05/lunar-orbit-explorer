/**
 * Separation between two returned states, split by direction (PRD 5.6).
 *
 * Geometry only: every input is a state the propagator produced. Separations
 * reach hundreds of km on an 1,838 km orbit, so a flat projection onto the
 * reference velocity would be distorted by curvature. The split is
 * curvilinear instead, matching propagator/examples/validate.rs:
 *   height      |r_b| − |r_a|  (difference in altitude above mean radius)
 *   along-track  angle from r_a to r_b in a's orbit plane, times |r_a|,
 *                signed positive in a's direction of motion
 *   cross-track  r_b's offset out of a's orbit plane
 */

const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a) => Math.sqrt(dot(a, a));
const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

/**
 * @param {ArrayLike<number>} a reference state [x,y,z km, vx,vy,vz km/s]
 * @param {ArrayLike<number>} b compared state, same layout
 * @returns {{total:number,height:number,along:number,cross:number,dominant:'height'|'along-track'|'cross-track'}}
 */
export function separation(a, b) {
  const ra = [a[0], a[1], a[2]];
  const va = [a[3], a[4], a[5]];
  const rb = [b[0], b[1], b[2]];
  const n = cross(ra, va);
  const nn = norm(n);
  const nh = [n[0] / nn, n[1] / nn, n[2] / nn];

  const total = norm([rb[0] - ra[0], rb[1] - ra[1], rb[2] - ra[2]]);
  const height = norm(rb) - norm(ra);
  const crossTrack = dot(rb, nh);
  const rbIn = [rb[0] - crossTrack * nh[0], rb[1] - crossTrack * nh[1], rb[2] - crossTrack * nh[2]];
  const angle = Math.atan2(dot(cross(ra, rbIn), nh), dot(ra, rbIn));
  const along = angle * norm(ra);

  const parts = [['height', Math.abs(height)], ['along-track', Math.abs(along)], ['cross-track', Math.abs(crossTrack)]];
  parts.sort((p, q) => q[1] - p[1]);
  return { total, height, along, cross: crossTrack, dominant: parts[0][0] };
}
