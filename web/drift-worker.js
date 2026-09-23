/**
 * Act 5 precompute (PRD 5.6).
 *
 * Two propagators from the same initial state, one point mass and one at the
 * working degree, stepped in lockstep in fixed chunks. Progress messages count
 * real completed chunks. The act replays the returned table; nothing here
 * runs in the frame loop.
 */

import init, { Propagator } from '../propagator/pkg/propagator.js';

self.onmessage = async (e) => {
  const { gm, elements: el, degree, days, chunkS } = e.data;
  await init();
  const started = performance.now();

  const make = (deg) => {
    const p = new Propagator();
    p.init(gm);
    p.init_from_keplerian(el.sma, el.ecc, el.inc, el.raan, el.argp, el.ta);
    p.set_gravity_degree(deg);
    return p;
  };
  const pointMass = make(0);
  const real = make(degree);

  const steps = Math.round((days * 86400) / chunkS);
  const times = new Float64Array(steps + 1);
  const a = new Float64Array((steps + 1) * 6);
  const b = new Float64Array((steps + 1) * 6);
  a.set(pointMass.get_state(), 0);
  b.set(real.get_state(), 0);

  const reportEvery = Math.max(1, Math.floor(steps / 100));
  for (let k = 1; k <= steps; k++) {
    pointMass.step(chunkS);
    real.step(chunkS);
    times[k] = real.get_time();
    a.set(pointMass.get_state(), k * 6);
    b.set(real.get_state(), k * 6);
    if (k % reportEvery === 0 || k === steps) {
      postMessage({ type: 'progress', done: k, total: steps });
    }
  }
  const loadedDegree = real.get_loaded_degree();
  pointMass.free();
  real.free();

  postMessage(
    { type: 'done', times, a, b, degree: loadedDegree >= degree ? degree : loadedDegree,
      days, chunkS, elapsedMs: performance.now() - started },
    [times.buffer, a.buffer, b.buffer],
  );
};
