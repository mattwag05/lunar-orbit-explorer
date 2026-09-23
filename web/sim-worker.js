/**
 * Live simulation worker.
 *
 * Owns the running propagator so the frame loop never calls it (PRD 9,
 * Truthfulness). Steps at the requested warp within a per-tick compute
 * budget, and reports the rate it actually achieved so the UI can disclose
 * it. Also measures revolutions and period from the simulated motion itself.
 */

import init, { Propagator } from '../propagator/pkg/propagator.js';

const TICK_MS = 16;
const BUDGET_MS = 13;
const TRAIL_EVERY_S = 60;
const FULL_DEGREE_REQUEST = 0xffffffff; // propagator clamps to what it holds

let prop = null;
let shadow = null;           // same start, third bodies off (Act 6 comparison)
let gm = 0;
let cfg = { warp: 1, degree: 0, earth: false, sun: false, shadow: false };
let running = false;
let impacted = false;
let lastWall = 0;
let pendingTrail = [];
let achieved = 0;
let rev = null;              // revolution counter state
let lastPeriod = NaN;
let nextSample = TRAIL_EVERY_S; // sim time of the next evenly spaced trail point
let gen = 0;                  // reset generation, echoed so stale frames are ignored

function applyForces(p, withThirdBody) {
  p.set_gravity_degree(cfg.degree);
  p.enable_third_body(withThirdBody && cfg.earth, withThirdBody && cfg.sun);
}

function makeProp(el) {
  const p = new Propagator();
  p.init(gm);
  p.init_from_keplerian(el.sma, el.ecc, el.inc, el.raan, el.argp, el.ta);
  return p;
}

function planeAngleTracker(state) {
  // Cumulative angle swept around the orbit normal since the tracker started.
  return { prev: [state[0], state[1], state[2]], tPrev: 0, swept: 0, lastCross: 0, count: 0 };
}

function trackRevolution(state, t) {
  const r = [state[0], state[1], state[2]];
  const v = [state[3], state[4], state[5]];
  const h = [r[1] * v[2] - r[2] * v[1], r[2] * v[0] - r[0] * v[2], r[0] * v[1] - r[1] * v[0]];
  const p = rev.prev;
  const c = [p[1] * r[2] - p[2] * r[1], p[2] * r[0] - p[0] * r[2], p[0] * r[1] - p[1] * r[0]];
  const hn = Math.hypot(h[0], h[1], h[2]);
  const d = Math.atan2((c[0] * h[0] + c[1] * h[1] + c[2] * h[2]) / hn,
                       p[0] * r[0] + p[1] * r[1] + p[2] * r[2]);
  const before = rev.swept;
  rev.swept += d;
  const turns = Math.floor(rev.swept / (2 * Math.PI));
  if (turns > rev.count) {
    const frac = d > 0 ? (turns * 2 * Math.PI - before) / d : 1;
    const tCross = rev.tPrev + frac * (t - rev.tPrev);
    lastPeriod = tCross - rev.lastCross;
    rev.lastCross = tCross;
    rev.count = turns;
  }
  rev.prev = r;
  rev.tPrev = t;
}

/** One revolution from epoch: measured period and mean altitude (Act 2). */
function measureEpoch(el) {
  const p = makeProp(el);
  p.set_gravity_degree(cfg.degree);
  const saved = [rev, lastPeriod];
  rev = planeAngleTracker(p.get_state());
  let altSum = 0;
  let n = 0;
  const dt = 20;
  while (rev.count < 1 && p.get_time() < 86400) {
    p.step(dt);
    trackRevolution(p.get_state(), p.get_time());
    altSum += p.get_altitude();
    n++;
  }
  const period = rev.lastCross;
  p.free();
  [rev, lastPeriod] = saved;
  return { periodS: period, meanAltKm: altSum / n };
}

function reset(el, measure) {
  prop?.free();
  shadow?.free();
  prop = makeProp(el);
  applyForces(prop, true);
  shadow = null;
  if (cfg.shadow) {
    shadow = makeProp(el);
    applyForces(shadow, false);
  }
  impacted = false;
  pendingTrail = [];
  const s = prop.get_state();
  rev = planeAngleTracker(s);
  lastPeriod = NaN;
  nextSample = TRAIL_EVERY_S;
  const e0 = prop.get_orbital_elements()[1];
  const msg = { type: 'reset', gen, ecc0: e0, state: Array.from(s) };
  if (measure) msg.epoch = measureEpoch(el);
  postMessage(msg);
}

function maxChunk() {
  return Math.min(300, Math.max(TRAIL_EVERY_S, cfg.warp / 20));
}

function tick() {
  const now = performance.now();
  const wall = Math.min((now - lastWall) / 1000, 0.25);
  lastWall = now;

  if (prop && running && !impacted && cfg.warp > 0) {
    let remaining = wall * cfg.warp;
    const want = remaining;
    const start = performance.now();
    const chunk = maxChunk();
    while (remaining > 1e-9) {
      if (performance.now() - start > BUDGET_MS) break;
      const t = prop.get_time();
      // Land exactly on trail sample times so the trail is evenly spaced.
      const dt = Math.min(remaining, chunk, Math.max(nextSample - t, 1e-3));
      if (!prop.step(dt)) { impacted = true; break; }
      shadow?.step(dt);
      remaining -= dt;
      const s = prop.get_state();
      const tn = prop.get_time();
      trackRevolution(s, tn);
      if (tn >= nextSample - 1e-6) {
        pendingTrail.push(s[0], s[1], s[2]);
        nextSample += TRAIL_EVERY_S;
      }
      if (prop.get_altitude() <= 0) {
        impacted = true;
        postMessage({ type: 'impact', gen, t: tn, revolutions: rev.count });
        break;
      }
    }
    const done = want - remaining;
    achieved = wall > 0 ? done / wall : cfg.warp;
  } else {
    achieved = 0;
  }

  if (prop) {
    const msg = {
      type: 'frame',
      gen,
      degree: cfg.degree,
      t: prop.get_time(),
      state: Array.from(prop.get_state()),
      alt: prop.get_altitude(),
      speed: prop.get_speed(),
      elements: Array.from(prop.get_orbital_elements()),
      trail: pendingTrail,
      achieved,
      period: lastPeriod,
      revolutions: rev ? rev.count : 0,
      impacted,
    };
    if (shadow) msg.shadowState = Array.from(shadow.get_state());
    postMessage(msg);
    pendingTrail = [];
  }
  setTimeout(tick, TICK_MS);
}

self.onmessage = async (e) => {
  const m = e.data;
  switch (m.type) {
    case 'init': {
      await init();
      gm = m.gm;
      const probe = new Propagator();
      probe.set_gravity_degree(FULL_DEGREE_REQUEST);
      postMessage({ type: 'ready', fullDegree: probe.get_loaded_degree() });
      probe.free();
      lastWall = performance.now();
      tick();
      break;
    }
    case 'config': {
      Object.assign(cfg, m.cfg);
      if (prop) applyForces(prop, true);
      if (shadow) applyForces(shadow, false);
      break;
    }
    case 'reset': {
      if (m.cfg) Object.assign(cfg, m.cfg);
      gen = m.gen;
      reset(m.elements, m.measure);
      break;
    }
    case 'run': {
      running = m.running;
      lastWall = performance.now();
      break;
    }
  }
};
