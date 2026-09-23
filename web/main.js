/**
 * Lunar Orbit Explorer: guided explainer (docs/prd-guided-explainer.md).
 *
 * Boots the workers and the scene, runs the act sequence, and renders. The
 * frame loop only draws what the workers produced; it never propagates.
 */

import '@fontsource/ibm-plex-sans/300.css';
import '@fontsource/ibm-plex-sans/400.css';
import '@fontsource/ibm-plex-sans/600.css';
import '@fontsource/ibm-plex-mono/400.css';

import init, { Propagator } from '../propagator/pkg/propagator.js';
import { createScene, GRID } from './scene.js';
import { Stage, buildRail } from './ui.js';
import { ACTS, COCKPIT } from './acts.js';
import { COPY, EAGLE, ACT_TITLES, APP_NAME } from './copy.js';
import { warpChip } from './warp.js';
import {
  LUNAR_GM_KM3_S2, MOON_RADIUS_KM, OMEGA_MOON_RAD_S, SURFACE_GRAVITY_KM_S2, MGAL_KM_S2,
} from './physics-constants.js';

const MAX_TRAIL_POINTS = 2000;
const FULL_DEGREE_REQUEST = 0xffffffff; // the propagator clamps to what it holds
/**
 * Act 5 working degree, chosen up front (PRD 5.6). Measured natively on the
 * default orbit: degree 20 lands within 2.2 km of degree 100 at day 7
 * (676.4 vs 674.4 km apart) at about 1/30 of the cost. docs/validation.md.
 */
const DRIFT_DEGREE = 20;
const DRIFT_DAYS = 7;
const DRIFT_CHUNK_S = 300;
const RELIEF_TARGET = 0.035; // tallest drawn bump as a fraction of the radius

const reducedMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;
const phoneQuery = matchMedia('(max-width: 767px)');

// ─── Simulation client ───────────────────────────────────────────────────

class Sim {
  constructor() {
    this.worker = new Worker(new URL('./sim-worker.js', import.meta.url), { type: 'module' });
    this.gen = 0;
    this.latest = null;
    this.waiters = new Map();
    this.onFrame = () => {};
    this.onImpact = () => {};
    this.worker.onmessage = (e) => this.receive(e.data);
  }

  receive(m) {
    if (m.type === 'ready') { this.waiters.get('ready')?.(m); return; }
    if (m.gen !== this.gen) return;
    if (m.type === 'reset') { this.latest = null; this.waiters.get(`reset${m.gen}`)?.(m); }
    else if (m.type === 'frame') { this.latest = m; this.onFrame(m); }
    else if (m.type === 'impact') this.onImpact(m);
  }

  wait(key) { return new Promise((resolve) => this.waiters.set(key, resolve)); }

  init(gm) {
    const p = this.wait('ready');
    this.worker.postMessage({ type: 'init', gm });
    return p;
  }

  reset(elements, cfg, { measure = false } = {}) {
    this.gen++;
    const p = this.wait(`reset${this.gen}`);
    this.worker.postMessage({ type: 'reset', gen: this.gen, elements, cfg, measure });
    return p;
  }

  config(cfg) { this.worker.postMessage({ type: 'config', cfg }); }
  run(running) { this.worker.postMessage({ type: 'run', running }); }
}

// ─── Boot ────────────────────────────────────────────────────────────────

async function boot() {
  const loading = document.getElementById('loading');
  document.title = APP_NAME;

  const sim = new Sim();
  const [, ready] = await Promise.all([init(), sim.init(LUNAR_GM_KM3_S2)]);

  // Main-thread probe: field inspection only, never stepped.
  const probe = new Propagator();
  probe.set_gravity_degree(FULL_DEGREE_REQUEST);
  const loadedDegree = probe.get_loaded_degree();
  const coefficientCount = probe.get_coefficient_count();

  // Act 5's table starts computing now, so it is usually ready on arrival.
  const drift = { done: 0, total: 0, table: null };
  const driftWorker = new Worker(new URL('./drift-worker.js', import.meta.url), { type: 'module' });
  driftWorker.onmessage = (e) => {
    const m = e.data;
    if (m.type === 'progress') { drift.done = m.done; drift.total = m.total; }
    if (m.type === 'done') {
      drift.table = m;
      console.info(`[LOE] Act 5 drift precompute: ${(m.elapsedMs / 1000).toFixed(2)} s at degree ${m.degree}`);
    }
  };
  driftWorker.postMessage({
    gm: LUNAR_GM_KM3_S2, elements: EAGLE.elements, degree: DRIFT_DEGREE, days: DRIFT_DAYS, chunkS: DRIFT_CHUNK_S,
  });

  const first = await sim.reset(EAGLE.elements, {
    warp: 1, degree: ready.fullDegree, earth: false, sun: false, shadow: false,
  }, { measure: true });
  sim.run(true);

  const view3d = createScene('scene', { reducedMotion });
  const { moon, primary, secondary, rig } = view3d;

  // ─── Moon field cache ──────────────────────────────────────────────────
  const fields = new Map();
  function field(degree) {
    if (!fields.has(degree)) {
      const grid = probe.gravity_anomaly_grid(GRID.nLat, GRID.nLon, degree);
      let min = Infinity, max = -Infinity;
      for (const v of grid) { if (v < min) min = v; if (v > max) max = v; }
      // Tint saturates at the 98th percentile so a few mascons don't wash out the rest.
      const abs = Array.from(grid, Math.abs).sort((a, b) => a - b);
      const p98 = abs[Math.floor(abs.length * 0.98)] || 0;
      fields.set(degree, { grid, min, max, p98, maxAbs: Math.max(Math.abs(min), Math.abs(max)) });
    }
    return fields.get(degree);
  }
  function niceFloor(x) {
    const e = 10 ** Math.floor(Math.log10(x));
    const f = x / e;
    return (f >= 5 ? 5 : f >= 2 ? 2 : 1) * e;
  }
  let moonKey = null;
  function moonMode(mode) {
    if (mode === 'smooth') {
      if (moonKey !== 'smooth') moon.set(null, 0, null);
      moonKey = 'smooth';
      return null;
    }
    if (mode === 'tint') {
      const f = field(loadedDegree);
      if (moonKey !== 'tint') moon.set(f.grid, f.p98, null);
      moonKey = 'tint';
      return f;
    }
    const f = field(mode.degree);
    const key = `relief${mode.degree}`;
    let factor = 0;
    if (f.maxAbs > 0) {
      // Anomaly as a fraction of surface gravity, drawn as the same fraction
      // of the radius times `factor` (stated on screen).
      const fraction = (f.maxAbs * MGAL_KM_S2) / SURFACE_GRAVITY_KM_S2;
      factor = niceFloor(RELIEF_TARGET / fraction);
    }
    if (moonKey !== key) {
      let radii = null;
      if (factor) {
        radii = new Float64Array(f.grid.length);
        for (let k = 0; k < radii.length; k++) {
          radii[k] = MOON_RADIUS_KM * (1 + (factor * f.grid[k] * MGAL_KM_S2) / SURFACE_GRAVITY_KM_S2);
        }
      }
      moon.set(f.maxAbs > 0 ? f.grid : null, f.p98, radii);
      moonKey = key;
    }
    return { ...f, factor };
  }

  // ─── Trails ────────────────────────────────────────────────────────────
  let trailMode = 'append';
  let trailStartT = 0;
  let replaying = false;
  sim.onFrame = (f) => {
    if (replaying) return;
    if (trailMode === 'one-revolution' && f.t - trailStartT >= first.epoch.periodS) trailMode = 'hold';
    if (trailMode === 'append' || trailMode === 'one-revolution') {
      for (let i = 0; i < f.trail.length; i += 3) primary.push(f.trail[i], f.trail[i + 1], f.trail[i + 2], MAX_TRAIL_POINTS);
    }
  };
  sim.onImpact = (m) => {
    if (current === ACTS[7]) ACTS[7].onImpact(ctx, m.t, m.revolutions);
  };

  // ─── Context the acts drive ────────────────────────────────────────────
  let warpStep = null;
  const ctx = {
    sim, rig, reducedMotion, drift, loadedDegree, coefficientCount,
    fullDegree: ready.fullDegree,
    epoch: first.epoch,
    autoplay: !reducedMotion,
    onFirstDrag: null,
    moonMode,
    setLight: (x) => view3d.setLight(x),
    setWarp(step) { warpStep = step; },
    clearTrail() { primary.clear(); trailStartT = sim.latest?.t ?? 0; },
    liveTrail(mode) {
      trailMode = mode;
      trailStartT = sim.latest?.t ?? 0;
      replaying = mode === 'off';
      if (replaying) primary.clear();
      primary.visible = true;
    },
    styleDriftTracks(selected) {
      const { signal, warn } = view3d.colors;
      primary.setStyle(signal, selected === 0 ? 1 : 0.28, selected === 0 ? 3 : 2);
      secondary.setStyle(warn, selected === 1 ? 1 : 0.28, selected === 1 ? 3 : 2);
      secondary.visible = true;
    },
    drawDrift(table, k, t) {
      const from = Math.max(0, k - MAX_TRAIL_POINTS);
      primary.clear();
      secondary.clear();
      for (let i = from; i <= k; i++) {
        primary.push(table.a[6 * i], table.a[6 * i + 1], table.a[6 * i + 2], MAX_TRAIL_POINTS);
        secondary.push(table.b[6 * i], table.b[6 * i + 1], table.b[6 * i + 2], MAX_TRAIL_POINTS);
      }
      // Tips interpolate between the table rows so motion stays smooth.
      const j = Math.min(table.times.length - 2, k);
      const u = Math.max(0, Math.min(1, (t - table.times[j]) / table.chunkS));
      const lerp = (arr, c) => arr[6 * j + c] + (arr[6 * (j + 1) + c] - arr[6 * j + c]) * u;
      primary.setTip(lerp(table.a, 0), lerp(table.a, 1), lerp(table.a, 2));
      secondary.setTip(lerp(table.b, 0), lerp(table.b, 1), lerp(table.b, 2));
      moon.setRotation(OMEGA_MOON_RAD_S * t);
    },
    endDrift() {
      secondary.visible = false;
      primary.setStyle(view3d.colors.signal, 1, 3);
      primary.clear();
      replaying = false;
    },
    setCameraMode(mode) {
      rig.mode = mode === 'chase' || mode === 'view' ? mode : 'rig';
      rig.follow = mode === 'orbit';
      if (mode === 'orbit') rig.aim({ range: 450, minRange: 60, el: 0.4 });
      if (mode === 'rig') rig.aim({ tx: 0, ty: 0, tz: 0, minRange: 2100, range: 7200 });
    },
    advance: () => step(1),
    enterCockpit() { show(COCKPIT, 7); },
    exitCockpit() { if (current === COCKPIT) goto(7); },
  };

  // ─── Act controller ────────────────────────────────────────────────────
  const stage = new Stage(document.getElementById('stage'));
  let current = null;
  let index = 0;
  let actTime = 0;
  let view = null;

  const rail = buildRail(document.getElementById('rail'), ACT_TITLES, {
    onJump(i) { setAutoplay(false); goto(i); },
    onToggleAutoplay() { setAutoplay(!ctx.autoplay); },
  });
  function setAutoplay(on) {
    ctx.autoplay = on;
    rail.setAutoplay(on);
    document.body.classList.toggle('autoplay', on);
  }

  function show(act, railIndex) {
    current?.leave?.(ctx);
    sim.run(true);
    current = act;
    index = railIndex;
    actTime = 0;
    view = act.enter(ctx);
    stage.show(view);
    rail.setCurrent(railIndex);
    rail.setProgress(0);
    document.body.dataset.act = act === COCKPIT ? 'cockpit' : String(railIndex);
  }
  function goto(i) { show(ACTS[Math.max(0, Math.min(ACTS.length - 1, i))], Math.max(0, Math.min(ACTS.length - 1, i))); }
  function step(d) {
    if (current === COCKPIT) { if (d < 0) goto(7); return; }
    if (index + d >= 0 && index + d < ACTS.length) goto(index + d);
  }

  document.getElementById('provenance').textContent = COPY.provenance(loadedDegree);

  // ─── Layout: keep the Moon clear of the text ───────────────────────────
  // Act camera ranges were tuned on a 16:9 desktop, where the Moon fits the
  // viewport height. Elsewhere, scale range so the same framing fits the
  // space actually left for the scene.
  const TAN_HALF_FOV = Math.tan(Math.PI / 6);  // Cesium default fov, larger dimension
  const DESIGN_FIT = TAN_HALF_FOV * (9 / 16);
  function layout() {
    const w = innerWidth, hgt = innerHeight;
    const tanW = w >= hgt ? TAN_HALF_FOV : TAN_HALF_FOV * (w / hgt);
    const tanH = w >= hgt ? TAN_HALF_FOV * (hgt / w) : TAN_HALF_FOV;
    const box = document.getElementById('column').getBoundingClientRect();
    if (phoneQuery.matches) {
      const visible = Math.max(1, hgt - box.height);
      rig.layout.shiftX = 0;
      rig.layout.shiftY = (hgt / 2 - visible / 2) / hgt;
      // 20% margin: on a phone the Moon should never touch the screen edge.
      rig.layout.rangeScale = 1.2 * Math.max(1, DESIGN_FIT / Math.min(tanW, tanH * (visible / hgt)));
    } else {
      const free = Math.max(1, w - box.right);
      rig.layout.shiftX = (box.right + free / 2 - w / 2) / w;
      rig.layout.shiftY = 0;
      rig.layout.rangeScale = Math.max(1, DESIGN_FIT / Math.min(tanW * (free / w), tanH));
    }
  }
  addEventListener('resize', layout);
  new ResizeObserver(layout).observe(document.getElementById('column'));

  // ─── Input: drag and pinch on the scene, taps advance during autoplay ──
  const sceneEl = document.getElementById('scene');
  const pointers = new Map();
  let moved = 0;
  let pinchStart = 0;
  let dragFired = false;
  sceneEl.addEventListener('pointerdown', (e) => {
    sceneEl.setPointerCapture(e.pointerId);
    pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
    if (pointers.size === 1) { moved = 0; dragFired = false; }
    if (pointers.size === 2) {
      const [a, b] = [...pointers.values()];
      pinchStart = Math.hypot(a.x - b.x, a.y - b.y);
    }
  });
  sceneEl.addEventListener('pointermove', (e) => {
    const p = pointers.get(e.pointerId);
    if (!p) return;
    const dx = e.clientX - p.x, dy = e.clientY - p.y;
    p.x = e.clientX; p.y = e.clientY;
    if (pointers.size === 2) {
      const [a, b] = [...pointers.values()];
      const d = Math.hypot(a.x - b.x, a.y - b.y);
      if (pinchStart > 0 && d > 0) rig.zoom(pinchStart / d);
      pinchStart = d;
      moved += 10;
      return;
    }
    moved += Math.abs(dx) + Math.abs(dy);
    if (rig.mode === 'rig') rig.drag(dx, dy);
    if (moved > 8 && !dragFired && ctx.onFirstDrag) { dragFired = true; const f = ctx.onFirstDrag; ctx.onFirstDrag = null; f(); }
  });
  const release = (e) => {
    if (!pointers.delete(e.pointerId)) return;
    if (pointers.size > 0) return;
    if (dragFired) return;
    if (ctx.autoplay && current?.duration && current !== COCKPIT) step(1);
  };
  sceneEl.addEventListener('pointerup', release);
  sceneEl.addEventListener('pointercancel', (e) => pointers.delete(e.pointerId));
  sceneEl.addEventListener('wheel', (e) => { e.preventDefault(); rig.zoom(Math.exp(e.deltaY * 0.001)); }, { passive: false });

  addEventListener('keydown', (e) => {
    if (e.defaultPrevented || e.metaKey || e.ctrlKey || e.altKey) return;
    const t = e.target;
    const inControl = t instanceof HTMLElement && t.closest('input, button, select, textarea, a');
    if (e.key === 'Escape') { if (current === COCKPIT) { e.preventDefault(); ctx.exitCockpit(); } return; }
    if (e.key === 'ArrowRight' && !(t instanceof HTMLInputElement)) { e.preventDefault(); step(1); return; }
    if (e.key === 'ArrowLeft' && !(t instanceof HTMLInputElement)) { e.preventDefault(); step(-1); return; }
    if (inControl || ['Tab', 'Shift', 'CapsLock'].includes(e.key)) return;
    if (ctx.autoplay && current?.duration && current !== COCKPIT) step(1);
  });

  // ─── Frame loop: render only ───────────────────────────────────────────
  let last = performance.now();
  let readoutClock = 0;
  function frame(now) {
    const dt = Math.min((now - last) / 1000, 0.1);
    last = now;
    const f = sim.latest;

    if (f && !replaying) {
      primary.setTip(f.state[0], f.state[1], f.state[2]);
      moon.setRotation(OMEGA_MOON_RAD_S * f.t);
      if (rig.follow) rig.aim({ tx: f.state[0], ty: f.state[1], tz: f.state[2] });
      rig.craft = { r: f.state.slice(0, 3), v: f.state.slice(3, 6) };
    }
    current.frame?.(ctx, dt);
    rig.update(dt);

    readoutClock += dt;
    if (readoutClock > 0.1) {
      readoutClock = 0;
      if (f || replaying) current.readouts?.(ctx, f, view);
      view.chip('warp', warpStep ? warpChip(warpStep, replaying ? undefined : f?.achieved) : null);
    }

    if (ctx.autoplay && current.duration && current !== COCKPIT) {
      actTime += dt * 1000;
      rail.setProgress(actTime / current.duration);
      const waiting = current.ready && !current.ready();
      if (actTime >= current.duration && !waiting) step(1);
    } else {
      rail.setProgress(0);
    }
    requestAnimationFrame(frame);
  }

  // Inspection handle for debugging in DevTools; not used by the app.
  window.__loe = { viewer: view3d.viewer, rig, sim, ctx };

  setAutoplay(ctx.autoplay);
  goto(0);
  layout();
  loading.hidden = true;
  requestAnimationFrame(frame);
}

boot().catch((err) => {
  console.error('[LOE] Fatal:', err);
  const loading = document.getElementById('loading');
  if (loading) loading.textContent = `Could not start: ${err.message}`;
});
