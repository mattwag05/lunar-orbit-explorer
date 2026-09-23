/**
 * The eight acts plus the cockpit (PRD section 5).
 *
 * Each act owns what it emphasises and the one control it explains. The
 * simulation keeps running underneath (in sim-worker.js), except in Act 5,
 * which replays the precomputed drift table instead.
 */

import { ActView, h } from './ui.js';
import { COPY, EAGLE, PFS2, pfs2Elements, FROZEN_INCLINATIONS_DEG } from './copy.js';
import { MOON_RADIUS_KM } from './physics-constants.js';
import { COCKPIT_WARP, ACT_WARP, warpChip } from './warp.js';
import { num, signed, hoursMinutes, simSpan } from './format.js';
import { separation } from './drift-geometry.js';

const DAY = 86400;
const RAD = 180 / Math.PI;

/** Degree steps for Act 4: the labelled stops in the copy block, then full. */
function act4Stops(fullDegree) {
  return [...Object.keys(COPY.act4.stepLabels).map(Number), fullDegree];
}

// ─── Act 0 ───────────────────────────────────────────────────────────────

const act0 = {
  duration: 4000,
  enter(ctx) {
    ctx.moonMode('smooth');
    ctx.setLight(0.7);
    ctx.liveTrail('append');
    ctx.sim.config({ warp: ACT_WARP.real.mult, degree: ctx.fullDegree, earth: false, sun: false });
    ctx.setWarp(ACT_WARP.real);
    ctx.rig.rate = 0.18;
    ctx.rig.aim({ range: 6200, el: 0.28, az: -1.7, shift: 1 });
    return new ActView(COPY.act0);
  },
};

// ─── Act 1 ───────────────────────────────────────────────────────────────

const act1 = {
  duration: 16000,
  enter(ctx) {
    ctx.moonMode('smooth');
    ctx.setLight(2.0);
    ctx.liveTrail('append');
    ctx.sim.config({ warp: ACT_WARP.real.mult, degree: ctx.fullDegree, earth: false, sun: false });
    ctx.setWarp(ACT_WARP.real);
    ctx.rig.rate = 0.6;
    ctx.rig.aim({ range: 13000, el: 0.2, az: -1.9, shift: 1 });
    const view = new ActView({
      eyebrow: COPY.act1.eyebrow,
      headline: COPY.act1.headline(),
      readouts: [
        { key: 'alt', label: 'Altitude', unit: 'km' },
        { key: 'speed', label: 'Speed', unit: 'km/s' },
      ],
    });
    view.chip('scale', COPY.act1.chip);
    return view;
  },
  readouts(ctx, f, view) {
    view.set('alt', num(f.alt, 1));
    view.set('speed', num(f.speed, 3));
  },
};

// ─── Act 2 ───────────────────────────────────────────────────────────────

const act2 = {
  duration: 20000,
  enter(ctx) {
    ctx.moonMode('smooth');
    ctx.setLight(2.0);
    ctx.clearTrail();
    ctx.liveTrail('one-revolution');
    ctx.sim.config({ warp: ACT_WARP.twentyMin.mult, degree: ctx.fullDegree, earth: false, sun: false });
    ctx.setWarp(ACT_WARP.twentyMin);
    ctx.rig.rate = 0.8;
    ctx.rig.aim({ range: 7200, el: 0.62, az: -2.3, shift: 1 });
    ctx.onFirstDrag = () => ctx.advance();
    return new ActView({
      eyebrow: COPY.act2.eyebrow,
      headline: COPY.act2.headline(ctx.epoch.periodS / 60, ctx.epoch.meanAltKm),
      body: COPY.act2.body,
      readouts: [
        { key: 'period', label: 'Period, last orbit', unit: 'min' },
        { key: 'alt', label: 'Altitude', unit: 'km' },
        { key: 'speed', label: 'Speed', unit: 'km/s' },
      ],
    });
  },
  readouts(ctx, f, view) {
    view.set('period', Number.isFinite(f.period) ? num(f.period / 60, 1) : 'measuring');
    view.set('alt', num(f.alt, 1));
    view.set('speed', num(f.speed, 3));
  },
  leave(ctx) { ctx.onFirstDrag = null; },
};

// ─── Act 3 ───────────────────────────────────────────────────────────────

const act3 = {
  duration: 22000,
  enter(ctx) {
    ctx.moonMode('smooth');
    ctx.setLight(2.0);
    ctx.clearTrail();
    ctx.liveTrail('append');
    ctx.sim.config({ warp: ACT_WARP.sixHour.mult, degree: 0, earth: false, sun: false });
    ctx.setWarp(ACT_WARP.sixHour);
    ctx.rig.rate = 0.8;
    ctx.rig.aim({ range: 7000, el: 1.05, az: -2.6, shift: 1 });
    this.start = null;
    return new ActView({
      eyebrow: COPY.act3.eyebrow,
      headline: COPY.act3.headline,
      readouts: [
        { key: 'ecc', label: 'Eccentricity', unit: '' },
        { key: 'de', label: 'Δe since this act began', unit: '' },
        { key: 'alt', label: 'Altitude', unit: 'km' },
        { key: 'elapsed', label: 'Simulated time', unit: '' },
      ],
    });
  },
  readouts(ctx, f, view) {
    // Wait for the first frame after the switch to point mass.
    if (!this.start && f.degree === 0) this.start = { t: f.t, e: f.elements[1] };
    const e = f.elements[1];
    view.set('ecc', num(e, 6));
    view.set('de', this.start ? signed(e - this.start.e, 6) : '—');
    view.set('alt', num(f.alt, 1));
    view.set('elapsed', this.start ? simSpan(f.t - this.start.t) : '—');
  },
};

// ─── Act 4 ───────────────────────────────────────────────────────────────

const act4 = {
  duration: 30000,
  enter(ctx) {
    ctx.setLight(2.2);
    ctx.liveTrail('append');
    ctx.setWarp(ACT_WARP.fourMin);
    ctx.rig.rate = 0.7;
    ctx.rig.aim({ range: 5600, el: 0.3, shift: 1 });
    ctx.rig.spin = ctx.reducedMotion ? 0 : 0.05;
    const stops = act4Stops(ctx.fullDegree);
    this.stops = stops;
    this.touched = false;
    this.stepTimer = 0;

    const label = h('p', { class: 'step-label', 'aria-live': 'polite' });
    const slider = h('input', {
      type: 'range', min: 0, max: stops.length - 1, step: 1, value: 0, id: 'degree-slider',
      'aria-label': 'Spherical harmonic degree',
    });
    const ticks = h('div', { class: 'slider-ticks', 'aria-hidden': 'true' },
      stops.map((d) => h('span', { text: String(d) })));
    slider.addEventListener('input', () => { this.touched = true; this.apply(ctx, Number(slider.value)); });
    this.slider = slider;
    this.label = label;

    const view = new ActView({
      eyebrow: COPY.act4.eyebrow,
      headline: COPY.act4.headline,
      controls: [h('label', { class: 'control-label', for: 'degree-slider', text: 'Degree' }), slider, ticks, label],
      readouts: [
        { key: 'range', label: 'Surface anomaly range', unit: 'mGal' },
        { key: 'degree', label: 'Degree driving this orbit', unit: '' },
      ],
    });
    this.view = view;
    this.apply(ctx, 0);
    return view;
  },
  apply(ctx, index) {
    const degree = this.stops[index];
    this.index = index;
    this.slider.value = String(index);
    this.slider.setAttribute('aria-valuetext', `degree ${degree}`);
    ctx.sim.config({ warp: ACT_WARP.fourMin.mult, degree, earth: false, sun: false });
    const full = index === this.stops.length - 1;
    this.label.textContent = full
      ? COPY.act4.fullLabel(ctx.loadedDegree, ctx.coefficientCount)
      : COPY.act4.stepLabels[degree];
    const field = ctx.moonMode({ degree, relief: true });
    this.view.set('degree', String(degree));
    this.view.set('range', field.max > field.min ? `${signed(field.min, 0)} to ${signed(field.max, 0)}` : '0');
    this.view.chip('exaggeration', field.factor ? COPY.act4.chip(field.factor) : null);
  },
  frame(ctx, dt) {
    if (this.touched || ctx.reducedMotion || !ctx.autoplay) return;
    this.stepTimer += dt;
    const next = Math.min(this.stops.length - 1, Math.floor(this.stepTimer / 6.5));
    if (next !== this.index) this.apply(ctx, next);
  },
  leave(ctx) { ctx.rig.spin = 0; },
};

// ─── Act 5 ───────────────────────────────────────────────────────────────

const REPLAY = ACT_WARP.halfDay;
const CHECKPOINT_DAYS = [1, 3, 7];

const act5 = {
  duration: 25000,
  enter(ctx) {
    ctx.setLight(2.0);
    ctx.moonMode('tint');
    ctx.liveTrail('off');
    ctx.sim.run(false);
    ctx.rig.rate = 0.7;
    ctx.rig.aim({ range: 6800, el: 1.2, az: -2.2, shift: 1 });
    ctx.setWarp(REPLAY);
    this.replayT = 0;
    this.nextCheckpoint = 0;
    this.finished = false;
    this.show = 1; // 0 point mass, 1 real gravity

    const [pmLabel, realLabel] = COPY.act5.toggle;
    const pm = h('button', { type: 'button', class: 'seg', 'aria-pressed': 'false', disabled: true, text: pmLabel });
    const real = h('button', { type: 'button', class: 'seg', 'aria-pressed': 'true', disabled: true, text: realLabel });
    const select = (i) => {
      this.show = i;
      pm.setAttribute('aria-pressed', String(i === 0));
      real.setAttribute('aria-pressed', String(i === 1));
      ctx.styleDriftTracks(i);
    };
    pm.addEventListener('click', () => select(0));
    real.addEventListener('click', () => select(1));
    this.buttons = [pm, real];

    this.bar = h('div', {
      class: 'progress', role: 'progressbar', 'aria-label': 'Computing a week of drift',
      'aria-valuemin': 0, 'aria-valuemax': 100, 'aria-valuenow': 0,
    }, h('span', { class: 'progress-fill' }));
    this.barText = h('p', { class: 'progress-text', 'aria-live': 'polite' });

    const view = new ActView({
      eyebrow: COPY.act5.eyebrow,
      headline: COPY.act5.headline,
      controls: [h('div', { class: 'segmented', role: 'group', 'aria-label': 'Gravity model' }, pm, real), this.bar, this.barText],
      readouts: [
        { key: 'elapsed', label: 'Simulated time', unit: '' },
        { key: 'sep', label: 'Apart now', unit: 'km' },
        { key: 'height', label: 'Height', unit: 'km' },
        { key: 'along', label: 'Along-track', unit: 'km' },
        { key: 'cross', label: 'Cross-track', unit: 'km' },
      ],
      lines: true,
    });
    this.view = view;
    ctx.styleDriftTracks(1);
    this.onProgress(ctx);
    return view;
  },
  onProgress(ctx) {
    const d = ctx.drift;
    const pct = d.total ? Math.round((100 * d.done) / d.total) : 0;
    this.bar.setAttribute('aria-valuenow', String(pct));
    this.bar.firstChild.style.transform = `scaleX(${pct / 100})`;
    if (d.table) {
      this.bar.hidden = true;
      this.barText.textContent = '';
      this.buttons.forEach((b) => { b.disabled = false; });
      this.view.chip('degree', COPY.act5.chipDegree(d.table.degree, d.table.days, d.table.chunkS));
    } else {
      this.barText.textContent = `Computing a week of drift with two propagators: ${pct}%`;
    }
  },
  frame(ctx, dt) {
    const table = ctx.drift.table;
    if (!table) { this.onProgress(ctx); return; }
    const end = table.times[table.times.length - 1];
    if (!this.finished) this.replayT = Math.min(end, this.replayT + dt * REPLAY.mult);
    const k = Math.min(table.times.length - 1, Math.floor(this.replayT / table.chunkS));
    ctx.drawDrift(table, k, this.replayT);

    while (this.nextCheckpoint < CHECKPOINT_DAYS.length && this.replayT >= CHECKPOINT_DAYS[this.nextCheckpoint] * DAY) {
      const days = CHECKPOINT_DAYS[this.nextCheckpoint];
      const idx = Math.round((days * DAY) / table.chunkS);
      const s = separation(table.a.subarray(idx * 6, idx * 6 + 6), table.b.subarray(idx * 6, idx * 6 + 6));
      this.view.addLine(COPY.act5.checkpoint(days, s.total, s.dominant), 'checkpoint');
      this.nextCheckpoint++;
    }
    if (!this.finished && this.replayT >= end) {
      this.finished = true;
      this.finishedAt = performance.now();
      this.view.addLine(COPY.act5.conclusion(), 'conclusion');
      ctx.setWarp(null);
    }
  },
  readouts(ctx, f, view) {
    const table = ctx.drift.table;
    if (!table) return;
    const k = Math.min(table.times.length - 1, Math.floor(this.replayT / table.chunkS));
    const s = separation(table.a.subarray(k * 6, k * 6 + 6), table.b.subarray(k * 6, k * 6 + 6));
    view.set('elapsed', simSpan(this.replayT));
    view.set('sep', num(s.total, 1));
    view.set('height', signed(s.height, 1));
    view.set('along', signed(s.along, 1));
    view.set('cross', signed(s.cross, 1));
  },
  /** Autoplay waits for the table, the whole replay, and time to read. */
  ready() {
    return this.finished && performance.now() - this.finishedAt > 6000;
  },
  leave(ctx) {
    ctx.endDrift();
    ctx.sim.run(true);
  },
};

// ─── Act 6 ───────────────────────────────────────────────────────────────

const act6 = {
  duration: 15000,
  enter(ctx) {
    ctx.setLight(2.0);
    ctx.moonMode('smooth');
    ctx.rig.rate = 0.7;
    ctx.rig.aim({ range: 8200, el: 0.5, az: -1.6, shift: 1 });
    ctx.setWarp(ACT_WARP.day);
    this.forces = { earth: false, sun: false };
    this.touched = false;
    this.t = 0;
    const box = (key, text) => {
      const input = h('input', { type: 'checkbox', id: `tb-${key}` });
      input.addEventListener('change', () => {
        this.touched = true;
        this.forces[key] = input.checked;
        this.restart(ctx);
      });
      this[`${key}Box`] = input;
      return h('label', { class: 'check', for: `tb-${key}` }, input, h('span', { text }));
    };
    const view = new ActView({
      eyebrow: COPY.act6.eyebrow,
      headline: COPY.act6.headline(),
      body: COPY.act6.body,
      controls: [h('div', { class: 'checks' }, box('earth', 'Earth'), box('sun', 'Sun'))],
      readouts: [
        { key: 'moved', label: 'Moved by Earth and Sun', unit: 'km' },
        { key: 'elapsed', label: 'Simulated time', unit: '' },
      ],
    });
    this.view = view;
    this.restart(ctx);
    return view;
  },
  restart(ctx) {
    ctx.clearTrail();
    ctx.liveTrail('append');
    ctx.sim.reset(EAGLE.elements, {
      warp: ACT_WARP.day.mult, degree: 0, earth: this.forces.earth, sun: this.forces.sun, shadow: true,
    });
  },
  frame(ctx, dt) {
    if (this.touched || !ctx.autoplay || ctx.reducedMotion) return;
    this.t += dt;
    if (this.t > 2 && !this.forces.earth) { this.earthBox.checked = true; this.forces.earth = true; this.restart(ctx); }
    if (this.t > 7 && !this.forces.sun) { this.sunBox.checked = true; this.forces.sun = true; this.restart(ctx); }
  },
  readouts(ctx, f, view) {
    if (f.shadowState) view.set('moved', num(separation(f.shadowState, f.state).total, 1));
    view.set('elapsed', simSpan(f.t));
  },
  leave(ctx) { ctx.sim.config({ shadow: false }); },
};

// ─── Act 7 ───────────────────────────────────────────────────────────────

const act7 = {
  duration: null,
  enter(ctx) {
    ctx.setLight(2.0);
    ctx.moonMode('tint');
    ctx.rig.rate = 0.7;
    ctx.rig.aim({ range: 7200, el: 0.45, az: -2.0, shift: 1 });
    const [moveLabel, breakLabel, rideLabel] = COPY.act7.actions;
    this.panel = h('div', { class: 'panel' });
    const view = new ActView({
      eyebrow: COPY.act7.eyebrow,
      headline: COPY.act7.headline,
      controls: [
        h('div', { class: 'actions' },
          h('button', { type: 'button', class: 'action', text: moveLabel, onclick: () => this.move(ctx) }),
          h('button', { type: 'button', class: 'action warn', text: breakLabel, onclick: () => this.breakIt(ctx) }),
          h('button', { type: 'button', class: 'action', text: rideLabel, onclick: () => ctx.enterCockpit() })),
        this.panel,
      ],
    });
    view.el.append(h('p', { class: 'footer-hint', text: COPY.act7.footer }));
    this.view = view;
    this.mode = null;
    this.move(ctx);
    return view;
  },

  move(ctx) {
    this.mode = 'move';
    const el = { ...EAGLE.elements };
    const sliders = [
      { key: 'sma', label: 'Semi-major axis', unit: 'km', min: MOON_RADIUS_KM + 30, max: MOON_RADIUS_KM + 1500, step: 1, value: el.sma, fmt: (v) => num(v, 0) },
      { key: 'ecc', label: 'Eccentricity', unit: '', min: 0, max: 0.3, step: 0.001, value: el.ecc, fmt: (v) => num(v, 3) },
      { key: 'inc', label: 'Inclination', unit: '°', min: 0, max: 180, step: 0.5, value: el.inc * RAD, fmt: (v) => num(v, 1) },
    ];
    let pending = 0;
    const restart = () => {
      clearTimeout(pending);
      pending = setTimeout(() => {
        ctx.clearTrail();
        ctx.sim.reset(el, { warp: ACT_WARP.twentyMin.mult, degree: ctx.fullDegree, earth: true, sun: true, shadow: false });
        this.impactLine.textContent = '';
      }, 150);
    };
    const rows = sliders.map((s) => {
      const out = h('output', { class: 'slider-value', for: `el-${s.key}`, text: `${s.fmt(s.value)} ${s.unit}` });
      const input = h('input', { type: 'range', id: `el-${s.key}`, min: s.min, max: s.max, step: s.step, value: s.value });
      input.addEventListener('input', () => {
        const v = Number(input.value);
        el[s.key] = s.key === 'inc' ? v / RAD : v;
        out.textContent = `${s.fmt(v)} ${s.unit}`;
        restart();
      });
      return h('div', { class: 'slider-row' },
        h('label', { for: `el-${s.key}`, class: 'control-label', text: s.label }), input, out);
    });
    this.impactLine = h('p', { class: 'warn-line', 'aria-live': 'polite' });
    this.panel.replaceChildren(...rows, this.impactLine);
    this.resetReadouts([
      { key: 'sma', label: 'Semi-major axis', unit: 'km' },
      { key: 'ecc', label: 'Eccentricity', unit: '' },
      { key: 'inc', label: 'Inclination', unit: '°' },
      { key: 'raan', label: 'Ascending node (RAAN)', unit: '°' },
      { key: 'argp', label: 'Argument of periapsis', unit: '°' },
      { key: 'ta', label: 'True anomaly', unit: '°' },
      { key: 'alt', label: 'Altitude', unit: 'km' },
    ]);
    ctx.setWarp(ACT_WARP.twentyMin);
    ctx.clearTrail();
    ctx.liveTrail('append');
    ctx.sim.reset(el, { warp: ACT_WARP.twentyMin.mult, degree: ctx.fullDegree, earth: true, sun: true, shadow: false });
  },

  breakIt(ctx) {
    this.mode = 'break';
    const eagleIncFromEquator = 180 - EAGLE.elements.inc * RAD;
    this.impactLine = h('p', { class: 'warn-line', 'aria-live': 'polite' });
    this.panel.replaceChildren(
      h('p', { class: 'body', text:
        `${PFS2.name}, ${PFS2.periseleneKm} x ${PFS2.aposeleneKm} km, ${PFS2.incPublishedDeg}° from the lunar equator. ` +
        `${EAGLE.name} flew ${num(eagleIncFromEquator, 1)}° from it. Stable low orbits sit near ` +
        `${FROZEN_INCLINATIONS_DEG.join('°, ')}°. Neither is near one.` }),
      h('p', { class: 'note', text:
        'The published record gives the altitude and tilt but not where the orbit\'s node and low point were. ' +
        `This run assumes a ${PFS2.assumedRaanDeg}° node and a ${PFS2.assumedArgpDeg}° argument of periapsis. ` +
        'Other choices fall sooner or later; the source repository\'s docs/validation.md has the sweep.' }),
      this.impactLine,
    );
    this.resetReadouts([
      { key: 'alt', label: 'Altitude', unit: 'km' },
      { key: 'peri', label: 'Lowest point this orbit', unit: 'km' },
      { key: 'elapsed', label: 'Simulated time', unit: '' },
      { key: 'revs', label: 'Revolutions', unit: '' },
    ]);
    ctx.setWarp(ACT_WARP.day);
    ctx.clearTrail();
    ctx.liveTrail('append');
    ctx.sim.reset(pfs2Elements(), { warp: ACT_WARP.day.mult, degree: ctx.fullDegree, earth: true, sun: true, shadow: false });
  },

  resetReadouts(rows) {
    this.view.readoutsEl.replaceChildren();
    this.view.values.clear();
    rows.forEach((r) => this.view.addReadout(r));
  },

  onImpact(ctx, t, revs) {
    const text = `Impact after ${num(t / DAY, 1)} days and ${num(revs, 0)} revolutions.` +
      (this.mode === 'break' ? ` The real ${PFS2.name} lasted ${PFS2.lifetimeDays} days.` : '');
    this.impactLine.textContent = text;
    ctx.setWarp(null);
  },

  readouts(ctx, f, view) {
    if (this.mode === 'move') {
      const [a, e, i, raan, argp, ta] = f.elements;
      view.set('sma', num(a, 1));
      view.set('ecc', num(e, 4));
      view.set('inc', num(i * RAD, 2));
      view.set('raan', num(raan * RAD, 1));
      view.set('argp', num(argp * RAD, 1));
      view.set('ta', num(ta * RAD, 1));
      view.set('alt', num(f.alt, 1));
    } else if (this.mode === 'break') {
      const [a, e] = f.elements;
      view.set('alt', num(f.alt, 1));
      view.set('peri', num(a * (1 - e) - MOON_RADIUS_KM, 1));
      view.set('elapsed', simSpan(f.t));
      view.set('revs', num(f.revolutions, 0));
    }
  },
};

// ─── Cockpit (Act 7b) ────────────────────────────────────────────────────

const cockpit = {
  duration: null,
  enter(ctx) {
    ctx.setLight(2.2);
    ctx.moonMode('tint');
    ctx.clearTrail();
    ctx.liveTrail('append');
    this.warpIndex = 3;
    ctx.sim.reset(EAGLE.elements, {
      warp: COCKPIT_WARP[this.warpIndex].mult, degree: ctx.fullDegree, earth: true, sun: true, shadow: false,
    });
    ctx.setWarp(COCKPIT_WARP[this.warpIndex]);

    const warpButtons = COCKPIT_WARP.map((step, i) => h('button', {
      type: 'button', class: 'seg', 'aria-pressed': String(i === this.warpIndex), text: step.label,
      onclick: () => {
        this.warpIndex = i;
        warpButtons.forEach((b, k) => b.setAttribute('aria-pressed', String(k === i)));
        ctx.sim.config({ warp: step.mult });
        ctx.setWarp(step);
      },
    }));
    const modes = ['chase', 'orbit', 'view'];
    const camButtons = modes.map((m) => h('button', {
      type: 'button', class: 'seg', 'aria-pressed': String(m === 'chase'), text: m,
      onclick: () => {
        camButtons.forEach((b) => b.setAttribute('aria-pressed', String(b.textContent === m)));
        ctx.setCameraMode(m);
      },
    }));
    ctx.setCameraMode('chase');

    const back = h('button', { type: 'button', class: 'back', text: '← Back', onclick: () => ctx.exitCockpit() });
    const view = new ActView({
      eyebrow: COPY.cockpit.eyebrow,
      headline: COPY.cockpit.headline,
      controls: [
        back,
        h('p', { class: 'control-label', text: 'Camera' }),
        h('div', { class: 'segmented', role: 'group', 'aria-label': 'Camera mode' }, camButtons),
        h('p', { class: 'control-label', text: 'Time' }),
        h('div', { class: 'segmented wrap', role: 'group', 'aria-label': 'Time warp' }, warpButtons),
      ],
      readouts: [
        { key: 'speed', label: 'Speed', unit: 'km/s' },
        { key: 'alt', label: 'Height above mean surface', unit: 'km' },
        { key: 'period', label: 'Orbit period', unit: '' },
        { key: 'elapsed', label: 'Simulated time', unit: '' },
      ],
    });
    view.el.append(h('p', { class: 'footer-hint', text: COPY.act7.footer }));
    return view;
  },
  readouts(ctx, f, view) {
    view.set('speed', num(f.speed, 3));
    view.set('alt', num(f.alt, 1));
    view.set('period', Number.isFinite(f.period) ? hoursMinutes(f.period / 60) : 'measuring');
    view.set('elapsed', simSpan(f.t));
  },
  leave(ctx) { ctx.setCameraMode('rig'); },
};

export const ACTS = [act0, act1, act2, act3, act4, act5, act6, act7];
export const COCKPIT = cockpit;
