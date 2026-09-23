// Grep-style acceptance checks from PRD section 9.
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync, readdirSync } from 'node:fs';

const dir = new URL('../', import.meta.url);
const files = readdirSync(dir).filter((f) => f.endsWith('.js'));
const read = (f) => readFileSync(new URL(f, dir), 'utf8');
const presentation = files.filter((f) => f !== 'physics-constants.js');
// Code only: comments may cite measured figures and their sources.
const code = (f) => read(f).replace(/\/\*[\s\S]*?\*\//g, '').replace(/(^|[^:])\/\/.*$/gm, '$1');

test('Class B literals live only in physics-constants.js', () => {
  for (const f of presentation) {
    const text = code(f);
    for (const lit of ['1737', '384400', '384_400', '149597870', '1.495', '4902.8', '2.6617']) {
      assert.ok(!text.includes(lit), `${f} contains ${lit}`);
    }
  }
});

test('coefficient count and full degree are never literals', () => {
  for (const f of presentation) {
    const text = code(f);
    assert.ok(!/5,?151/.test(text), `${f} contains 5151`);
    assert.ok(!/degree[^\n]{0,24}\b100\b/i.test(text), `${f} hardcodes degree 100`);
  }
});

test('Act 5 UI never mentions RAAN', () => {
  const acts = read('acts.js');
  const act5 = acts.slice(acts.indexOf('// ─── Act 5'), acts.indexOf('// ─── Act 6'));
  assert.ok(act5.length > 100);
  assert.ok(!/RAAN/i.test(act5));
});

test('no drift figures hardcoded in the presentation layer', () => {
  for (const f of presentation) {
    const text = code(f);
    for (const lit of ['645', '674', '676', '89.1', '31.8']) {
      assert.ok(!text.includes(lit), `${f} contains drift figure ${lit} outside a comment`);
    }
  }
});

test('frame loop never calls the propagator', () => {
  const main = read('main.js');
  const loop = main.slice(main.indexOf('function frame('), main.indexOf('setAutoplay(ctx.autoplay);'));
  assert.ok(loop.length > 100);
  assert.ok(!/\.step\(|probe\.|gravity_anomaly_grid/.test(loop));
});

test('act copy matches PRD section 5', () => {
  const copy = read('copy.js');
  for (const s of [
    'NO GPS AT THE MOON', 'The Moon has no satellites to guide you. So how did Eagle find its way home?',
    'WHERE YOU ARE', 'THE ORBIT', 'IF THE MOON WERE SMOOTH',
    'Treat the Moon as a single point of mass and this orbit never changes. Not in a day. Not in a year.',
    'THE MOON IS NOT SMOOTH', 'The Moon is lumpy. There is extra mass buried under the maria, and it pulls.',
    'LET A WEEK PASS', 'Now let the real gravity act, and watch what a week does to the same orbit.',
    'AND TWO MORE THINGS PULL', 'NOW TRY IT YOURSELF', 'Move the orbit. Break it. Ride along.',
    'Distances shown to scale. Orbit shown at 1x.', 'Esc to return.',
  ]) assert.ok(copy.includes(s), `missing: ${s}`);
});
