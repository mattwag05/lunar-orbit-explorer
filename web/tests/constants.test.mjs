// Class B constants must match what the Rust propagator actually uses (PRD 5.0).
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import * as C from '../physics-constants.js';

const src = (f) => readFileSync(new URL(`../../propagator/src/${f}`, import.meta.url), 'utf8');
const rustConst = (text, name) => {
  const m = text.match(new RegExp(`const ${name}: f64 = ([0-9_.e+-]+)`));
  assert.ok(m, `${name} not found`);
  return Number(m[1].replaceAll('_', ''));
};

test('GM matches lib.rs DEFAULT_GM', () => {
  assert.equal(C.LUNAR_GM_KM3_S2, rustConst(src('lib.rs'), 'DEFAULT_GM'));
});

test('radius matches the get_altitude constant', () => {
  const m = src('lib.rs').match(/fn get_altitude[\s\S]*?\.sqrt\(\) - ([0-9.]+)/);
  assert.ok(m);
  assert.equal(C.MOON_RADIUS_KM, Number(m[1]));
});

test('third-body distances match third_body.rs', () => {
  const t = src('third_body.rs');
  assert.equal(C.EARTH_MOON_SMA_KM, rustConst(t, 'EARTH_MOON_SMA'));
  assert.equal(C.EARTH_SUN_SMA_KM, rustConst(t, 'EARTH_SUN_SMA'));
});

test('rotation rate matches frames.rs', () => {
  assert.equal(C.OMEGA_MOON_RAD_S, rustConst(src('frames.rs'), 'OMEGA_MOON'));
});
