import test from 'node:test';
import assert from 'node:assert/strict';
import { separation } from '../drift-geometry.js';

const R = 1838;
const V = 1.633;
// Reference on the x axis moving +y in the xy plane (normal +z).
const a = [R, 0, 0, 0, V, 0];
const at = (angle, r = R, z = 0) => [r * Math.cos(angle), r * Math.sin(angle), z, 0, 0, 0];

test('pure along-track lead is positive along-track', () => {
  const s = separation(a, at(0.1));
  assert.ok(Math.abs(s.along - R * 0.1) < 1e-9);
  assert.ok(Math.abs(s.height) < 1e-9);
  assert.equal(s.dominant, 'along-track');
});

test('lag is negative along-track', () => {
  assert.ok(separation(a, at(-0.05)).along < 0);
});

test('height and cross-track', () => {
  const s = separation(a, at(0, R - 30));
  assert.ok(Math.abs(s.height + 30) < 1e-9);
  assert.equal(s.dominant, 'height');
  const c = separation(a, at(0, R, 12));
  assert.ok(Math.abs(c.cross - 12) < 1e-9);
  assert.equal(c.dominant, 'cross-track');
});

test('total is the straight-line distance', () => {
  const b = at(0.2, R + 5, 3);
  const s = separation(a, b);
  assert.ok(Math.abs(s.total - Math.hypot(b[0] - R, b[1], b[2])) < 1e-9);
});
