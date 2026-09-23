import test from 'node:test';
import assert from 'node:assert/strict';
import { num, signed, words, hoursMinutes, simSpan } from '../format.js';

test('num groups thousands and never shows negative zero', () => {
  assert.equal(num(1737.4, 0), '1,737');
  assert.equal(num(-0.0000001, 6), '0.000000');
  assert.equal(num(NaN), '—');
});

test('signed uses a true minus and no sign on zero', () => {
  assert.equal(signed(4.25, 1), '+4.3');
  assert.equal(signed(-31.8, 1), '−31.8');
  assert.equal(signed(1e-9, 6), '0.000000');
});

test('words renders large distances', () => {
  assert.equal(words(149597870.7), '150 million');
  assert.equal(words(384400), '384,400');
});

test('durations', () => {
  assert.equal(hoursMinutes(118.4), '1 h 58 min');
  assert.equal(hoursMinutes(45), '45 min');
  assert.equal(simSpan(3 * 86400), '3.0 days');
  assert.equal(simSpan(7200), '2.0 h');
});
