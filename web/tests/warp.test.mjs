import test from 'node:test';
import assert from 'node:assert/strict';
import { COCKPIT_WARP, warpChip } from '../warp.js';

test('cockpit ladder matches PRD 5.9', () => {
  assert.deepEqual(COCKPIT_WARP.map((s) => [s.label, s.mult]), [
    ['Pause', 0], ['Real', 1], ['1 min/s', 60], ['4 min/s', 240], ['20 min/s', 1200], ['1 h/s', 3600],
  ]);
});

test('chip states multiplier and label, and only when time is distorted', () => {
  assert.equal(warpChip(COCKPIT_WARP[3]), 'Time shown 240x faster than real (4 min/s).');
  assert.equal(warpChip(COCKPIT_WARP[0]), null);
  assert.equal(warpChip(COCKPIT_WARP[1]), null);
});

test('chip discloses a shortfall with the achieved rate', () => {
  assert.match(warpChip(COCKPIT_WARP[5], 1500), /about 1,500x/);
});
