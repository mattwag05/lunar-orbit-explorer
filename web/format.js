/** Render-time number formatting. No physics lives here. */

const nf = (digits) => new Intl.NumberFormat('en-US', {
  minimumFractionDigits: digits, maximumFractionDigits: digits,
});
const cache = new Map();
function fmt(digits) {
  if (!cache.has(digits)) cache.set(digits, nf(digits));
  return cache.get(digits);
}

/** Fixed decimals with thousands separators; never renders "-0". */
export function num(value, digits = 0) {
  if (!Number.isFinite(value)) return '—';
  const s = fmt(digits).format(value);
  return /^-0(\.0+)?$/.test(s) ? s.slice(1) : s;
}

/** Signed fixed decimals ("+4.2", "−31.8"), using a true minus sign. */
export function signed(value, digits = 1) {
  const s = num(value, digits);
  if (s === '—') return s;
  if (/^0(\.0+)?$/.test(s)) return s;
  return value < 0 ? `−${s.slice(1)}` : `+${s}`;
}

/** Large distances in words: 149597870.7 → "150 million". */
export function words(value) {
  if (value >= 1e9) return `${num(value / 1e9, 0)} billion`;
  if (value >= 1e6) return `${num(value / 1e6, 0)} million`;
  return num(value, 0);
}

/** Minutes → "1 h 58 min" (or "58 min"). */
export function hoursMinutes(minutes) {
  if (!Number.isFinite(minutes)) return '—';
  const total = Math.round(minutes);
  const h = Math.floor(total / 60);
  const m = total % 60;
  return h > 0 ? `${h} h ${m} min` : `${m} min`;
}

/** Seconds of simulated time → "3.4 days" / "5.2 h" / "12 min". */
export function simSpan(seconds) {
  if (seconds >= 86400 * 0.95) return `${num(seconds / 86400, 1)} days`;
  if (seconds >= 3600) return `${num(seconds / 3600, 1)} h`;
  return `${num(seconds / 60, 0)} min`;
}
