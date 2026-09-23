/**
 * Text column: cross-fading act views, readouts, honesty chips, act rail.
 */

export function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === undefined || v === null || v === false) continue;
    if (k === 'class') el.className = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else if (k === 'text') el.textContent = v;
    else el.setAttribute(k, v === true ? '' : v);
  }
  for (const c of children.flat()) {
    if (c === null || c === undefined || c === false) continue;
    el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return el;
}

/**
 * A rendered act. `readouts` rows are right-aligned tabular numerals with a
 * unit column, so digits never shift as values change.
 */
export class ActView {
  constructor({ eyebrow, headline, body, controls = [], readouts = [], lines = false }) {
    this.values = new Map();
    this.el = h('section', { class: 'act-view' });
    this.eyebrowEl = h('p', { class: 'eyebrow', text: eyebrow });
    this.headlineEl = h('h1', { class: 'headline', text: headline ?? '' });
    this.bodyEl = h('p', { class: 'body', text: body ?? '' });
    if (!body) this.bodyEl.hidden = true;
    this.controlsEl = h('div', { class: 'controls' }, controls);
    this.readoutsEl = h('div', { class: 'readouts', role: 'group', 'aria-label': 'Live readouts' });
    for (const r of readouts) this.addReadout(r);
    this.linesEl = h('ol', { class: 'lines', 'aria-live': 'polite' });
    if (!lines) this.linesEl.hidden = true;
    this.chipsEl = h('div', { class: 'chips' });
    this.el.append(this.eyebrowEl, this.headlineEl, this.bodyEl, this.controlsEl,
      this.readoutsEl, this.linesEl, this.chipsEl);
  }

  addReadout({ key, label, unit = '' }) {
    const v = h('span', { class: 'v', text: '—' });
    const row = h('div', { class: 'row' },
      h('span', { class: 'k', text: label }), v, h('span', { class: 'u', text: unit }));
    this.readoutsEl.append(row);
    this.values.set(key, v);
  }

  set(key, text) {
    const el = this.values.get(key);
    if (el && el.textContent !== text) el.textContent = text;
  }

  setHeadline(text) { if (this.headlineEl.textContent !== text) this.headlineEl.textContent = text; }

  setBody(text) {
    this.bodyEl.hidden = !text;
    if (this.bodyEl.textContent !== text) this.bodyEl.textContent = text ?? '';
  }

  addLine(text, cls) {
    this.linesEl.hidden = false;
    this.linesEl.append(h('li', { class: cls, text }));
  }

  clearLines() { this.linesEl.replaceChildren(); this.linesEl.hidden = true; }

  /** Chips are keyed so each distortion is stated once and updated in place. */
  chip(key, text) {
    let el = this.chipsEl.querySelector(`[data-chip="${key}"]`);
    if (!text) { el?.remove(); return; }
    if (!el) {
      el = h('p', { class: 'chip', 'data-chip': key });
      this.chipsEl.append(el);
    }
    if (el.textContent !== text) el.textContent = text;
  }
}

/** Two stacked slots so the outgoing and incoming act truly cross-fade. */
export class Stage {
  constructor(root) {
    this.root = root;
    this.current = null;
  }

  show(view) {
    const old = this.current;
    view.el.classList.add('entering');
    this.root.append(view.el);
    // Force a style flush so the transition runs from the entering state.
    void view.el.offsetWidth;
    view.el.classList.remove('entering');
    if (old) {
      old.el.classList.add('leaving');
      old.el.setAttribute('aria-hidden', 'true');
      old.el.inert = true;
      const done = () => old.el.remove();
      old.el.addEventListener('transitionend', done, { once: true });
      setTimeout(done, 700);
    }
    this.current = view;
  }
}

export function buildRail(root, titles, { onJump, onToggleAutoplay }) {
  const buttons = titles.map((t, i) => h('button', {
    class: 'rail-step', type: 'button', 'aria-label': `Act ${i}: ${t}`, title: t,
    onclick: () => onJump(i),
  }, h('span', { class: 'rail-num', text: String(i) })));
  const auto = h('button', { class: 'rail-auto', type: 'button', onclick: onToggleAutoplay });
  const progress = h('span', { class: 'rail-progress', 'aria-hidden': 'true' });
  root.replaceChildren(h('div', { class: 'rail-steps' }, buttons), auto, progress);
  return {
    setCurrent(i) {
      buttons.forEach((b, k) => {
        if (k === i) b.setAttribute('aria-current', 'step');
        else b.removeAttribute('aria-current');
      });
    },
    setAutoplay(on) {
      auto.textContent = on ? 'Autoplay on' : 'Autoplay off';
      auto.setAttribute('aria-pressed', String(on));
    },
    setProgress(fraction) {
      progress.style.transform = `scaleX(${Math.max(0, Math.min(1, fraction))})`;
    },
  };
}
