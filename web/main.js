/**
 * Lunar Orbit Explorer — main.js (Phase 2)
 *
 * Wires up the WASM propagator to a CesiumJS lunar visualization.
 * Phase 2 additions: force model selection, orbital elements readout,
 * eccentricity & altitude time-series plots.
 */

import * as Cesium from 'cesium';
import 'cesium/Build/Cesium/Widgets/widgets.css';
import init, { Propagator } from '../propagator/pkg/propagator.js';

window.Cesium = Cesium;

// ─── Constants ────────────────────────────────────────────────────────────

const LUNAR_GM      = 4902.800066;
const LUNAR_RADIUS  = 1737.4;
const METERS_PER_KM = 1000.0;
const RAD2DEG       = 180.0 / Math.PI;

const ORBIT = {
  sma:  1838.13,
  ecc:  0.0076,
  inc:  179.07 * Math.PI / 180.0,
  raan: 183.41 * Math.PI / 180.0,
  argp: 179.86 * Math.PI / 180.0,
  ta:   0.0,
};

const MAX_TRAIL_POINTS = 2000;
const WARP_LEVELS = [1, 10, 100, 500, 1000, 5000];
const MAX_PLOT_POINTS = 600;

const ionToken = import.meta.env.VITE_CESIUM_TOKEN ?? '';
if (ionToken) Cesium.Ion.defaultAccessToken = ionToken;

// ─── Viewer ───────────────────────────────────────────────────────────────

function createViewer() {
  const viewer = new Cesium.Viewer('cesiumContainer', {
    baseLayerPicker: false, geocoder: false, homeButton: false,
    sceneModePicker: false, navigationHelpButton: false,
    animation: false, timeline: false, fullscreenButton: false,
    infoBox: false, selectionIndicator: false,
    globe: false, baseLayer: false,
  });
  viewer.scene.backgroundColor = Cesium.Color.BLACK;
  viewer.scene.skyBox.show = false;
  viewer.scene.sun.show = false;
  viewer.scene.moon.show = false;
  viewer.scene.skyAtmosphere.show = false;
  window.viewer = viewer;
  return viewer;
}

// ─── Moon sphere ──────────────────────────────────────────────────────────

function createMoonSphere(viewer) {
  return viewer.scene.primitives.add(
    new Cesium.Primitive({
      geometryInstances: new Cesium.GeometryInstance({
        geometry: new Cesium.EllipsoidGeometry({
          radii: new Cesium.Cartesian3(
            LUNAR_RADIUS * METERS_PER_KM,
            LUNAR_RADIUS * METERS_PER_KM,
            LUNAR_RADIUS * METERS_PER_KM
          ),
          vertexFormat: Cesium.MaterialAppearance.VERTEX_FORMAT,
        }),
      }),
      appearance: new Cesium.MaterialAppearance({
        material: Cesium.Material.fromType('Color', {
          color: new Cesium.Color(0.58, 0.56, 0.52, 1.0),
        }),
        faceForward: true, flat: true, closed: true, translucent: false,
      }),
      asynchronous: false,
    })
  );
}

// ─── Coordinate conversion ────────────────────────────────────────────────

function mciToCartesian3(x, y, z) {
  return new Cesium.Cartesian3(x * METERS_PER_KM, y * METERS_PER_KM, z * METERS_PER_KM);
}

// ─── Orbit trail & spacecraft ─────────────────────────────────────────────

function createTrailEntity(viewer) {
  return viewer.entities.add({
    polyline: {
      positions: new Cesium.CallbackProperty(() => trailPositions, false),
      width: 3.0,
      material: new Cesium.PolylineGlowMaterialProperty({
        glowPower: 0.4, color: Cesium.Color.CYAN.withAlpha(1.0),
      }),
      arcType: Cesium.ArcType.NONE,
      depthFailMaterial: new Cesium.PolylineGlowMaterialProperty({
        glowPower: 0.25, color: Cesium.Color.CYAN.withAlpha(0.25),
      }),
    },
  });
}

function createSpacecraftEntity(viewer) {
  return viewer.entities.add({
    position: new Cesium.CallbackProperty(() => spacecraftPosition, false),
    point: {
      pixelSize: 7, color: Cesium.Color.WHITE,
      outlineColor: Cesium.Color.CYAN, outlineWidth: 2.0,
      heightReference: Cesium.HeightReference.NONE,
      disableDepthTestDistance: Number.POSITIVE_INFINITY,
    },
  });
}

// ─── Module-level mutable state ───────────────────────────────────────────

let trailPositions     = [];
let spacecraftPosition = Cesium.Cartesian3.ZERO;
let propagator         = null;
let isPlaying          = true;
let warpIndex          = 2;
let eccHistory         = [];
let altHistory         = [];

// ─── DOM references ───────────────────────────────────────────────────────

const hudAlt     = document.getElementById('hud-alt');
const hudVel     = document.getElementById('hud-vel');
const hudTime    = document.getElementById('hud-time');
const warpLabel  = document.getElementById('warp-label');
const statusBar  = document.getElementById('status-bar');
const warpSlider = document.getElementById('warp-slider');
const hudPeriod  = document.getElementById('hud-period');

const oeSma  = document.getElementById('oe-sma');
const oeEcc  = document.getElementById('oe-ecc');
const oeInc  = document.getElementById('oe-inc');
const oeRaan = document.getElementById('oe-raan');
const oeArgp = document.getElementById('oe-argp');
const oeTa   = document.getElementById('oe-ta');

const plotEccCanvas = document.getElementById('plot-ecc');
const plotAltCanvas = document.getElementById('plot-alt');

// ─── HUD + orbital elements update ────────────────────────────────────────

function formatTime(s) {
  const d = Math.floor(s / 86400);
  const h = Math.floor((s % 86400) / 3600);
  const m = Math.floor((s % 3600) / 60);
  const sec = Math.floor(s % 60);
  if (d > 0) return `${d}d ${h}h ${m}m`;
  if (h > 0) return `${h}h ${m}m ${sec}s`;
  return `${m}m ${sec}s`;
}

function updateHUD() {
  if (!propagator) return;
  const alt = propagator.get_altitude();
  hudAlt.textContent  = alt.toFixed(1);
  hudVel.textContent  = propagator.get_speed().toFixed(3);
  hudTime.textContent = formatTime(propagator.get_time());

  // Orbital elements
  const oe = propagator.get_orbital_elements();
  const sma = oe[0], ecc = oe[1], inc = oe[2];
  const raan = oe[3], argp = oe[4], ta = oe[5];
  oeSma.textContent  = sma.toFixed(2);
  oeEcc.textContent  = ecc.toFixed(6);
  oeInc.textContent  = (inc * RAD2DEG).toFixed(2);
  oeRaan.textContent = (raan * RAD2DEG).toFixed(2);
  oeArgp.textContent = (argp * RAD2DEG).toFixed(2);
  oeTa.textContent   = (ta * RAD2DEG).toFixed(2);

  // Dynamic period from current SMA
  const period_s = 2.0 * Math.PI * Math.sqrt(sma * sma * sma / LUNAR_GM);
  hudPeriod.textContent = '~' + (period_s / 60.0).toFixed(1);

  // Accumulate time-series data
  const t = propagator.get_time();
  eccHistory.push({ t, val: ecc });
  altHistory.push({ t, val: alt });
  if (eccHistory.length > MAX_PLOT_POINTS) eccHistory.shift();
  if (altHistory.length > MAX_PLOT_POINTS) altHistory.shift();

  drawPlot(plotEccCanvas, eccHistory, '#38bdf8', 'ecc');
  drawPlot(plotAltCanvas, altHistory, '#4ade80', 'alt');
}

// ─── Canvas time-series plotter ───────────────────────────────────────────

function drawPlot(canvas, data, color, label) {
  if (!canvas || data.length < 2) return;
  const ctx = canvas.getContext('2d');
  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  let minV = Infinity, maxV = -Infinity;
  for (const d of data) {
    if (d.val < minV) minV = d.val;
    if (d.val > maxV) maxV = d.val;
  }
  // Add 10% padding
  const range = maxV - minV || 1e-6;
  minV -= range * 0.1;
  maxV += range * 0.1;
  const tMin = data[0].t;
  const tMax = data[data.length - 1].t;
  const tRange = tMax - tMin || 1;

  // Grid lines
  ctx.strokeStyle = 'rgba(80,180,255,0.12)';
  ctx.lineWidth = 1;
  for (let i = 1; i < 4; i++) {
    const y = (H * i) / 4;
    ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke();
  }

  // Data line
  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  for (let i = 0; i < data.length; i++) {
    const x = ((data[i].t - tMin) / tRange) * W;
    const y = H - ((data[i].val - minV) / (maxV - minV)) * H;
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  }
  ctx.stroke();

  // Labels: min/max on right edge
  ctx.fillStyle = 'rgba(200,230,255,0.6)';
  ctx.font = '9px Courier New';
  ctx.textAlign = 'right';
  const fmt = label === 'ecc' ? (v => v.toFixed(5)) : (v => v.toFixed(1));
  ctx.fillText(fmt(maxV + range * 0.1), W - 2, 10);
  ctx.fillText(fmt(minV + range * 0.1), W - 2, H - 3);
}

// ─── Force model controls ─────────────────────────────────────────────────

function applyForceModel() {
  if (!propagator) return;
  const mode = document.querySelector('input[name="gravity"]:checked').value;
  const shControls = document.getElementById('sh-controls');

  if (mode === 'pointmass') {
    propagator.set_gravity_degree(0);
    shControls.style.display = 'none';
  } else {
    const deg = parseInt(document.getElementById('sh-degree').value, 10);
    propagator.set_gravity_degree(deg);
    shControls.style.display = 'flex';
  }

  const earthOn = document.getElementById('chk-earth').checked;
  const sunOn   = document.getElementById('chk-sun').checked;
  propagator.enable_third_body(earthOn, sunOn);
}

function bindForceControls() {
  document.querySelectorAll('input[name="gravity"]').forEach(el => {
    el.addEventListener('change', applyForceModel);
  });
  const shSlider = document.getElementById('sh-degree');
  const shLabel  = document.getElementById('sh-degree-label');
  shSlider.addEventListener('input', () => {
    shLabel.textContent = shSlider.value;
    applyForceModel();
  });
  document.getElementById('chk-earth').addEventListener('change', applyForceModel);
  document.getElementById('chk-sun').addEventListener('change', applyForceModel);
}

// ─── Propagation loop ─────────────────────────────────────────────────────

function resetOrbit() {
  if (!propagator) return;
  propagator.init(LUNAR_GM);
  propagator.init_from_keplerian(
    ORBIT.sma, ORBIT.ecc, ORBIT.inc, ORBIT.raan, ORBIT.argp, ORBIT.ta
  );
  trailPositions = [];
  eccHistory = [];
  altHistory = [];
  const s = propagator.get_state();
  spacecraftPosition = mciToCartesian3(s[0], s[1], s[2]);
  applyForceModel();
  updateHUD();
}

let lastTimestamp = null;

function tick(timestamp) {
  if (!propagator) { requestAnimationFrame(tick); return; }
  if (lastTimestamp === null) lastTimestamp = timestamp;
  const wallDt = Math.min((timestamp - lastTimestamp) / 1000.0, 0.1);
  lastTimestamp = timestamp;

  if (isPlaying) {
    const warp   = WARP_LEVELS[warpIndex];
    const simDt  = wallDt * warp;
    const SUB_STEP = 60.0;
    let remaining  = simDt;

    while (remaining > 0) {
      const dt = Math.min(remaining, SUB_STEP);
      propagator.step(dt);
      remaining -= dt;

      const s   = propagator.get_state();
      const pos = mciToCartesian3(s[0], s[1], s[2]);
      trailPositions.push(pos);
      if (trailPositions.length > MAX_TRAIL_POINTS) trailPositions.shift();
      spacecraftPosition = pos;
    }

    updateHUD();
  }
  requestAnimationFrame(tick);
}

// ─── Camera ───────────────────────────────────────────────────────────────

function positionCamera(viewer) {
  const offset = new Cesium.Cartesian3(
    0,
    -(LUNAR_RADIUS + 4000) * METERS_PER_KM,
     (LUNAR_RADIUS + 2000) * METERS_PER_KM
  );
  viewer.camera.lookAt(Cesium.Cartesian3.ZERO, offset);
  viewer.camera.lookAtTransform(Cesium.Matrix4.IDENTITY);

  viewer.scene.preRender.addEventListener(() => {
    viewer.camera.frustum.near = 1.0;
    viewer.camera.frustum.far  = 1.0e8;
  });
  viewer.camera.frustum.near = 1.0;
  viewer.camera.frustum.far  = 1.0e8;
}

// ─── Playback controls ───────────────────────────────────────────────────

function bindControls() {
  document.getElementById('btn-play').addEventListener('click', () => {
    isPlaying = true; setActiveBtn('btn-play');
  });
  document.getElementById('btn-pause').addEventListener('click', () => {
    isPlaying = false; setActiveBtn('btn-pause');
  });
  document.getElementById('btn-reset').addEventListener('click', () => {
    resetOrbit(); isPlaying = true; setActiveBtn('btn-play');
  });
  warpSlider.addEventListener('input', () => {
    warpIndex = parseInt(warpSlider.value, 10);
    warpLabel.textContent = `${WARP_LEVELS[warpIndex]}×`;
  });
  warpLabel.textContent = `${WARP_LEVELS[warpIndex]}×`;
}

function setActiveBtn(id) {
  ['btn-play', 'btn-pause', 'btn-reset'].forEach(bid => {
    document.getElementById(bid).classList.toggle('active', bid === id);
  });
}

// ─── Entry point ──────────────────────────────────────────────────────────

async function main() {
  console.log('[LOE] main() — start');
  statusBar.textContent = 'Loading WASM propagator…';

  await init();
  propagator = new Propagator();
  statusBar.textContent = 'WASM ready — initialising orbit…';

  resetOrbit();
  console.log('[LOE] Orbit initialised', {
    alt_km: propagator.get_altitude().toFixed(1),
    spd_kms: propagator.get_speed().toFixed(3),
  });

  statusBar.textContent = 'Building 3D scene…';
  let viewer;
  try { viewer = createViewer(); }
  catch (err) {
    statusBar.textContent = `Cesium error: ${err.message}`;
    console.error('[LOE] createViewer() threw:', err);
    return;
  }

  createMoonSphere(viewer);
  createTrailEntity(viewer);
  createSpacecraftEntity(viewer);
  positionCamera(viewer);

  bindControls();
  bindForceControls();

  statusBar.textContent = 'Running';
  console.log('[LOE] main() — running');
  requestAnimationFrame(tick);
}

main().catch(err => {
  console.error('[LOE] Fatal:', err);
  if (statusBar) statusBar.textContent = `Error: ${err.message}`;
});
