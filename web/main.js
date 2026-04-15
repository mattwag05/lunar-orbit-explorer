/**
 * Lunar Orbit Explorer — main.js
 *
 * Wires up the WASM two-body propagator to a CesiumJS lunar visualization.
 * Eagle-like demo orbit: sma=1838 km, ecc=0.01, inc=179°
 *
 * Rendering strategy (Phase 1):
 *   globe: false  — disables CesiumJS's Earth Globe/tile pipeline entirely.
 *   The Moon is an EllipsoidGraphics entity at the scene origin. This avoids
 *   a CesiumJS 1.115 issue where custom-ellipsoid Globes never issue tile
 *   requests (the LOD geometric-error thresholds are calibrated for Earth's
 *   radius and silently suppress all loads for smaller bodies).
 */

import * as Cesium from 'cesium';
import 'cesium/Build/Cesium/Widgets/widgets.css';
import init, { Propagator } from '../propagator/pkg/propagator.js';

// Expose Cesium on window so we can inspect from the browser console.
window.Cesium = Cesium;

// ─── Constants ────────────────────────────────────────────────────────────

const LUNAR_GM     = 4902.800066;   // km³/s²
const LUNAR_RADIUS = 1737.4;        // km (mean)
const METERS_PER_KM = 1000.0;

// Eagle-like initial orbit (Apollo 11 LM ascent stage approximation)
// ecc=0.01 → perigee ~82.6 km, apogee ~118.8 km, mean ~100 km altitude
const ORBIT = {
  sma:  1838.0,                         // km
  ecc:  0.01,                           // near-circular; range ~83–119 km alt
  inc:  179.0 * Math.PI / 180.0,        // near-retrograde
  raan: 0.0,
  argp: 0.0,
  ta:   0.0,
};

// Trail: keep last N positions (~one full orbit at 100× warp)
const MAX_TRAIL_POINTS = 2000;

// Time warp levels (slider index → sim multiplier)
const WARP_LEVELS = [1, 10, 100, 500, 1000, 5000];

// ─── Ion token ────────────────────────────────────────────────────────────

const ionToken = import.meta.env.VITE_CESIUM_TOKEN ?? '';
if (ionToken) {
  Cesium.Ion.defaultAccessToken = ionToken;
}

// ─── Viewer ───────────────────────────────────────────────────────────────

/**
 * Build the CesiumJS viewer.
 *
 * globe: false — disables the Globe terrain/tile system entirely.
 * The Moon is rendered as a separate EllipsoidGraphics entity (see
 * createMoonSphere). This approach is reliable across CesiumJS versions
 * because it avoids the tile LOD system altogether.
 */
function createViewer() {
  console.log('[LOE] createViewer() — start');

  const viewer = new Cesium.Viewer('cesiumContainer', {
    // Disable all Earth-centric UI
    baseLayerPicker: false,
    geocoder: false,
    homeButton: false,
    sceneModePicker: false,
    navigationHelpButton: false,
    animation: false,
    timeline: false,
    fullscreenButton: false,
    infoBox: false,
    selectionIndicator: false,

    // Disable Globe/terrain entirely — Moon rendered as entity sphere.
    // With a custom-ellipsoid Globe, CesiumJS 1.115's tile LOD scheduler
    // (designed for Earth's radius) issues zero tile requests, producing a
    // solid-black globe even when an imagery provider is registered.
    globe: false,

    // No base imagery layer (no Globe → no imageryLayers collection)
    baseLayer: false,
  });

  // Pure black space background
  viewer.scene.backgroundColor = Cesium.Color.BLACK;
  viewer.scene.skyBox.show = false;
  viewer.scene.sun.show = false;
  viewer.scene.moon.show = false;  // hide Cesium's built-in Moon billboard
  viewer.scene.skyAtmosphere.show = false;

  // Note: scene.light is intentionally NOT set here.
  // The Moon sphere uses MaterialAppearance({ flat: true }) which compiles
  // the unlit shader variant — czm_lightDirectionEC is never read, so the
  // scene light has zero effect on its appearance.

  // Expose on window for interactive browser console debugging.
  window.viewer = viewer;

  console.log('[LOE] createViewer() — done', {
    hasGlobe:   !!viewer.scene.globe,
    sceneMode:  viewer.scene.mode,
    entityCount: viewer.entities.values.length,
  });

  return viewer;
}

// ─── Moon sphere ──────────────────────────────────────────────────────────

/**
 * Add the Moon as a scene Primitive (EllipsoidGeometry + MaterialAppearance).
 *
 * WHY NOT EllipsoidGraphics entity:
 *   EllipsoidGraphics uses PerInstanceColorAppearance, which runs a Phong
 *   shader keyed on czm_lightDirectionEC.  With globe:false, CesiumJS does
 *   not guarantee that uniform is populated for the entity pass — the sphere
 *   renders jet-black regardless of scene.light settings.
 *
 * WHY flat:true:
 *   MaterialAppearance({ flat: true }) compiles the "unlit" fragment shader
 *   variant that never reads czm_lightDirectionEC.  The color renders at
 *   100 % intensity no matter what any light source does.
 *
 * Phase 2: swap Cesium.Material.fromType('Color') for a
 *   Cesium.Material.fromType('Image') with an LROC equirectangular texture.
 */
function createMoonSphere(viewer) {
  console.log('[LOE] createMoonSphere() — adding Primitive (flat/unlit)', {
    radii_m: LUNAR_RADIUS * METERS_PER_KM,
  });

  const moonPrimitive = viewer.scene.primitives.add(
    new Cesium.Primitive({
      geometryInstances: new Cesium.GeometryInstance({
        geometry: new Cesium.EllipsoidGeometry({
          radii: new Cesium.Cartesian3(
            LUNAR_RADIUS * METERS_PER_KM,
            LUNAR_RADIUS * METERS_PER_KM,
            LUNAR_RADIUS * METERS_PER_KM
          ),
          // flat:true skips normals in fragment shader, but the geometry
          // still needs them for the vertex shader; VERTEX_FORMAT is safe.
          vertexFormat: Cesium.MaterialAppearance.VERTEX_FORMAT,
        }),
      }),
      appearance: new Cesium.MaterialAppearance({
        material: Cesium.Material.fromType('Color', {
          // Warm lunar-highland grey.
          // Phase 2: replace with Material.fromType('Image', { image: lrocUrl })
          color: new Cesium.Color(0.58, 0.56, 0.52, 1.0),
        }),
        faceForward: true,
        flat:        true,   // ← UNLIT: no czm_lightDirectionEC dependency
        closed:      true,
        translucent: false,
      }),
      asynchronous: false,   // compile synchronously — no loading-state frame
    })
  );

  console.log('[LOE] createMoonSphere() — primitive added', {
    sceneCount: viewer.scene.primitives.length,
  });

  return moonPrimitive;
}


// ─── Coordinate conversion ────────────────────────────────────────────────

/**
 * Moon-centered inertial Cartesian [km] → Cesium Cartesian3 [m].
 * We treat MCI ≈ MCMF (no rotation) for the Phase 1 time window.
 * The scene origin (0,0,0) = Moon center.
 */
function mciToCartesian3(x, y, z) {
  return new Cesium.Cartesian3(
    x * METERS_PER_KM,
    y * METERS_PER_KM,
    z * METERS_PER_KM
  );
}

// ─── Orbit trail ──────────────────────────────────────────────────────────

/**
 * Polyline entity that traces the spacecraft's recent history.
 * depthFailMaterial renders the behind-Moon portion as dim translucent
 * cyan so you can see the full orbital geometry at a glance.
 */
function createTrailEntity(viewer) {
  return viewer.entities.add({
    polyline: {
      positions: new Cesium.CallbackProperty(() => trailPositions, false),
      width: 3.0,
      material: new Cesium.PolylineGlowMaterialProperty({
        glowPower: 0.4,
        color: Cesium.Color.CYAN.withAlpha(1.0),
      }),
      arcType: Cesium.ArcType.NONE,
      // Show occluded segments as dim ghost trail
      depthFailMaterial: new Cesium.PolylineGlowMaterialProperty({
        glowPower: 0.25,
        color: Cesium.Color.CYAN.withAlpha(0.25),
      }),
    },
  });
}

/**
 * Bright white point that marks the spacecraft's current position.
 */
function createSpacecraftEntity(viewer) {
  return viewer.entities.add({
    position: new Cesium.CallbackProperty(() => spacecraftPosition, false),
    point: {
      pixelSize: 7,
      color: Cesium.Color.WHITE,
      outlineColor: Cesium.Color.CYAN,
      outlineWidth: 2.0,
      heightReference: Cesium.HeightReference.NONE,
      disableDepthTestDistance: Number.POSITIVE_INFINITY,
    },
  });
}

// ─── Module-level mutable state ───────────────────────────────────────────

let trailPositions   = [];
let spacecraftPosition = Cesium.Cartesian3.ZERO;
let propagator       = null;
let isPlaying        = true;
let warpIndex        = 2;  // 100×

// ─── HUD helpers ──────────────────────────────────────────────────────────

const hudAlt     = document.getElementById('hud-alt');
const hudVel     = document.getElementById('hud-vel');
const hudTime    = document.getElementById('hud-time');
const warpLabel  = document.getElementById('warp-label');
const statusBar  = document.getElementById('status-bar');
const warpSlider = document.getElementById('warp-slider');

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
  hudAlt.textContent  = propagator.get_altitude().toFixed(1);
  hudVel.textContent  = propagator.get_speed().toFixed(3);
  hudTime.textContent = formatTime(propagator.get_time());
}

// ─── Propagation loop ─────────────────────────────────────────────────────

function resetOrbit() {
  if (!propagator) return;
  propagator.init(LUNAR_GM);
  propagator.init_from_keplerian(
    ORBIT.sma, ORBIT.ecc, ORBIT.inc, ORBIT.raan, ORBIT.argp, ORBIT.ta
  );
  trailPositions = [];
  const s = propagator.get_state();
  spacecraftPosition = mciToCartesian3(s[0], s[1], s[2]);
  updateHUD();
}

let lastTimestamp = null;

function tick(timestamp) {
  if (!propagator) {
    requestAnimationFrame(tick);
    return;
  }
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

/**
 * Position the camera ~6800 km from Moon center — far enough to see a full
 * orbital loop. Camera is off-axis so the retrograde orbit appears as an
 * inclined ellipse rather than a head-on circle.
 *
 * setView() sets position + orientation independently, so a heading/pitch
 * that looks correct in Earth-ECEF often points away from (0,0,0) in our
 * Moon-centric scene.  lookAt() is explicit: it computes the direction from
 * the camera offset to the target, guaranteeing the Moon is centred on
 * startup.  lookAtTransform(IDENTITY) immediately unlocks the camera so the
 * user can freely orbit/zoom.
 *
 * FRUSTUM:
 * With globe:false, CesiumJS still derives far from WGS84 altitude.
 * Our camera is only ~477 km "above Earth" in ECEF, giving far ≈ 3962 km.
 * The Moon sphere's nearest surface is 5111 km away — beyond that far plane.
 * We override and pin near/far every preRender frame.
 */
function positionCamera(viewer) {
  // Offset from target (Moon centre) to camera, expressed as a Cartesian3.
  // lookAt places the camera at (target + offset) and aims it at target.
  const offset = new Cesium.Cartesian3(
    0,
    -(LUNAR_RADIUS + 4000) * METERS_PER_KM,   // −5737 km  (South)
     (LUNAR_RADIUS + 2000) * METERS_PER_KM    // +3737 km  (Up)
  );
  viewer.camera.lookAt(Cesium.Cartesian3.ZERO, offset);

  // Release the look-at lock so the user can orbit/zoom freely.
  viewer.camera.lookAtTransform(Cesium.Matrix4.IDENTITY);

  // Pin frustum near/far — CesiumJS recomputes from WGS84 altitude each frame.
  viewer.scene.preRender.addEventListener(() => {
    viewer.camera.frustum.near = 1.0;       // 1 m  (Moon surface is 5111 km away)
    viewer.camera.frustum.far  = 1.0e8;     // 100 000 km
  });
  // Also set immediately so the very first frame is correct.
  viewer.camera.frustum.near = 1.0;
  viewer.camera.frustum.far  = 1.0e8;

  const f = viewer.camera.frustum;
  console.log('[LOE] Camera aimed at origin via lookAt', {
    offset_y_km: (offset.y / 1000).toFixed(0),
    offset_z_km: (offset.z / 1000).toFixed(0),
    frustum_near: f.near,
    frustum_far_km: (f.far / 1000).toFixed(0),
  });
}

// ─── Controls ─────────────────────────────────────────────────────────────

function bindControls() {
  document.getElementById('btn-play').addEventListener('click', () => {
    isPlaying = true;
    setActiveBtn('btn-play');
  });
  document.getElementById('btn-pause').addEventListener('click', () => {
    isPlaying = false;
    setActiveBtn('btn-pause');
  });
  document.getElementById('btn-reset').addEventListener('click', () => {
    resetOrbit();
    isPlaying = true;
    setActiveBtn('btn-play');
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

  // 1. Load WASM propagator
  await init();
  propagator = new Propagator();
  console.log('[LOE] WASM loaded — Propagator ready');
  statusBar.textContent = 'WASM ready — initialising orbit…';

  // 2. Set initial state
  resetOrbit();
  console.log('[LOE] Orbit initialised', {
    alt_km: propagator.get_altitude().toFixed(1),
    spd_kms: propagator.get_speed().toFixed(3),
  });

  // 3. Build Cesium viewer
  statusBar.textContent = 'Building 3D scene…';
  let viewer;
  try {
    viewer = createViewer();
  } catch (err) {
    statusBar.textContent = `Cesium error: ${err.message}`;
    console.error('[LOE] createViewer() threw:', err);
    return;
  }

  // 4. Add Moon sphere (Primitive, flat/unlit — no lighting dependency)
  createMoonSphere(viewer);

  // 5. Add orbit trail + spacecraft
  createTrailEntity(viewer);
  createSpacecraftEntity(viewer);

  console.log('[LOE] Entities:', viewer.entities.values.length,
              '| Primitives:', viewer.scene.primitives.length);

  // 6. Position camera
  positionCamera(viewer);
  const cp  = viewer.camera.position;
  const mag = Cesium.Cartesian3.magnitude(cp);
  console.log('[LOE] Camera position (metres → km for readability)', {
    x_km:              (cp.x  / 1000).toFixed(1),
    y_km:              (cp.y  / 1000).toFixed(1),
    z_km:              (cp.z  / 1000).toFixed(1),
    dist_from_origin:  (mag   / 1000).toFixed(1) + ' km',
    moon_radius:       LUNAR_RADIUS + ' km',
    camera_above_surf: ((mag / 1000) - LUNAR_RADIUS).toFixed(1) + ' km',
  });

  // 7. Bind HUD controls
  bindControls();

  statusBar.textContent = 'Running ✓';
  console.log('[LOE] main() — running');

  // 8. Start propagation loop
  requestAnimationFrame(tick);
}

main().catch(err => {
  console.error('[LOE] Fatal:', err);
  if (statusBar) statusBar.textContent = `Error: ${err.message}`;
});
