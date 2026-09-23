/**
 * 3D scene: Moon mesh, trails, spacecraft markers, and an eased camera rig.
 *
 * Nothing here propagates. Positions arrive from the workers; the Moon's
 * tint and relief come from `gravity_anomaly_grid` via main.js.
 */

import * as Cesium from 'cesium';
import 'cesium/Build/Cesium/Widgets/widgets.css';
import { MOON_RADIUS_KM } from './physics-constants.js';

const M = 1000; // metres per km
export const GRID = { nLat: 91, nLon: 181 };

function cssColor(name, alpha = 1) {
  const v = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  return Cesium.Color.fromCssColorString(v).withAlpha(alpha);
}

// ─── Moon ────────────────────────────────────────────────────────────────

function moonGeometry(radiiKm) {
  const { nLat, nLon } = GRID;
  const n = nLat * nLon;
  const positions = new Float64Array(n * 3);
  const st = new Float32Array(n * 2);
  for (let i = 0; i < nLat; i++) {
    const lat = (90 - (i * 180) / (nLat - 1)) * (Math.PI / 180);
    for (let j = 0; j < nLon; j++) {
      const lon = (-180 + (j * 360) / (nLon - 1)) * (Math.PI / 180);
      const k = i * nLon + j;
      const r = (radiiKm ? radiiKm[k] : MOON_RADIUS_KM) * M;
      positions[3 * k] = r * Math.cos(lat) * Math.cos(lon);
      positions[3 * k + 1] = r * Math.cos(lat) * Math.sin(lon);
      positions[3 * k + 2] = r * Math.sin(lat);
      st[2 * k] = j / (nLon - 1);
      st[2 * k + 1] = 1 - i / (nLat - 1);
    }
  }
  const indices = new Uint32Array((nLat - 1) * (nLon - 1) * 6);
  let q = 0;
  for (let i = 0; i < nLat - 1; i++) {
    for (let j = 0; j < nLon - 1; j++) {
      const a = i * nLon + j, b = a + 1, c = a + nLon, d = c + 1;
      indices.set([a, c, b, b, c, d], q);
      q += 6;
    }
  }
  const geometry = new Cesium.Geometry({
    attributes: {
      position: new Cesium.GeometryAttribute({
        componentDatatype: Cesium.ComponentDatatype.DOUBLE, componentsPerAttribute: 3, values: positions,
      }),
      st: new Cesium.GeometryAttribute({
        componentDatatype: Cesium.ComponentDatatype.FLOAT, componentsPerAttribute: 2, values: st,
      }),
    },
    indices,
    primitiveType: Cesium.PrimitiveType.TRIANGLES,
    boundingSphere: Cesium.BoundingSphere.fromVertices(positions),
  });
  return Cesium.GeometryPipeline.computeNormal(geometry);
}

/** Equirectangular tint texture; row 0 is +90° latitude, like the grid. */
function tintTexture(grid, scale) {
  const { nLat, nLon } = GRID;
  const canvas = document.createElement('canvas');
  canvas.width = nLon;
  canvas.height = nLat;
  const ctx = canvas.getContext('2d');
  const img = ctx.createImageData(nLon, nLat);
  const base = [66, 63, 61];
  const warm = [201, 138, 75];   // --accent-warm
  const cool = [30, 38, 52];
  for (let k = 0; k < nLat * nLon; k++) {
    const t = grid && scale > 0 ? Math.max(-1, Math.min(1, grid[k] / scale)) : 0;
    const to = t >= 0 ? warm : cool;
    const w = Math.abs(t);
    img.data[4 * k] = base[0] + (to[0] - base[0]) * w;
    img.data[4 * k + 1] = base[1] + (to[1] - base[1]) * w;
    img.data[4 * k + 2] = base[2] + (to[2] - base[2]) * w;
    img.data[4 * k + 3] = 255;
  }
  ctx.putImageData(img, 0, 0);
  return canvas.toDataURL();
}

class Moon {
  constructor(scene) {
    this.scene = scene;
    this.primitive = null;
    this.rotation = 0;
    this.key = null;
    this.set(null, 0, null);
  }

  /**
   * @param grid   anomaly grid [mGal] or null for a smooth, untinted Moon
   * @param scale  |mGal| mapped to full tint
   * @param radii  per-vertex radii [km] for exaggerated relief, or null
   */
  set(grid, scale, radii) {
    const primitive = new Cesium.Primitive({
      geometryInstances: new Cesium.GeometryInstance({ geometry: moonGeometry(radii) }),
      appearance: new Cesium.MaterialAppearance({
        material: Cesium.Material.fromType('Image', { image: tintTexture(grid, scale) }),
        materialSupport: Cesium.MaterialAppearance.MaterialSupport.TEXTURED,
        faceForward: false, closed: true, translucent: false, flat: false,
      }),
      asynchronous: false,
      modelMatrix: this.matrix(),
    });
    this.scene.primitives.add(primitive);
    if (this.primitive) this.scene.primitives.remove(this.primitive);
    this.primitive = primitive;
  }

  matrix() {
    return Cesium.Matrix4.fromRotationTranslation(Cesium.Matrix3.fromRotationZ(this.rotation));
  }

  /** Body-fixed → inertial is R_z(+ω t) (propagator/src/frames.rs). */
  setRotation(theta) {
    this.rotation = theta;
    if (this.primitive) this.primitive.modelMatrix = this.matrix();
  }
}

// ─── Trails and markers ──────────────────────────────────────────────────

class Track {
  constructor(viewer, color, width, alpha) {
    this.points = [];
    this.tip = Cesium.Cartesian3.ZERO;
    this.entity = viewer.entities.add({
      polyline: {
        positions: new Cesium.CallbackProperty(() => this.points, false),
        width,
        material: new Cesium.PolylineGlowMaterialProperty({ glowPower: 0.35, color: color.withAlpha(alpha) }),
        arcType: Cesium.ArcType.NONE,
        depthFailMaterial: new Cesium.PolylineGlowMaterialProperty({ glowPower: 0.2, color: color.withAlpha(alpha * 0.22) }),
      },
    });
    this.marker = viewer.entities.add({
      position: new Cesium.CallbackProperty(() => this.tip, false),
      point: {
        pixelSize: 7, color: Cesium.Color.WHITE.withAlpha(Math.max(alpha, 0.5)),
        outlineColor: color, outlineWidth: 2,
        disableDepthTestDistance: Number.POSITIVE_INFINITY,
      },
    });
  }

  // Entity.show is a plain boolean, not a Property.
  get visible() { return this.entity.show; }
  set visible(on) { this.entity.show = on; this.marker.show = on; }

  clear() { this.points = []; }

  push(xKm, yKm, zKm, max) {
    this.points.push(new Cesium.Cartesian3(xKm * M, yKm * M, zKm * M));
    if (this.points.length > max) this.points.splice(0, this.points.length - max);
  }

  setTip(xKm, yKm, zKm) {
    this.tip = new Cesium.Cartesian3(xKm * M, yKm * M, zKm * M);
  }

  setStyle(color, alpha, width) {
    this.entity.polyline.material = new Cesium.PolylineGlowMaterialProperty({ glowPower: 0.35, color: color.withAlpha(alpha) });
    this.entity.polyline.width = width;
    this.marker.point.color = Cesium.Color.WHITE.withAlpha(Math.max(alpha, 0.5));
  }
}

// ─── Camera rig ──────────────────────────────────────────────────────────

class Rig {
  constructor(camera, canvas, reducedMotion) {
    this.camera = camera;
    this.canvas = canvas;
    this.reduced = reducedMotion;
    this.cur = { az: -1.9, el: 0.35, range: 11000, tx: 0, ty: 0, tz: 0, shift: 1 };
    this.goal = { ...this.cur };
    this.spin = 0;          // rad/s of slow automatic orbit
    this.rate = 1.6;        // easing rate [1/s]; acts slow it for gentle moves
    this.follow = false;    // cockpit 'orbit' mode: target tracks the spacecraft
    this.mode = 'rig';      // 'rig' | 'chase' | 'view'
    this.craft = null;      // {r:[km], v:[km/s]} for cockpit modes
    // Screen-space offset of the target, as fractions of viewport width and
    // height, so the Moon clears the text column (right of it on desktop,
    // above the bottom sheet on phones). rangeScale keeps the whole Moon in
    // frame on narrow portrait screens.
    this.layout = { shiftX: 0, shiftY: 0, rangeScale: 1 };
  }

  aim(goal) { Object.assign(this.goal, goal); }

  drag(dxPx, dyPx) {
    this.goal.az -= dxPx * 0.005;
    this.goal.el = Math.max(-1.45, Math.min(1.45, this.goal.el + dyPx * 0.005));
    if (this.reduced) { this.cur.az = this.goal.az; this.cur.el = this.goal.el; }
  }

  zoom(factor) {
    this.goal.range = Math.max(this.goal.minRange ?? 2100, Math.min(40000, this.goal.range * factor));
  }

  update(dt) {
    this.goal.az += this.spin * dt;
    const k = this.reduced ? 1 : 1 - Math.exp(-dt * this.rate);
    if (this.follow) {
      this.cur.tx = this.goal.tx; this.cur.ty = this.goal.ty; this.cur.tz = this.goal.tz;
    }
    for (const key of ['az', 'el', 'range', 'tx', 'ty', 'tz', 'shift']) {
      this.cur[key] += (this.goal[key] - this.cur[key]) * k;
    }
    if (this.mode === 'rig' || !this.craft) this.applyRig();
    else this.applyCraft();
  }

  applyRig() {
    const c = this.cur;
    const t = new Cesium.Cartesian3(c.tx * M, c.ty * M, c.tz * M);
    const off = new Cesium.Cartesian3(
      Math.cos(c.el) * Math.cos(c.az), Math.cos(c.el) * Math.sin(c.az), Math.sin(c.el));
    const range = c.range * this.layout.rangeScale;
    const pos = Cesium.Cartesian3.add(t, Cesium.Cartesian3.multiplyByScalar(off, range * M, new Cesium.Cartesian3()), new Cesium.Cartesian3());
    const fwd = Cesium.Cartesian3.normalize(Cesium.Cartesian3.negate(off, new Cesium.Cartesian3()), new Cesium.Cartesian3());
    const z = Cesium.Cartesian3.UNIT_Z;
    let right = Cesium.Cartesian3.cross(fwd, z, new Cesium.Cartesian3());
    if (Cesium.Cartesian3.magnitude(right) < 1e-6) right = Cesium.Cartesian3.clone(Cesium.Cartesian3.UNIT_X);
    Cesium.Cartesian3.normalize(right, right);
    const up = Cesium.Cartesian3.cross(right, fwd, new Cesium.Cartesian3());
    // Cesium's fov spans the larger viewport dimension.
    const canvas = this.canvas;
    const w = canvas?.clientWidth || 16;
    const hgt = canvas?.clientHeight || 9;
    const tanMax = Math.tan((this.camera.frustum.fov ?? Math.PI / 3) / 2);
    const tanW = w >= hgt ? tanMax : tanMax * (w / hgt);
    const tanH = w >= hgt ? tanMax * (hgt / w) : tanMax;
    const lx = 2 * this.layout.shiftX * c.shift * tanW;
    const ly = 2 * this.layout.shiftY * c.shift * tanH;
    const S = Cesium.Cartesian3;
    const dir = S.normalize(
      S.subtract(S.subtract(fwd, S.multiplyByScalar(right, lx, new S()), new S()), S.multiplyByScalar(up, ly, new S()), new S()),
      new S());
    // Cesium uses `up` as given, so it must stay orthogonal to the tilted direction.
    const upOrtho = S.normalize(S.cross(S.cross(dir, up, new S()), dir, new S()), new S());
    this.camera.setView({ destination: pos, orientation: { direction: dir, up: upOrtho } });
  }

  applyCraft() {
    const r = new Cesium.Cartesian3(this.craft.r[0] * M, this.craft.r[1] * M, this.craft.r[2] * M);
    const rh = Cesium.Cartesian3.normalize(r, new Cesium.Cartesian3());
    const vh = Cesium.Cartesian3.normalize(new Cesium.Cartesian3(...this.craft.v), new Cesium.Cartesian3());
    const S = Cesium.Cartesian3;
    let pos, dir;
    if (this.mode === 'chase') {
      pos = S.add(r, S.add(S.multiplyByScalar(vh, -90 * M, new S()), S.multiplyByScalar(rh, 30 * M, new S()), new S()), new S());
      const look = S.add(r, S.multiplyByScalar(vh, 120 * M, new S()), new S());
      dir = S.normalize(S.subtract(look, pos, new S()), new S());
    } else {
      // 'view': from the spacecraft, looking ahead and down toward the horizon.
      pos = S.clone(r);
      dir = S.normalize(S.add(vh, S.multiplyByScalar(rh, -0.35, new S()), new S()), new S());
    }
    let right = S.normalize(S.cross(dir, rh, new S()), new S());
    let up = S.cross(right, dir, new S());
    // On phones, look slightly lower so the view clears the bottom sheet.
    if (this.layout.shiftY) {
      const w = this.canvas?.clientWidth || 9;
      const hgt = this.canvas?.clientHeight || 16;
      const tanMax = Math.tan((this.camera.frustum.fov ?? Math.PI / 3) / 2);
      const tanH = w >= hgt ? tanMax * (hgt / w) : tanMax;
      dir = S.normalize(S.subtract(dir, S.multiplyByScalar(up, 2 * this.layout.shiftY * tanH, new S()), new S()), new S());
      right = S.normalize(S.cross(dir, rh, new S()), new S());
      up = S.cross(right, dir, new S());
    }
    this.camera.setView({ destination: pos, orientation: { direction: dir, up } });
  }
}

// ─── Scene facade ────────────────────────────────────────────────────────

export function createScene(container, { reducedMotion }) {
  const viewer = new Cesium.Viewer(container, {
    baseLayerPicker: false, geocoder: false, homeButton: false,
    sceneModePicker: false, navigationHelpButton: false,
    animation: false, timeline: false, fullscreenButton: false,
    infoBox: false, selectionIndicator: false,
    globe: false, baseLayer: false, skyBox: false, skyAtmosphere: false,
    // 3D only: no 2D projection, so Cesium never splits the Moon mesh at ±180°.
    scene3DOnly: true,
  });
  const scene = viewer.scene;
  scene.backgroundColor = cssColor('--bg');
  if (scene.sun) scene.sun.show = false;
  if (scene.moon) scene.moon.show = false;
  scene.screenSpaceCameraController.enableInputs = false;
  scene.light = new Cesium.DirectionalLight({
    direction: Cesium.Cartesian3.normalize(new Cesium.Cartesian3(0.55, 0.8, -0.25), new Cesium.Cartesian3()),
    intensity: 2.0,
  });
  const setFrustum = () => {
    viewer.camera.frustum.near = 1.0;
    viewer.camera.frustum.far = 1.0e9;
  };
  scene.preRender.addEventListener(setFrustum);
  setFrustum();

  const signal = cssColor('--accent-signal');
  const warn = cssColor('--accent-warn');
  const ink = cssColor('--ink');

  const moon = new Moon(scene);
  const primary = new Track(viewer, signal, 3, 1);
  const secondary = new Track(viewer, warn, 2.5, 0.9);
  secondary.visible = false;
  const rig = new Rig(viewer.camera, scene.canvas, reducedMotion);

  return {
    viewer, moon, primary, secondary, rig,
    colors: { signal, warn, ink },
    setLight(intensity) {
      scene.light.intensity = intensity;
    },
  };
}
