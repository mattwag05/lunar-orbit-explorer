# Lunar Orbit Explorer

A guided explainer for one idea: the Moon's lumpy gravity destroys low lunar
orbits. Eight short acts walk from "where you are" to "watch a week of drift",
with a live orbit propagator running behind the text the whole time. The final
act hands over control: move the orbit, break it, or ride along.

Every number on screen comes from the propagator or from a sourced constant.
Anything distorted (sped-up time, exaggerated relief, a reduced gravity degree
for speed) says so on screen.

The spec is [docs/prd-guided-explainer.md](docs/prd-guided-explainer.md). The
physics checks behind Act 5 are in [docs/validation.md](docs/validation.md).

## What is inside

| Part | What it does |
|---|---|
| `propagator/` | Rust crate compiled to WebAssembly. DOP853 adaptive integrator, GRGM1200A spherical-harmonic gravity to degree 100 (5,151 coefficient pairs from NASA GSFC), Earth and Sun third-body terms, and a surface free-air anomaly grid. |
| `web/sim-worker.js` | Runs the live propagator in a Web Worker at the requested time warp and reports the rate it actually achieved. |
| `web/drift-worker.js` | Act 5's precompute: two propagators (point mass and degree 20) stepped for 7 days in 300 s chunks. |
| `web/acts.js`, `web/copy.js` | The act sequence and its copy. `copy.js` is the only place mission facts live, each with a source. |
| `web/physics-constants.js` | Constants mirrored from the Rust sources; a test fails if they drift. |
| `web/scene.js` | CesiumJS scene: Moon mesh tinted and (in Act 4) raised by the gravity anomaly, orbit trails, an eased camera rig. |

The layout works on desktop and on phones, where the text becomes a bottom
sheet and the scene takes drag and pinch.

## Build and run

Prerequisites: Rust with the `wasm32-unknown-unknown` target, `wasm-pack`, and
Node 20 or newer.

```bash
npm install
npm run build:wasm
npm run dev
```

Production build (WASM plus front end, output in `dist/`):

```bash
npm run build:all
```

## Tests

```bash
npm test
```

```bash
cd propagator && cargo test --release
```

`npm test` covers formatting, the drift geometry, the warp ladder, the
constants cross-check against the Rust sources, and the grep-style acceptance
checks from the PRD. The Rust suite covers the integrator, frames, third-body
terms, coefficient loading, ALF orthonormality, gravity gradients against
independent potentials, and the anomaly grid.

## Validation harness

```bash
cd propagator
cargo run --release --example validate -- eagle 20
cargo run --release --example validate -- pfs2 100 180 180
```

See [docs/validation.md](docs/validation.md) for what these reproduce.

## Propagator API (WASM)

```typescript
class Propagator {
  constructor();
  init(gm: number): void;                                   // km³/s²
  set_state(x, y, z, vx, vy, vz: number): void;             // km, km/s; resets clock
  init_from_keplerian(sma, ecc, inc, raan, argp, ta: number): void; // radians; resets clock
  set_gravity_degree(degree: number): void;                 // 0 or 1 = point mass, up to 100
  load_coefficients(data: Uint8Array): void;
  enable_third_body(earth: boolean, sun: boolean): void;
  step(dt: number): boolean;                                // seconds
  get_state(): Float64Array;                                // [x,y,z,vx,vy,vz]
  get_time(): number;
  get_altitude(): number;                                   // km above 1737.4 km
  get_speed(): number;
  get_orbital_elements(): Float64Array;                     // [sma, ecc, inc, raan, argp, ta]
  get_loaded_degree(): number;
  get_coefficient_count(): number;
  gravity_anomaly_grid(n_lat: number, n_lon: number, degree: number): Float64Array; // mGal
}
```

## Data

- Gravity: GRGM1200A, NASA Goddard Space Flight Center, truncated to 100x100.
- Mission facts: NASA NSSDCA catalog entries for the Apollo 15 and Apollo 16
  subsatellites and the NASA Science Apollo 16 Subsatellite page. Sources are
  cited inline in `web/copy.js`.

This is a simulation, not a live tracking feed.
