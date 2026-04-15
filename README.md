# Lunar Orbit Explorer — Phase 1

A proof-of-concept that retires two scary integration risks:

1. **Rust → WASM propagation** running in the browser (Dormand-Prince 45 two-body integrator)
2. **CesiumJS** rendering a 3D lunar globe with a live orbit trail

The default demo orbit approximates the Apollo 11 LM ascent stage (Eagle):
near-retrograde, ~100–115 km altitude, near-circular.

---

## Architecture

```
lunar-orbit-explorer/
├── propagator/          # Rust crate → compiled to WebAssembly
│   ├── Cargo.toml
│   └── src/lib.rs       # DP45 two-body integrator + wasm-bindgen API
├── web/
│   ├── index.html       # Entry point
│   ├── main.js          # WASM init + CesiumJS wiring + animation loop
│   └── style.css        # HUD + minimal dark theme
├── package.json
├── vite.config.js
└── README.md
```

---

## Prerequisites

### 1. Rust + wasm-pack

```bash
# Install Rust (if not already installed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source "$HOME/.cargo/env"

# Add WebAssembly target
rustup target add wasm32-unknown-unknown

# Install wasm-pack
cargo install wasm-pack
```

### 2. Node.js 18+

```bash
# macOS
brew install node
```

---

## Build & Run

### Step 1 — Compile the WASM propagator

```bash
cd propagator
wasm-pack build --target web
# Output: propagator/pkg/  (JS bindings + .wasm)
```

### Step 2 — Install JS dependencies

```bash
cd ..        # back to project root
npm install
```

### Step 3 — Start the dev server

```bash
npm run dev
# Opens http://localhost:3000
```

### One-shot (build WASM + start dev server)

```bash
npm run dev:full
```

### Production build

```bash
npm run build
# Output: dist/
```

---

## Cesium Ion Token

Basic Moon rendering works without a token. For full ion asset access:

```bash
# .env.local (gitignored)
VITE_CESIUM_TOKEN=your_token_here
```

---

## HUD Controls

| Control       | Description                                    |
|---------------|------------------------------------------------|
| ▶ Play        | Resume propagation                             |
| ⏸ Pause       | Freeze propagation (camera still controllable) |
| ↺ Reset       | Return to epoch state                          |
| Time Warp     | 1× / 10× / 100× / 500× / 1000× / 5000×        |
| Altitude      | Current altitude above mean lunar surface (km) |
| Speed         | Current orbital speed (km/s)                   |
| Elapsed       | Simulated mission elapsed time                 |

---

## Propagator API (WASM)

```typescript
class Propagator {
  constructor();
  init(gm: f64): void;                                          // set GM (km³/s²)
  set_state(x, y, z, vx, vy, vz: f64): void;                   // Cartesian (km, km/s)
  init_from_keplerian(sma, ecc, inc, raan, argp, ta: f64): void; // angles in radians
  step(dt: f64): boolean;                                       // propagate dt seconds
  get_state(): Float64Array;                                    // [x,y,z,vx,vy,vz]
  get_time(): f64;                                              // elapsed seconds
  get_altitude(): f64;                                          // km above 1737.4 km
  get_speed(): f64;                                             // km/s
}
```

---

## Orbit Parameters (Default)

| Parameter | Value | Notes |
|-----------|-------|-------|
| SMA | 1838 km | ~100 km altitude |
| Eccentricity | 0.04 | Near-circular |
| Inclination | 179° | Near-retrograde |
| RAAN | 0° | |
| Arg. Perigee | 0° | |
| True Anomaly | 0° (epoch) | |
| Period | ~118.7 min | |

---

## Running Rust Tests

```bash
cd propagator
cargo test
```

Tests validate:
- One-orbit position closure (circular orbit returns to start)
- Specific orbital energy conservation (< 1e-8 relative error)
- Angular momentum conservation
- Eagle-like orbit period (~118 min)
- Keplerian element → Cartesian conversion round-trip
- Zero-dt propagation identity

---

## Phase 2 Roadmap

- Spherical harmonics (LP165P gravity model) for high-fidelity lunar perturbations
- Dormand-Prince 78 for higher accuracy
- STL terrain mesh for real surface elevation
- Multiple spacecraft / mission phases
- TLE import
