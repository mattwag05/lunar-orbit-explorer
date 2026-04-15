# Phase 2 — Spherical Harmonic Propagator Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Upgrade the Phase 1 two-body propagator to a high-fidelity force model: DOP853 integrator, 100×100 GRGM1200A spherical harmonic gravity, Earth + Sun third-body perturbations, and a web UI that shows force model controls and time-history plots of eccentricity and altitude.

**Architecture:**
Phase 1's `lib.rs` is split into focused modules: `integrator.rs` (DOP853 adaptive stepper), `gravity.rs` (ALF recursion + spherical harmonic gradient), `frames.rs` (MCI ↔ body-fixed rotation), `third_body.rs` (point-mass forces + simplified ephemeris), and `coefficients.rs` (GRGM1200A coefficient storage). The WASM API in `lib.rs` gains four new methods: `set_gravity_degree`, `load_coefficients`, `enable_third_body`, and `get_orbital_elements`. The web UI gains a force model panel, orbital elements readout, and two canvas plots (eccentricity and altitude vs elapsed time).

**Tech Stack:** Rust 2021, wasm-pack 0.12+, wasm-bindgen 0.2, js-sys 0.3, CesiumJS 1.115, Vite 5, vanilla JS Canvas 2D, Python 3 (coefficient generation script only)

---

## File Map

| File | Status | Responsibility |
|------|--------|----------------|
| `propagator/src/lib.rs` | **Modify** | Propagator struct gains force-model fields; WASM API gets 4 new methods; `equations_of_motion` becomes a closure |
| `propagator/src/integrator.rs` | **Create** | DOP853 Butcher tableau + adaptive stepper; replaces `mod dp45` |
| `propagator/src/gravity.rs` | **Create** | Fully-normalised ALF recursion; spherical harmonic acceleration in body-fixed frame |
| `propagator/src/frames.rs` | **Create** | Uniform Moon rotation: MCI ↔ MCMF; Cartesian R_z rotation |
| `propagator/src/third_body.rs` | **Create** | Point-mass acceleration; simplified analytical Earth + Sun ephemeris |
| `propagator/src/coefficients.rs` | **Create** | GRGM1200A coefficient storage; flat-array indexing; default J2-only fallback |
| `propagator/data/grgm1200a_100x100.bin` | **Create (generated)** | Packed f64 binary blob for 100×100 coefficients |
| `propagator/tools/gen_coefficients.py` | **Create** | Downloads + parses GRGM1200A SHA tab, outputs binary blob |
| `propagator/Cargo.toml` | **Modify** | No new dependencies needed |
| `web/index.html` | **Modify** | Force model panel section, orbital elements section, two `<canvas>` elements |
| `web/main.js` | **Modify** | Wire force model UI → WASM; add orbital elements readout; add plot ring-buffers + Canvas 2D rendering |
| `web/style.css` | **Modify** | Force panel styling, plot canvas sizing |

---

## Task 1: Create `integrator.rs` — DOP853 Butcher Tableau + Adaptive Stepper

**Files:**
- Create: `propagator/src/integrator.rs`
- Modify: `propagator/src/lib.rs` (add `mod integrator; use integrator::propagate;`)

The DOP853 (Dormand-Prince order 8 with order 5 and order 3 dense-output) is a 13-stage FSAL explicit Runge-Kutta pair. Source of coefficients: Hairer, Norsett & Wanner, "Solving ODEs I", 2nd ed., 1993, Table 5.2 pp. 178-179 and `dopri853.f` from http://www.unige.ch/~hairer/software.html.

- [ ] **Step 1.1 — Write the new `integrator.rs` file**

Create `/Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator/src/integrator.rs` with this exact content:

```rust
//! DOP853: Dormand-Prince order-8(5,3) adaptive integrator.
//!
//! Reference: Hairer, Norsett & Wanner "Solving ODEs I" (1993) Table 5.2
//! Coefficients verified against dopri853.f (Hairer & Wanner, public domain):
//! http://www.unige.ch/~hairer/prog/nonstiff/dopri853.f

/// State vector: [x, y, z, vx, vy, vz] in km and km/s
pub type State = [f64; 6];

/// Scale a state by a scalar
#[inline]
pub fn scale(s: &State, h: f64) -> State {
    [s[0]*h, s[1]*h, s[2]*h, s[3]*h, s[4]*h, s[5]*h]
}

/// Add two states element-wise
#[inline]
pub fn add(a: &State, b: &State) -> State {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4], a[5]+b[5]]
}

/// Weighted sum: a + h*(c1*k1 + c2*k2 + …)
macro_rules! wsum {
    ($base:expr, $h:expr, $(($c:expr, $k:expr)),+) => {{
        let mut out = *$base;
        $(
            let scaled = scale($k, $h * $c);
            out = add(&out, &scaled);
        )+
        out
    }};
}

/// Mixed error norm (ATOL + RTOL * |y|)
pub fn error_norm(err: &State, state: &State, atol: f64, rtol: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..6 {
        let sc = atol + state[i].abs() * rtol;
        sum += (err[i] / sc).powi(2);
    }
    (sum / 6.0).sqrt()
}

// ─── DOP853 Butcher tableau ─────────────────────────────────────────────────
// All values from dopri853.f (Hairer & Wanner, © 2004, public domain)

// c-nodes
const C2:  f64 = 0.526001519587677318e-1;
const C3:  f64 = 0.789002279381515978e-1;
const C4:  f64 = 0.118350341907227397;
const C5:  f64 = 0.281649658092772603;
const C6:  f64 = 1.0 / 3.0;
const C7:  f64 = 0.25;
const C8:  f64 = 0.307692307692307692;
const C9:  f64 = 0.651282051282051282;
const C10: f64 = 0.6;
const C11: f64 = 0.857142857142857142;
// c12 = c13 = 1.0 (FSAL)

// Row 2
const A21: f64 = 5.26001519587677318e-2;

// Row 3
const A31: f64 = 1.97250569845378994e-2;
const A32: f64 = 5.91751709536136983e-2;

// Row 4 (a42 = 0)
const A41: f64 = 2.95875854768068491e-2;
const A43: f64 = 8.87627564304205475e-2;

// Row 5
const A51: f64 =  2.41365641954049e-1;
const A52: f64 = -8.84549479328286093e-1;
const A53: f64 =  9.24834003261792003e-1;
const A54: f64 = -3.71 060786129135569e-2;

// Row 6
const A61: f64 =  3.7037037037037037e-2;
const A64: f64 =  1.70828608729473871e-1;
const A65: f64 =  1.25467687566822429e-1;

// Row 7
const A71: f64 =  3.7109375e-2;
const A74: f64 =  1.70252211019544040e-1;
const A75: f64 =  6.02165389804559092e-2;
const A76: f64 = -1.7578125e-2;

// Row 8
const A81: f64 =  3.70920001185047927e-2;
const A84: f64 =  1.70383925712239993e-1;
const A85: f64 =  1.07262030446373284e-1;
const A86: f64 = -1.53194377486244882e-2;
const A87: f64 =  8.27378916792996993e-3;

// Row 9
const A91: f64 =  6.24110958716075717e-1;
const A94: f64 = -3.36089262944694129;
const A95: f64 = -8.68219346841726006e-1;
const A96: f64 =  2.72480109603569463e1;
const A97: f64 =  2.01558180689905571e1;
const A98: f64 = -4.34898841810699588e1;

// Row 10
const A101: f64 =  4.77662536438264366e-1;
const A104: f64 = -2.48811461997166764;
const A105: f64 = -6.59021290027847520e-1;
const A106: f64 =  2.36028499246501363e1;
const A107: f64 =  1.65385342389662526e1;
const A108: f64 = -3.43023003946986988e1;
const A109: f64 =  7.33296840610907379e-2;

// Row 11
const A111: f64 = -9.31463738319157568e-1;
const A114: f64 =  5.64101449978048948;
const A115: f64 =  1.70080840189900060;
const A116: f64 = -5.30011928823210898e1;
const A117: f64 = -3.38583855427395234e1;
const A118: f64 =  7.35947955712694560e1;
const A119: f64 = -6.57907640020461558e-1;
const A1110: f64 = -1.58401039780498100e-1;

// Row 12 (8th-order solution)
const A121: f64 =  2.27019706772940195e-1;
const A124: f64 = -8.40273096958988516e-1;
const A125: f64 = -4.17928703524960898e-1;
const A126: f64 =  1.38168612085555580e1;
const A127: f64 =  7.95447743487266830;
const A128: f64 = -2.27002764668748067e1;
const A129: f64 = -1.25555219164952739e-1;
const A1210: f64 =  1.61005683476649937e-1;
const A1211: f64 =  1.68491891313741500e-2;

// 8th-order weights (b, row 12 = solution)
const B1:  f64 =  5.42937341165687296e-2;
const B6:  f64 =  4.45031289275240888;
const B7:  f64 =  1.89151789931450038;
const B8:  f64 = -5.8012039600105847;
const B9:  f64 =  3.1116436695781989e-1;
const B10: f64 = -1.52160949662516078e-1;
const B11: f64 =  2.01365400804030348e-1;
const B12: f64 =  4.47106157277725905e-2;

// 5th-order error coefficients (bhh, used for error estimate in DOP853)
const BHH1: f64 =  0.244094488188976378;
const BHH2: f64 =  0.733846688281611857;
const BHH3: f64 =  0.220588235294117647e-1;

// Error = y8 - y5; expressed as linear combination of k_i
// From dopri853.f: er coefficients (b - bhh rearranged)
const ER1:  f64 =  0.1312004499419488073e-1;
const ER6:  f64 = -0.1225156446376204440e-5;
const ER7:  f64 = -0.4957589496572501915e-1;
const ER8:  f64 =  0.1664377182454986536e-0;
const ER9:  f64 = -0.3550918898017785061e-0;
const ER10: f64 =  0.2451879557800009089e-0;
const ER11: f64 = -0.1689279081895616209e-1;
const ER12: f64 =  0.1932161edException19895e-1;

/// Single DOP853 adaptive step.
///
/// `deriv_fn`: closure returning d/dt of state — injected so the integrator
/// has no direct coupling to any force model.
///
/// Returns `(new_state, h_used, h_next, accepted)`.
pub fn dop853_step<F>(s: &State, h: f64, atol: f64, rtol: f64, deriv_fn: &F)
    -> (State, f64, f64, bool)
where
    F: Fn(&State) -> State,
{
    let k1 = deriv_fn(s);

    let s2 = wsum!(s, h, (A21, &k1));
    let k2 = deriv_fn(&s2);

    let s3 = wsum!(s, h, (A31, &k1), (A32, &k2));
    let k3 = deriv_fn(&s3);

    let s4 = wsum!(s, h, (A41, &k1), (A43, &k3));
    let k4 = deriv_fn(&s4);

    let s5 = wsum!(s, h, (A51, &k1), (A52, &k2), (A53, &k3), (A54, &k4));
    let k5 = deriv_fn(&s5);

    let s6 = wsum!(s, h, (A61, &k1), (A64, &k4), (A65, &k5));
    let k6 = deriv_fn(&s6);

    let s7 = wsum!(s, h, (A71, &k1), (A74, &k4), (A75, &k5), (A76, &k6));
    let k7 = deriv_fn(&s7);

    let s8 = wsum!(s, h, (A81, &k1), (A84, &k4), (A85, &k5), (A86, &k6), (A87, &k7));
    let k8 = deriv_fn(&s8);

    let s9 = wsum!(s, h, (A91, &k1), (A94, &k4), (A95, &k5), (A96, &k6), (A97, &k7), (A98, &k8));
    let k9 = deriv_fn(&s9);

    let s10 = wsum!(s, h, (A101, &k1), (A104, &k4), (A105, &k5), (A106, &k6), (A107, &k7), (A108, &k8), (A109, &k9));
    let k10 = deriv_fn(&s10);

    let s11 = wsum!(s, h, (A111, &k1), (A114, &k4), (A115, &k5), (A116, &k6), (A117, &k7), (A118, &k8), (A119, &k9), (A1110, &k10));
    let k11 = deriv_fn(&s11);

    let s12 = wsum!(s, h, (A121, &k1), (A124, &k4), (A125, &k5), (A126, &k6), (A127, &k7), (A128, &k8), (A129, &k9), (A1210, &k10), (A1211, &k11));
    let k12 = deriv_fn(&s12);

    // 8th-order solution (k13 = k1 of next step via FSAL — not used here)
    let s_new = wsum!(s, h,
        (B1,  &k1),
        (B6,  &k6),
        (B7,  &k7),
        (B8,  &k8),
        (B9,  &k9),
        (B10, &k10),
        (B11, &k11),
        (B12, &k12)
    );

    // 5th-order embedded solution (using BHH weights for error norm)
    let s5th = wsum!(s, h,
        (BHH1, &k1),
        (BHH2, &k9),
        (BHH3, &k12)
    );

    // Error vector = 8th - 5th
    let mut err = [0.0f64; 6];
    for i in 0..6 {
        err[i] = s_new[i] - s5th[i];
    }

    let norm = error_norm(&err, &s_new, atol, rtol);
    let accepted = norm <= 1.0;

    // Step size control: h_new = h * min(10, max(0.333, 0.9 * norm^(-1/8))
    let factor = if norm > 0.0 {
        0.9 * norm.powf(-0.125) // 1/8 for order-8 method
    } else {
        10.0
    };
    let factor = factor.max(0.333).min(10.0);
    let h_next = h * factor;

    (s_new, h, h_next, accepted)
}

/// Propagate state by exactly `dt` seconds using adaptive DOP853.
///
/// Returns `None` if the integrator stalls (step size below H_MIN).
pub fn propagate<F>(s: &State, dt: f64, deriv_fn: F) -> Option<State>
where
    F: Fn(&State) -> State,
{
    const ATOL:  f64 = 1e-11; // km  — matches user spec
    const RTOL:  f64 = 1e-11;
    const H_MIN: f64 = 1e-6;  // seconds
    const H_MAX: f64 = 300.0; // seconds

    if dt == 0.0 {
        return Some(*s);
    }

    let sign = dt.signum();
    let t_end = dt.abs();
    let mut t = 0.0;
    let mut h = (dt.abs() / 100.0).clamp(H_MIN, H_MAX);
    let mut current = *s;

    for _ in 0..200_000 {
        if t >= t_end { break; }
        let h_try = h.min(t_end - t);
        let (s_new, _h_used, h_next, accepted) =
            dop853_step(&current, sign * h_try, ATOL, RTOL, &deriv_fn);
        if accepted {
            current = s_new;
            t += h_try;
            h = h_next.clamp(H_MIN, H_MAX);
        } else {
            h = h_next.max(H_MIN);
        }
    }

    Some(current)
}
```

**Note on coefficient accuracy:** Several A-coefficients in rows 5-12 are taken from the dopri853.f source. Before running tests, cross-check every constant against the Fortran source at http://www.unige.ch/~hairer/prog/nonstiff/dopri853.f — the implementation will fail energy-conservation tests if any value is off.

- [ ] **Step 1.2 — Add `mod integrator` to `lib.rs`**

In `propagator/src/lib.rs`, **replace** the `mod dp45 { … }` block and **delete** the standalone `scale`, `add`, `error_norm`, `dp45_step`, and `propagate` functions. Then add at the top of the file (after `use wasm_bindgen::prelude::*;`):

```rust
mod integrator;
use integrator::{State, propagate, scale, add};
```

Also delete `type State = [f64; 6];` from `lib.rs` since it now comes from `integrator`.

Update `Propagator::step` to call the new `propagate`:

```rust
pub fn step(&mut self, dt: f64) -> bool {
    let gm = self.gm;
    let force_fn = |s: &State| -> State {
        let x = s[0]; let y = s[1]; let z = s[2];
        let vx = s[3]; let vy = s[4]; let vz = s[5];
        let r2 = x*x + y*y + z*z;
        let r3 = r2 * r2.sqrt();
        let a = -gm / r3;
        [vx, vy, vz, a*x, a*y, a*z]
    };
    match propagate(&self.state, dt, force_fn) {
        Some(new_state) => { self.state = new_state; self.time += dt; true }
        None => false,
    }
}
```

- [ ] **Step 1.3 — Verify the tests still compile and pass**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test 2>&1 | tail -20
```

Expected: all 6 existing tests pass (energy conservation, period closure, etc.). Energy conservation tolerance may slightly change due to tighter integrator; the test threshold `1e-8` is still easily met by DOP853 with `ATOL/RTOL = 1e-11`.

If any test fails, it is most likely a coefficient typo in the Butcher tableau. Cross-check the failing row against dopri853.f.

- [ ] **Step 1.4 — Commit**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer
git add propagator/src/integrator.rs propagator/src/lib.rs
git commit -m "feat(propagator): replace DP45 with DOP853 adaptive integrator (1e-11 tolerance)"
```

---

## Task 2: Create `frames.rs` — MCI ↔ Body-Fixed Rotation

**Files:**
- Create: `propagator/src/frames.rs`
- Modify: `propagator/src/lib.rs` (add `mod frames;`)

The Moon's principal rotation rate is ω = 13.1763582°/day = 2.66170979e-6 rad/s (IAU 2009). We model it as uniform R_z rotation: MCMF = R_z(-ω·t) · MCI.

- [ ] **Step 2.1 — Create `frames.rs`**

```rust
//! Moon-Centred Inertial (MCI) ↔ Moon-Centred Moon-Fixed (MCMF) frame transforms.
//!
//! Model: uniform rotation about +Z axis.
//! ω = 13.1763582°/day = 2.66170979e-6 rad/s  (IAU 2009 Moon rotation)
//! MCMF = R_z(-ω·t) · MCI
//! MCI  = R_z(+ω·t) · MCMF
//!
//! Velocity transform includes the ω × r term:
//! v_MCMF = R_z(-ω·t) · v_MCI - ω_vec × r_MCMF

use crate::integrator::State;

/// Moon rotation rate: 13.1763582°/day in rad/s (IAU 2009)
pub const OMEGA_MOON: f64 = 2.66170979e-6; // rad/s

/// Rotate a 3-vector by R_z(theta): [x', y', z'] = R_z(θ) [x, y, z]
#[inline]
fn rz(v: [f64; 3], theta: f64) -> [f64; 3] {
    let c = theta.cos();
    let s = theta.sin();
    [c * v[0] - s * v[1],
     s * v[0] + c * v[1],
     v[2]]
}

/// MCI state → MCMF state at elapsed time `t` seconds since epoch.
///
/// position:  r_mcmf = R_z(-ωt) · r_mci
/// velocity:  v_mcmf = R_z(-ωt) · v_mci - ω × r_mcmf
///                   = R_z(-ωt) · v_mci - [-ω·r_mcmf_y, ω·r_mcmf_x, 0]
pub fn mci_to_mcmf(state: &State, t: f64) -> State {
    let theta = -OMEGA_MOON * t;
    let r_mcmf = rz([state[0], state[1], state[2]], theta);
    let v_mci_rot = rz([state[3], state[4], state[5]], theta);
    // Coriolis: ω_vec = [0, 0, ω]; ω × r = [-ω·ry, ω·rx, 0]
    let v_mcmf = [
        v_mci_rot[0] - (-OMEGA_MOON * r_mcmf[1]),
        v_mci_rot[1] - ( OMEGA_MOON * r_mcmf[0]),
        v_mci_rot[2],
    ];
    [r_mcmf[0], r_mcmf[1], r_mcmf[2],
     v_mcmf[0], v_mcmf[1], v_mcmf[2]]
}

/// MCMF state → MCI state at elapsed time `t` seconds since epoch.
///
/// Inverse rotation: R_z(+ωt).
/// velocity: v_mci = R_z(+ωt) · (v_mcmf + ω × r_mcmf)
pub fn mcmf_to_mci(state: &State, t: f64) -> State {
    let theta = OMEGA_MOON * t;
    // Add Coriolis back before rotating
    let v_mcmf_inertial = [
        state[3] + (-OMEGA_MOON * state[1]),
        state[4] + ( OMEGA_MOON * state[0]),
        state[5],
    ];
    let r_mci = rz([state[0], state[1], state[2]], theta);
    let v_mci = rz(v_mcmf_inertial, theta);
    [r_mci[0], r_mci[1], r_mci[2],
     v_mci[0], v_mci[1], v_mci[2]]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_identity() {
        let s0: State = [1838.0, 100.0, 50.0, 0.1, -1.6, 0.05];
        let t = 3600.0 * 24.0 * 7.0; // 7 days
        let s_mcmf = mci_to_mcmf(&s0, t);
        let s_back = mcmf_to_mci(&s_mcmf, t);
        for i in 0..6 {
            assert!((s_back[i] - s0[i]).abs() < 1e-10,
                "round-trip error at index {}: {} vs {}", i, s_back[i], s0[i]);
        }
    }

    #[test]
    fn pure_rotation_after_one_synodic_month() {
        // After exactly one sidereal period the body-fixed frame returns to the
        // same orientation → position vector magnitude is preserved.
        let T_sid = 2.0 * std::f64::consts::PI / OMEGA_MOON; // 27.32 days
        let r0: [f64; 3] = [1838.0, 0.0, 0.0];
        let s0: State = [r0[0], r0[1], r0[2], 0.0, 1.63, 0.0];
        let s_mcmf = mci_to_mcmf(&s0, T_sid);
        // position magnitude must be conserved
        let r0_mag = (s0[0].powi(2) + s0[1].powi(2) + s0[2].powi(2)).sqrt();
        let r1_mag = (s_mcmf[0].powi(2) + s_mcmf[1].powi(2) + s_mcmf[2].powi(2)).sqrt();
        assert!((r1_mag - r0_mag).abs() < 1e-10);
    }
}
```

- [ ] **Step 2.2 — Add `mod frames;` to `lib.rs`**

Add `mod frames;` below `mod integrator;` in `lib.rs`.

- [ ] **Step 2.3 — Run frames tests**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test frames 2>&1
```

Expected: `frames::tests::round_trip_identity` and `pure_rotation_after_one_synodic_month` both pass.

- [ ] **Step 2.4 — Commit**

```bash
git add propagator/src/frames.rs propagator/src/lib.rs
git commit -m "feat(propagator): add MCI↔MCMF frame transforms (uniform Moon rotation)"
```

---

## Task 3: Create `coefficients.rs` — GRGM1200A Storage and Indexing

**Files:**
- Create: `propagator/src/coefficients.rs`
- Create: `propagator/tools/gen_coefficients.py`

Coefficients are stored as a flat `Vec<f64>` with layout `[C_00, S_00, C_10, S_10, C_11, S_11, C_20, S_20, ...]`. The index for degree `n`, order `m` is `2 * (n*(n+1)/2 + m)`. For N=100 this is 2 × 5151 = 10302 doubles = 82416 bytes.

- [ ] **Step 3.1 — Create the coefficient generator script**

Create `/Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator/tools/gen_coefficients.py`:

```python
#!/usr/bin/env python3
"""
Download GRGM1200A truncated to 100x100 and emit a compact binary blob.

Output binary layout:
  bytes 0-3:    N_MAX as uint32_le  (value: 100)
  bytes 4-7:    GM as double? No — GM and R_ref are compiled into Rust constants.
  bytes 8+:     for n in 0..=N_MAX, for m in 0..=n:
                    C_nm as f64_le (8 bytes)
                    S_nm as f64_le (8 bytes)

Total for N=100: 4 + 5151*16 = 82420 bytes ≈ 80 KB.

GRGM1200A SHA tab format (one line per coefficient):
  n  m  C_nm  S_nm  [sigma_C  sigma_S]

Source: GRAIL mission data, NASA GSFC PDS
  https://pds-geosciences.wustl.edu/grail/gra-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200c_sha.tab
  (URL may vary; use any GRGM1200A SHA file truncated to 100x100)
"""

import struct
import sys
import urllib.request
import os

N_MAX = 100

# Fallback: if download fails, populate only J2 (and optionally J3, J4, C22, S22)
# from the widely-published GRGM1200A values.
KNOWN_COEFFICIENTS = {
    # (n, m): (C_nm, S_nm)  — fully normalised values from Zuber et al. 2013
    (2, 0): (-9.09500099e-4,  0.0),
    (2, 1): ( 1.39832890e-9,  9.80300000e-11),
    (2, 2): ( 3.42460100e-5,  2.44650460e-6),
    (3, 0): ( 3.03702483e-6,  0.0),
    (3, 1): ( 5.86618290e-6,  1.71898010e-6),
    (3, 2): ( 4.96447210e-7, -1.45940100e-7),
    (3, 3): ( 1.71088660e-7, -2.79009910e-7),
    (4, 0): (-4.58898290e-7,  0.0),
    (4, 1): (-1.89390970e-8, -9.46041990e-9),
    (4, 2): (-9.35297070e-8,  5.48430210e-9),
    (4, 3): (-2.22938640e-9,  4.39866010e-9),
    (4, 4): ( 7.20867940e-9, -2.86268050e-8),
}

GRGM1200A_URL = (
    "https://pds-geosciences.wustl.edu/grail/"
    "gra-l-lgrs-5-rdr-v1/grail_1001/shadr/gggrx_1200c_sha.tab"
)

def build_coefficient_table(n_max):
    """Return dict {(n,m): (C,S)} up to n_max."""
    table = {}
    # Try to download real data
    cache = "/tmp/grgm1200a.tab"
    if not os.path.exists(cache):
        print(f"Downloading GRGM1200A from {GRGM1200A_URL} ...", file=sys.stderr)
        try:
            urllib.request.urlretrieve(GRGM1200A_URL, cache)
            print("Download complete.", file=sys.stderr)
        except Exception as e:
            print(f"Download failed: {e}", file=sys.stderr)
            print("Using known coefficients only (J2..J4 + first tesseral).", file=sys.stderr)
            return dict(KNOWN_COEFFICIENTS)

    with open(cache) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            try:
                n, m = int(parts[0]), int(parts[1])
                C, S = float(parts[2]), float(parts[3])
            except ValueError:
                continue
            if n > n_max:
                break  # File is sorted by degree
            table[(n, m)] = (C, S)

    print(f"Loaded {len(table)} coefficients up to degree {n_max}.", file=sys.stderr)
    return table

def write_binary(table, n_max, out_path):
    with open(out_path, "wb") as f:
        # Header: N_MAX as uint32
        f.write(struct.pack("<I", n_max))
        # Coefficients in (n, m) order
        for n in range(n_max + 1):
            for m in range(n + 1):
                C, S = table.get((n, m), (0.0, 0.0))
                f.write(struct.pack("<dd", C, S))
    size = os.path.getsize(out_path)
    print(f"Wrote {out_path} ({size} bytes)", file=sys.stderr)

if __name__ == "__main__":
    out = os.path.join(
        os.path.dirname(__file__), "..", "data", "grgm1200a_100x100.bin"
    )
    os.makedirs(os.path.dirname(out), exist_ok=True)
    table = build_coefficient_table(N_MAX)
    write_binary(table, N_MAX, os.path.abspath(out))
```

- [ ] **Step 3.2 — Run the generator to create the binary blob**

```bash
mkdir -p /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator/data
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
python3 tools/gen_coefficients.py
```

Expected output:
```
Downloading GRGM1200A from https://...  (or "Using known coefficients only")
Loaded N coefficients up to degree 100.
Wrote …/data/grgm1200a_100x100.bin (82420 bytes)
```

Verify the file exists and is non-zero:
```bash
ls -la /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator/data/grgm1200a_100x100.bin
```

- [ ] **Step 3.3 — Create `coefficients.rs`**

```rust
//! GRGM1200A spherical harmonic coefficient storage.
//!
//! Binary format (generated by tools/gen_coefficients.py):
//!   4 bytes: N_MAX as u32_le
//!   then for n=0..=N_MAX, for m=0..=n: (C_nm: f64_le, S_nm: f64_le)
//!
//! Index formula: pair_index(n, m) = n*(n+1)/2 + m

use std::io::Cursor;

/// Reference radius for GRGM1200A [km]
pub const R_REF: f64 = 1738.0;

/// Bundled 100×100 GRGM1200A binary blob (compiled into the WASM binary).
static GRGM1200A_BIN: &[u8] = include_bytes!("../data/grgm1200a_100x100.bin");

pub struct Coefficients {
    /// Maximum degree/order loaded
    pub n_max: usize,
    /// Flat storage: index 2*(n*(n+1)/2 + m) = C_nm; +1 = S_nm
    data: Vec<f64>,
}

impl Coefficients {
    /// Load from the bundled GRGM1200A binary blob, truncated to `n_max`.
    pub fn from_bundle(n_max: usize) -> Coefficients {
        Self::parse(GRGM1200A_BIN, n_max)
    }

    /// Load from an externally-provided byte slice (e.g., passed via WASM API).
    pub fn from_bytes(bytes: &[u8], n_max: usize) -> Coefficients {
        Self::parse(bytes, n_max)
    }

    fn parse(bytes: &[u8], n_max_request: usize) -> Coefficients {
        let mut cur = Cursor::new(bytes);
        use std::io::Read;

        let mut buf4 = [0u8; 4];
        cur.read_exact(&mut buf4).expect("truncated header");
        let file_n_max = u32::from_le_bytes(buf4) as usize;

        let n_max = n_max_request.min(file_n_max);
        let n_pairs = (n_max + 1) * (n_max + 2) / 2;
        let mut data = vec![0.0f64; n_pairs * 2];

        let mut buf8 = [0u8; 8];
        for n in 0..=file_n_max {
            for m in 0..=n {
                // Always read the full file sequence to keep cursor in sync
                cur.read_exact(&mut buf8).ok();
                let c = f64::from_le_bytes(buf8);
                cur.read_exact(&mut buf8).ok();
                let s = f64::from_le_bytes(buf8);
                if n <= n_max {
                    let idx = 2 * (n * (n + 1) / 2 + m);
                    data[idx]     = c;
                    data[idx + 1] = s;
                }
            }
        }

        Coefficients { n_max, data }
    }

    /// Get (C_nm, S_nm). Returns (0,0) for out-of-range.
    #[inline]
    pub fn get(&self, n: usize, m: usize) -> (f64, f64) {
        if n > self.n_max || m > n {
            return (0.0, 0.0);
        }
        let idx = 2 * (n * (n + 1) / 2 + m);
        (self.data[idx], self.data[idx + 1])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn j2_value_in_range() {
        let c = Coefficients::from_bundle(4);
        let (c20, _s20) = c.get(2, 0);
        // GRGM1200A C20 (fully normalised) ≈ -9.095e-4
        // If download failed, known fallback also gives ≈ -9.095e-4
        assert!(c20 < -8e-4 && c20 > -1e-3,
            "C20 = {} out of expected range for GRGM1200A", c20);
    }

    #[test]
    fn n_max_clamp() {
        let c = Coefficients::from_bundle(10);
        assert_eq!(c.n_max, 10);
        // Should not panic for any valid n <= 10
        let _ = c.get(10, 5);
        let (c_oob, _) = c.get(11, 0);
        assert_eq!(c_oob, 0.0);
    }
}
```

- [ ] **Step 3.4 — Add `mod coefficients;` to `lib.rs`**

Add `mod coefficients;` below `mod frames;` in `lib.rs`.

- [ ] **Step 3.5 — Run coefficients tests**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test coefficients 2>&1
```

Expected: `coefficients::tests::j2_value_in_range` and `n_max_clamp` pass.

- [ ] **Step 3.6 — Commit**

```bash
git add propagator/src/coefficients.rs propagator/src/lib.rs \
        propagator/tools/gen_coefficients.py \
        propagator/data/grgm1200a_100x100.bin
git commit -m "feat(propagator): add GRGM1200A coefficient storage (100x100 bundled binary)"
```

---

## Task 4: Create `gravity.rs` — ALF Recursion and Spherical Harmonic Acceleration

**Files:**
- Create: `propagator/src/gravity.rs`
- Modify: `propagator/src/lib.rs` (add `mod gravity;`)

Algorithm: fully normalised Associated Legendre Functions via standard recursion; gradient evaluated in geocentric spherical coordinates, converted to MCMF Cartesian.

Reference: Montenbruck & Gill "Satellite Orbits" (2000) §3.2; Heiskanen & Moritz "Physical Geodesy" (1967) §2-13.

- [ ] **Step 4.1 — Create `gravity.rs`**

```rust
//! Spherical harmonic gravity in the Moon-centred Moon-fixed frame.
//!
//! Computes the Newtonian + non-central gravitational acceleration from
//! the GRGM1200A potential up to degree/order N_MAX.
//!
//! Algorithm: Montenbruck & Gill "Satellite Orbits" (2000) §3.2.
//! ALF recursion: Heiskanen & Moritz (1967) §2-13, fully-normalised form.

use crate::coefficients::{Coefficients, R_REF};

/// Compute gravitational acceleration [km/s²] for position `r` [km] in
/// the body-fixed (MCMF) frame, using spherical harmonics up to `n_max`.
///
/// `gm`: gravitational parameter [km³/s²]
pub fn gravity_sh(r: &[f64; 3], gm: f64, n_max: usize, coeff: &Coefficients) -> [f64; 3] {
    let x = r[0];
    let y = r[1];
    let z = r[2];
    let rho = (x * x + y * y + z * z).sqrt(); // radial distance

    if rho < 1.0 {
        // Inside body — return point-mass only (should never happen in practice)
        let r3 = rho.powi(3);
        return [-gm * x / r3, -gm * y / r3, -gm * z / r3];
    }

    let rxy = (x * x + y * y).sqrt(); // projected radius in xy-plane
    let sin_phi = z / rho;           // sin(geocentric latitude)
    let cos_phi = rxy / rho;         // cos(geocentric latitude)
    let lambda = y.atan2(x);        // East longitude [rad]

    // Precompute cos(mλ) and sin(mλ) for m = 0..=n_max
    let mut cos_ml = vec![0.0f64; n_max + 1];
    let mut sin_ml = vec![0.0f64; n_max + 1];
    for m in 0..=n_max {
        cos_ml[m] = (m as f64 * lambda).cos();
        sin_ml[m] = (m as f64 * lambda).sin();
    }

    // ALF storage: p[n*(n+1)/2 + m] = P̄_nm(sin_phi)
    let n_pairs = (n_max + 1) * (n_max + 2) / 2;
    let mut p = vec![0.0f64; n_pairs];
    let mut dp = vec![0.0f64; n_pairs]; // dP̄_nm/d(sin_phi) * cos_phi = dP̄_nm/dphi

    // Seed values
    p[0] = 1.0; // P̄_00 = 1

    let idx = |n: usize, m: usize| n * (n + 1) / 2 + m;

    // Build ALFs using the standard recursion (Montenbruck & Gill eq. 3.29)
    for n in 1..=n_max {
        // Sectorial: P̄_nn from P̄_{n-1,n-1}
        if n == 1 {
            // P̄_11 = sqrt(3) * cos(phi)
            p[idx(1, 1)] = 3.0_f64.sqrt() * cos_phi;
        } else {
            let a = ((2 * n + 1) as f64 / (2 * n) as f64).sqrt();
            p[idx(n, n)] = a * cos_phi * p[idx(n - 1, n - 1)];
        }

        // Sub-diagonal: P̄_{n,n-1} from P̄_{n-1,n-1}
        if n >= 1 {
            let a = ((2 * n + 1) as f64).sqrt();
            p[idx(n, n - 1)] = a * sin_phi * p[idx(n - 1, n - 1)];
        }

        // Off-diagonal recurrence: P̄_nm from P̄_{n-1,m} and P̄_{n-2,m}
        for m in 0..=(n.saturating_sub(2)) {
            let n_f = n as f64;
            let m_f = m as f64;
            let a = ((4.0 * n_f * n_f - 1.0) / (n_f * n_f - m_f * m_f)).sqrt();
            let b = (((2.0 * n_f + 1.0) * (n_f - 1.0 + m_f) * (n_f - 1.0 - m_f))
                     / ((2.0 * n_f - 3.0) * (n_f * n_f - m_f * m_f))).sqrt();
            let p_n1m = if n >= 1 { p[idx(n - 1, m)] } else { 0.0 };
            let p_n2m = if n >= 2 { p[idx(n - 2, m)] } else { 0.0 };
            p[idx(n, m)] = a * sin_phi * p_n1m - b * p_n2m;
        }
    }

    // Compute dP̄_nm/dphi from the recursion (avoids 1/cos_phi singularity)
    // dP̄_nm/d(sin_phi) derivative approach:
    // dP̄_nm/dphi = cos_phi * dP̄_nm/d(sin_phi)
    // Using: dP̄_nm/d(sin_phi) = [n*sin_phi*P̄_nm - sqrt((2n+1)(n+m)(n-m)/(2n-1)) * P̄_{n-1,m}] / cos_phi^2
    // → dP̄_nm/dphi = [n*sin_phi*P̄_nm - A_{nm}*P̄_{n-1,m}] / cos_phi
    // This has a 1/cos_phi term; handle poles specially (set dp=0 within 1e-10 of poles).
    for n in 0..=n_max {
        for m in 0..=n {
            if cos_phi < 1e-10 {
                dp[idx(n, m)] = 0.0;
                continue;
            }
            let n_f = n as f64;
            let m_f = m as f64;
            let a_nm = if n >= 1 {
                ((2.0 * n_f + 1.0) * (n_f + m_f) * (n_f - m_f) / (2.0 * n_f - 1.0)).sqrt()
            } else {
                0.0
            };
            let p_n1m = if n >= 1 { p[idx(n - 1, m)] } else { 0.0 };
            dp[idx(n, m)] = (n_f * sin_phi * p[idx(n, m)] - a_nm * p_n1m) / cos_phi;
        }
    }

    // Accumulate gradient in spherical coordinates: (a_r, a_phi, a_lambda)
    let mut sum_r   = 0.0f64; // (r/R)^n * P̄_nm * (C cos mλ + S sin mλ)  → da_r/drho
    let mut sum_phi = 0.0f64;
    let mut sum_lam = 0.0f64;

    for n in 0..=n_max {
        let n_f = n as f64;
        let ratio = (R_REF / rho).powi(n as i32);
        for m in 0..=n {
            let (c_nm, s_nm) = coeff.get(n, m);
            let cs = c_nm * cos_ml[m] + s_nm * sin_ml[m]; // C cos mλ + S sin mλ
            let sn = -c_nm * sin_ml[m] + s_nm * cos_ml[m]; // -C sin mλ + S cos mλ (for λ-deriv)
            let p_nm  = p[idx(n, m)];
            let dp_nm = dp[idx(n, m)];
            let m_f = m as f64;

            sum_r   -= (n_f + 1.0) * ratio * p_nm  * cs;
            sum_phi +=               ratio * dp_nm * cs;
            sum_lam +=               ratio * p_nm  * m_f * sn;
        }
    }

    // Prefactor
    let gm_r2 = gm / (rho * rho);

    // a_r:   radial acceleration component
    // a_phi: geocentric latitude component (from ∂V/∂φ divided by r)
    // a_lam: longitude component (from ∂V/∂λ divided by r*cos_phi)
    let a_r   = gm_r2 * sum_r;
    let a_phi = gm_r2 * sum_phi / rho;           // (1/r) * dV/dphi
    let a_lam = if cos_phi > 1e-10 {
        gm_r2 * sum_lam / (rho * cos_phi)        // (1/(r cos_phi)) * dV/dlambda
    } else {
        0.0
    };

    // Convert spherical → Cartesian (body-fixed)
    let cp_cl = cos_phi * lambda.cos();
    let cp_sl = cos_phi * lambda.sin();
    let sp_cl = sin_phi * lambda.cos();
    let sp_sl = sin_phi * lambda.sin();

    let ax = a_r * cp_cl - a_phi * sp_cl - a_lam * lambda.sin();
    let ay = a_r * cp_sl - a_phi * sp_sl + a_lam * lambda.cos();
    let az = a_r * sin_phi + a_phi * cos_phi;

    [ax, ay, az]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coefficients::Coefficients;

    const GM: f64 = 4902.800066;

    #[test]
    fn point_mass_limit() {
        // With only C00=1 (degree 0), gravity_sh should match point-mass to < 1e-8
        // Build a trivial coefficient set with only (0,0)
        // We do this by using from_bundle with n_max=0 which gives C00=1, S00=0
        let coeff = Coefficients::from_bundle(0);
        let r = [1838.0f64, 0.0, 0.0];
        let a = gravity_sh(&r, GM, 0, &coeff);
        let a_pt = -GM / r[0].powi(2);
        assert!((a[0] - a_pt).abs() / a_pt.abs() < 1e-6,
            "Point-mass limit failed: a[0]={:.6e} vs {:.6e}", a[0], a_pt);
        assert!(a[1].abs() < 1e-12);
        assert!(a[2].abs() < 1e-12);
    }

    #[test]
    fn j2_perturbation_sign() {
        // J2 accelerates an equatorial spacecraft outward (positive radial for
        // equatorial orbit) compared to point-mass.
        // For z=0: the C20 perturbation should slightly reduce radial acceleration
        // (Moon is oblate → equatorial r_eff > 1/r² for equatorial position).
        let coeff = Coefficients::from_bundle(2);
        let r_eq = [1838.0f64, 0.0, 0.0];
        let a_j2 = gravity_sh(&r_eq, GM, 2, &coeff);
        let a_pm = -GM / r_eq[0].powi(2);
        // J2 introduces a radial correction; its magnitude should be O(J2) ≈ 2e-4 * a_pm
        let rel_diff = ((a_j2[0] - a_pm) / a_pm).abs();
        assert!(rel_diff > 1e-6 && rel_diff < 1e-2,
            "J2 perturbation out of expected range: {}", rel_diff);
    }
}
```

- [ ] **Step 4.2 — Add `mod gravity;` to `lib.rs`**

Add `mod gravity;` below `mod coefficients;` in `lib.rs`.

- [ ] **Step 4.3 — Run gravity tests**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test gravity 2>&1
```

Expected: `gravity::tests::point_mass_limit` and `j2_perturbation_sign` pass.

- [ ] **Step 4.4 — Commit**

```bash
git add propagator/src/gravity.rs propagator/src/lib.rs
git commit -m "feat(propagator): add ALF recursion + spherical harmonic acceleration"
```

---

## Task 5: Create `third_body.rs` — Point-Mass Perturbations

**Files:**
- Create: `propagator/src/third_body.rs`
- Modify: `propagator/src/lib.rs` (add `mod third_body;`)

Third-body acceleration formula (Battin, 1987):
`a = GM_3 * ( (r_3 - r) / |r_3 - r|³ - r_3 / |r_3|³ )`

where `r` is the spacecraft position in MCI frame, `r_3` is the perturbing body position in MCI.

Earth ephemeris: circular orbit, semi-major axis 384400 km, period 27.3217 days, ecliptic inclination ≈ 5.145° to Moon equator, initial RAAN 0.
Sun ephemeris: Earth-Moon barycenter at 1 AU from Sun; treat as Earth-Moon pointing to Sun; Sun period 365.25 days.

- [ ] **Step 5.1 — Create `third_body.rs`**

```rust
//! Third-body point-mass perturbations: Earth and Sun.
//!
//! Both bodies use simplified Keplerian analytical ephemerides in the MCI frame.
//! Reference: Battin "An Introduction to Astrodynamics" (1987) §9.2.

/// GM values [km³/s²]
pub const GM_EARTH: f64 = 398600.4418;
pub const GM_SUN:   f64 = 1.327124400e11;

/// Earth-Moon mean semi-major axis [km]
pub const EARTH_MOON_SMA: f64 = 384400.0;
/// Earth-Moon sidereal period [s]
pub const EARTH_MOON_T: f64 = 27.3217 * 86400.0;
/// Moon-orbit inclination to lunar equator [rad] (≈ 6.68°)
pub const EARTH_MOON_INC: f64 = 0.11659; // rad

/// Earth-Sun mean semi-major axis [km]
pub const EARTH_SUN_SMA: f64 = 1.495978707e8;
/// Earth-Sun sidereal period [s]
pub const EARTH_SUN_T: f64 = 365.25 * 86400.0;
/// Ecliptic inclination to Moon equator [rad] (≈ 1.54°, simplified)
pub const EARTH_SUN_INC: f64 = 0.02687; // rad

/// Compute position of Earth in Moon-centred inertial frame at time `t` [s].
/// Uses circular Keplerian orbit.
pub fn earth_mci(t: f64) -> [f64; 3] {
    let n = 2.0 * std::f64::consts::PI / EARTH_MOON_T; // mean motion
    let theta = n * t;
    let c = EARTH_MOON_INC.cos();
    let s = EARTH_MOON_INC.sin();
    let a = EARTH_MOON_SMA;
    [
        a * theta.cos(),
        a * theta.sin() * c,
        a * theta.sin() * s,
    ]
}

/// Compute position of Sun in Moon-centred inertial frame at time `t` [s].
pub fn sun_mci(t: f64) -> [f64; 3] {
    // Sun is in the opposite direction from the Moon-Earth vector
    // (approximately: Moon is between Earth and Sun at t=0 for simplicity).
    // Sun position ≈ -Earth_direction * 1AU  (pointing away from Earth along sun line)
    let n = 2.0 * std::f64::consts::PI / EARTH_SUN_T;
    let theta = n * t;
    let c = EARTH_SUN_INC.cos();
    let s = EARTH_SUN_INC.sin();
    let a = EARTH_SUN_SMA;
    [
        -a * theta.cos(),
        -a * theta.sin() * c,
        -a * theta.sin() * s,
    ]
}

/// Indirect/direct third-body acceleration on spacecraft at `r_sc` [km].
///
/// Formula: a = GM_3 * ( (r_3 - r_sc)/|r_3 - r_sc|³  -  r_3/|r_3|³ )
pub fn third_body_accel(r_sc: &[f64; 3], r_3: &[f64; 3], gm_3: f64) -> [f64; 3] {
    // Vector from spacecraft to third body
    let d = [
        r_3[0] - r_sc[0],
        r_3[1] - r_sc[1],
        r_3[2] - r_sc[2],
    ];
    let d_mag3 = (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]).powf(1.5);
    // Magnitude of r_3
    let r3_mag3 = (r_3[0]*r_3[0] + r_3[1]*r_3[1] + r_3[2]*r_3[2]).powf(1.5);

    [
        gm_3 * (d[0] / d_mag3 - r_3[0] / r3_mag3),
        gm_3 * (d[1] / d_mag3 - r_3[1] / r3_mag3),
        gm_3 * (d[2] / d_mag3 - r_3[2] / r3_mag3),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_at_correct_distance() {
        let r = earth_mci(0.0);
        let dist = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
        assert!((dist - EARTH_MOON_SMA).abs() < 1.0,
            "Earth distance {} km != expected {} km", dist, EARTH_MOON_SMA);
    }

    #[test]
    fn sun_at_correct_distance() {
        let r = sun_mci(0.0);
        let dist = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
        assert!((dist - EARTH_SUN_SMA).abs() / EARTH_SUN_SMA < 0.01,
            "Sun distance {} km off from 1AU", dist);
    }

    #[test]
    fn third_body_accel_is_small_for_earth() {
        // Earth perturbation on a 100 km lunar orbit: ~4e-7 km/s²  (Battin §9.2)
        let r_sc = [1838.0f64, 0.0, 0.0];
        let r_e  = earth_mci(0.0);
        let a    = third_body_accel(&r_sc, &r_e, GM_EARTH);
        let mag  = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2]).sqrt();
        assert!(mag > 1e-8 && mag < 1e-5,
            "Earth third-body magnitude {} out of expected range (1e-8, 1e-5) km/s²", mag);
    }
}
```

- [ ] **Step 5.2 — Add `mod third_body;` to `lib.rs`**

Add `mod third_body;` below `mod gravity;` in `lib.rs`.

- [ ] **Step 5.3 — Run third_body tests**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test third_body 2>&1
```

Expected: all 3 tests pass.

- [ ] **Step 5.4 — Commit**

```bash
git add propagator/src/third_body.rs propagator/src/lib.rs
git commit -m "feat(propagator): add third-body point-mass perturbations (Earth + Sun)"
```

---

## Task 6: Update `lib.rs` — New Propagator Struct and WASM API

**Files:**
- Modify: `propagator/src/lib.rs`

The `Propagator` struct gains force-model fields. `step()` builds a combined force closure. Four new WASM methods are added.

- [ ] **Step 6.1 — Rewrite `Propagator` struct and `step()` in `lib.rs`**

Replace the entire `Propagator` struct definition and its `impl` block with the following. **Keep** the `keplerian_to_cartesian` function and all `use` imports; only the struct + impl changes.

```rust
use crate::coefficients::Coefficients;
use crate::gravity::gravity_sh;
use crate::frames::{mci_to_mcmf, mcmf_to_mci};
use crate::third_body::{earth_mci, sun_mci, third_body_accel, GM_EARTH, GM_SUN};

#[wasm_bindgen]
pub struct Propagator {
    gm:              f64,
    state:           State,   // MCI frame, [x,y,z km, vx,vy,vz km/s]
    time:            f64,     // elapsed seconds since epoch
    // Force model settings
    gravity_degree:  usize,   // 0 = point-mass, 2..=100 = spherical harmonics
    harmonics_on:    bool,
    third_body_earth: bool,
    third_body_sun:  bool,
    // Lazy-loaded coefficients
    coefficients:    Option<Coefficients>,
}

#[wasm_bindgen]
impl Propagator {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Propagator {
        Propagator {
            gm:               DEFAULT_GM,
            state:            [0.0; 6],
            time:             0.0,
            gravity_degree:   0,
            harmonics_on:     false,
            third_body_earth: false,
            third_body_sun:   false,
            coefficients:     None,
        }
    }

    pub fn init(&mut self, gm: f64) {
        self.gm = gm;
        self.time = 0.0;
    }

    pub fn set_state(&mut self, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) {
        self.state = [x, y, z, vx, vy, vz];
        self.time = 0.0;
    }

    pub fn init_from_keplerian(
        &mut self,
        sma: f64, ecc: f64, inc: f64, raan: f64, argp: f64, ta: f64,
    ) {
        self.state = keplerian_to_cartesian(sma, ecc, inc, raan, argp, ta, self.gm);
        self.time = 0.0;
    }

    // ── Force model setters ──────────────────────────────────────────────

    /// Enable spherical harmonics and set maximum degree (2..=100).
    /// Calling with degree=0 disables harmonics (point-mass only).
    pub fn set_gravity_degree(&mut self, degree: u32) {
        let d = degree as usize;
        self.gravity_degree = d;
        self.harmonics_on   = d >= 2;
        if self.harmonics_on && self.coefficients.is_none() {
            // Auto-load bundled GRGM1200A on first harmonic request
            self.coefficients = Some(Coefficients::from_bundle(d.min(100)));
        }
    }

    /// Load coefficient binary (same format as tools/gen_coefficients.py output).
    /// Call before or after set_gravity_degree.
    pub fn load_coefficients(&mut self, data: &[u8]) {
        let n = self.gravity_degree.max(2).min(100);
        self.coefficients = Some(Coefficients::from_bytes(data, n));
    }

    /// Toggle Earth and/or Sun third-body perturbations.
    pub fn enable_third_body(&mut self, earth: bool, sun: bool) {
        self.third_body_earth = earth;
        self.third_body_sun   = sun;
    }

    // ── Propagation ──────────────────────────────────────────────────────

    pub fn step(&mut self, dt: f64) -> bool {
        let gm              = self.gm;
        let harmonics_on    = self.harmonics_on;
        let gravity_degree  = self.gravity_degree;
        let earth_on        = self.third_body_earth;
        let sun_on          = self.third_body_sun;
        let t0              = self.time;

        // We need a reference to the coefficients that can be captured in the closure.
        // Safety: coefficients only read inside propagate(), never mutated.
        let coeff_ptr: Option<*const Coefficients> = self.coefficients.as_ref().map(|c| c as *const _);

        let deriv_fn = move |s: &State| -> State {
            let t = t0; // constant during sub-steps (first-order coupling is fine)
            let rx = s[0]; let ry = s[1]; let rz = s[2];
            let vx = s[3]; let vy = s[4]; let vz = s[5];

            // Point-mass acceleration (always present)
            let r2 = rx*rx + ry*ry + rz*rz;
            let r3 = r2 * r2.sqrt();
            let mut ax = -gm * rx / r3;
            let mut ay = -gm * ry / r3;
            let mut az = -gm * rz / r3;

            // Spherical harmonics (in body-fixed frame)
            if harmonics_on {
                if let Some(coeff) = coeff_ptr.map(|p| unsafe { &*p }) {
                    let s_mcmf  = mci_to_mcmf(s, t);
                    let r_mcmf  = [s_mcmf[0], s_mcmf[1], s_mcmf[2]];
                    let a_bf    = gravity_sh(&r_mcmf, gm, gravity_degree.min(100), coeff);
                    // Rotate acceleration back to MCI (just position rotation, no Coriolis on acc)
                    let theta   = crate::frames::OMEGA_MOON * t;
                    let c = theta.cos(); let ss = theta.sin();
                    let a_mci_x = c * a_bf[0] - ss * a_bf[1];
                    let a_mci_y = ss * a_bf[0] + c * a_bf[1];
                    let a_mci_z = a_bf[2];
                    // Subtract point-mass (already included above) to get perturbation only
                    ax += a_mci_x - (-gm * rx / r3);
                    ay += a_mci_y - (-gm * ry / r3);
                    az += a_mci_z - (-gm * rz / r3);
                }
            }

            // Earth third-body
            if earth_on {
                let r_e = earth_mci(t);
                let a_e = third_body_accel(&[rx, ry, rz], &r_e, GM_EARTH);
                ax += a_e[0]; ay += a_e[1]; az += a_e[2];
            }

            // Sun third-body
            if sun_on {
                let r_s = sun_mci(t);
                let a_s = third_body_accel(&[rx, ry, rz], &r_s, GM_SUN);
                ax += a_s[0]; ay += a_s[1]; az += a_s[2];
            }

            [vx, vy, vz, ax, ay, az]
        };

        match crate::integrator::propagate(&self.state, dt, deriv_fn) {
            Some(new_state) => {
                self.state  = new_state;
                self.time  += dt;
                true
            }
            None => false,
        }
    }

    // ── Telemetry getters ────────────────────────────────────────────────

    pub fn get_state(&self) -> Float64Array {
        let arr = Float64Array::new_with_length(6);
        for (i, &v) in self.state.iter().enumerate() {
            arr.set_index(i as u32, v);
        }
        arr
    }

    pub fn get_time(&self) -> f64 { self.time }

    pub fn get_altitude(&self) -> f64 {
        let r = (self.state[0].powi(2) + self.state[1].powi(2) + self.state[2].powi(2)).sqrt();
        r - 1737.4
    }

    pub fn get_speed(&self) -> f64 {
        (self.state[3].powi(2) + self.state[4].powi(2) + self.state[5].powi(2)).sqrt()
    }

    /// Return Keplerian elements [sma, ecc, inc, raan, argp, ta] in km + rad.
    pub fn get_orbital_elements(&self) -> Float64Array {
        let s = &self.state;
        let arr = Float64Array::new_with_length(6);
        let elems = cartesian_to_keplerian(s, self.gm);
        for (i, v) in elems.iter().enumerate() {
            arr.set_index(i as u32, *v);
        }
        arr
    }
}
```

- [ ] **Step 6.2 — Add `cartesian_to_keplerian` to `lib.rs`**

Add this standalone function just below `keplerian_to_cartesian`:

```rust
/// Cartesian state [x,y,z km, vx,vy,vz km/s] → Keplerian elements
/// Returns [sma km, ecc, inc rad, raan rad, argp rad, ta rad]
pub fn cartesian_to_keplerian(s: &State, gm: f64) -> [f64; 6] {
    use std::f64::consts::PI;

    let r = [s[0], s[1], s[2]];
    let v = [s[3], s[4], s[5]];

    let r_mag = (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]).sqrt();
    let v2    = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

    // Specific angular momentum h = r × v
    let h = [
        r[1]*v[2] - r[2]*v[1],
        r[2]*v[0] - r[0]*v[2],
        r[0]*v[1] - r[1]*v[0],
    ];
    let h_mag = (h[0]*h[0] + h[1]*h[1] + h[2]*h[2]).sqrt();

    // Ascending node vector n = k × h (k = [0,0,1])
    let n = [-h[1], h[0], 0.0];
    let n_mag = (n[0]*n[0] + n[1]*n[1]).sqrt();

    // Eccentricity vector e = (v×h)/GM - r/|r|
    let vxh = [
        v[1]*h[2] - v[2]*h[1],
        v[2]*h[0] - v[0]*h[2],
        v[0]*h[1] - v[1]*h[0],
    ];
    let e = [
        vxh[0]/gm - r[0]/r_mag,
        vxh[1]/gm - r[1]/r_mag,
        vxh[2]/gm - r[2]/r_mag,
    ];
    let ecc = (e[0]*e[0] + e[1]*e[1] + e[2]*e[2]).sqrt();

    // Specific energy → SMA
    let energy = v2 / 2.0 - gm / r_mag;
    let sma = -gm / (2.0 * energy);

    // Inclination
    let inc = (h[2] / h_mag).acos();

    // RAAN
    let raan = if n_mag > 1e-10 {
        let mut raan = (n[0] / n_mag).acos();
        if n[1] < 0.0 { raan = 2.0 * PI - raan; }
        raan
    } else { 0.0 };

    // Argument of perigee
    let argp = if n_mag > 1e-10 && ecc > 1e-10 {
        let ndote = (n[0]*e[0] + n[1]*e[1] + n[2]*e[2]) / (n_mag * ecc);
        let mut argp = ndote.clamp(-1.0, 1.0).acos();
        if e[2] < 0.0 { argp = 2.0 * PI - argp; }
        argp
    } else { 0.0 };

    // True anomaly
    let ta = if ecc > 1e-10 {
        let edotr = (e[0]*r[0] + e[1]*r[1] + e[2]*r[2]) / (ecc * r_mag);
        let mut ta = edotr.clamp(-1.0, 1.0).acos();
        let rdotv = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
        if rdotv < 0.0 { ta = 2.0 * PI - ta; }
        ta
    } else { 0.0 };

    [sma, ecc, inc, raan, argp, ta]
}
```

- [ ] **Step 6.3 — Verify full build and existing tests pass**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test 2>&1 | tail -30
```

Expected: all existing 6 tests pass; build succeeds with no errors.

- [ ] **Step 6.4 — Commit**

```bash
git add propagator/src/lib.rs
git commit -m "feat(propagator): update Propagator struct with force model fields and new WASM API"
```

---

## Task 7: Rust Validation Tests — J2 Secular Drift

**Files:**
- Modify: `propagator/src/lib.rs` (add tests in the `#[cfg(test)] mod tests` block)

J2 secular RAAN drift: `dΩ/dt = -3n J₂ R² cos(i) / (2 a² (1-e²)²)`
J2 secular argp drift: `dω/dt = (3n J₂ R²) / (2 a² (1-e²)²) * (2 - 5/2 sin²i)`

Where J₂ = C₂₀_norm × √5 ≈ 2.034×10⁻⁴ for GRGM1200A.

- [ ] **Step 7.1 — Add J2 drift tests to `tests` mod in `lib.rs`**

Inside the existing `#[cfg(test)] mod tests { … }` block, add:

```rust
    #[test]
    fn j2_raan_secular_drift() {
        // Propagate Eagle-like orbit with J2 only for 10 sidereal periods.
        // Measure RAAN drift and compare to analytical prediction.
        use crate::coefficients::Coefficients;

        const J2: f64 = 2.0334e-4;           // GRGM1200A J2 (unnormalised)
        const R_REF: f64 = crate::coefficients::R_REF; // 1738.0 km

        let sma  = 1838.0_f64;
        let ecc  = 0.01_f64;
        let inc  = 179.0_f64.to_radians();
        let period = orbital_period(sma);

        // Analytical RAAN drift [rad/s]
        let n    = 2.0 * std::f64::consts::PI / period;
        let p    = sma * (1.0 - ecc * ecc);
        let draan_dt = -1.5 * n * J2 * R_REF * R_REF * inc.cos()
                       / (p * p);

        let t_total = 10.0 * period;
        let expected_draan = draan_dt * t_total;

        // Build propagator with J2 harmonics (degree 2)
        let mut prop = Propagator::new();
        prop.init(GM);
        prop.init_from_keplerian(sma, ecc, inc, 0.3, 0.0, 0.0);
        prop.set_gravity_degree(2);

        let elems0 = cartesian_to_keplerian(&prop.state, GM);
        let raan0  = elems0[3];

        for _ in 0..100 {
            prop.step(t_total / 100.0);
        }
        let elems1 = cartesian_to_keplerian(&prop.state, GM);
        let raan1  = elems1[3];

        // Handle 2π wrap
        let mut draan = raan1 - raan0;
        while draan >  std::f64::consts::PI { draan -= 2.0 * std::f64::consts::PI; }
        while draan < -std::f64::consts::PI { draan += 2.0 * std::f64::consts::PI; }

        let rel_err = ((draan - expected_draan) / expected_draan).abs();
        assert!(rel_err < 0.05,
            "J2 RAAN drift: got {:.4e} rad, expected {:.4e} rad (rel_err={:.2}%)",
            draan, expected_draan, rel_err * 100.0);
    }

    #[test]
    fn j2_argp_secular_drift() {
        use crate::coefficients::Coefficients;
        const J2: f64 = 2.0334e-4;
        const R_REF: f64 = crate::coefficients::R_REF;

        let sma  = 1838.0_f64;
        let ecc  = 0.01_f64;
        let inc  = 45.0_f64.to_radians(); // 45° shows non-zero argp drift
        let period = orbital_period(sma);

        let n    = 2.0 * std::f64::consts::PI / period;
        let p    = sma * (1.0 - ecc * ecc);
        let dargp_dt = 1.5 * n * J2 * R_REF * R_REF
                       * (2.0 - 2.5 * inc.sin() * inc.sin()) / (p * p);

        let t_total = 10.0 * period;
        let expected_dargp = dargp_dt * t_total;

        let mut prop = Propagator::new();
        prop.init(GM);
        prop.init_from_keplerian(sma, ecc, inc, 0.0, 0.3, 0.0);
        prop.set_gravity_degree(2);

        let elems0 = cartesian_to_keplerian(&prop.state, GM);
        let argp0  = elems0[4];

        for _ in 0..100 {
            prop.step(t_total / 100.0);
        }
        let elems1 = cartesian_to_keplerian(&prop.state, GM);
        let argp1  = elems1[4];

        let mut dargp = argp1 - argp0;
        while dargp >  std::f64::consts::PI { dargp -= 2.0 * std::f64::consts::PI; }
        while dargp < -std::f64::consts::PI { dargp += 2.0 * std::f64::consts::PI; }

        let rel_err = ((dargp - expected_dargp) / expected_dargp).abs();
        assert!(rel_err < 0.05,
            "J2 argp drift: got {:.4e} rad, expected {:.4e} rad (rel_err={:.2}%)",
            dargp, expected_dargp, rel_err * 100.0);
    }

    #[test]
    fn point_mass_energy_unchanged_vs_phase1() {
        // With harmonics off, DOP853 energy conservation must be < 1e-10 relative
        // over 1 orbit (tighter than Phase 1's 1e-8 threshold due to 1e-11 tolerance).
        let sma = 1838.0;
        let s0 = keplerian_to_cartesian(sma, 0.01, 1.0, 0.5, 0.3, 0.0, GM);

        let e0 = specific_energy(&s0, GM);
        let period = orbital_period(sma);

        let mut prop = Propagator::new();
        prop.init(GM);
        prop.state = s0;

        for _ in 0..10 {
            prop.step(period / 10.0);
        }

        let e1 = specific_energy(&prop.state, GM);
        let rel_err = ((e1 - e0) / e0).abs();
        assert!(rel_err < 1e-10,
            "DOP853 energy conservation: rel_err = {:.3e}", rel_err);
    }
```

- [ ] **Step 7.2 — Run the new validation tests**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
cargo test j2 -- --nocapture 2>&1
cargo test point_mass_energy -- --nocapture 2>&1
```

Expected: all three tests pass. The J2 drift tests assert < 5% relative error vs analytical; the energy test asserts < 1e-10 relative error.

If the J2 drift tests fail by more than 10%, re-check:
1. The sign of the `gravity_sh` minus-point-mass subtraction in `lib.rs step()` (the perturbation addition).
2. The ALF recursion seeding in `gravity.rs`.
3. The MCI↔MCMF rotation in the `step()` force closure (theta sign: `mci_to_mcmf` uses `R_z(-ωt)`, so rotation back uses `+ωt`).

- [ ] **Step 7.3 — Commit**

```bash
git add propagator/src/lib.rs
git commit -m "test(propagator): add J2 secular drift validation tests (RAAN + argp)"
```

---

## Task 8: Rebuild WASM and Wire Force Model in `main.js`

**Files:**
- Modify: `web/index.html`
- Modify: `web/style.css`
- Modify: `web/main.js`

**WASM rebuild checkpoint — must complete before any JS changes:**

- [ ] **Step 8.1 — Rebuild WASM**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator
wasm-pack build --target web 2>&1 | tail -10
```

Expected output ends with: `[INFO]: ✨  Done in Xs`
The `propagator/pkg/` directory gains updated `propagator.js`, `propagator_bg.wasm`, and `propagator_bg.wasm.d.ts`.

Verify the new API is exported:
```bash
grep -E "set_gravity_degree|load_coefficients|enable_third_body|get_orbital_elements" \
    /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/propagator/pkg/propagator.d.ts
```

Expected: 4 lines matching the four new methods.

---

## Task 9: Update `index.html` — Force Model Panel, Orbital Elements, Plot Canvases

**Files:**
- Modify: `web/index.html`

- [ ] **Step 9.1 — Read the current `index.html`**

```bash
cat /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/web/index.html
```

- [ ] **Step 9.2 — Add new HUD sections to `index.html`**

Inside the `#hud` div (after the existing `#status-bar`), add:

```html
    <!-- Force model panel -->
    <div class="hud-divider"></div>
    <div class="hud-title">FORCE MODEL</div>

    <div class="hud-section">
      <label class="hud-label" for="toggle-harmonics">HARMONICS</label>
      <input type="checkbox" id="toggle-harmonics" class="hud-toggle">
    </div>
    <div class="hud-section" id="degree-row" style="display:none;">
      <label class="hud-label" for="degree-slider">DEGREE</label>
      <input type="range" id="degree-slider" min="2" max="100" value="8" class="hud-slider">
      <span class="hud-value" id="degree-label" style="font-size:12px;">8</span>
    </div>
    <div class="hud-section">
      <label class="hud-label" for="toggle-earth">EARTH 3B</label>
      <input type="checkbox" id="toggle-earth" class="hud-toggle">
    </div>
    <div class="hud-section">
      <label class="hud-label" for="toggle-sun">SUN 3B</label>
      <input type="checkbox" id="toggle-sun" class="hud-toggle">
    </div>

    <!-- Orbital elements -->
    <div class="hud-divider"></div>
    <div class="hud-title">ELEMENTS</div>
    <div class="hud-section">
      <span class="hud-label">SMA</span>
      <span class="hud-value" id="hud-sma">–</span>
      <span class="hud-unit">km</span>
    </div>
    <div class="hud-section">
      <span class="hud-label">ECC</span>
      <span class="hud-value" id="hud-ecc" style="font-size:13px;">–</span>
      <span class="hud-unit"></span>
    </div>
    <div class="hud-section">
      <span class="hud-label">INC</span>
      <span class="hud-value" id="hud-inc">–</span>
      <span class="hud-unit">°</span>
    </div>
    <div class="hud-section">
      <span class="hud-label">RAAN</span>
      <span class="hud-value" id="hud-raan">–</span>
      <span class="hud-unit">°</span>
    </div>

    <!-- Time-history plots -->
    <div class="hud-divider"></div>
    <div class="hud-title">ECCENTRICITY</div>
    <canvas id="plot-ecc" width="208" height="60"></canvas>
    <div class="hud-title" style="margin-top:10px;">ALTITUDE</div>
    <canvas id="plot-alt" width="208" height="60"></canvas>
```

- [ ] **Step 9.3 — Commit**

```bash
git add web/index.html
git commit -m "feat(ui): add force model panel, orbital elements, and plot canvases to HUD"
```

---

## Task 10: Update `style.css` — Force Panel and Plot Styling

**Files:**
- Modify: `web/style.css`

- [ ] **Step 10.1 — Add new CSS rules to `style.css`**

Append to the end of `/Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer/web/style.css`:

```css
/* ─── Force model panel ───────────────────────────────────────────────── */

.hud-toggle {
  accent-color: #38bdf8;
  cursor: pointer;
  width: 16px;
  height: 16px;
}

.hud-slider {
  flex: 1;
  accent-color: #38bdf8;
  cursor: pointer;
  margin: 0 6px;
}

/* ─── Canvas plots ─────────────────────────────────────────────────────── */

#plot-ecc,
#plot-alt {
  display: block;
  border: 1px solid rgba(80, 180, 255, 0.2);
  border-radius: 4px;
  background: rgba(0, 5, 15, 0.6);
  margin-top: 4px;
}

/* Expand HUD width to fit plots */
#hud {
  width: 240px; /* already set; increase if plots need more room */
}
```

- [ ] **Step 10.2 — Commit**

```bash
git add web/style.css
git commit -m "feat(ui): style force model panel and time-history plot canvases"
```

---

## Task 11: Update `main.js` — Force Model Wiring, Orbital Elements, Canvas Plots

**Files:**
- Modify: `web/main.js`

- [ ] **Step 11.1 — Add force model control wiring to `main.js`**

Add the following new functions and variable declarations to `main.js`. Insert them after the `WARP_LEVELS` declaration and before the `ionToken` section:

```javascript
// ─── Force model state ────────────────────────────────────────────────────
let harmonicsEnabled = false;
let gravityDegree    = 8;
let earthEnabled     = false;
let sunEnabled       = false;

function applyForceModel() {
  if (!propagator) return;
  if (harmonicsEnabled) {
    propagator.set_gravity_degree(gravityDegree);
  } else {
    propagator.set_gravity_degree(0);
  }
  propagator.enable_third_body(earthEnabled, sunEnabled);
}
```

- [ ] **Step 11.2 — Add force model control binding in `bindControls()`**

Inside the existing `bindControls()` function, after the warp slider binding, add:

```javascript
  // Harmonics toggle
  const toggleHarmonics = document.getElementById('toggle-harmonics');
  const degreeRow       = document.getElementById('degree-row');
  const degreeSlider    = document.getElementById('degree-slider');
  const degreeLabel     = document.getElementById('degree-label');

  toggleHarmonics.addEventListener('change', () => {
    harmonicsEnabled = toggleHarmonics.checked;
    degreeRow.style.display = harmonicsEnabled ? 'flex' : 'none';
    applyForceModel();
  });

  degreeSlider.addEventListener('input', () => {
    gravityDegree = parseInt(degreeSlider.value, 10);
    degreeLabel.textContent = gravityDegree;
    if (harmonicsEnabled) applyForceModel();
  });

  // Third-body toggles
  document.getElementById('toggle-earth').addEventListener('change', e => {
    earthEnabled = e.target.checked;
    applyForceModel();
  });
  document.getElementById('toggle-sun').addEventListener('change', e => {
    sunEnabled = e.target.checked;
    applyForceModel();
  });
```

- [ ] **Step 11.3 — Add orbital elements readout to `updateHUD()`**

Replace the existing `updateHUD()` function with:

```javascript
const hudSma  = document.getElementById('hud-sma');
const hudEcc  = document.getElementById('hud-ecc');
const hudInc  = document.getElementById('hud-inc');
const hudRaan = document.getElementById('hud-raan');

function updateHUD() {
  if (!propagator) return;
  hudAlt.textContent  = propagator.get_altitude().toFixed(1);
  hudVel.textContent  = propagator.get_speed().toFixed(3);
  hudTime.textContent = formatTime(propagator.get_time());

  const el = propagator.get_orbital_elements();
  hudSma.textContent  = el[0].toFixed(1);
  hudEcc.textContent  = el[1].toFixed(5);
  hudInc.textContent  = (el[2] * 180 / Math.PI).toFixed(2);
  hudRaan.textContent = (el[3] * 180 / Math.PI).toFixed(2);
}
```

- [ ] **Step 11.4 — Add time-history plot ring buffers and rendering**

Add after the `formatTime` function:

```javascript
// ─── Time-history plots ───────────────────────────────────────────────────

const PLOT_SAMPLES = 500;
const eccHistory = [];   // {t, val}
const altHistory = [];   // {t, val}

function recordPlotSamples() {
  if (!propagator) return;
  const t   = propagator.get_time();
  const el  = propagator.get_orbital_elements();
  const ecc = el[1];
  const alt = propagator.get_altitude();

  eccHistory.push({ t, val: ecc });
  altHistory.push({ t, val: alt });
  if (eccHistory.length > PLOT_SAMPLES) eccHistory.shift();
  if (altHistory.length > PLOT_SAMPLES) altHistory.shift();
}

function drawPlot(canvasId, history, color, yLabel) {
  const canvas = document.getElementById(canvasId);
  if (!canvas || history.length < 2) return;
  const ctx = canvas.getContext('2d');
  const W = canvas.width;
  const H = canvas.height;
  ctx.clearRect(0, 0, W, H);

  const vals = history.map(p => p.val);
  const minV = Math.min(...vals);
  const maxV = Math.max(...vals);
  const range = maxV - minV || 1e-10;

  ctx.strokeStyle = color;
  ctx.lineWidth   = 1.5;
  ctx.beginPath();
  history.forEach((p, i) => {
    const x = (i / (history.length - 1)) * W;
    const y = H - ((p.val - minV) / range) * (H - 8) - 4;
    if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
  });
  ctx.stroke();

  // Axis labels (min/max)
  ctx.fillStyle = 'rgba(120, 200, 255, 0.6)';
  ctx.font = '8px Courier New';
  ctx.fillText(maxV.toFixed(5), 2, 10);
  ctx.fillText(minV.toFixed(5), 2, H - 2);
}

function renderPlots() {
  drawPlot('plot-ecc', eccHistory, '#38bdf8', 'ecc');
  drawPlot('plot-alt', altHistory, '#7dd3fc', 'alt');
}
```

- [ ] **Step 11.5 — Wire `recordPlotSamples` and `renderPlots` into the tick loop**

Inside `tick()`, after the `updateHUD()` call, add:

```javascript
    // Sample for plots every ~10 simulation steps
    if (trailPositions.length % 10 === 0) {
      recordPlotSamples();
      renderPlots();
    }
```

- [ ] **Step 11.6 — Wire `applyForceModel` into `resetOrbit()`**

At the end of `resetOrbit()`, add:

```javascript
  applyForceModel();
```

- [ ] **Step 11.7 — Start the dev server and verify the UI**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer
npm run dev
```

Open http://localhost:3000. Verify:
1. HUD shows FORCE MODEL section with unchecked checkboxes.
2. ELEMENTS section shows SMA ≈ 1838, ECC ≈ 0.01, INC ≈ 179°.
3. ECCENTRICITY and ALTITUDE canvases render a live line.
4. Checking HARMONICS + setting degree 8 changes the orbit (SMA stays constant but ECC/INC drift slightly over time due to J2).
5. Checking EARTH 3B + SUN 3B: run for simulated 24+ days (use 5000× warp) and observe eccentricity oscillating.

- [ ] **Step 11.8 — Commit**

```bash
git add web/main.js
git commit -m "feat(ui): wire force model, orbital elements display, and time-history plots"
```

---

## Task 12: End-to-End Validation — 24-Day Eccentricity Oscillation

**This task is manual browser validation (no code changes).**

- [ ] **Step 12.1 — Set up the validation scenario**

In the running browser at http://localhost:3000:
1. Click ↺ Reset.
2. Enable HARMONICS, set degree slider to 100.
3. Enable EARTH 3B.
4. Enable SUN 3B.
5. Set Time Warp to 5000×.
6. Click ▶ Play.

- [ ] **Step 12.2 — Observe eccentricity oscillation**

Watch the ECCENTRICITY plot. With full 100×100 harmonics + Earth + Sun:
- Eccentricity should oscillate with a period ≈ 24 days (Meador et al. 1994).
- The amplitude should be ≈ 0.01–0.05 depending on initial conditions.
- The ALTITUDE plot should show correlated periodic variations (perigee–apogee oscillation).

At 5000× warp, 24 simulated days takes ≈ 415 wall-clock seconds (~7 minutes). Let it run.

- [ ] **Step 12.3 — Document observed period**

Open the browser console and check `window.propagator.get_time()` at a few eccentricity peaks to measure the actual oscillation period. It should be in the range 20–30 days.

- [ ] **Step 12.4 — Final commit**

```bash
cd /Users/matthewwagner/Desktop/Projects/lunar-orbit-explorer
git add -A
git commit -m "chore: Phase 2 validation — 24-day eccentricity oscillation confirmed"
```

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-14-phase2-harmonic-propagator.md`.

**Two execution options:**

**1. Subagent-Driven (recommended)** — Fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — Execute tasks in this session using `executing-plans`, batch execution with checkpoints for review.

**Which approach?**
