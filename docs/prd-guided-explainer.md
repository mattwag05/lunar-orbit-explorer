# PRD: Lunar Orbit Explorer as a guided explainer

**Status:** ready for implementation
**Revision:** 2 (fixes warp ladder, Act 5 timing, missing getters, headline sourcing; answers section 12)
**Owner:** Iris (product and design)
**Implementer:** Claude (coding agent)
**Target repo:** `mattwag05/lunar-orbit-explorer`

---

## 0. What changed in revision 2

Four contradictions from revision 1 are fixed, the open questions are answered, and one misleading measurement is corrected. Section 9 is the definition of done.

| Was | Is now |
|---|---|
| Cockpit warp "maps onto existing `WARP_LEVELS`" (which lacks 240x/2400x/3600x) | A labeled ladder that replaces `WARP_LEVELS` entirely, every step stating its multiplier (5.9) |
| Act 5 runs 30 days live in a 25-second act | 7-day horizon, precomputed in a Web Worker and replayed (5.6) |
| Act 4 needs unlisted Rust additions | Two named getters plus an anomaly grid, each with a test (5.5, 8) |
| Rule 4 bans literals but headlines contain them | Simulation values from getters, constants templated from physics constants, mission facts in one sourced copy block (5.0, section 7 rule 4) |
| Act 5 reported ΔRAAN as evidence of drift | Separation split by direction; RAAN dropped as an artifact of a near-equatorial orbit (5.6) |

---

## 1. Why this exists

The repo is currently an instrument panel. The reference video is a story that teaches one idea at a time, with a live simulation running behind it the whole way.

The goal in one line: make Lunar Orbit Explorer feel like the video, not like a telemetry dashboard.

What that means concretely. A visitor who knows nothing about orbital mechanics should be able to watch the app for two minutes and walk away understanding why the Moon's lumpy gravity destroys low lunar orbits. Right now they would walk away knowing that a number labeled ECC reads `0.007600`.

---

## 2. What the reference does

Source: "GPS, explained by Claude Opus 5.5", `https://youtu.be/K-pgPNFcAj4`. 1920x1080, 60 fps, 168.4 seconds.

Measured facts, not impressions. These came from frame sampling plus OCR and scene analysis over the full runtime.

**There is no narration.** The video has no spoken track and only about 2 KB of subtitle data against a music bed. Every idea is carried by on-screen text and animation. That constraint is why the text density is high and why nearly every claim arrives with a number attached.

### 2.1 Measured palette

| Role | Value | Notes |
|---|---|---|
| Background | `#000208` to `#010309`, and pure `#000000` | Deep navy-black, not neutral gray |
| Primary text | `#e9cbca` | Warm ivory. 14.2 percent of pixels across sampled frames, the single dominant non-background color |
| Signal accent | cyan family near 210 degrees hue | Light, radio, travel |
| Warm accent | near 30 degrees hue | Satellite hardware, clock faces |
| Warning accent | near 345 to 0 degrees hue | Error and divergence states |

Font names were not extractable from video. What is measurable is the treatment: a small letterspaced uppercase eyebrow, a large light-weight sentence beneath it, and monospace alignment for numeric readouts.

### 2.2 The repeating device

Every act follows the same shape:

1. A small uppercase eyebrow that names the concept. `THE SATELLITES`. `THE SIGNAL`. `DISTANCE`. `THE CLOCK PROBLEM`. `UP THERE`.
2. One plain sentence underneath. "The four signals race toward you at the speed of light." "But your phone has a cheap clock."
3. Live instrumentation with units. `77.294 ms`. `19,944 km`. `+45.8 µs/day`. `10,229,999.99543 ticks`.
4. A persistent honesty chip when time is distorted. "Signal shown 1,000x slower than light (1 ms takes 1 s)." "Light slowed 58x so you can watch it (80 ms takes 4.6 s)." "Time shown 120x faster."

### 2.3 The two devices worth stealing outright

**Progressive geometry build.** Distance from one satellite gives a sphere: "you are somewhere on this sphere." Two spheres intersect in a circle. Three leave two points. The fourth resolves the clock. The video adds one shape per act and labels each step. The viewer watches the answer assemble instead of being told the answer.

**The counterfactual toggle.** Partway through, the video offers `RELATIVITY CORRECTION ON / OFF` and says "let a day pass." With correction on, both dials read the same. Switched off, they separate, and the error accumulates on screen: `0.40 days on: 6.9 km off`, then `1 day on: 17.2 km off, mostly in height`, then the conclusion, "every distance comes out 11.6 km too short." The viewer does not have to trust the narrator. They watch the wrong answer get worse.

### 2.4 The finale

The video ends by handing over control: "Now try it yourself. Move your pin. Break it. Ride a satellite (solve it yourself)." Riding a satellite opens a live cockpit with camera modes (chase, orbit, view), a time warp labeled in human units (pause, real, 4 min/s, 40 min/s, 1 h/s), and telemetry that updates as it flies.

Our warp ladder is deliberately not a copy of the video's steps. The video is a fixed edit and picks rates that suit its own shots. Section 5.9 defines ours.

### 2.5 Provenance

A footer states the data source plainly: "True simulation from the US Coast Guard GPS almanac downloaded 20 September 2026, not a live radio capture." The honest framing is part of the appeal. It says: this is real math, here is where the numbers came from, and here is what it is not.

---

## 3. Where the repo is today

The engineering behind this repo is genuinely good, and this PRD adds a presentation layer without disturbing it.

| Item | State |
|---|---|
| `propagator/src/` | 1721 lines across 6 modules: `lib.rs` 638, `integrator.rs` 379, `gravity.rs` 264, `third_body.rs` 169, `coefficients.rs` 153, `frames.rs` 118 |
| Physics | DOP853 adaptive integrator, GRGM1200A spherical harmonics at 100x100 (5151 coefficients), Earth and Sun third-body |
| Tests | 32 Rust tests (12 in `lib.rs`, 5 `coefficients`, 5 `third_body`, 4 `gravity`, 3 `frames`, 3 `integrator`) |
| WASM API | `init`, `set_state`, `init_from_keplerian`, `set_gravity_degree`, `load_coefficients`, `enable_third_body`, `step`, `get_state`, `get_time`, `get_altitude`, `get_speed`, `get_orbital_elements` |
| Default orbit | Apollo 11 LM ascent stage (Eagle). SMA 1838.13 km, ecc 0.0076, inc 179.07 degrees, period about 118.7 minutes |
| Front end | CesiumJS 1.115 with a flat gray ellipsoid, cyan glowing polyline trail, Vite 5, vanilla Canvas 2D |
| Aesthetic | Cyberpunk HUD. Courier New, `#7dd3fc` and `#38bdf8` cyan, `rgba(0,8,20,0.82)` panels, 240 px right HUD and 250 px left force panel |

The repo is one asset short of the reference and it is not a physics asset. It has real gravity data and a real integrator. What it lacks is a narrative layer that decides what the viewer looks at first.

---

## 4. The gap

| Reference | Repo today |
|---|---|
| One idea per screen | Two dense panels fighting for attention, always both visible |
| Eyebrow plus plain sentence | Labels only (`SMA`, `ECC`, `INC`, `RAAN`, `AoP`, `TA`) |
| Numbers that change because something happened | Numbers that tick continuously with no explanation of why they matter |
| Distortion disclosed ("slowed 58x") | Time warp slider with no statement of what it implies |
| Counterfactual toggle that shows failure | Force model radios that silently change the answer |
| Progressive build, 1 then 2 then 3 | All elements present at once from frame one |
| Ends with "try it yourself" | No handover, no instruction, no arc |
| Provenance footer | None |

The repo answers "what is the state of this orbit." The reference answers "why does this matter." Different jobs.

---

## 5. The act structure

This is the core of the spec. Implement these acts in order. Each act is one screen state.

Design rule that applies to all acts: the 3D view is always live and always behind the text. Nothing pauses for a slide. Text sits in a reserved column and the simulation keeps running. This is what the reference does and it is the reason the video feels alive rather than like a slideshow.

### 5.0 Where the numbers in the copy come from

This resolves the contradiction between the headline copy and section 7 rule 4. Three classes, three rules.

**Class A: simulation values.** Anything describing the current run. Period, altitude, speed, drift, separation. Read from propagator getters. Never literals in the presentation layer.

**Class B: defined constants.** The Moon's radius, the Earth distance, the Sun distance. These get templated from the same constants the physics uses, so copy and code cannot disagree:

| Copy figure | Source of truth |
|---|---|
| `1,737 km` radius | the `1737.4` constant in `lib.rs` `get_altitude` |
| `384,400 km` Earth distance | `EARTH_MOON_SMA` in `third_body.rs` |
| `150 million km` Sun distance | `EARTH_SUN_SMA` in `third_body.rs` |

Each constant is exported once into a single module the presentation layer imports. Formatting (thousands separators, rounding) happens at render.

**Class C: mission facts.** Apollo 11, the LM ascent stage, "Eagle", PFS-2. These are the only true literals. They live in one copy block, each with a source comment.

### 5.1 Act 0, the hook (0:00)

- Eyebrow: `NO GPS AT THE MOON`
- Headline: "The Moon has no satellites to guide you. So how did Eagle find its way home?"
- Visual: the Moon, dark, one small point orbiting it. Camera slowly closing.
- Interaction: after 4 seconds, advance on its own, or advance immediately on any input.

### 5.2 Act 1, where you are (0:12)

- Eyebrow: `WHERE YOU ARE`
- Headline: `The Moon is {MOON_RADIUS_KM} km in radius and has no air. Nothing up here can hear you.` (Class B, renders as "1,737 km")
- Live numbers: current altitude, current speed. (Class A)
- Visual: pull back until the Moon is a disc against black. No grid, no axes.
- Honesty chip: `Distances shown to scale. Orbit shown at 1x.`

### 5.3 Act 2, the orbit (0:28)

- Eyebrow: `THE ORBIT`
- Headline: `Eagle circled the Moon every {epoch_period_min} minutes, {epoch_alt_km} km above the ground.` (Class A, sampled from the propagator at epoch, renders as "118 minutes" and "100 km")
- Live numbers: current period, altitude, speed.
- Visual: draw the full trail once, then hold it. The trail is the answer to "what shape is this."
- Interaction: drag to orbit the camera. The first time the user drags, advance.

### 5.4 Act 3, the assumption (0:48)

- Eyebrow: `IF THE MOON WERE SMOOTH`
- Headline: "Treat the Moon as a single point of mass and this orbit never changes. Not in a day. Not in a year."
- Visual: switch the propagator to `set_gravity_degree(0)`. The trail becomes a closed ellipse that retraces itself exactly.
- Live numbers: eccentricity and altitude, both steady. Show `Δe = 0.000000` accumulating over simulated days.
- This act establishes the baseline the next act destroys.

### 5.5 Act 4, the lumps (1:10)

- Eyebrow: `THE MOON IS NOT SMOOTH`
- Headline: "The Moon is lumpy. There is extra mass buried under the maria, and it pulls."
- Visual: this is the progressive build. Tint the Moon by gravity anomaly and raise the spherical harmonic degree in steps, driven by the existing `set_gravity_degree`:
  - degree 0: smooth sphere. Label `A point of mass.`
  - degree 2: the Moon flattens slightly at the poles. Label `Two: the Moon is not round.`
  - degree 20: the maria appear as warm anomalies. Label `Twenty: the buried mass shows up.`
  - full degree: complete GRGM1200A detail. Label `{loaded_degree}: {coefficient_count} measured coefficients.` (Class A, renders as "100: 5,151 measured coefficients")
- Interaction: a single slider drives the degree. The tint and the label both follow it.
- Honesty chip: `Gravity field exaggerated {factor}x so you can see it.` Real lunar free-air anomalies are far too small to see at true scale, so the exaggeration factor is required and must be stated. A legend shows the actual data range in mGal, derived from the returned grid.

**Rust additions Act 4 requires, named:**

| Method | Returns | Purpose |
|---|---|---|
| `gravity_anomaly_grid(n_lat: u32, n_lon: u32, degree: u32) -> Float64Array` | surface gravity anomaly in mGal, row-major over a regular lat/lon grid | the tint. Regenerated when the degree slider moves. |
| `get_loaded_degree() -> u32` | the degree actually loaded from the coefficient blob (`Coefficients.n_max`) | the label |
| `get_coefficient_count() -> u32` | `(n+1)(n+2)/2` pairs for the loaded degree | the label. Renders as 5151 at degree 100. |

Each needs a test. `get_loaded_degree` and `get_coefficient_count` are trivially checkable against `from_bundle(100)`. The grid needs a test that degree 0 returns a flat field (no anomaly) and that the grid's extent matches the requested resolution.

Computing min/max or a range from a returned array in JavaScript is allowed. That is arithmetic on data the propagator produced, not synthesizing a physical value.

### 5.6 Act 5, the drift (1:40)

This is the counterfactual toggle and it is the most important act in the app. It is also the act with a real performance constraint, so read this section carefully.

- Eyebrow: `LET A WEEK PASS`
- Headline: "Now let the real gravity act, and watch what a week does to the same orbit."
- Control: `POINT MASS / REAL GRAVITY`, mirroring the reference's correction toggle.

**Horizon is 7 days, not 30.** Measured on the default orbit, one day already separates the two orbits by 89 km, and by seven days the separation is 645 km with 32 km of altitude change. Thirty days costs roughly four times the runtime for no additional insight. Checkpoints at 1, 3 and 7 days.

**Readout: separation, split by direction. Do not show RAAN for this orbit.** RAAN is not a usable number for a near-equatorial orbit, and the default orbit is 0.93 degrees from equatorial. At that tilt the ascending node is barely defined, so a small real plane change produces an enormous apparent RAAN swing. The measured ΔRAAN of 130 degrees in seven days back-solves to only about 2.1 degrees of actual plane change, amplified roughly 62x by 1/sin(i). J2 nodal regression on this orbit is about 1.2 degrees per day, so roughly 8 degrees over the horizon; the other 122 degrees is the coordinate singularity, not physics.

The readout is therefore the separation between the two spacecraft, decomposed into three components, which also fills the reference's `mostly in <direction>` slot directly:

| Component | Definition |
|---|---|
| `height` | difference in altitude above mean lunar radius |
| `along-track` | separation projected on the velocity direction of the reference run |
| `cross-track` | remainder, perpendicular to both |

State the dominant component, for example `645 km off, mostly along-track`. If an orbital-element readout is wanted instead, use the angle to closest approach measured from a fixed direction (`Ω + ω`), which stays well-defined when the orbit is nearly flat. Do not put RAAN on screen for a near-equatorial orbit, and do not use it as a validation quantity.

**Precompute, do not propagate live.** A 7-day two-propagator run measured 34.2 seconds of native ARM CPU at degree 100 with 300-second chunks. That cannot happen in a frame loop.

Requirements:

1. The drift table is produced in a **Web Worker**, using two propagator instances from the same initial state: one at degree 0, one at full degree. Diff their states at each checkpoint.
2. The act shows a **real progress indicator** while it computes. The toggle is disabled until the table is ready. Do not fake progress.
3. Budget: **10 seconds** on target hardware. If full degree cannot make the budget, compute at a reduced degree and say so in the chip: `Drift computed at degree {n} for speed.` Never hide a reduced-degree run.
4. **No propagation in the frame loop.** The act replays the precomputed table. The two trails and the readouts scrub along it.
5. The chip must also state the horizon and chunk size used.

**Measured evidence, so Claude does not have to rediscover it:**

| Configuration | Wall time | Separation | Δaltitude |
|---|---|---|---|
| 1 day, both props, 60s chunks | 14.5s | 89.1 km | +4.3 km |
| 7 days, both props, 60s chunks | 76.7s | 645.1 km | −31.8 km |
| 7 days, both props, 300s chunks | 34.2s | 645.1 km | −31.8 km |
| 30 days, single prop, 60s chunks | 323.6s | n/a | n/a |

ΔRAAN was measured alongside these runs and is deliberately not reported here. See the readout note above: it is dominated by the coordinate singularity, not by the drift.

Two findings worth acting on:

- **Chunk size is free.** 300-second chunks give results identical to 60-second chunks at half the cost, because `H_MAX` in `integrator.rs` is already 300 seconds. Requesting 60-second chunks makes the adaptive stepper restart five times per 300 seconds for no accuracy gain. Use 300.
- **Where the time goes is unconfirmed.** `gravity_sh` allocates four `Vec<f64>` per call (`cos_ml`, `sin_ml`, `p`, `dp`), and a 7-day pair run is on the order of a million RHS evaluations. That allocation pattern is a plausible dominant cost, but profile before optimizing.

**Plan on a reduced degree from the start.** The 34.2-second figure is a native ARM release build. The browser WASM build will be slower, so the 10-second budget is unlikely to be met at full degree 100. Choose the working degree up front rather than discovering the shortfall after wiring the act, and name it in the chip.

**Physics validation gate.** Before this act ships, validate the drift against an independent published result. Two rules for the check:

1. **Compare separation and altitude decay, not RAAN.** RAAN is unusable for a near-equatorial orbit, for the reason given above.
2. **Use the anchor mission's own starting orbit.** PFS-2 decayed after roughly 34 days (425 revolutions), but it flew a different orbit than Eagle. PFS-2's inclination was about 10 degrees (some sources say 11); Eagle's is 179.07 degrees, or under 1 degree from the equatorial plane. Both are low-inclination, so they are not interchangeable, and validating Eagle's altitude decay against PFS-2's lifetime is not a like-for-like comparison. Either propagate PFS-2's actual orbital elements (sourced and cited, see 5.8) and compare the resulting lifetime against ~34 days, or use a published result computed for a near-equatorial low lunar orbit.

Worth knowing for the gate: stable low lunar orbits cluster at the frozen inclinations of 27, 50, 76 and 86 degrees, where the mascon perturbations balance. Neither orbit here sits at one. PFS-2's 10 degrees is well clear of all four, and Eagle's near-equatorial plane is likewise not frozen, so both are expected to decay. That is the physics the act is teaching, and it is also why the anchor has to match the orbit being validated.

If the simulation disagrees with the published result, the simulation is wrong, not the lesson.

The 645 km separation at 7 days is a plausible size rather than an obviously wrong one. Both runs start from the same instantaneous state, so the lumpy field changes the effective orbital period, and a timing gap of roughly 4.6 seconds per orbit accumulates to that order over about 12 orbits. The gate exists because plausible is not verified.

**Conclusion line**, Class C plus a Class A number: `Apollo 16 released PFS-2 into a low lunar orbit in 1972. It was expected to last a year and a half. It fell after {days} days.` Renders with the historical figure from the sourced copy block. One source must be chosen: NASA Science says 34 days, other sources say tracked for 35. Pick one and cite it in the copy block.

### 5.7 Act 6, the other pullers (2:05)

- Eyebrow: `AND TWO MORE THINGS PULL`
- Headline: `Earth is {EARTH_MOON_SMA} km away and the Sun is {EARTH_SUN_SMA} km. Both still change this orbit.` (Class B)
- Control: checkboxes for `Earth` and `Sun`, wired to the existing `enable_third_body`.
- Visual: enable each and let the trail respond over simulated weeks.

### 5.8 Act 7, try it yourself (2:20)

- Eyebrow: `NOW TRY IT YOURSELF`
- Headline: "Move the orbit. Break it. Ride along."
- Three named actions, matching the reference's "move your pin / break it / ride a satellite":
  1. `Move the orbit` opens the Keplerian element controls as sliders (SMA, ECC, INC).
  2. `Break it` loads PFS-2's own published orbital elements (inclination from the chosen source set in 5.8) and starts the run at full degree. This is also the configuration the validation gate compares against PFS-2's recorded lifetime, so the elements must be the real ones, sourced and cited, not a generic low orbit. Show the inclination against Eagle's explicitly: neither orbit is at a frozen inclination, and both fail, over different timescales.
     Inclination convention: Eagle's 179.07 degrees is retrograde in the propagator's frame, while published Apollo figures (Apollo 11 about 1.25 degrees, PFS-2 10 or 11 degrees) are often given unsigned. Confirm each source's convention and convert before calling `init_from_keplerian`. PFS-2 is expected to land near 170 degrees, not 10, if it followed the same retrograde pattern; verify this, do not assume it.
  3. `Ride along` enters the cockpit, described next.
- Footer: `Esc to return.`

### 5.9 Act 7b, the cockpit

Reached from `Ride along`.

**Warp ladder.** This replaces `WARP_LEVELS` in `web/main.js`. The old array `[1, 10, 100, 500, 1000, 5000]` goes away. Every step labels what the viewer is actually watching:

| Label | Multiplier |
|---|---|
| `Pause` | 0 |
| `Real` | 1x |
| `1 min/s` | 60x |
| `4 min/s` | 240x |
| `20 min/s` | 1200x |
| `1 h/s` | 3600x |

Define this as a `COCKPIT_WARP` array of `{label, mult}` so no label can drift from its multiplier again.

- Camera modes: `chase`, `orbit`, `view`.
- Live telemetry in the reference's format, stacked with units: speed in km/s, height above mean surface in km, orbit period as hours and minutes.
- Honesty chip whenever warp is not `Pause` or `Real`: `Time shown {mult}x faster than real ({label}).` For example, `Time shown 240x faster than real (4 min/s).`

### 5.10 Provenance footer

Present in every act, bottom of the text column, small and low contrast:

`True simulation. Gravity from GRGM1200A ({loaded_degree}x{loaded_degree}, NASA GSFC). Propagated with DOP853. Not a live tracking feed.`

Degree comes from `get_loaded_degree()`.

---

## 6. Design system

### 6.1 Color

Measured from the reference, carried over by role.

| Token | Value | Role |
|---|---|---|
| `--bg` | `#010309` | Page background |
| `--bg-deep` | `#000000` | Behind the globe, deepest layer |
| `--ink` | `#e9cbca` | Headlines and body |
| `--ink-dim` | `rgba(233,203,202,0.62)` | Secondary body, provenance |
| `--eyebrow` | `rgba(233,203,202,0.72)` | Eyebrow text |
| `--accent-signal` | `#4fc3d9` | Orbit trail, signal, motion |
| `--accent-warm` | `#c98a4b` | Gravity anomaly, satellite hardware |
| `--accent-warn` | `#d4564a` | Divergence, decay, breaking the orbit |
| `--rule` | `rgba(233,203,202,0.14)` | Dividers |
| `--panel` | `rgba(1,3,9,0.72)` | Panel fill over the 3D scene |

The cyan signal accent is close to the repo's existing `#38bdf8` family, so the trail does not have to change much. The bigger change is the ink color. Warm ivory on navy-black replaces cool cyan on black. That single swap does most of the work of moving the app away from cyberpunk.

### 6.2 Type

- Eyebrow: 11 px, weight 600, letter-spacing 0.18em, uppercase, `--eyebrow`.
- Headline: 30 to 36 px, weight 300, line-height 1.25, `--ink`. Light weight is important. The reference's headlines are not bold.
- Body: 15 px, weight 400, line-height 1.6, `--ink-dim`.
- Numerals: monospace, tabular, right-aligned in readout rows, so digits do not shift as values change.
- No font is currently loaded. Add one display family and one monospace family. Self-host them, do not add a runtime CDN dependency.

### 6.3 Layout

Two zones, not three.

- Text column: fixed 360 px wide, left edge, vertically centered in the viewport. Holds eyebrow, headline, body, live readouts, honesty chip, provenance.
- Scene: everything else, full bleed behind.

The left force panel and right HUD both go away as persistent furniture. Their contents survive, but each control appears inside the act that needs it. The orbital elements readout moves into Act 7's `Move the orbit` panel.

Transition between acts: the text column cross-fades over 400 ms. The scene never cuts.

### 6.4 What to remove

- The 240 px right HUD and 250 px left force panel as always-visible boxes.
- The 220x60 sparkline canvases. They are too small to read and they duplicate the live readouts. Act 5's divergence readout replaces them and does the job better.
- Visible panel borders. The reference has no boxed panels at all. Text sits on the scene with a soft scrim behind it for legibility.

---

## 7. Interaction rules

1. **Never block the simulation for text.** The integrator runs continuously. Acts change what is emphasized, not whether time advances.
2. **Autoplay with escape.** Each act advances on a timer. Any click, keypress, or drag advances immediately. A visible act rail (thin, bottom of the text column) lets the user jump anywhere and return to free exploration.
3. **Disclose every distortion.** Any time warping, any exaggerated field, any reduced-degree computation, any scale change gets a chip on screen stating the factor. No exceptions.
4. **Single-source every number.** No physics value is hardcoded in the presentation layer. Simulation values come from propagator getters. Defined constants are templated from the constants the physics uses. Only documented mission facts may be literals, and they live in one copy block with a source comment. See 5.0 for the three classes.
5. **One control per act.** A control appears where it is explained and leaves when the act ends. `set_gravity_degree` belongs to Act 4, `enable_third_body` to Act 6, warps to the cockpit.
6. **Keyboard.** Left and right arrows step acts. `Escape` exits the cockpit and returns to the act rail. All controls reachable by tab.
7. **Reduce motion.** Under `prefers-reduced-motion: reduce`, disable camera easing and act autoplay. Keep cross-fades short and non-directional.
8. **Long work reports its own progress.** Anything over one second shows a real progress indicator. The Act 5 precompute is the only such job in this spec.

---

## 8. Technical constraints

- Keep the Rust physics as it is. This PRD adds a presentation layer and changes no physics.
- **Rust additions are named, not implied.** Act 4 requires `gravity_anomaly_grid`, `get_loaded_degree` and `get_coefficient_count` (5.5). Each gets a test. No other Rust change is authorized by this PRD.
- Act 5 needs two propagator instances. The current `main.js` holds one module-level `propagator`. Change it to hold two, driven from the worker.
- Act 5's precompute runs in a Web Worker. Vite supports this directly. Do not propagate on the main thread.
- Use 300-second chunks for Act 5's precompute. See 5.6 for the measurement.
- Vite 5, vanilla JS, no framework. Do not introduce React or a state library.
- CesiumJS stays. Replace the flat gray `MaterialAppearance` in `createMoonSphere` with something that can carry the Act 4 anomaly tint. Keep the existing trail entity and its glow material; only the color moves to `--accent-signal`.
- Performance targets are split, because Act 5 has two different phases:
  - Act 5 **replay**: 60 fps, measured in Chrome DevTools with the trail at `MAX_TRAIL_POINTS` of 2000.
  - Act 5 **precompute**: 10 seconds on target hardware, with a progress indicator.
- Mobile is out of scope for this pass, but the text column must not overflow below 1280 px wide. Below that, reduce the column to 300 px and the headline to 26 px.

---

## 9. Acceptance criteria

Claude, treat this list as the definition of done. Each item is checkable.

**Narrative**

- [ ] All eight acts exist and advance in order, with the eyebrow and headline copy specified in section 5, verbatim after templating.
- [ ] An act rail is visible and lets the user jump to any act.
- [ ] The cockpit is reachable from Act 7 and exits on `Escape`.
- [ ] Camera moves continuously through act transitions. There is no cut to a static frame anywhere.

**Truthfulness**

- [ ] No propagator call exists in the frame loop. Act 5 is table replay.
- [ ] Act 5's drift figures come from the worker's table, which came from two propagator instances at degree 0 and full degree. Grep the presentation layer for the drift numbers: zero hardcoded occurrences.
- [ ] Act 5's precompute completes within 10 seconds on target hardware, or the chip names the reduced degree that made the budget.
- [ ] Act 5's readout shows separation split into height, along-track and cross-track, and names the dominant component. Grep the Act 5 UI for `RAAN`: zero occurrences.
- [ ] Act 5's drift has been validated against one independent published result using separation or altitude decay, on a matching starting orbit, and the comparison is recorded (in a comment or a note next to the code). See 5.6.
- [ ] Act 4's degree and coefficient-count labels read from `get_loaded_degree()` and `get_coefficient_count()`. Grep for `5151` and `100` as literals in the presentation layer: zero occurrences.
- [ ] The Act 4 anomaly tint comes from `gravity_anomaly_grid`. No anomaly math in JavaScript.
- [ ] The copy block is the only file containing physics literals. Grep the presentation layer for `1737`, `384400`, `149597870.7`: they appear only as imports from the constants module.
- [ ] Every readout carries its unit, including the new mGal legend in Act 4.
- [ ] Every honesty chip in section 5 exists: Act 1 scale, Act 4 exaggeration factor, Act 5 degree and chunk, cockpit warp.
- [ ] The provenance footer names GRGM1200A, the loaded degree, NASA GSFC and DOP853, and states it is not a live feed.

**Physics unchanged**

- [ ] `cargo test` in `propagator/` still passes all 32 existing tests.
- [ ] Existing `lib.rs` public method signatures are unchanged. The only additions are the three named in 5.5, each with a test.

**Visual**

- [ ] Palette measured against section 6.1. Sample a rendered frame and confirm the background sits at `#010309` and the primary text at `#e9cbca`.
- [ ] No panel has a visible border. Legibility comes from the scrim.
- [ ] Headlines render at light weight. A screenshot of any headline shows a stroke noticeably thinner than the eyebrow.
- [ ] Numerals are tabular and right-aligned. Watch a readout for 10 seconds and confirm no digit shifts horizontally.

**Performance**

- [ ] 60 fps during Act 5 replay, measured in Chrome DevTools on the trail at 2000 points.
- [ ] Act 5 precompute within 10 seconds, with a progress indicator that reflects real work.
- [ ] No layout overflow at 1280 px wide.

**Accessibility**

- [ ] Every control reachable by keyboard, with a visible focus ring.
- [ ] `prefers-reduced-motion: reduce` disables camera easing and autoplay.
- [ ] Text over the scene holds at least 4.5:1 contrast against the scrim.
- [ ] The Act 5 progress indicator announces itself to assistive tech and does not trap focus.

---

## 10. Out of scope

- Mobile layout and touch gestures.
- Any change to the physics, the integrator, or the gravity data.
- Multi-spacecraft or mission-phase features from the existing Phase 2 roadmap.
- TLE import.
- Narration or audio. The reference proves text and motion are enough, and audio adds a whole production track.

## 11. Version 2

- Real STL terrain mesh so Act 4's anomalies sit on actual elevation.
- A second orbit flown side by side, a stable one against a decaying one, extending Act 5's comparison into a lasting mode.
- Linkable act URLs, so a specific act can be shared and lands the viewer where it matters.
- Higher-precision Act 5: raise the degree once the precompute is fast enough.

## 12. Decisions on the open questions

### 12.1 Act 5's conclusion names a real mission

**Apollo 16's PFS-2.** Released April 1972 into a low lunar orbit, expected to last about a year and a half, destroyed by mascon-driven orbital decay after roughly 34 days. Apollo 15's PFS-1 lasted far longer, and the contrast is not altitude: PFS-1 was at 102 x 139 km and PFS-2 at 90 x 130 km, close enough to be the same class of orbit. What differed was inclination, PFS-1 at 28.5 degrees and PFS-2 at 10 degrees (some sources say 11). PFS-1 happened to sit near the 27-degree frozen inclination where mascon perturbations balance; PFS-2 did not. That contrast is the lesson, stated as fact rather than as a claim.

Source conflicts to resolve when the copy block is written: inclination is given as 10 degrees (Gunter's Space Page, with periselene 90 km and aposelene 130 km) or 11 degrees (Wikipedia); tracked lifetime is given as 34 days or 35 days. Pick one source set, cite it, and use it consistently in both the copy and the validation gate.

### 12.2 No wordmark

The reference ends with a brand mark. This app's credibility comes from its provenance footer, not a logo. The name appears in the act rail and the footer, and nowhere else. Adding a wordmark would undercut the honesty posture the rest of this spec is built on.

### 12.3 Eagle stays the default; PFS-2 becomes the Act 7 target

The concern was that a stable default would leave Act 5 with nothing to show. The measurement answers it: Eagle already drifts hard, 645 km of separation and 32 km of altitude change in seven days. The default does not need to change.

(An earlier draft cited 130 degrees of RAAN change here. That figure is real but misleading for this orbit: it is a coordinate artifact, not drift. See 5.6.)

Instead, PFS-2 becomes the vehicle for Act 7's `Break it` control. That gives the button a real historical target, and it ties the interactive act back to the Act 5 conclusion.
