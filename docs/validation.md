# Validation record

This is the physics check the PRD requires before Act 5 ships (PRD 5.6,
"Physics validation gate"), plus the measurements behind the choices the app
makes. Every run below is reproducible with
`propagator/examples/validate.rs`:

```bash
cd propagator
cargo build --release --example validate
./target/release/examples/validate eagle <degree>
./target/release/examples/validate pfs2 <degree> <raan_deg> <argp_deg> [days]
./target/release/examples/validate orbit <degree> <peri_km> <apo_km> <inc_deg> <raan_deg> <argp_deg> <days>
```

All runs: native ARM64 release build on pironman (Raspberry Pi 5 class, 4
cores), 23 September 2026, 300 s chunks. Impact means altitude above the 1,737.4 km
mean radius reached zero.

## 1. A normalisation bug found and fixed first

Before any of the numbers below were taken, `gravity.rs` seeded the Legendre
recursion with P̄₁₁ = √(3/2)·cos φ. Under the 4π normalisation GRGM1200A uses,
the correct value is √3·cos φ, so every m ≥ 1 term (C₂₂, all tesserals, the
mascon signal) was scaled by 1/√2. No existing test exercised an m ≥ 1 term
against a known value. Commit `4264088` fixes it and adds three tests that fail
on the old code:

- ALF orthonormality through degree 30 (∫P̄² cos φ dφ = 2 or 4).
- The degree-2 gradient against closed-form P̄₂₀, P̄₂₁, P̄₂₂.
- The degree-30 gradient against the potential summed from the ALF table.

The drift figures in PRD revision 2 (645 km, −31.8 km at day 7) were measured
on the old code. The figures below supersede them.

## 2. Act 5: Eagle, point mass against real gravity

Default orbit (the `EAGLE` block in `web/copy.js`), no third bodies, 7 days.
Separation is split the same way as `web/drift-geometry.js`: height is the
difference in altitude, along-track is arc length in the point-mass run's
orbit plane, cross-track is the out-of-plane offset.

| Degree | Day | Apart | Height | Along-track | Cross-track | Wall time (7 days, both runs) |
|---|---|---|---|---|---|---|
| 20 | 1 | 89.2 km | +5.75 km | +88.8 km | −4.73 km | 1.1 s |
| 20 | 3 | 303.4 km | −20.53 km | +304.5 km | +12.05 km | |
| 20 | 7 | 676.4 km | −42.73 km | +686.6 km | +19.69 km | |
| 100 | 1 | 93.0 km | +6.23 km | +92.6 km | −4.60 km | 35.4 s |
| 100 | 3 | 301.1 km | −21.87 km | +302.3 km | +11.17 km | |
| 100 | 7 | 674.4 km | −44.92 km | +684.9 km | +17.30 km | |

**Working degree for Act 5: 20.** At every checkpoint it is within 4 km of
degree 100 (under 1 percent at day 7) at about 1/30 of the cost. In the browser
the whole 7-day, two-propagator table computes in 0.55 s on an M5 MacBook Pro,
well inside the 10 s budget, and the app's day 1/3/7 readouts match the degree-20
row above to the kilometre (89, 303, 676 km, all mostly along-track). The chip
names the degree.

RAAN is not reported: the orbit is 0.93° from equatorial, where RAAN is
dominated by the coordinate singularity (PRD 5.6).

## 3. Independent check: PFS-2 and PFS-1

Sources (also cited in `web/copy.js`):

- NASA NSSDCA, Apollo 16 Subsatellite, 1972-031D: periselene 90 km,
  aposelene 130 km, 10° to the lunar equator, "clockwise as viewed from north",
  impacted 29 May 1972 "after 34 days (425 revolutions)".
- NASA Science, Apollo 16 Subsatellite: "34 days in orbit rather than the
  planned one year".
- NASA NSSDCA, Apollo 15 Subsatellite, 1971-063D: perilune 102 km, apolune
  139 km, 28.5° to the lunar equator, clockwise as viewed from north; returned
  data until January 1973.
- Science@NASA, "Bizarre Lunar Orbits" (2006): PFS-2's closest approach fell
  to about 10 km within 2.5 weeks, rose again to about 30 miles, then the
  satellite crashed.

"Clockwise as viewed from north" is retrograde, so the propagator
inclinations are 180° − 10° = 170° for PFS-2 and 180° − 28.5° = 151.5° for
PFS-1. Neither source gives the node or the argument of periapsis, so both
were swept.

### PFS-2: degree 100, Earth and Sun on, 45 days

| RAAN | argp 0° | argp 90° | argp 180° | argp 270° |
|---|---|---|---|---|
| 0° | survived (min 10.4 km) | survived (min 36.7 km) | survived (min 42.5 km) | survived (min 18.7 km) |
| 90° | impact day 28.9 | impact day 34.3 | impact day 35.1 | impact day 30.0 |
| 180° | impact day 12.7 | impact day 32.3 | **impact day 34.0** | impact day 31.5 |
| 270° | impact day 35.5 | impact day 38.2 | impact day 43.5 | impact day 39.5 |

12 of 16 geometries hit the Moon within 45 days. Their median is about 34
days, against the recorded 34. The four survivors all dipped to between 10 and
43 km. Several cases also reproduce the reported shape of the decay: the low
point falls to around 10 km in the second to third week, rises again, then
falls to impact (for example RAAN 270°, argp 0°: periapsis 10.9 km at day 15,
1.4 km at day 20, 35.9 km at day 30, impact at day 35.5).

The app's `Break it` uses RAAN 180°, argp 180° (impact at 34.0 days here) and
says on screen that the node and periapsis are assumed.

### PFS-1 control: degree 100, Earth and Sun on, 90 days, argp 180°

| RAAN | Result | Lowest altitude |
|---|---|---|
| 0° | survived 90 days | 42.0 km |
| 90° | survived 90 days | 62.2 km |
| 180° | survived 90 days | 40.6 km |
| 270° | survived 90 days | 45.5 km |

All four survive, with the low point oscillating between about 40 and 115 km
and never trending toward the surface over the 90 days.

Same model, same forces, same code. The orbit that the record says kept
returning data for about 17 months stays up in this model; the one that lasted 34 days
falls in about that time.

## 4. Verdict

The gate passes. The model reproduces the one independent, published lifetime
available for a low lunar orbit of this class, reproduces the qualitative
periapsis history described for it, and keeps the control case alive. Act 5's
drift uses the same force model and code path at degree 20, which is within
1 percent of degree 100 for the Eagle orbit.

Limits: the published record does not fix PFS-2's node or periapsis, so this
is a distribution check, not a reconstruction. Impact is judged against the
mean radius, not local terrain. Third-body positions use the simplified
circular ephemerides in `third_body.rs`.
