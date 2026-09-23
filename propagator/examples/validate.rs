//! Native validation harness for the Act 5 drift and the PFS-2 gate.
//!
//!   cargo run --release --example validate -- eagle <degree>
//!   cargo run --release --example validate -- pfs2 <degree> <raan_deg> <argp_deg> [days]
//!   cargo run --release --example validate -- orbit <degree> <peri_km> <apo_km> <inc_deg> <raan_deg> <argp_deg> <days>
//!
//! `eagle` diffs a point-mass run against a degree-N run from the default
//! orbit for 7 days in 300 s chunks, exactly as the Act 5 worker does.
//! `pfs2` flies the NSSDC PFS-2 orbit (90 x 130 km, 10° retrograde, so 170°
//! in this frame) with Earth and Sun on and reports periapsis history and
//! the impact time, if any.

use propagator::{cartesian_to_keplerian, Propagator};
use std::time::Instant;

const GM: f64 = 4902.800066;
const R_MOON: f64 = 1737.4;
const CHUNK: f64 = 300.0;
const DAY: f64 = 86400.0;

fn sub(a: &[f64], b: &[f64]) -> [f64; 3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
fn dot(a: &[f64], b: &[f64]) -> f64 { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
fn norm(a: &[f64]) -> f64 { dot(a, a).sqrt() }
fn cross(a: &[f64], b: &[f64]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}

/// Same decomposition as web/drift-geometry.js.
fn decompose(a: &[f64; 6], b: &[f64; 6]) -> (f64, f64, f64, f64) {
    let (ra, va, rb) = (&a[0..3], &a[3..6], &b[0..3]);
    let n = cross(ra, va);
    let nn = norm(&n);
    let nh = [n[0]/nn, n[1]/nn, n[2]/nn];
    let total = norm(&sub(rb, ra));
    let height = norm(rb) - norm(ra);
    let cross_t = dot(rb, &nh);
    let rb_in = [rb[0]-cross_t*nh[0], rb[1]-cross_t*nh[1], rb[2]-cross_t*nh[2]];
    let ang = dot(&cross(ra, &rb_in), &nh).atan2(dot(ra, &rb_in));
    let along = ang * norm(ra);
    (total, height, along, cross_t)
}

fn eagle(degree: u32) {
    let (sma, ecc, inc, raan, argp) = (1838.13, 0.0076, 179.07f64, 183.41f64, 179.86f64);
    let mut a = Propagator::new();
    let mut b = Propagator::new();
    for p in [&mut a, &mut b] {
        p.init(GM);
        p.init_from_keplerian(sma, ecc, inc.to_radians(), raan.to_radians(), argp.to_radians(), 0.0);
    }
    b.set_gravity_degree(degree);
    let t0 = Instant::now();
    let steps = (7.0 * DAY / CHUNK) as usize;
    for k in 1..=steps {
        a.step(CHUNK);
        b.step(CHUNK);
        let t = k as f64 * CHUNK;
        if [1.0, 3.0, 7.0].iter().any(|d| (t - d * DAY).abs() < 1.0) {
            let (sa, sb) = (a.state_vec(), b.state_vec());
            let (tot, h, al, cr) = decompose(&sa, &sb);
            println!("deg {degree:>3} day {:>3.0}: sep {tot:8.1} km  height {h:+8.2}  along {al:+9.1}  cross {cr:+8.2}  alt_real {:7.2} alt_pm {:7.2}",
                t / DAY, b.get_altitude(), a.get_altitude());
        }
    }
    println!("deg {degree:>3} wall {:.1} s", t0.elapsed().as_secs_f64());
}

fn pfs2(degree: u32, raan: f64, argp: f64, days: f64) {
    // NSSDC 1972-031D: periselene 90 km, aposelene 130 km, 10° to the lunar
    // equator, clockwise as viewed from north (retrograde) → 170°.
    orbit("pfs2", degree, 90.0, 130.0, 170.0, raan, argp, days);
}

fn orbit(tag: &str, degree: u32, peri: f64, apo: f64, inc: f64, raan: f64, argp: f64, days: f64) {
    let rp = R_MOON + peri;
    let ra = R_MOON + apo;
    let sma = (rp + ra) / 2.0;
    let ecc = (ra - rp) / (ra + rp);
    let mut p = Propagator::new();
    p.init(GM);
    p.init_from_keplerian(sma, ecc, inc.to_radians(), raan.to_radians(), argp.to_radians(), 0.0);
    p.set_gravity_degree(degree);
    p.enable_third_body(true, true);
    let t0 = Instant::now();
    let steps = (days * DAY / CHUNK) as usize;
    let mut min_day = f64::MAX;
    let mut peri_log = Vec::new();
    for k in 1..=steps {
        p.step(CHUNK);
        let t = k as f64 * CHUNK;
        let alt = p.get_altitude();
        if alt <= 0.0 {
            println!("{tag} deg {degree} raan {raan:>5.1} argp {argp:>5.1}: IMPACT day {:.2}  wall {:.1}s  peri(d5,10,15,20,25,30) {:?}",
                t / DAY, t0.elapsed().as_secs_f64(), peri_log);
            return;
        }
        if (t % (5.0 * DAY)).abs() < 1.0 {
            let el = cartesian_to_keplerian(&p.state_vec(), GM);
            peri_log.push(((el[0] * (1.0 - el[1]) - R_MOON) * 10.0).round() / 10.0);
        }
        min_day = min_day.min(alt);
    }
    println!("{tag} deg {degree} raan {raan:>5.1} argp {argp:>5.1}: survived {days} d, min alt {min_day:.1} km  wall {:.1}s  peri(d5..) {:?}",
        t0.elapsed().as_secs_f64(), peri_log);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("eagle") => eagle(args[2].parse().unwrap()),
        Some("pfs2") => pfs2(
            args[2].parse().unwrap(),
            args[3].parse().unwrap(),
            args[4].parse().unwrap(),
            args.get(5).map_or(45.0, |d| d.parse().unwrap()),
        ),
        Some("orbit") => {
            let f = |i: usize| args[i].parse::<f64>().unwrap();
            orbit("orbit", args[2].parse().unwrap(), f(3), f(4), f(5), f(6), f(7), f(8));
        }
        _ => eprintln!("usage: validate eagle <deg> | pfs2 <deg> <raan> <argp> [days] | orbit <deg> <peri> <apo> <inc> <raan> <argp> <days>"),
    }
}
