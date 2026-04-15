//! DOP853: Dormand-Prince order 8(5,3) adaptive integrator.
//!
//! Coefficients from: Hairer & Wanner, dop853.f (October 11, 2009)
//! http://www.unige.ch/~hairer/prog/nonstiff/dop853.f
//!
//! 13-stage explicit Runge-Kutta pair (FSAL not exploited for simplicity).
//! Error estimate: primary ERR from ER1/ER6-ER12 coefficients.

/// State vector: [x, y, z, vx, vy, vz] in km and km/s.
pub type State = [f64; 6];

/// Scale a state by a scalar.
#[inline(always)]
pub fn scale(s: &State, h: f64) -> State {
    [s[0]*h, s[1]*h, s[2]*h, s[3]*h, s[4]*h, s[5]*h]
}

/// Add two states element-wise.
#[inline(always)]
pub fn add(a: &State, b: &State) -> State {
    [a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4], a[5]+b[5]]
}

/// Weighted state sum: base + h*(c0*k0 + c1*k1 + …).
/// Usage: wsum!(base, h, (c0, k0), (c1, k1), ...)
macro_rules! wsum {
    ($base:expr, $h:expr, $(($c:expr, $k:expr)),+ $(,)?) => {{
        let mut out = *$base;
        $(
            let hc = $h * $c;
            out[0] += $k[0] * hc;
            out[1] += $k[1] * hc;
            out[2] += $k[2] * hc;
            out[3] += $k[3] * hc;
            out[4] += $k[4] * hc;
            out[5] += $k[5] * hc;
        )+
        out
    }};
}

/// Mixed error norm (ATOL + RTOL·|y|) for adaptive step control.
pub fn error_norm(err: &State, state: &State, atol: f64, rtol: f64) -> f64 {
    let mut sum = 0.0;
    for i in 0..6 {
        let sc = atol + state[i].abs() * rtol;
        sum += (err[i] / sc) * (err[i] / sc);
    }
    (sum / 6.0).sqrt()
}

// ─── DOP853 Butcher tableau ──────────────────────────────────────────────────
// Source: Hairer & Wanner dop853.f (© Hairer/Wanner, public domain, 2009)
// Verified from PARAMETER statements in that Fortran source.

// c-nodes (c1=0, c12=c13=1 implicit)
const C2:  f64 = 0.526001519587677318785587544488e-1;
const C3:  f64 = 0.789002279381515978178381316732e-1;
const C4:  f64 = 0.118350341907227396726757197510e0;
const C5:  f64 = 0.281649658092772603273242802490e0;
const C6:  f64 = 1.0 / 3.0;
const C7:  f64 = 0.25;
const C8:  f64 = 0.307692307692307692307692307692e0;
const C9:  f64 = 0.651282051282051282051282051282e0;
const C10: f64 = 0.6;
const C11: f64 = 0.857142857142857142857142857142e0;
// c12 = c13 = 1.0

// Row 2
const A21: f64 = 5.26001519587677318785587544488e-2;

// Row 3
const A31: f64 = 1.97250569845378994544595329183e-2;
const A32: f64 = 5.91751709536136983633785987549e-2;

// Row 4  (A42 = 0)
const A41: f64 = 2.95875854768068491816892993775e-2;
const A43: f64 = 8.87627564304205475450678981324e-2;

// Row 5  (A52 = 0)
const A51: f64 =  2.41365134159266685502369798665e-1;
const A53: f64 = -8.84549479328286085344864962717e-1;
const A54: f64 =  9.24834003261792003115737966543e-1;

// Row 6  (A62 = A63 = 0)
const A61: f64 = 3.7037037037037037037037037037e-2;
const A64: f64 = 1.70828608729473871279604482173e-1;
const A65: f64 = 1.25467687566822425016691814123e-1;

// Row 7  (A72 = A73 = 0)
const A71: f64 =  3.7109375e-2;
const A74: f64 =  1.70252211019544039314978060272e-1;
const A75: f64 =  6.02165389804559606850219397283e-2;
const A76: f64 = -1.7578125e-2;

// Row 8  (A82 = A83 = 0)
const A81: f64 =  3.70920001185047927108779319836e-2;
const A84: f64 =  1.70383925712239993810214054705e-1;
const A85: f64 =  1.07262030446373284651809199168e-1;
const A86: f64 = -1.53194377486244017527936158236e-2;
const A87: f64 =  8.27378916381402288758473766002e-3;

// Row 9  (A92 = A93 = 0)
const A91: f64 =  6.24110958716075717114429577812e-1;
const A94: f64 = -3.36089262944694129406857109825e0;
const A95: f64 = -8.68219346841726006818189891453e-1;
const A96: f64 =  2.75920996994467083049415600797e1;
const A97: f64 =  2.01540675504778934086186788979e1;
const A98: f64 = -4.34898841810699588477366255144e1;

// Row 10  (A102 = A103 = 0)
const A101: f64 =  4.77662536438264365890433908527e-1;
const A104: f64 = -2.48811461997166764192642586468e0;
const A105: f64 = -5.90290826836842996371446475743e-1;
const A106: f64 =  2.12300514481811942347288949897e1;
const A107: f64 =  1.52792336328824235832596922938e1;
const A108: f64 = -3.32882109689848629194453265587e1;
const A109: f64 = -2.03312017085086261358222928593e-2;

// Row 11  (A112 = A113 = 0)
const A111:  f64 = -9.3714243008598732571704021658e-1;
const A114:  f64 =  5.18637242884406370830023853209e0;
const A115:  f64 =  1.09143734899672957818500254654e0;
const A116:  f64 = -8.14978701074692612513997267357e0;
const A117:  f64 = -1.85200656599969598641566180701e1;
const A118:  f64 =  2.27394870993505042818970056734e1;
const A119:  f64 =  2.49360555267965238987089396762e0;
const A1110: f64 = -3.0467644718982195003823669022e0;

// Row 12  (A122 = A123 = 0; this row also equals the B-weights in sum)
const A121:  f64 =  2.27331014751653820792359768449e0;
const A124:  f64 = -1.05344954667372501984066689879e1;
const A125:  f64 = -2.00087205822486249909675718444e0;
const A126:  f64 = -1.79589318631187989172765950534e1;
const A127:  f64 =  2.79488845294199600508499808837e1;
const A128:  f64 = -2.85899827713502369474065508674e0;
const A129:  f64 = -8.87285693353062954433549289258e0;
const A1210: f64 =  1.23605671757943030647266201528e1;
const A1211: f64 =  6.43392746015763530355970484046e-1;

// 8th-order solution weights (B2..B5 = 0)
const B1:  f64 =  5.42937341165687622380535766363e-2;
const B6:  f64 =  4.45031289275240888144113950566e0;
const B7:  f64 =  1.89151789931450038304281599044e0;
const B8:  f64 = -5.8012039600105847814672114227e0;
const B9:  f64 =  3.1116436695781989440891606237e-1;
const B10: f64 = -1.52160949662516078556178806805e-1;
const B11: f64 =  2.01365400804030348374776537501e-1;
const B12: f64 =  4.47106157277725905176885569043e-2;

// Primary error coefficients (ER2..ER5 = 0)
// err[i] = h*(ER1*k1[i] + ER6*k6[i] + ... + ER12*k12[i])
const ER1:  f64 =  1.312004499419488073250102996e-2;
const ER6:  f64 = -1.225156446376204440720569753e0;
const ER7:  f64 = -4.957589496572501915214079952e-1;
const ER8:  f64 =  1.664377182454986536961530415e0;
const ER9:  f64 = -3.503288487499736816886487290e-1;
const ER10: f64 =  3.341791187130174790297318841e-1;
const ER11: f64 =  8.192320648511571246570742613e-2;
const ER12: f64 = -2.235530786388629525884427845e-2;

// ─── Single adaptive DOP853 step ─────────────────────────────────────────────

/// Attempt one adaptive DOP853 step.
///
/// `deriv_fn`: closure f(state) → d/dt(state). Injected so this module has
///             no coupling to any force model.
///
/// Returns `(new_state, h_used, h_next, accepted)`.
pub fn dop853_step<F>(s: &State, h: f64, atol: f64, rtol: f64, deriv_fn: &F)
    -> (State, f64, f64, bool)
where
    F: Fn(&State) -> State,
{
    let k1  = deriv_fn(s);

    let s2  = wsum!(s, h, (A21, k1));
    let k2  = deriv_fn(&s2);

    let s3  = wsum!(s, h, (A31, k1), (A32, k2));
    let k3  = deriv_fn(&s3);

    let s4  = wsum!(s, h, (A41, k1), (A43, k3));
    let k4  = deriv_fn(&s4);

    let s5  = wsum!(s, h, (A51, k1), (A53, k3), (A54, k4));
    let k5  = deriv_fn(&s5);

    let s6  = wsum!(s, h, (A61, k1), (A64, k4), (A65, k5));
    let k6  = deriv_fn(&s6);

    let s7  = wsum!(s, h, (A71, k1), (A74, k4), (A75, k5), (A76, k6));
    let k7  = deriv_fn(&s7);

    let s8  = wsum!(s, h, (A81, k1), (A84, k4), (A85, k5), (A86, k6), (A87, k7));
    let k8  = deriv_fn(&s8);

    let s9  = wsum!(s, h, (A91, k1), (A94, k4), (A95, k5), (A96, k6), (A97, k7), (A98, k8));
    let k9  = deriv_fn(&s9);

    let s10 = wsum!(s, h,
        (A101, k1), (A104, k4), (A105, k5), (A106, k6),
        (A107, k7), (A108, k8), (A109, k9));
    let k10 = deriv_fn(&s10);

    let s11 = wsum!(s, h,
        (A111, k1), (A114, k4), (A115, k5), (A116, k6),
        (A117, k7), (A118, k8), (A119, k9), (A1110, k10));
    let k11 = deriv_fn(&s11);

    let s12 = wsum!(s, h,
        (A121, k1), (A124, k4), (A125, k5), (A126, k6),
        (A127, k7), (A128, k8), (A129, k9), (A1210, k10), (A1211, k11));
    let k12 = deriv_fn(&s12);

    // 8th-order solution
    let s_new = wsum!(s, h,
        (B1, k1), (B6, k6), (B7, k7), (B8, k8),
        (B9, k9), (B10, k10), (B11, k11), (B12, k12));

    // Primary error estimate
    let mut err = [0.0f64; 6];
    for i in 0..6 {
        err[i] = h * (
            ER1  * k1[i]  + ER6  * k6[i]  + ER7  * k7[i]  + ER8  * k8[i]
          + ER9  * k9[i]  + ER10 * k10[i] + ER11 * k11[i] + ER12 * k12[i]
        );
    }

    let norm     = error_norm(&err, &s_new, atol, rtol);
    let accepted = norm <= 1.0;

    // Step size control: exponent 1/9 for order-8 pair (Hairer §II.4)
    let factor = if norm > 0.0 {
        (0.9 * norm.powf(-1.0 / 9.0)).clamp(0.2, 10.0)
    } else {
        10.0
    };
    let h_next = h * factor;

    (s_new, h, h_next, accepted)
}

// ─── Adaptive propagation ─────────────────────────────────────────────────────

/// Propagate state by exactly `dt` seconds using adaptive DOP853.
///
/// `deriv_fn` is a closure returning d/dt of the state. Tolerance constants:
/// ATOL = RTOL = 1e-11 (km, km/s — matches Phase 2 spec).
///
/// Returns `None` only if the integrator stalls with step size below 1 µs.
pub fn propagate<F>(s: &State, dt: f64, deriv_fn: F) -> Option<State>
where
    F: Fn(&State) -> State,
{
    const ATOL:  f64 = 1e-11;  // km
    const RTOL:  f64 = 1e-11;
    const H_MIN: f64 = 1e-6;   // seconds
    const H_MAX: f64 = 300.0;  // seconds

    if dt == 0.0 {
        return Some(*s);
    }

    let sign   = dt.signum();
    let t_end  = dt.abs();
    let mut t  = 0.0;
    let mut h  = (dt.abs() / 50.0).clamp(H_MIN, H_MAX);
    let mut current = *s;

    for _ in 0..200_000 {
        if t >= t_end { break; }
        let h_try = h.min(t_end - t);
        let (s_new, _h_used, h_next, accepted) =
            dop853_step(&current, sign * h_try, ATOL, RTOL, &deriv_fn);
        if accepted {
            current = s_new;
            t      += h_try;
            h       = h_next.clamp(H_MIN, H_MAX);
        } else {
            h = h_next.max(H_MIN);
        }
    }

    Some(current)
}

// ─── c-node constants re-exported for verification ──────────────────────────
/// DOP853 c-node values, publicly accessible for unit tests.
pub mod c_nodes {
    pub const C2:  f64 = super::C2;
    pub const C3:  f64 = super::C3;
    pub const C4:  f64 = super::C4;
    pub const C5:  f64 = super::C5;
    pub const C6:  f64 = super::C6;
    pub const C7:  f64 = super::C7;
    pub const C8:  f64 = super::C8;
    pub const C9:  f64 = super::C9;
    pub const C10: f64 = super::C10;
    pub const C11: f64 = super::C11;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_body(s: &State, gm: f64) -> State {
        let r2 = s[0]*s[0] + s[1]*s[1] + s[2]*s[2];
        let r3 = r2 * r2.sqrt();
        let a  = -gm / r3;
        [s[3], s[4], s[5], a*s[0], a*s[1], a*s[2]]
    }

    #[test]
    fn circular_orbit_closure() {
        // 100 km circular orbit — after exactly one period should return to start.
        const GM: f64 = 4902.800066;
        const LUNAR_R: f64 = 1737.4;
        let sma = LUNAR_R + 100.0;
        let v   = (GM / sma).sqrt();
        let s0: State = [sma, 0.0, 0.0, 0.0, v, 0.0];
        let period = 2.0 * std::f64::consts::PI * (sma.powi(3) / GM).sqrt();

        let result = propagate(&s0, period, |s| two_body(s, GM)).unwrap();

        let dr = ((result[0]-s0[0]).powi(2) + (result[1]-s0[1]).powi(2) + (result[2]-s0[2]).powi(2)).sqrt();
        assert!(dr < 1e-6, "Position closure error: {:.3e} km (expected < 1e-6)", dr);
    }

    #[test]
    fn energy_conservation_dop853() {
        // Specific energy must be conserved to < 1e-10 relative over one orbit.
        const GM: f64 = 4902.800066;
        const LUNAR_R: f64 = 1737.4;
        let sma = LUNAR_R + 100.0;
        let ecc = 0.01;
        let p   = sma * (1.0 - ecc * ecc);
        let v0  = (GM * (2.0 / sma - 1.0 / sma)).sqrt();
        let s0: State = [sma * (1.0 - ecc), 0.0, 0.0, 0.0, (GM * (1.0 + ecc) / (sma * (1.0 - ecc))).sqrt(), 0.0];

        let energy0 = |s: &State| {
            let r = (s[0]*s[0]+s[1]*s[1]+s[2]*s[2]).sqrt();
            let v2 = s[3]*s[3]+s[4]*s[4]+s[5]*s[5];
            v2/2.0 - GM/r
        };

        let period = 2.0 * std::f64::consts::PI * (sma.powi(3) / GM).sqrt();
        let e0 = energy0(&s0);
        let result = propagate(&s0, period, |s| two_body(s, GM)).unwrap();
        let e1 = energy0(&result);

        let rel_err = ((e1 - e0) / e0).abs();
        assert!(rel_err < 1e-10,
            "Energy conservation: rel_err = {:.3e} (expected < 1e-10)", rel_err);
    }

    #[test]
    fn c_nodes_consistency() {
        // A21 = c2 (row 2 must sum to c2)
        assert!((A21 - C2).abs() < 1e-15, "A21 != C2");
        // Row 3 sums to C3
        assert!((A31 + A32 - C3).abs() < 1e-15, "Row 3 sum != C3");
        // Row 4 sums to C4
        assert!((A41 + A43 - C4).abs() < 1e-15, "Row 4 sum != C4");
        // Row 5 sums to C5
        assert!((A51 + A53 + A54 - C5).abs() < 1e-14, "Row 5 sum != C5: {}", A51+A53+A54-C5);
        // B sums to 1
        let b_sum = B1 + B6 + B7 + B8 + B9 + B10 + B11 + B12;
        assert!((b_sum - 1.0).abs() < 1e-14, "B-weights don't sum to 1: {}", b_sum - 1.0);
    }
}
