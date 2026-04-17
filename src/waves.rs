use crate::params::{Derived, DEPTH_ATT_FACTOR, DRAFT, G, PITCH_FORCE_COEFF, RHO};
use ndarray::Array1;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use std::f64::consts::PI;

pub fn jonswap_spectrum(omega: &Array1<f64>, hs: f64, tp: f64, gamma: f64) -> Array1<f64> {
    let wp = 2.0 * PI / tp;
    let alpha = 0.0624 * hs * hs * wp.powi(4) / (G * G) * (1.0 - 0.185 * (1.9 + gamma).recip());
    let mut s = Array1::<f64>::zeros(omega.len());
    for (i, &w) in omega.iter().enumerate() {
        let sigma = if w <= wp { 0.07 } else { 0.09 };
        let r = (-0.5 * ((w - wp) / (sigma * wp)).powi(2)).exp();
        let val = alpha * G * G / w.powi(5)
            * (-1.25 * (wp / w).powi(4)).exp()
            * gamma.powf(r);
        s[i] = val;
    }
    s
}

#[allow(dead_code)]
pub struct WaveSeries {
    pub t: Array1<f64>,
    pub eta: Array1<f64>,
    pub f3: Array1<f64>,
    pub f5: Array1<f64>,
}

pub fn generate_wave_time_series(
    t: Array1<f64>,
    hs: f64,
    tp: f64,
    gamma: f64,
    n: usize,
    seed: u64,
    d: &Derived,
) -> WaveSeries {
    let dw = (d.omega_max - d.omega_min) / (n as f64 - 1.0);
    let omega = Array1::linspace(d.omega_min, d.omega_max, n);
    let s = jonswap_spectrum(&omega, hs, tp, gamma);
    let mut amp = Array1::<f64>::zeros(n);
    for i in 0..n {
        amp[i] = (2.0 * s[i] * dw).sqrt();
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut phi = Array1::<f64>::zeros(n);
    for i in 0..n {
        phi[i] = rng.gen_range(0.0..(2.0 * PI));
    }
    let mut k_arr = Array1::<f64>::zeros(n);
    let mut dk = Array1::<f64>::zeros(n);
    for i in 0..n {
        let k = omega[i] * omega[i] / G;
        k_arr[i] = k;
        dk[i] = (-k * DRAFT * DEPTH_ATT_FACTOR).exp();
    }
    let c_f3 = RHO * G * d.aw;
    let c_f5 = PITCH_FORCE_COEFF * RHO * G * d.v * d.gm;
    let nt = t.len();
    let mut eta = Array1::<f64>::zeros(nt);
    let mut f3 = Array1::<f64>::zeros(nt);
    let mut f5 = Array1::<f64>::zeros(nt);
    for j in 0..nt {
        let tj = t[j];
        let mut e = 0.0;
        let mut ff3 = 0.0;
        let mut ff5 = 0.0;
        for i in 0..n {
            let ang = omega[i] * tj + phi[i];
            let c = ang.cos();
            let sn = ang.sin();
            let ac = amp[i] * c;
            e += ac;
            ff3 += c_f3 * amp[i] * dk[i] * c;
            ff5 += c_f5 * k_arr[i] * amp[i] * dk[i] * sn;
        }
        eta[j] = e;
        f3[j] = ff3;
        f5[j] = ff5;
    }
    WaveSeries { t, eta, f3, f5 }
}
