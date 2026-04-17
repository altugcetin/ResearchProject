use std::f64::consts::PI;

pub const RHO: f64 = 1025.0;
pub const G: f64 = 9.81;

pub const D: f64 = 32.31;
pub const DRAFT: f64 = 153.92;
pub const KG: f64 = 90.39;
pub const KB: f64 = 107.69;
pub const M33_KG: f64 = 56_401_450.0;
#[allow(dead_code)]
pub const GYRADIUS: f64 = 60.96;
pub const TN3: f64 = 24.40;
pub const TN5: f64 = 39.10;
pub const A12: f64 = 0.0797;
pub const B12: f64 = 0.0696;

pub const HS: f64 = 12.20;
pub const TP: f64 = 14.0;
pub const GAMMA: f64 = 3.3;

pub const T_SIM: f64 = 10800.0;
pub const DT: f64 = 0.5;
pub const T_TRANSIENT: f64 = 2000.0;
pub const N_WAVE: usize = 200;
pub const RNG_SEED: u64 = 7;
pub const RTOL: f64 = 1e-5;
pub const ATOL: f64 = 1e-7;

pub const NPERSEG: usize = 2048;

pub const DEPTH_ATT_FACTOR: f64 = 0.35;
pub const PITCH_FORCE_COEFF: f64 = 0.08;

#[allow(dead_code)]
#[derive(Clone, Copy)]
pub struct Derived {
    pub aw: f64,
    pub hg: f64,
    pub iw: f64,
    pub v: f64,
    pub gm: f64,
    pub c33: f64,
    pub c55: f64,
    pub wn3: f64,
    pub wn5: f64,
    pub m33a33: f64,
    pub m55a55: f64,
    pub d33: f64,
    pub d55: f64,
    pub kb_coupling: f64,
    pub wp: f64,
    pub omega_min: f64,
    pub omega_max: f64,
}

pub fn derived() -> Derived {
    let aw = PI / 4.0 * D * D;
    let hg = DRAFT - KG;
    let iw = PI / 64.0 * D.powi(4);
    let v = M33_KG / RHO;
    let gm = KB - KG + iw / v;
    let c33 = RHO * G * aw;
    let c55 = RHO * G * v * gm;
    let wn3 = 2.0 * PI / TN3;
    let wn5 = 2.0 * PI / TN5;
    let m33a33 = c33 / (wn3 * wn3);
    let m55a55 = c55 / (wn5 * wn5);
    let d33 = A12 * m33a33;
    let d55 = B12 * m55a55;
    let kb_coupling = 0.5 * RHO * G * (v + 2.0 * aw * gm);
    let wp = 2.0 * PI / TP;
    let omega_min = 2.0 * PI / (3.0 * TP);
    let omega_max = 2.0 * PI / 0.8;
    Derived {
        aw,
        hg,
        iw,
        v,
        gm,
        c33,
        c55,
        wn3,
        wn5,
        m33a33,
        m55a55,
        d33,
        d55,
        kb_coupling,
        wp,
        omega_min,
        omega_max,
    }
}
