mod dynamics;
mod params;
mod plotting;
mod psd;
mod waves;

use anyhow::Result;
use ndarray::Array1;

fn main() -> Result<()> {
    let d = params::derived();
    let n_t = (params::T_SIM / params::DT) as usize + 1;
    let t = Array1::linspace(0.0, params::T_SIM, n_t);
    let ws = waves::generate_wave_time_series(
        t.clone(),
        params::HS,
        params::TP,
        params::GAMMA,
        params::N_WAVE,
        params::RNG_SEED,
        &d,
    );
    let eta_i = dynamics::LinearInterp::new(&ws.eta, params::DT);
    let f3_i = dynamics::LinearInterp::new(&ws.f3, params::DT);
    let f5_i = dynamics::LinearInterp::new(&ws.f5, params::DT);
    let problem = dynamics::Problem {
        d,
        eta: eta_i,
        f3: f3_i,
        f5: f5_i,
    };
    let sol = dynamics::integrate(problem)?;
    let nt_sol = sol.xi3.len();
    let t_sol = Array1::linspace(0.0, params::DT * (nt_sol as f64 - 1.0), nt_sol);
    let mask: Vec<usize> = t_sol
        .iter()
        .enumerate()
        .filter(|(_, &v)| v > params::T_TRANSIENT)
        .map(|(i, _)| i)
        .collect();
    let ts: Array1<f64> = Array1::from(mask.iter().map(|&i| t_sol[i]).collect::<Vec<_>>());
    let h: Array1<f64> = Array1::from(mask.iter().map(|&i| sol.xi3[i]).collect::<Vec<_>>());
    let p5: Array1<f64> = Array1::from(mask.iter().map(|&i| sol.xi5[i]).collect::<Vec<_>>());
    let (hmean, hstd, _, hmax) = stats(&h);
    let (pmean, pstd, _, pmax) = stats(&p5);
    println!(
        "Heave: mu={:.3} m  sigma={:.3} m  max={:.2} m",
        hmean, hstd, hmax
    );
    println!(
        "Pitch: mu={:.4} rad  sigma={:.4} rad  max={:.4} rad",
        pmean, pstd, pmax
    );
    let fs = 1.0 / params::DT;
    let (om3, s3) = psd::welch(&h, fs, params::NPERSEG);
    let (om5, s5) = psd::welch(&p5, fs, params::NPERSEG);
    let t_minutes = ts.mapv(|v| v / 60.0);
    plotting::generate_figure(
        "sonuc_grafikleri.png",
        &t_minutes,
        &h,
        &p5,
        &om3,
        &s3,
        &om5,
        &s5,
        &d,
    )?;
    Ok(())
}

fn stats(x: &Array1<f64>) -> (f64, f64, f64, f64) {
    let n = x.len() as f64;
    let mean = x.iter().sum::<f64>() / n;
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n;
    let std = var.sqrt();
    let max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let min = x.iter().cloned().fold(f64::INFINITY, f64::min);
    (mean, std, min, max)
}
