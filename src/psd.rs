use ndarray::Array1;
use num_complex::Complex;
use rustfft::FftPlanner;
use std::f64::consts::PI;

pub fn welch(signal: &Array1<f64>, fs: f64, nperseg: usize) -> (Array1<f64>, Array1<f64>) {
    let n = signal.len();
    let noverlap = nperseg / 2;
    let step = nperseg - noverlap;
    let mut window = vec![0.0_f64; nperseg];
    for i in 0..nperseg {
        window[i] = 0.5 * (1.0 - (2.0 * PI * i as f64 / (nperseg as f64 - 1.0)).cos());
    }
    let w2_sum: f64 = window.iter().map(|x| x * x).sum();
    let nfreq = nperseg / 2 + 1;
    let mut acc = vec![0.0_f64; nfreq];
    let mut n_seg = 0usize;
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(nperseg);
    let mut start = 0usize;
    while start + nperseg <= n {
        let mut buf: Vec<Complex<f64>> = (0..nperseg)
            .map(|i| Complex::new(signal[start + i] * window[i], 0.0))
            .collect();
        fft.process(&mut buf);
        for i in 0..nfreq {
            let mag2 = buf[i].re * buf[i].re + buf[i].im * buf[i].im;
            let mut v = mag2;
            if i > 0 && i < nfreq - 1 {
                v *= 2.0;
            }
            acc[i] += v;
        }
        n_seg += 1;
        start += step;
    }
    let denom = fs * w2_sum * n_seg as f64;
    let mut pxx = Array1::<f64>::zeros(nfreq);
    let mut freqs = Array1::<f64>::zeros(nfreq);
    for i in 0..nfreq {
        pxx[i] = acc[i] / denom;
        freqs[i] = 2.0 * PI * (i as f64) * fs / nperseg as f64;
    }
    (freqs, pxx)
}
