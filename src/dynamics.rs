use crate::params::{Derived, ATOL, DT, RTOL, T_SIM};
use ndarray::Array1;
use ode_solvers::dopri5::Dopri5;
use ode_solvers::{SVector, System};

pub type State = SVector<f64, 4>;

pub struct LinearInterp {
    pub y: Vec<f64>,
    pub dt: f64,
    pub n: usize,
}

impl LinearInterp {
    pub fn new(y: &Array1<f64>, dt: f64) -> Self {
        Self {
            y: y.to_vec(),
            dt,
            n: y.len(),
        }
    }

    pub fn eval(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return self.y[0];
        }
        let tn = (self.n as f64 - 1.0) * self.dt;
        if t >= tn {
            return self.y[self.n - 1];
        }
        let idx_f = t / self.dt;
        let i = idx_f.floor() as usize;
        let frac = idx_f - i as f64;
        let a = self.y[i];
        let b = self.y[i + 1];
        a + frac * (b - a)
    }
}

pub struct Problem {
    pub d: Derived,
    pub eta: LinearInterp,
    pub f3: LinearInterp,
    pub f5: LinearInterp,
}

impl System<f64, State> for Problem {
    fn system(&self, t: f64, y: &State, dy: &mut State) {
        let x3 = y[0];
        let dx3 = y[1];
        let x5 = y[2];
        let dx5 = y[3];
        let eta = self.eta.eval(t);
        let f3 = self.f3.eval(t);
        let f5 = self.f5.eval(t);
        let d = &self.d;
        let ddx3 = (f3
            - d.d33 * dx3.abs() * dx3
            - d.c33 * (x3 - eta - 0.5 * d.hg * x5 * x5))
            / d.m33a33;
        let restoring5 = crate::params::RHO * crate::params::G * d.v * d.gm * x5
            - d.kb_coupling * x3 * x5
            + d.kb_coupling * eta * x5;
        let ddx5 = (f5 - d.d55 * dx5.abs() * dx5 - restoring5) / d.m55a55;
        dy[0] = dx3;
        dy[1] = ddx3;
        dy[2] = dx5;
        dy[3] = ddx5;
    }
}

pub struct Solution {
    pub xi3: Array1<f64>,
    pub xi5: Array1<f64>,
}

pub fn integrate(problem: Problem) -> anyhow::Result<Solution> {
    let y0 = State::new(0.0, 0.0, 0.0, 0.0);
    let mut stepper = Dopri5::new(problem, 0.0, T_SIM, DT, y0, RTOL, ATOL);
    stepper
        .integrate()
        .map_err(|e| anyhow::anyhow!("ODE integration failed: {:?}", e))?;
    let y_out = stepper.y_out();
    let n = y_out.len();
    let mut xi3 = Array1::<f64>::zeros(n);
    let mut xi5 = Array1::<f64>::zeros(n);
    for (i, y) in y_out.iter().enumerate() {
        xi3[i] = y[0];
        xi5[i] = y[2];
    }
    Ok(Solution { xi3, xi5 })
}
