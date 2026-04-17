# Spar Heave-Pitch

Time-domain simulation of the coupled heave-pitch response of a truss spar platform subjected to random, long-crested sea states. The numerical model follows Liu et al. (2016), *Journal of Marine Science and Application*, 15: 166-174, and is implemented entirely in Rust.

## Overview

The program generates an irregular sea surface from a JONSWAP spectrum, computes the resulting hydrodynamic excitation on the platform, integrates the coupled nonlinear equations of motion for heave and pitch, and produces a single summary figure containing the time histories, the phase portrait, and the power spectral densities of both degrees of freedom.

The reference sea state corresponds to an extreme Gulf of Mexico condition with a significant wave height of 12.2 m, a peak period of 14.0 s, and a peak enhancement factor of 3.3. A three-hour simulation window is used and the first 2000 seconds are discarded as transient before any statistics or spectra are computed.

## Physical Model

- Truss spar with a circular hard tank of diameter 32.31 m and a draft of 153.92 m.
- Quadratic viscous damping in both heave and pitch, calibrated through the non-dimensional damping ratios reported in the reference paper.
- Nonlinear coupling between heave and pitch through the metacentric stiffness and through the quadratic term in the heave restoring force.
- Irregular wave excitation synthesized from 200 discrete components using the Longuet-Higgins superposition with uniformly distributed random phases and a fixed random seed for reproducibility.
- Depth attenuation applied to the pressure-induced loads at the centre of the hard tank.

The equations of motion are integrated with an adaptive Runge-Kutta 4(5) scheme (Dormand-Prince) via the `ode_solvers` crate. Wave elevation and excitation forces are precomputed on a uniform 0.5 s grid and evaluated during integration through an O(1) linear-interpolation lookup.

## Project Layout

```
Cargo.toml
src/
  main.rs        entry point, simulation pipeline, statistics
  params.rs      physical constants and derived quantities
  waves.rs       JONSWAP spectrum and wave/force time-series generation
  dynamics.rs    ODE system, linear interpolator, DOPRI5 integration
  psd.rs         Welch power spectral density (Hann window, 50% overlap)
  plotting.rs    figure composition with the plotters crate
```

## Dependencies

- `ndarray` for numerical arrays
- `ode_solvers` for Dormand-Prince RK45 integration
- `rand`, `rand_distr` for reproducible random phases
- `rustfft`, `num-complex` for the Welch PSD estimator
- `plotters` (with `ab_glyph`) for PNG figure generation
- `anyhow` for error propagation

## Building and Running

The project targets the stable Rust toolchain, Rust 2021 edition. A release build is recommended because the simulation involves a three-hour time series with 200 wave components.

```
cargo run --release
```

On successful completion the program writes `sonuc_grafikleri.png` in the project root and prints the heave and pitch statistics (mean, standard deviation, maximum) to standard output.

## Output

The generated figure contains five panels arranged on a three-row grid:

1. Heave time history over the full post-transient window, with the mean and one-sigma band overlaid.
2. Pitch time history in degrees, with the mean line overlaid.
3. Heave-pitch phase portrait, coloured by time.
4. Heave power spectral density on a logarithmic ordinate, with the heave natural frequency and the wave peak frequency marked.
5. Pitch power spectral density on a logarithmic ordinate, with the pitch natural frequency and the wave peak frequency marked.

Representative statistics produced by the reference run (seed 7):

```
Heave: mu=0.002 m  sigma=1.189 m  max=2.42 m
Pitch: mu=-0.0000 rad  sigma=0.0003 rad  max=0.0007 rad
```

## Reproducibility

All stochastic components are driven by a single `StdRng` seeded with a fixed integer, so repeated executions on the same platform yield bit-identical time series and therefore identical figures.

## Reference Hardware

The project has been validated on the following workstation:

- CPU: AMD Ryzen 9 9950X
- GPU: NVIDIA GeForce RTX 5090
- Memory: 96 GB DDR5-6000
- OS: Arch Linux

Under this configuration a full release-mode run, including figure rendering, completes in a few seconds.

## Reference

Liu, C., Li, H., Zhou, Z., & Zhang, Q. (2016). Coupled heave and pitch motions of a truss spar platform in random seas. *Journal of Marine Science and Application*, 15(2), 166-174.
