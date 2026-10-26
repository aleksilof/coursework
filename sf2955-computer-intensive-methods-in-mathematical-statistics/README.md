# SF2955 – Computer Intensive Methods in Mathematical Statistics

Coursework from **SF2955** at KTH Royal Institute of Technology (Spring 2026).
The course covers simulation-based methods for statistical inference: Monte Carlo
integration, importance sampling, sequential Monte Carlo (particle filters),
Markov chain Monte Carlo and Hamiltonian Monte Carlo.

All code is written in **MATLAB**.

## Assignments

| # | Assignment | Topics |
|---|---|---|
| 1 | [SMC-based mobility tracking](ha1-smc-mobility-tracking/) | Hidden Markov models, sequential importance sampling (SIS), particle filtering with resampling (SISR), likelihood-based parameter calibration |
| 2 | [MCMC for coal mine disasters and HMC](ha2-mcmc-coal-mine-and-hmc/) | Bayesian change-point model, hybrid Gibbs / Metropolis–Hastings sampler, Hamiltonian Monte Carlo vs. random-walk MH |

## Highlights

- **Particle filter tracking (HA1):** Tracked a moving target in a cellular network
  from noisy signal-strength measurements using a vectorized particle filter
  (N = 10,000 particles), then estimated an unknown noise parameter by
  maximizing the SMC-estimated log-likelihood.
- **Change-point detection (HA2):** Inferred when the rate of British coal mine
  disasters changed (1851–1963), using a hybrid MCMC sampler with a varying
  number of breakpoints.
- **Hamiltonian Monte Carlo (HA2):** Sampled a ring-shaped posterior with HMC and
  showed how it explores the support far better than random-walk Metropolis–Hastings.

## Structure

```
sf2955-computer-intensive-methods-in-mathematical-statistics/
├── ha1-smc-mobility-tracking/
└── ha2-mcmc-coal-mine-and-hmc/
```

Each folder contains:
- `README.md`: a short summary of the assignment
- `hax_report.pdf`: the full written report for home assignment number x
- `src/`: MATLAB code, with a main script that reproduces the analysis

## Running the code

Requires MATLAB. Open MATLAB in an assignment's `src/` folder
and run its main script. 
