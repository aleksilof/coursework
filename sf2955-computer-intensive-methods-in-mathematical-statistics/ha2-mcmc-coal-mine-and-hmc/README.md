\# HA2 – MCMC for Coal Mine Disasters and Hamiltonian Monte Carlo



Two Bayesian sampling problems where the posterior can't be sampled directly.



\*\*Part 1: Coal mine disasters (1851–1963).\*\* Disaster times are modelled as a Poisson

process whose rate changes at unknown breakpoints. A hybrid MCMC sampler draws from

the posterior, using Gibbs updates for the rates and hyperparameter and

Metropolis–Hastings for the breakpoints.



\*\*Part 2: Ring-shaped posterior.\*\* HMC with a leapfrog integrator is compared with

random-walk Metropolis–Hastings on a posterior concentrated on a circle.



\*\*What I did\*\*

\- Derived the full conditional posteriors and identified their distributions

\- Ran the sampler with 1–4 breakpoints and checked sensitivity to the hyperprior

&#x20; and the proposal tuning parameter ρ

\- Derived the log-posterior gradient and implemented HMC

\- Compared acceptance rates, autocorrelation and coverage of the ring for different

&#x20; step sizes ε, leapfrog steps L and random-walk scales ζ



\*\*Result:\*\* <!-- e.g. "Main rate change around 18XX; HMC covers the full ring while RW-MH gets stuck" -->



📄 \[Report](report.pdf) · 💻 \[Code](src/)

