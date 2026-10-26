\# HA1 – SMC-based Mobility Tracking in Cellular Networks



Estimating the position of a moving target from the signal strengths it receives

from six base stations, using sequential Monte Carlo.



\*\*Model:\*\* A hidden Markov model. The target's position, velocity and acceleration

are driven by a Markov chain of steering commands plus Gaussian noise. Observations

follow a log-distance path-loss model with Gaussian noise.



\*\*What I did\*\*

\- Simulated trajectories from the motion model

\- Implemented SIS and showed weight degeneracy via the effective sample size

\- Implemented SISR (particle filter with resampling), fully vectorized for N = 10,000

&#x20; particles, and estimated the target's path and most probable driving commands

\- Estimated the unknown observation noise ς by maximizing the SMC log-likelihood

&#x20; over a grid



\*\*Result:\*\* <!-- e.g. "SISR tracks the target closely; estimated ς ≈ X.XX" -->



📄 \[Report](report.pdf) · 💻 \[Code](src/)

