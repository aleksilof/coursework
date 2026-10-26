# Coursework

Selected assignments from my MSc in Applied and Computational Mathematics
(Financial Mathematics) at **KTH Royal Institute of Technology**, Stockholm.

Each course has its own folder with a README describing the assignments and methods.
Reports are included as PDFs alongside the code used in the project.

## Courses

| Course | Term | Language | Topics |
|---|---|---|---|
| [SF2524 – Matrix Computations for Large-scale Systems](sf2524-matrix-computations-for-large-scale-systems/) | HT 2025 | MATLAB | Krylov methods (Arnoldi, GMRES, CG), preconditioning, QR method, matrix functions |
| [SF2930 – Regression Analysis](sf2930-regression-analysis/) | VT 2026 | R, Python | Linear regression, diagnostics, variable selection, ridge/Lasso, GLMs, insurance pricing |
| [SF2955 – Computer Intensive Methods in Mathematical Statistics](sf2955-computer-intensive-methods-in-mathematical-statistics/) | VT 2026 | MATLAB | Sequential Monte Carlo, particle filters, MCMC, Hamiltonian Monte Carlo |


<!--
To add a course, add one row to the table above:
| [CODE – Course Name](folder-name/) | HT/VT YYYY | Language | Topic, topic, topic |

Each course folder only needs a README.md. Everything else
(assignment subfolders, report.pdf, src/, figures/) is optional.
-->

## Highlights

- **Particle filter tracking** (SF2955): tracked a moving target from noisy
  cellular signal strengths with a vectorized particle filter, and calibrated
  model noise by SMC-based maximum likelihood.
- **Hamiltonian Monte Carlo** (SF2955): sampled a ring-shaped posterior where
  random-walk Metropolis–Hastings struggles.
- **Eigenfaces** (SF2524): computed the dominant eigenvectors of a 4096 × 4096
  covariance matrix with the Arnoldi method, without ever forming the matrix.
- **Preconditioned Krylov solvers** (SF2524): GMRES and PCG on a FEM heatsink
  model, compared with a direct solver.

## Repository structure

```
coursework/
├── README.md                     ← this file
├── .gitignore
├── sf2524-matrix-computations-for-large-scale-systems/
│   ├── README.md
│   └── hw1-…/ hw2-…/ hw3-…/
├── sf2930-regression-analysis/
│   ├── README.md
│   ├── project1-regression-model-building/
│   └── project2-glm-insurance-pricing/
└── sf2955-computer-intensive-methods-in-mathematical-statistics/
    ├── README.md
    └── ha1-…/ ha2-…/
```

Folder layout varies by course. See each course README for how to run its code.

## Notes

- Course-provided data files and template code are not always included.
- Assignments were designed by the respective course staff at KTH. The
  solutions, code and reports are my own (or my group's).
