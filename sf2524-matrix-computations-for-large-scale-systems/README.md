# SF2524 – Matrix Computations for Large-scale Systems

Coursework from **SF2524** at KTH Royal Institute of Technology (Autumn 2025).
The course covers iterative methods for large, sparse matrices: eigenvalue
problems, linear systems and matrix functions. The focus is on understanding
convergence theory and comparing it with what happens in practice.

All code is written in **MATLAB**.

## Assignments

| # | Assignment | Topics |
|---|---|---|
| 1 | [Iterative eigenvalue methods](hw1-krylov-eigenvalue-methods/) | Power method, inverse iteration, Rayleigh quotient iteration, Arnoldi with Gram–Schmidt variants, shift-and-invert Arnoldi, eigenfaces (PCA) |
| 2 | [Krylov methods for linear systems](hw2-krylov-linear-solvers/) | GMRES, CG, CGN, convergence bounds, preconditioning on a heatsink FEM model |
| 3 | [QR method and matrix functions](hw3-qr-method-and-matrix-functions/) | Basic and shifted QR, Hessenberg reduction, Householder reflectors, Schur–Parlett, scaling-and-squaring, Estrada index centrality |

## Highlights

- **Eigenfaces (HW1):** Computed the dominant eigenvectors of a 4096 × 4096
  covariance matrix from the Olivetti faces dataset. Using the factored product
  `Xc*(Xc'*b)` avoids forming the matrix at all.
- **Preconditioned solvers (HW2):** Compared GMRES and PCG with diagonal and
  incomplete-factorization preconditioners on a FEM heatsink model, and found
  the accuracy level at which iterative methods beat a direct solve.
- **Graph centrality (HW3):** Ranked a set of pixel-art images by Estrada index,
  using the matrix exponential, and compared it with a cheap truncated Taylor
  approximation.

## Structure

```
sf2524-matrix-computations-for-large-scale-systems/
├── hw1-krylov-eigenvalue-methods/
├── hw2-krylov-linear-solvers/
└── hw3-qr-method-and-matrix-functions/
```

Each folder contains:
- `README.md`: the problem, method and key results
- `report.pdf`: the full written report
- `src/`: MATLAB code
- `figures/`: plots used in the report

## Running the code

Requires MATLAB (tested on R20XX). Open MATLAB in an assignment's `src/` folder
and run the main script for each problem, e.g.:

```matlab
run('problem2_power_method.m')
```

Some problems use data files provided by the course (e.g. `Bwedge.mat`,
`olivetti_faces.mat`, heatsink matrices). These are not included in this repo.
See each assignment's README for details.

## Acknowledgements

Assignments designed by Elias Jarlebring, KTH. Template code provided by the
course (e.g. `arnoldi.m`, `pcg_v1.m`, `schur_parlett.m`) is not included.
