# SF2930 – Regression Analysis

Coursework from **SF2930** at KTH Royal Institute of Technology (Spring 2026).
The course covers linear regression and its extensions: model diagnostics,
variable selection, regularization, high-dimensional regression and generalized
linear models (GLMs).

Code is written in **R** (Project 1) and **Python** (Project 2).

## Projects

| # | Project | Topics |
|---|---|---|
| 1 | [Regression model building](project1-regression-model-building/) | Residual analysis, outliers and influence, multicollinearity, variable selection, cross-validation, bootstrap / PCR, PLS, ridge, Lasso, post-selection inference |
| 2 | [GLM-based insurance pricing](project2-glm-insurance-pricing/) | Poisson and Gamma GLMs, claim frequency and severity, risk-factor grouping, likelihood ratio tests, premium leveling |

## Highlights

- **Model building (Project 1):** Built and validated a regression model from
  start to finish, from diagnostics and transformations to variable selection
  and bootstrap-based assessment.
- **Insurance pricing (Project 2):** Built a multiplicative tariff for business
  travel insurance on real claims data from If P&C, modelling claim frequency
  and severity with separate GLMs and leveling the base premium to a 90% loss
  ratio target.

## Structure

```
sf2930-regression-analysis/
├── project1-regression-model-building/
└── project2-glm-insurance-pricing/
```

Each folder contains a short `README.md`, the report (`projectx_report.pdf`) and the code.

## Acknowledgements

Projects designed by the SF2930 course staff, KTH. Project 2 was set in
collaboration with If P&C Insurance. The datasets are not included in this repo.
Completed in a group of two.
