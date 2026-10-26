# Project 2 – GLM-based Insurance Pricing

Built a price model for business travel insurance (cancellation cover) using
claims data from If P&C, covering 2018–2023.

**Model:** A multiplicative tariff, where the price is a base premium times one
risk factor per rating variable. Claim frequency (Poisson GLM) and claim
severity (Gamma GLM) are modelled separately, and the final risk factors are
their product.

**What I did**
- Grouped rating variables (e.g. number of persons, financial rating, activity
  code) into risk-homogeneous groups with enough data, and handled missing and
  odd values
- Compared candidate models with likelihood ratio tests and goodness-of-fit
  measures
- Leveled the base premium so the portfolio meets a 90% target ratio of claims
  to premium

**Result:** <!-- e.g. final rating variables and base premium γ0 -->

📄 [Report](project2_report.pdf) · 💻 [Notebook](src/)
