# Project 1 – Regression Model Building

<!-- Keep the scenario you chose and delete the other. -->

**Scenario I: Predicting body fat (p < n).** Built a multiple linear regression
model that predicts body density from age and anthropometric measurements, as a
cheap alternative to underwater weighing.

- Residual analysis with scaled residuals and diagnostic plots
- Outlier, leverage and influence diagnostics (Cook's distance, COVRATIO)
- Variable transformations and multicollinearity checks (VIF)
- Variable selection with AIC, BIC, Mallows' Cp and adjusted R², evaluated by
  cross-validation
- Bootstrap confidence intervals for the coefficients

**Scenario II: Riboflavin production (p ≫ n).** Predicted riboflavin production
from the expression levels of 4,088 genes using only 71 samples.

- PCR and PLS regression with cross-validated number of components
- Ridge and Lasso: coefficient paths, CV-selected penalty, test-set MSE
- Post-selection inference for the Lasso with multi-sample splitting (`hdi`)

**Result:** <!-- e.g. final model, or "Lasso gave the lowest test MSE (X.XX)" -->

📄 [Report](project1_report.pdf) · 💻 [Code](src/)
