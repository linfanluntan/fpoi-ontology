# OpenKBP controlled validation

Contains derived artifacts for the OpenKBP perturbation and classification-rule experiments. OpenKBP source data are not redistributed.

## C4 published-rule comparison — completed

The published-rule experiment uses the same 60 Experiment-A controlled lesions (20 high-dose-seeded, 20 gradient-seeded, 20 low-dose-seeded) and the same 10-mm spherical construction used in the manuscript.

Two published methods were executed without forcing their native semantics into a common label vocabulary:

1. **Ferreira et al. 2015, method TV** (doi:10.1186/s13014-015-0345-4). The paper defines Treated Volume as the 95% prescription isodose. Controlled lesions are InField when >95% of their volume lies in TV, Marginal for 20–95%, and OutOfField below 20%.
2. **Mohamed et al. 2016, combined centroid/fD95 typology** (doi:10.1186/s13014-016-0678-7). OpenKBP PTV70 is used as the high-dose target and PTV63/PTV56 as elective-dose targets; native A–E labels are retained.

**Li et al. 2014** (doi:10.1186/1748-717X-9-87) is explicitly marked `NOT_EXECUTABLE_FAITHFULLY`: its published rule requires primary GTV, CTV1, and CTV2, whereas OpenKBP provides PTV70/PTV63/PTV56. Substituting PTVs would change the rule semantics.

### Results

The Experiment-A reconstruction exactly reproduces the manuscript's illustrative dominant-vs-strict disagreement counts: 0/20 gradient-seeded, 19/20 high-dose-seeded, and 2/20 low-dose-seeded.

Against the published Ferreira-TV method, `dominant_fraction_v1` disagrees in 19/20 gradient-seeded, 19/20 high-dose-seeded, and 18/20 low-dose-seeded lesions: **56/60 overall (93.3%; exact 95% CI 83.8–98.2%)**.

Ferreira-TV labels all 20 high-dose-seeded lesions Marginal, 19/20 gradient-seeded lesions OutOfField, and all 20 low-dose-seeded lesions OutOfField. Mohamed's native typology labels all 20 high-dose-seeded lesions Type B (peripheral high dose) and all gradient/low-dose lesions Type E (extraneous).

These are controlled semantic stress tests, not claims that one published scheme is clinically superior.
