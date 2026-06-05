<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## Introduction

This repository is intended as an online supplement to the working paper:

Lautier, J. P., Pozdnyakov, V., Yan, J. (2026) "Discrete Time-to-Event Regression Analysis Under Left-Truncation with Applications to Consumer Finance."
(see [https://jacksonlautier.com/publications](https://jacksonlautier.com/publications)
for current working papers)

Please attribute any citations of this repository to the original
manuscript.

This repository includes:

- **raw-data** Scraped loan demographic and performance data from the ABS bond
AART 2017-3.

- **data-clean** Cleaned raw data into files used within the manuscript.  These
files are identical to the files created by `data-processing.R' in the **code**
folder.

- **code** Replication code files.  First run `aart-data-processing.R` to create the
clean data files in a new folder, **processed-data** (alternatively, rename the
**data-clean** folder as **processed-data**).  Second, all results in the
manuscript can be replicated using the files listed below.
All results will either print in the R console or be
stored in a new folder, **results**.

| Coding File | Description |
| --- | --- |
| `aart-data-analysis-36mo.R` | Replication code to reproduce AART 2017-3 application results |
| `aart-data-analysis-60mo.R` | Replication code to reproduce AART 2017-3 application results |
| `aart-data-processing.R` | Replication code to process **raw-data** into **data-clean**|
| `binom.p.est.aoas.R` | Supporting file left-truncation results |
| `censoring.h.star.param.est.R` | Supporting file left-truncation plus right-censoring results |
| `censoring.hessian.std.error.R` | Supporting file left-truncation plus right-censoring results |
| `default_time.R` | Supporting file for data processing |
| `h.star.param.est.R` | Supporting file left-truncation results |
| `h-dist-cens-simulation-function-beta-parameters.R` | Supporting file left-truncation plus right-censoring results |
| `hessian.std.error.R` | Supporting file left-truncation results |
| `h-star-simulation-function-beta-parameters.R` | Supporting file left-truncation results |
| `proportional-hazards-assumption.R` | Replication code to illustrate prop hazards assumption unfit |
| `sim-study-lt-chi-sq.R` | Replication code for simulation study, left-truncation |
| `sim-study-lt-point-est-se.R` | Replication code for simulation study, left-truncation |
| `sim-study-rc-chi-sq.R` | Replication code for simulation study, left-truncation plus right-censoring |
| `sim-study-rc-point-est-se.R` | Replication code for simulation study, left-truncation plus right-censoring |

## Screenshots

![Chi-Square Sim Study](/illustrative-figures/chi-sq-sim-cens.pdf)

![AART Application](/illustrative-figures/aart-g-est.pdf)

![Proportional Hazards Assumption Review](/illustrative-figures/aart-prop-haz.pdf)

## Lead, Corresponding Author

**Jackson P. Lautier**

- [Website](https://jacksonlautier.com/)

## Complete Authors

**Vladimir Pozdnyakov**

- [Website](https://vladimir-pozdnyakov.github.io/)

**Jun Yan**

- [Website](http://merlot.stat.uconn.edu/~jyan/)