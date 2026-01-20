<h1 align="center"><project-name></h1>

<p align="center"><project-description></p>

## Introduction

This repository is intended as an online supplement to the working paper:

Jackson P. Lautier. Vladimir Pozdnyakov. Jun Yan. "Discrete Time-to-Event Regression Analysis Under Left-Truncation."
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
manuscript can be replicated with `aart-data-analysis.R`, `sim-study-thm-2-1.R`,
and `sim-study-thm-2-2.R`.
All results will either print in the R console or be
stored in a new folder, **results**.


## Screenshots

![Chi-Square Test](/illustrative-figures/chi-sq-sim.pdf)

![AART Application](/illustrative-figures/aart-g-est.pdf)

## Lead, Corresponding Author

**Jackson P. Lautier**

- [Website](https://jacksonlautier.com/)

## Complete Authors

**Vladimir Pozdnyakov**

- [Website](https://vladimir-pozdnyakov.github.io/)

**Jun Yan**

- [Website](http://merlot.stat.uconn.edu/~jyan/)