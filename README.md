This repository contains the R code and Stan files used for the simulation study in the manuscript:

> **"A comparative analysis of Phase I dose-finding designs incorporating pharmacokinetics information."**

The study evaluates the performance of various model-based Bayesian dose-finding designs, including methods that explicitly integrate pharmacokinetic (PK) data, under a variety of simulation scenarios.

---

## 📦 Contents

- Scripts to run all simulation scenarios
- Functions to implement Bayesian PK dose-finding models
- Tools for generating figures and summarizing results for the manuscript

---

## 📚 Sources for code implementation

- **TITE-PK implementation** adapted from the code by Günhan et al. (2021): https://github.com/gunhanb/TITEPK_code
- **PKLOGIT implementation** inspired by the `dfpk` package by Toumazi et al. (2018), based on the work of Ursino et al. (2017): https://github.com/artemis-toumazi/dfpk

---

## 🛠️ Environment and Dependencies

The code was developed and tested in the following environment:

- **R version**: 4.2.3 (2023-06-18)  
- **Platform**: aarch64-apple-darwin20 (64-bit)

### Required R packages:

```r
assertthat (>= 0.2.1)
loo (>= 2.8.0)
functional (>= 0.6)
dplyr (>= 1.1.4)
saemix (>= 3.3)
npde (>= 3.5)
rstan (>= 2.32.6)
StanHeaders (>= 2.32.10)
deSolve (>= 1.40)
cubature (>= 2.1.1)
PK (>= 1.3-6)
patchwork (>= 1.3.0)
ggridges (>= 0.5.6)
gridExtra (>= 2.3)
ggplot2 (>= 3.5.1)
```
