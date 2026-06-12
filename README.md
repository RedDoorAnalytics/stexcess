# stexcess

Stata command for **modelled excess hazard models** — a relative-survival /
excess-hazard model in which the expected (reference) rate is itself a fitted
hazard model estimated from a control cohort, rather than taken from external
population life tables:

    h_total(t | x) = h_ref(t | x) + h_exc(t | x)

Both components are flexible parametric models (restricted cubic splines on
the log-hazard scale), estimated jointly so that the uncertainty in the
reference model propagates into every excess and net quantity.

Version 2 is a ground-up reimplementation with a pure Stata/Mata core:
designs are built once in Mata and maximised with `optimize()` using an exact
analytic log likelihood, gradient and Hessian; predictions, contrasts and
g-formula standardisation use analytic delta-method Jacobians. There are no
dependencies outside Stata, and estimation is roughly two orders of magnitude
faster than version 1 (~5 seconds for 200,000 records on a laptop).

## Features

- v1-compatible syntax: per-component `df()`/`knots()`, `tvc()`/`dftvc()`
  time-varying effects, multiple timescales (`time2()`–`time5()` with
  `offset()`/`moffset()`, e.g. attained age), `time`/`tvctime`
  identity-scale splines, `offset()` on the baseline, `noconstant`
- factor variables and interactions (`i.stage`, `c.age##i.sex`), with
  `at()` operating on the underlying variables so indicator and interaction
  terms are recomputed automatically
- joint maximum likelihood or `twostage` estimation (reference fitted to
  controls only, with a stacked M-estimation sandwich variance)
- stset weights: fweights (exactly equivalent to expanded data), pweights
  (robust sandwich variance), iweights; delayed entry (left truncation)
- `vce(robust)` and `vce(cluster clustvar)` sandwich variances
- predictions: hazard, cumulative hazard, survival, CIF, RMST, time lost,
  net (excess-only) survival/hazard/RMST, differences and ratios — all with
  analytic delta-method confidence intervals; conditional versions
  (`ltruncated()`, e.g. S(t | t0))
- regression-standardised (g-formula) predictions over the estimation
  sample, including counterfactual `at()` overrides

## Installation

Requires Stata 19.5 or later.

```stata
net install stexcess, from("https://raw.githubusercontent.com/RedDoorAnalytics/stexcess/main/")
```

## Getting started

Using the package's simulated example data (a cancer cohort with matched
population controls; `net get stexcess` downloads a local copy):

```stata
use https://raw.githubusercontent.com/RedDoorAnalytics/stexcess/main/data/stexcess_example.dta, clear
stset stime, failure(died)
stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)

predict h,  hazard ci                       // observed covariates + arm
predict sn, netsurvival ci                  // net survival, observed covariates

range tt 0 5 100
predict ms, survival standardise timevar(tt) ci at(patient 1)
predict rsr, sratio at1(age 70 female 0 patient 1) ///
    at2(age 70 female 0 patient 0) timevar(tt) ci   // relative survival ratio
```

See `help stexcess` and `help stexcess postestimation` for the full syntax.

## Performance

Wall-clock times on an Apple Silicon laptop (StataNow/MP 19.5), 200,000
records, df(3) on both components:

| task | stexcess v2 | v1 (merlin-based) |
|---|---:|---:|
| fit, 1 covariate per equation | ~5 s | ~500 s |
| fit, 20 covariates per equation (48 parameters) | ~18 s | — |
| standardised survival + CI, 50-point grid | ~11 s | — |

Standardised predictions cost is proportional to the number of distinct
evaluation times multiplied by the population size — use a `timevar()`
grid rather than the default (every observation's `_t`).

## Validation

The Mata core is certified against the likelihood of the original
(merlin-based) stexcess v1.1.1 evaluated at the same parameters — exact
agreement across multi-timescale, offset and time-varying-effect
configurations — and reproduces a v1 multiple-timescale fit
coefficient-for-coefficient when given v1's knots (`testing/`).

## Author

Michael J. Crowther, Red Door Analytics AB, Stockholm
(michael.crowther@reddooranalytics.se)
