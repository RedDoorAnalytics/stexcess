# stexcess

Flexible parametric excess hazard models in Stata, with a modelled expected rate.

`stexcess` fits an excess hazard model in which the expected (reference) hazard is not read from a life table but modelled, jointly with the excess hazard, from a reference population in the same data. It is built on [`merlin`](https://github.com/RedDoorAnalytics/merlin).

## What it fits

After `stset`, with an indicator variable coded 0 for the reference population and 1 for the population with an excess hazard, `stexcess` takes two models, `(reference_model)(excess_model)`:

- each has its own covariates and its own baseline, a restricted cubic spline of log time (or of time) on the log hazard scale;
- **time-dependent effects** in either model, through restricted cubic splines of log time or time (`tvc()`, `dftvc()`);
- **further timescales** in either model (`time2()` to `time5()`);
- delayed entry, taken from the `stset`.

Each linear predictor accepts `merlin`'s extended syntax, so splines and fractional polynomials of continuous covariates can be written directly in the model.

After estimation, `predict` gives survival, hazard, cumulative hazard and cumulative incidence functions, restricted mean survival time and time lost, and differences and ratios of these between covariate patterns, with confidence intervals.

## Requirements

stexcess 1.1.2 supports Stata 17 or later, merlin 2.5.0 or 2.6.0, and stmerlin 1.1.2. It does not support merlin 3.0.0.

```stata
net install merlin, from("https://reddooranalytics.se/install/stata/merlin/2.5.0/")
net install stmerlin, from("https://reddooranalytics.se/install/stata/stmerlin/1.1.2/")
```

## Installation

```stata
net install stexcess, from("https://raw.githubusercontent.com/RedDoorAnalytics/stexcess/main/")
```

## Example

```stata
// a reference population and a cancer population in one dataset: everyone has the
// reference hazard, and those with cancer==1 an excess hazard as well
clear
set seed 725
set obs 3000
generate byte cancer = runiform()>0.5
generate double age = rnormal(0,10)
generate byte sex = runiform()>0.5
generate double t1 = (-ln(runiform())/(0.01*exp(0.01*age + 0.2*sex)))^(1/1.5)
generate double t2 = (-ln(runiform())/(0.02*exp(0.02*age + 0.3*sex)))^(1/1.2) if cancer
generate double stime = min(t1, t2, 10)
generate byte died = stime<10
stset stime, failure(died)

// a flexible parametric excess hazard model
stexcess (age sex, df(3))(age sex, df(2)), indicator(cancer)

// a time-dependent effect of age in the reference model
stexcess (age sex, df(3) tvc(age) dftvc(1))(age sex, df(2)), indicator(cancer)
```

Further detail is in the help files: `help stexcess` and `help stexcess postestimation`.

## Checked with

The example above, predictions after it, a spline and a fractional polynomial of a covariate, a model with a further timescale and a model with delayed entry ran with merlin 2.5.0 and with merlin 2.6.0, each with stmerlin 1.1.2, on Stata 19.5.

## Version

Version 1.1.2.

**1.1.2.** The time-dependent effects option, `tvc()`, works: in 1.1.1 it stopped with an error in every call. The examples in the help file are corrected.

**1.1.1.** Predictions re-synced.

**1.1.0.** A new evaluator calculating the full log likelihood, with analytic scores.

## Author and licence

Michael J. Crowther (Red Door Analytics, Stockholm) — `michael.crowther@reddooranalytics.se`.

Licence: **GPL-3.0-or-later** — see [`LICENSE`](LICENSE).

Copyright (C) 2022-2026 Red Door Analytics.

> **Licence, 2026-10-09.** stexcess was published without a licence file from its first commit (2022-10-18) until this date. It is **GPL-3** from 9 October 2026, as is the rest of the family.
