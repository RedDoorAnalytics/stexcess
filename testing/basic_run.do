// basic_run.do -- end-to-end smoke test of stexcess v2 (pure Mata core).
// Needs only Stata 19.5+ with the repo's ado/ on the adopath:
//   adopath ++ "`c(pwd)'/ado"
//   mata: mata mlib index     // refresh library search after adopath changes
// then:  do testing/basic_run.do
//
// Syntax follows v1 (merlin-based stexcess): per-component options inside
// each (...) group; prediction statistics and at()/zeros semantics as in v1
// (covariates and the excess indicator default to their observed values).

clear all
set seed 12345

// ---- simulate an additive ref + exc model (controls vs patients) ----
local N = 20000
set obs `=2*`N''
gen byte excess = (_n > `N')              // 1 = patient/excess record
gen double age  = rnormal()

// crude additive-hazard simulation (Weibull ref + Weibull exc)
gen double _E = -ln(runiform())
gen double _lr = 0.30*age
gen double _le = 0.45*age
// invert H_ref + excess*H_exc = E by a short bisection on log t
gen double _lo = -30
gen double _hi =  30
forvalues it = 1/70 {
    gen double _mid = 0.5*(_lo+_hi)
    gen double _t0v = exp(_mid)
    gen double _H = 0.30*exp(_lr)*_t0v^1.20 + excess*0.15*exp(_le)*_t0v^0.80
    replace _hi = _mid if _H > _E
    replace _lo = _mid if _H <= _E
    drop _mid _t0v _H
}
gen double survtime = exp(0.5*(_lo+_hi))
gen byte died = survtime <= 5
replace survtime = min(survtime, 5)
drop _E _lr _le _lo _hi

stset survtime, failure(died)

// ---- fit ----
stexcess (age, df(3))(age, df(3)), indicator(excess)

// truth-ish: age effect on excess ~ 0.45, on reference ~ 0.30
matrix list e(b)

// ---- covariate-pattern predictions with CIs (v1 semantics) ----
predict h_tot,  hazard           ci          // observed age + indicator
predict s_net,  netsurvival      ci at(age 0)
predict s_all,  survival         ci at(age 1 excess 1)

list survtime h_tot s_net s_all in 1/5

// hazard = h_ref + ind*h_exc: at(excess 1) - at(excess 0) == excesshazard
predict h_t1, hazard at(excess 1)
predict h_t0, hazard at(excess 0)
predict h_e,  excesshazard
assert reldif(h_t1 - h_t0, h_e) < 1e-10 | missing(h_e)

// zeros: all covariates (and the indicator) to 0 unless overridden by at()
predict h_z, hazard zeros at(excess 1)
predict h_a, hazard at(age 0 excess 1)
assert reldif(h_z, h_a) < 1e-12 | missing(h_a)

// ---- statistic identities ----
range tt 0 5 50
predict s1,  survival timevar(tt) ci at(age 1 excess 1)
predict c1,  cif      timevar(tt) ci at(age 1 excess 1)
assert reldif(c1, 1 - s1) < 1e-10 if !missing(c1)
assert reldif(c1_lci, 1 - s1_uci) < 1e-10 if !missing(c1_lci)
assert reldif(c1_uci, 1 - s1_lci) < 1e-10 if !missing(c1_uci)

predict ch1, chazard    timevar(tt) at(age 1 excess 1)
predict lc1, logchazard timevar(tt) at(age 1 excess 1)
assert reldif(lc1, ln(ch1)) < 1e-10 if !missing(lc1) & tt > 0

predict rm1, rmst     timevar(tt) ci at(age 0 excess 1)
predict tl1, timelost timevar(tt)    at(age 0 excess 1)
assert reldif(tl1, tt - rm1) < 1e-10 if !missing(tl1)

// ---- standardised (g-formula) predictions over the estimation sample ----
predict ms, survival standardise timevar(tt) ci at(excess 1)
list tt ms ms_lci ms_uci if tt < . in 1/10

// standardised value == mean of the per-row predictions at the same time
qui su tt if _n == 25
local t25 = r(mean)
tempvar tfix
gen double `tfix' = `t25'
predict s_row, survival timevar(`tfix') at(excess 1)
qui su s_row if e(sample)
assert reldif(ms[25], r(mean)) < 1e-8

// ---- RMST (overall and net) at a covariate pattern ----
predict rm, rmst    timevar(tt) ci at(age 0 excess 1)
predict rn, rmstnet timevar(tt) ci at(age 0)
list tt rm rn if inlist(tt, 1, 3, 5)

// ---- contrasts: v1 difference/ratio statistics ----
predict sd, sdifference at1(age 1 excess 1) at2(age 0 excess 1) ///
    timevar(tt) ci
list tt sd sd_lci sd_uci if tt < . in 1/10
predict sr, sratio at1(age 1 excess 1) at2(age 0 excess 1) timevar(tt) ci
predict s0, survival timevar(tt) at(age 0 excess 1)
assert reldif(sd, s1 - s0) < 1e-10 if !missing(sd)
assert reldif(sr, s1 / s0) < 1e-10 if !missing(sr) & tt > 0

// total-vs-reference survival ratio via the indicator (v1 cert idiom)
predict relsr, sratio at1(age 0 excess 1) at2(age 0 excess 0) timevar(tt) ci

// ---- a time-varying excess effect of age (per-component options) ----
stexcess (age, df(3))(age, df(3) tvc(age) dftvc(2)), indicator(excess)
predict hr_tvc, excesshazard timevar(tt) ci at(age 1)

// ---- two-stage estimation (reference fitted to controls only) ----
stexcess (age, df(3))(age, df(3)), indicator(excess) twostage
assert "`e(method)'" == "twostage"
predict s_net2, netsurvival ci at(age 0)
list survtime s_net2 s_net2_lci s_net2_uci in 1/5

// ---- explicit knot placement (log-time scale; how v1 fits are reproduced) --
stexcess (age, df(3))(age, df(3)), indicator(excess)
local kref "`e(knotsref)'"
local kexc "`e(knotsexc)'"
local ll_df3 = e(ll)
di as txt "knots used (log t): ref = `kref'" _n _col(21) "exc = `kexc'"
// refitting with the stored knots must reproduce the df-based fit exactly
stexcess (age, knots(`kref'))(age, knots(`kexc')), indicator(excess)
assert e(dfref) == 3 & e(dfexc) == 3
assert reldif(e(ll), `ll_df3') < 1e-8
// numlist parsing may reformat digits, so compare knot values numerically
forvalues i = 1/4 {
    assert reldif(real(word("`e(knotsref)'", `i')), real(word("`kref'", `i'))) < 1e-10
    assert reldif(real(word("`e(knotsexc)'", `i')), real(word("`kexc'", `i'))) < 1e-10
}

// ---- from(): restarting at the optimum reproduces the fit ----
stexcess (age, df(3))(age, df(3)), indicator(excess)
matrix b0 = e(b)
stexcess (age, df(3))(age, df(3)), indicator(excess) from(b0)
assert reldif(e(ll), `ll_df3') < 1e-8

// ---- delayed entry (left truncation) ----
preserve
gen double entry = runiform()
keep if survtime > entry
stset survtime, failure(died) enter(entry)
stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
assert e(converged) == 1
predict sLT, survival ci at(age 0 excess 1)
assert !missing(sLT[1])
restore
stset survtime, failure(died)

// ---- eform display ----
stexcess (age, df(3))(age, df(3)), indicator(excess) eform
stexcess, eform              // replay

di as txt "basic_run.do completed."
