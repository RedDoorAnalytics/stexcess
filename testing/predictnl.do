// predictnl.do -- stexcess predictions must work inside predictnl, i.e. when
// predict is re-invoked with a perturbed e(b) (and, for predictnl's CI check,
// a temporarily reposted e(V)). This certifies two things:
//   (1) predict reads the CURRENTLY active coefficients, so the point estimate
//       from predictnl is identical to a plain predict; and
//   (2) predictnl's numerical delta-method CI matches our analytic Jacobian CI
//       to numerical-derivative precision, compared on the SAME transform scale
//       predict uses (identity for logchazard, log for hazards/rmst, cloglog
//       for survival) -- an independent check of the analytic Jacobians.
// The staleness guard against estimates restore of a different fit is covered
// in errors.do; here we only confirm predictnl no longer trips it.
// Self-contained; needs only the package on the adopath.

clear all
set seed 20260723

local N = 2500
set obs `=2*`N''
gen byte excess = (_n > `N')
gen double age  = rnormal()

gen double _E = -ln(runiform())
gen double _lr = 0.30*age
gen double _le = 0.45*age
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

stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
range tt 0.5 5 10

// predictnl differentiates numerically (step ~1e-2 of |b|), so CI endpoints
// agree with the analytic delta method to ~1e-4 relative, point estimates ~exact
local ptol 1e-7
local ctol 1e-4

//----------------------------------------------------------------------
// logchazard: predict CI is on the identity scale == predictnl's default
//----------------------------------------------------------------------
predict lc, logchazard timevar(tt) at(age 1 excess 1) ci
predictnl double lcp = predict(logchazard timevar(tt) at(age 1 excess 1)), ///
    ci(lcp_lci lcp_uci)
assert reldif(lc,     lcp)     < `ptol' if tt<.
assert reldif(lc_lci, lcp_lci) < `ctol' if tt<.
assert reldif(lc_uci, lcp_uci) < `ctol' if tt<.

//----------------------------------------------------------------------
// hazard: predict CI is on the log scale -> compare to predictnl of log(h)
//----------------------------------------------------------------------
predict h, hazard timevar(tt) at(age 1 excess 1) ci
predictnl double lh = log(predict(hazard timevar(tt) at(age 1 excess 1))), ///
    ci(lh_lci lh_uci)
assert reldif(h,     exp(lh))     < `ptol' if tt<.
assert reldif(h_lci, exp(lh_lci)) < `ctol' if tt<.
assert reldif(h_uci, exp(lh_uci)) < `ctol' if tt<.

//----------------------------------------------------------------------
// survival: predict CI is on the cloglog scale; S = exp(-exp(cll)) is
// DECREASING in cll, so cll's lower CI maps to S's upper CI and vice versa
//----------------------------------------------------------------------
predict s, survival timevar(tt) at(age 1 excess 1) ci
predictnl double cll = log(-log(predict(survival timevar(tt) at(age 1 excess 1)))), ///
    ci(cll_lci cll_uci)
assert reldif(s,     exp(-exp(cll)))     < `ptol' if tt<.
assert reldif(s_lci, exp(-exp(cll_uci))) < `ctol' if tt<.
assert reldif(s_uci, exp(-exp(cll_lci))) < `ctol' if tt<.

//----------------------------------------------------------------------
// rmst: log scale (as for hazards)
//----------------------------------------------------------------------
predict rm, rmst timevar(tt) at(age 0 excess 1) ci
predictnl double lrm = log(predict(rmst timevar(tt) at(age 0 excess 1))), ///
    ci(lrm_lci lrm_uci)
assert reldif(rm,     exp(lrm))     < `ptol' if tt<. & tt>0
assert reldif(rm_lci, exp(lrm_lci)) < `ctol' if tt<. & tt>0
assert reldif(rm_uci, exp(lrm_uci)) < `ctol' if tt<. & tt>0

//----------------------------------------------------------------------
// standardised (g-formula) survival inside predictnl
//----------------------------------------------------------------------
predict ms, survival standardise timevar(tt) at(excess 1) ci
predictnl double mcll = ///
    log(-log(predict(survival standardise timevar(tt) at(excess 1)))), ///
    ci(mcll_lci mcll_uci)
assert reldif(ms,     exp(-exp(mcll)))      < `ptol' if tt<.
assert reldif(ms_lci, exp(-exp(mcll_uci)))  < 2e-4   if tt<.
assert reldif(ms_uci, exp(-exp(mcll_lci)))  < 2e-4   if tt<.

//----------------------------------------------------------------------
// a difference contrast inside predictnl (identity scale both sides)
//----------------------------------------------------------------------
predict sd, sdifference at1(age 1 excess 1) at2(age 0 excess 1) timevar(tt) ci
predictnl double sdp = ///
    predict(sdifference at1(age 1 excess 1) at2(age 0 excess 1) timevar(tt)), ///
    ci(sdp_lci sdp_uci)
assert reldif(sd,     sdp)     < `ptol' if tt<.
assert reldif(sd_lci, sdp_lci) < `ctol' if tt<.
assert reldif(sd_uci, sdp_uci) < `ctol' if tt<.

di as txt _n "predictnl.do completed."
