// factor_vars.do -- certify factor-variable support against hand-coded
// dummy-variable fits: identical likelihoods, coefficients, variances and
// predictions, with at()/zeros recomputing indicators from the underlying
// variables. Self-contained; needs only the package on the adopath.

clear all
set seed 31415

local N = 4000
set obs `=2*`N''
gen byte excess = (_n > `N')
gen double age  = rnormal()
gen byte group  = 1 + floor(3*runiform())     // 1/2/3
gen byte female = runiform() > 0.5

gen double _E = -ln(runiform())
gen double _lr = 0.20*age + 0.30*(group==2) + 0.50*(group==3)
gen double _le = 0.45*age - 0.25*female + 0.10*age*female
gen double _lo = -30
gen double _hi =  30
forvalues it = 1/70 {
    gen double _mid = 0.5*(_lo+_hi)
    gen double _t0v = exp(_mid)
    gen double _H = 0.20*exp(_lr)*_t0v^1.20 + excess*0.10*exp(_le)*_t0v^0.80
    replace _hi = _mid if _H > _E
    replace _lo = _mid if _H <= _E
    drop _mid _t0v _H
}
gen double survtime = exp(0.5*(_lo+_hi))
gen byte died = survtime <= 5
replace survtime = min(survtime, 5)
drop _E _lr _le _lo _hi

// hand-coded dummies / interactions for the oracle fits
gen byte g2 = group == 2
gen byte g3 = group == 3
gen double agefem = age*female

stset survtime, failure(died)

// ======================================================================== //
// A: i.group == hand-coded dummies (coefficients, ll, V, predictions)
// ======================================================================== //
stexcess (age i.group, df(2))(age, df(2)), indicator(excess) nolog
local ll_fv = e(ll)
matrix BF = e(b)
matrix VF = e(V)
assert e(df_m) == 4                  // ref: age + 2 non-base levels; exc: age
// base level present in e(b) with a zero coefficient
local names : colnames e(b)
assert strpos("`names'", "1b.group")
assert _b[ref:1b.group] == 0
assert reldif(_b[ref:2.group], _b[ref:2.group]) == 0   // accessible by name

stexcess (age g2 g3, df(2))(age, df(2)), indicator(excess) nolog
assert reldif(e(ll), `ll_fv') < 1e-8
assert reldif(_b[ref:g2], BF[1, colnumb(BF, "ref:2.group")]) < 1e-6
assert reldif(_b[ref:g3], BF[1, colnumb(BF, "ref:3.group")]) < 1e-6
assert reldif(_se[ref:g2], sqrt(VF[colnumb(BF, "ref:2.group"), ///
    colnumb(BF, "ref:2.group")])) < 1e-6

// predictions at observed values agree
predict h_d, hazard ci
qui stexcess (age i.group, df(2))(age, df(2)), indicator(excess) nolog
predict h_f, hazard ci
assert reldif(h_f, h_d) < 1e-6 & reldif(h_f_lci, h_d_lci) < 1e-5 ///
    if !missing(h_f)

// at() on the underlying variable recomputes the indicators
range tt 0.5 4.5 20
predict s_f3, survival timevar(tt) ci at(group 3 age 0 excess 1)
qui stexcess (age g2 g3, df(2))(age, df(2)), indicator(excess) nolog
predict s_d3, survival timevar(tt) ci at(g2 0 g3 1 age 0 excess 1)
assert reldif(s_f3, s_d3) < 1e-6 if !missing(s_f3)
assert reldif(s_f3_lci, s_d3_lci) < 1e-5 if !missing(s_f3_lci)

// zeros: all underlying covariates to 0 -> base level of group
predict h_dz, hazard zeros at(excess 1)
qui stexcess (age i.group, df(2))(age, df(2)), indicator(excess) nolog
predict h_fz, hazard zeros at(excess 1)
// NB: zeros sets group = 0 (out of range), so all indicators are 0 -- the
// same linear predictor as the dummy model with g2 = g3 = 0
assert reldif(h_fz, h_dz) < 1e-6 if !missing(h_fz)

// standardised survival agrees, including under a counterfactual at()
predict ms_f, survival standardise timevar(tt) ci at(excess 1)
predict ms_f3, survival standardise timevar(tt) ci at(excess 1 group 3)
qui stexcess (age g2 g3, df(2))(age, df(2)), indicator(excess) nolog
predict ms_d, survival standardise timevar(tt) ci at(excess 1)
predict ms_d3, survival standardise timevar(tt) ci at(excess 1 g2 0 g3 1)
assert reldif(ms_f, ms_d) < 1e-6 if !missing(ms_f)
assert reldif(ms_f_lci, ms_d_lci) < 1e-5 if !missing(ms_f_lci)
assert reldif(ms_f3, ms_d3) < 1e-6 if !missing(ms_f3)

// contrasts across factor levels
qui stexcess (age i.group, df(2))(age, df(2)), indicator(excess) nolog
predict sd_f, sdifference at1(group 3 excess 1) at2(group 1 excess 1) ///
    timevar(tt) ci
qui stexcess (age g2 g3, df(2))(age, df(2)), indicator(excess) nolog
predict sd_d, sdifference at1(g2 0 g3 1 excess 1) at2(g2 0 g3 0 excess 1) ///
    timevar(tt) ci
assert reldif(sd_f, sd_d) < 1e-6 if !missing(sd_f) & abs(sd_d) > 1e-10
assert reldif(sd_f_lci, sd_d_lci) < 1e-5 if !missing(sd_f_lci) ///
    & abs(sd_d_lci) > 1e-10

// ======================================================================== //
// B: interactions and polynomials in the excess equation
// ======================================================================== //
stexcess (age, df(2))(c.age##i.female, df(2)), indicator(excess) nolog
local ll_fv = e(ll)
matrix BF = e(b)
stexcess (age, df(2))(age female agefem, df(2)), indicator(excess) nolog
assert reldif(e(ll), `ll_fv') < 1e-8
assert reldif(_b[exc:agefem], BF[1, colnumb(BF, "exc:1.female#c.age")]) < 1e-5

stexcess (age, df(2))(c.age##c.age, df(2)), indicator(excess) nolog
local ll_fv = e(ll)
gen double agesq = age*age
stexcess (age, df(2))(age agesq, df(2)), indicator(excess) nolog
assert reldif(e(ll), `ll_fv') < 1e-8

// ======================================================================== //
// C: housekeeping -- from() roundtrip, replay, twostage, errors
// ======================================================================== //
stexcess (age i.group, df(2))(c.age##i.female, df(2)), indicator(excess) nolog
local ll = e(ll)
matrix b0 = e(b)                               // full layout incl. base terms
stexcess (age i.group, df(2))(c.age##i.female, df(2)), indicator(excess) ///
    from(b0) nolog
assert reldif(e(ll), `ll') < 1e-8
stexcess                                       // replay
stexcess (age i.group, df(2))(age, df(2)), indicator(excess) twostage nolog
predict s2s, netsurvival ci at(age 0)
assert !missing(s2s[1])

// factor variables are not allowed in tvc()
rcof "stexcess (age, df(2))(age, df(2) tvc(i.female) dftvc(1)), indicator(excess)" == 101
// at() must name an underlying variable, not an expanded term
rcof "predict zz, survival at(2.group 1)" == 198

di as txt _n "factor_vars.do completed."
