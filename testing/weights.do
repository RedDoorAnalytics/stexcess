// weights.do -- certify stset weights. Frequency weights are validated
// against the expanded dataset (exact equality: knots, ll, b, V, e(N),
// predictions, standardisation); iweights/pweights share point estimates;
// pweights get a robust sandwich and reduce to the unweighted fit when the
// weights are constant. Conventions follow streg (confirmed empirically):
// e(N) = sum of weights under fw, physical count otherwise; record counts
// weighted everywhere; vcetype Robust under pw.

clear all
set seed 5150

local N = 2500
set obs `=2*`N''
gen byte excess = (_n > `N')
gen double age  = rnormal()
gen byte fw = 1 + floor(3*runiform())
gen double pw = 0.5 + 2*runiform()
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
tempfile base
save `base'

// ======================================================================== //
// A: fweights == expanded data, exactly
// ======================================================================== //
stset survtime [fw=fw], failure(died)
qui su fw
local sumw = r(sum)
stexcess (age, df(2))(age, df(2)), indicator(excess) nolog
assert e(N) == `sumw'
qui su fw if died
assert reldif(e(N_fail), r(sum)) < 1e-12
local ll_w = e(ll)
local kr_w "`e(knotsref)'"
matrix BW = e(b)
matrix VW = e(V)
range tt 0.5 4.5 15
predict snw, netsurvival at(age 1 excess 1) timevar(tt) ci
predict msw, survival standardise timevar(tt) ci at(excess 1)
mkmat snw snw_lci msw msw_lci in 1/15, matrix(PW1)

preserve
expand fw
stset survtime, failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) nolog
assert e(N) == `sumw'
assert reldif(e(ll), `ll_w') < 1e-10
assert "`e(knotsref)'" == "`kr_w'"                  // weighted siting exact
matrix BE = e(b)
matrix VE = e(V)
assert mreldif(BW, BE) < 1e-8 & mreldif(VW, VE) < 1e-6
range tt2 0.5 4.5 15
predict sne, netsurvival at(age 1 excess 1) timevar(tt2) ci
predict mse, survival standardise timevar(tt2) ci at(excess 1)
mkmat sne sne_lci mse mse_lci in 1/15, matrix(PE1)
assert mreldif(PW1, PE1) < 1e-8
restore

// twostage with fweights == twostage expanded (sandwich included)
stexcess (age, df(2))(age, df(2)), indicator(excess) twostage nolog
local ll2_w = e(ll)
matrix B2W = e(b)
matrix V2W = e(V)
preserve
expand fw
stset survtime, failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) twostage nolog
assert reldif(e(ll), `ll2_w') < 1e-10
assert mreldif(e(b), B2W) < 1e-8
matrix V2E = e(V)
assert mreldif(V2W, V2E) < 1e-6
restore

// ======================================================================== //
// B: iweights and pweights -- shared (pseudo-)likelihood
// ======================================================================== //
use `base', clear
stset survtime [iw=pw], failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) nolog
assert e(N) == 2*`N'                                 // physical count
assert "`e(vce)'" == "oim"
local ll_iw = e(ll)
matrix BI = e(b)
matrix VI = e(V)

stset survtime [pw=pw], failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) nolog
assert e(N) == 2*`N'
assert "`e(vce)'" == "robust" & "`e(vcetype)'" == "Robust"
assert reldif(e(ll), `ll_iw') < 1e-10                // same pseudo-ll
assert mreldif(e(b), BI) < 1e-8                      // same point estimates
matrix VP = e(V)
// robust and oim differ, but not absurdly
local r = sqrt(VP[1,1]/VI[1,1])
assert `r' > 0.3 & `r' < 3
qui su pw if died
assert reldif(e(N_fail), r(sum)) < 1e-10             // weighted failures

// pweights with constant weights: same fit as unweighted, given the same
// knots (default knot siting uses expanded-order-statistic percentiles, so
// it is exactly fw==expand-consistent but not invariant to scaling all
// weights -- hence knots are matched explicitly here) and ll scaled by w
gen double cw = 2
stset survtime [pw=cw], failure(died)
stexcess (age, df(2) noorthog)(age, df(2) noorthog), indicator(excess) nolog
matrix BC = e(b)
local kr "`e(knotsref)'"
local ke "`e(knotsexc)'"
local llc = e(ll)
stset survtime, failure(died)
stexcess (age, knots(`kr') noorthog)(age, knots(`ke') noorthog), indicator(excess) nolog
assert mreldif(e(b), BC) < 1e-8
assert reldif(2*e(ll), `llc') < 1e-10

// twostage under pweights runs and is Robust
stset survtime [pw=pw], failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) twostage nolog
assert "`e(vcetype)'" == "Robust"
predict s2, netsurvival ci at(age 0)
assert !missing(s2[1])

di as txt _n "weights.do completed."
