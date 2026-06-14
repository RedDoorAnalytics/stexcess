// multiple_timescales.do -- smoke test of an additional (attained-age)
// timescale: reference hazard modelled on time since entry AND attained age
// (time2 with offset), v1 syntax. Self-contained; needs only the package.

clear all
set seed 725

local N = 2500
set obs `=2*`N''
gen byte cancer = (_n > `N')
gen double age    = rnormal(0, 5)
gen double agecat = 60 + 10*runiform()        // age at time 0 (offset)

// additive hazards with an attained-age effect in the reference component
gen double _E = -ln(runiform())
gen double _lr = 0.01*age + 0.02*(agecat - 65)
gen double _le = 0.02*age
gen double _lo = -30
gen double _hi =  30
forvalues it = 1/70 {
    gen double _mid = 0.5*(_lo+_hi)
    gen double _t0v = exp(_mid)
    gen double _H = 0.05*exp(_lr)*_t0v^1.10 + cancer*0.04*exp(_le)*_t0v^1.30
    replace _hi = _mid if _H > _E
    replace _lo = _mid if _H <= _E
    drop _mid _t0v _H
}
gen double stime = exp(0.5*(_lo+_hi))
gen byte died = stime <= 10
replace stime = min(stime, 10)
drop _E _lr _le _lo _hi

stset stime, failure(died)

// reference: baseline in log time + attained-age timescale; excess: age
stexcess (age, df(2) time2(df(2) offset(agecat)))(age, df(2)), ///
    indicator(cancer)
assert e(converged) == 1
assert "`e(knotsref_t2)'" != ""

predict s1, survival ci
assert !missing(s1[1])

// total-vs-reference survival ratio at a covariate pattern (v1 cert idiom)
predict s3, sratio at1(cancer 1 age 0) at2(cancer 0 age 0) ci
assert !missing(s3[1])
assert s3[1] > 0 & s3[1] <= 1.0001

// refit with the stored knots reproduces the fit
local ll = e(ll)
stexcess (age, knots(`e(knotsref)') ///
               time2(knots(`e(knotsref_t2)') offset(agecat))) ///
         (age, knots(`e(knotsexc)')), indicator(cancer) nolog
assert reldif(e(ll), `ll') < 1e-8

// ---- three timescales (baseline + time2 + time3) parse, fit and predict ----
// guards time3()-time5(), which the syntax accepts but only time2() exercised
gen double agecat2 = 50 + 8*runiform()        // a second offset timescale
stexcess (age, df(2) time2(df(2) offset(agecat)) time3(df(2) offset(agecat2))) ///
         (age, df(2)), indicator(cancer) nolog
assert e(converged) == 1
assert "`e(knotsref_t2)'" != "" & "`e(knotsref_t3)'" != ""
predict s4, survival ci at(age 0 cancer 1)
assert !missing(s4[1]) & s4[1] > 0 & s4[1] <= 1.0001

// ---- baseline offset()/moffset() (shift _t before the baseline transform) --
// identity-scale baseline so the shift cannot drive the argument non-positive;
// exercises the baseline offset/moffset path outside merlin_oracle
gen double shift0 = 0.25
stexcess (age, df(2) time offset(agecat) moffset(shift0))(age, df(2)), ///
    indicator(cancer) nolog
assert e(converged) == 1
predict s5, survival ci at(age 0 cancer 1)
assert !missing(s5[1]) & s5[1] > 0 & s5[1] <= 1.0001

di as txt _n "multiple_timescales.do completed."
