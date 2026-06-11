// errors.do -- exercise stexcess error paths and validation messages.
// Every check asserts the documented return code, fitting in seconds on tiny
// data. Run via testing/run_live.do (adopath set).

clear all
set seed 2468

set obs 4000
gen byte excess = (_n > 2000)
gen double age  = rnormal()
gen double age2 = 2*age                    // collinear with age
gen byte one    = 1                        // constant covariate
gen double w    = ceil(3*runiform())
gen double t    = 0.1 + 4.9*runiform()
gen byte d      = runiform() < 0.6
gen double bigneg = -10                    // offset that breaks log(t + off)

// ---- not stset ----
rcof "stexcess (age, df(2))(age, df(2)), indicator(excess)" == 119

stset t, failure(d)

// ---- model-spec validation ----
rcof "stexcess (age)(age age)(age), indicator(excess)" == 198  // 3 specs
rcof "stexcess (age, df(0))(age, df(2)), indicator(excess)" == 125
rcof "stexcess (age, df(2) knots(0 1))(age, df(2)), indicator(excess)" == 198
rcof "stexcess (age, df(2) tvc(age))(age, df(2)), indicator(excess)" == 198
rcof "stexcess (age, df(2) tvc(age) dftvc(1 2))(age, df(2)), indicator(excess)" == 198
rcof "stexcess (age, df(2) knots(-1 -1 1))(age, df(2)), indicator(excess)" == 124
rcof "stexcess (age, df(2) time2(offset(age)))(age, df(2)), indicator(excess)" == 198
rcof "stexcess (age, df(2) time2(df(1) knots(0 1) offset(age)))(age, df(2)), indicator(excess)" == 198
rcof "stexcess (age, df(2))(age, df(2)), indicator(excess) chintpoints(0)" == 198

// log of a non-positive timescale: t + offset <= 0 somewhere
rcof "stexcess (age, df(2) time2(df(1) offset(bigneg)))(age, df(2)), indicator(excess)" == 459

// ---- sample validation ----
gen byte ind3 = excess
replace ind3 = 2 in 1
rcof "stexcess (age, df(2))(age, df(2)), indicator(ind3)" == 450
gen byte noexc = 0
rcof "stexcess (age, df(2))(age, df(2)), indicator(noexc)" == 2000
gen byte allexc = 1
rcof "stexcess (age, df(2))(age, df(2)), indicator(allexc) twostage" == 2000

// ---- collinearity: dropped with a note; constant tvc errors ----
stexcess (age age2, df(2))(age, df(2)), indicator(excess)
assert `: word count `e(refvars)'' == 1              // one of age/age2 dropped
assert e(N_ref) == 2000 & e(N_exc) == 2000
stexcess (age one, df(2))(age, df(2)), indicator(excess)  // constant covariate
local rv `e(refvars)'
assert !`: list posof "one" in rv'
rcof "stexcess (age, df(2) tvc(one) dftvc(1))(age, df(2)), indicator(excess)" == 459

// ---- predict validation ----
stexcess (age, df(2))(age, df(2)), indicator(excess)
rcof "predict z1, hazard survival" == 198            // two statistics
rcof "predict z2, survival at(bmi 25)" == 198        // unknown covariate
rcof "predict z3, survival at(age)" == 198           // missing value
rcof "predict z4, survival at(age x)" == 198         // non-numeric value
rcof "predict z5, survival at1(age 1) at2(age 0)" == 198  // needs ?diff/?ratio
rcof "predict z6, sdifference at1(age 1) at2(age 0) standardise" == 198
rcof "predict z7, sdifference at(age 1) at1(age 1) at2(age 0)" == 198
rcof "predict z8, logchazard standardise" == 198
gen double z9_lci = .
rcof "predict z9, survival ci" == 110                // CI name collision

// ---- stale / missing Mata model guards ----
stexcess (age, df(2))(age, df(2)), indicator(excess)
estimates store first
stexcess (age, df(3))(age, df(3)), indicator(excess)
estimates restore first
rcof "predict z10, hazard" == 301                   // store holds the df(3) fit
estimates drop first

stexcess (age, df(2))(age, df(2)), indicator(excess)
clear mata                                            // wipes the model store
mata: mata mlib index
rcof "predict z11, hazard" == 301
stexcess (age, df(2))(age, df(2)), indicator(excess)  // refit -> predict works
predict z12, hazard ci
assert !missing(z12[1])

di as txt _n "errors.do completed."
