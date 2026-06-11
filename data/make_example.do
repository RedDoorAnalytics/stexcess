// make_example.do -- simulate the packaged example dataset
// (data/stexcess_example.dta): a cancer cohort with matched population
// controls, used by the runnable examples in the help files.
//   stata-mp -q -b do data/make_example.do        (run from the repo root)

clear all
set seed 280461

local N = 3000                                // per arm
set obs `=2*`N''
gen byte patient = (_n > `N')
label var patient "1 = cancer patient, 0 = population control"

gen double age = min(max(rnormal(70, 8.5), 45), 90)
label var age "Age at diagnosis (years)"
gen byte female = runiform() < 0.5
label var female "Female"
label define female 0 "Male" 1 "Female"
label values female female

// stage of the matched patient (kept for controls so the variable is
// complete; it enters the excess equation only)
gen byte stage = 1 + (runiform() > 0.5) + (runiform() > 0.8)
label var stage "Stage at diagnosis"
label define stage 1 "Stage I" 2 "Stage II" 3 "Stage III"
label values stage stage

// additive hazards: Weibull reference (age, sex) + Weibull excess
// (age, sex, stage); inverted by bisection on log t
gen double _E = -ln(runiform())
gen double _lr = 0.080*(age - 70) - 0.30*female
gen double _le = 0.025*(age - 70) + 0.10*female + ///
    0.70*(stage == 2) + 1.40*(stage == 3)
gen double _lo = -30
gen double _hi =  30
forvalues it = 1/70 {
    gen double _mid = 0.5*(_lo + _hi)
    gen double _t0v = exp(_mid)
    gen double _H = 0.040*exp(_lr)*_t0v^1.30 + ///
        patient*0.085*exp(_le)*_t0v^0.90
    replace _hi = _mid if _H > _E
    replace _lo = _mid if _H <= _E
    drop _mid _t0v _H
}
gen double stime = exp(0.5*(_lo + _hi))
gen byte died = stime <= 5                     // 5-year administrative censoring
replace stime = min(stime, 5)
label var stime "Follow-up time (years)"
label var died "Died"
drop _E _lr _le _lo _hi

compress
label data "Simulated cancer cohort with population controls (stexcess example)"
save "data/stexcess_example.dta", replace

stset stime, failure(died)
tab patient died, row
