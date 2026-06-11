// weights_mc.do -- Monte-Carlo validation of the pweight robust sandwich:
// age-dependent (informative) sampling from a simulated population, with
// weights = 1/Pr(sampled). Checks that the mean reported robust SE matches
// the empirical SD of the estimates and that 95% CI coverage is nominal.
// Slower than the regular suite (~200 weighted fits); run on demand:
//   stata-mp -q -b do testing/weights_mc.do   (with ado/ on the adopath)

clear all
set seed 271828
local R = 200
local truth = 0.45                    // excess log-hazard ratio for age

tempname res
postfile `res' b se cover using "`c(tmpdir)'/stx_wmc.dta", replace

forvalues r = 1/`R' {
    qui {
        clear
        set obs 10000
        gen byte excess = (_n > 5000)
        gen double age = rnormal()
        gen double _E = -ln(runiform())
        gen double _lr = 0.30*age
        gen double _le = `truth'*age
        gen double _lo = -30
        gen double _hi =  30
        forvalues it = 1/70 {
            gen double _mid = 0.5*(_lo+_hi)
            gen double _t0v = exp(_mid)
            gen double _H = 0.30*exp(_lr)*_t0v^1.20 + ///
                excess*0.15*exp(_le)*_t0v^0.80
            replace _hi = _mid if _H > _E
            replace _lo = _mid if _H <= _E
            drop _mid _t0v _H
        }
        gen double survtime = exp(0.5*(_lo+_hi))
        gen byte died = survtime <= 5
        replace survtime = min(survtime, 5)

        // informative sampling: inclusion probability depends on age
        gen double p = 0.15 + 0.70*invlogit(1.2*age)
        keep if runiform() < p
        gen double w = 1/p

        stset survtime [pw=w], failure(died)
        capture stexcess (age, df(2))(age, df(2)), indicator(excess) nolog
        if _rc | e(converged) != 1 continue
        local bb = _b[exc:age]
        local ss = _se[exc:age]
        post `res' (`bb') (`ss') ///
            (abs(`bb' - `truth') < invnormal(0.975)*`ss')
    }
}
postclose `res'

use "`c(tmpdir)'/stx_wmc.dta", clear
qui su b
local bias = r(mean) - `truth'
local sd = r(sd)
qui su se
local sebar = r(mean)
qui su cover
local cov = r(mean)
di as txt _n "pweight sandwich MC (" _N " replications):"
di as txt "  mean estimate = " %8.5f `bias' + `truth' "   (truth `truth')"
di as txt "  empirical SD  = " %8.5f `sd'
di as txt "  mean robust SE = " %7.5f `sebar' ///
    "   ratio SE/SD = " %5.3f `sebar'/`sd'
di as txt "  95% CI coverage = " %5.3f `cov'
assert abs(`bias') < 0.02
assert inrange(`sebar'/`sd', 0.85, 1.18)
assert inrange(`cov', 0.91, 0.98)
di as txt _n "weights_mc.do completed."
