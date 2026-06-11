// merlin_oracle.do -- certify the v1-feature extensions (multiple timescales,
// offsets, tvctime, identity-scale splines) against merlin's stexcess
// likelihood (the engine beneath the original v1 command).
//
// Strategy: fit each config with the Mata core, then hand-build the merlin
// model v1 would have constructed -- explicit knots everywhere so knot-siting
// conventions cannot differ, noorthog so parameters are directly comparable --
// and evaluate merlin's log likelihood AT the Mata-core optimum (iterate(0),
// from(b)). Agreement certifies the design construction and quadrature.
// Requires the merlin package; runs WITH the repo ado/ on the adopath.

clear all
set seed 9182

// ---- simulate: age covariate, entry age for attained-age timescales ----
local N = 2000
set obs `=2*`N''
gen byte excess = (_n > `N')
gen double age    = rnormal()
gen double ageoff = 40 + 20*runiform()       // age at time 0 (offset)
gen double mvar   = 2*runiform()             // arbitrary moffset variable

gen double _E = -ln(runiform())
gen double _lr = 0.30*age + 0.01*ageoff
gen double _le = 0.45*age
gen double _lo = -30
gen double _hi =  30
forvalues it = 1/70 {
    gen double _mid = 0.5*(_lo+_hi)
    gen double _t0v = exp(_mid)
    gen double _H = 0.05*exp(_lr)*_t0v^1.20 + excess*0.15*exp(_le)*_t0v^0.80
    replace _hi = _mid if _H > _E
    replace _lo = _mid if _H <= _E
    drop _mid _t0v _H
}
gen double survtime = exp(0.5*(_lo+_hi))
gen byte died = survtime <= 5
replace survtime = min(survtime, 5)
drop _E _lr _le _lo _hi

stset survtime, failure(died)
global ind excess                            // v1 sets this for the evaluator

capture which merlin
if _rc {
    di as err "merlin not installed; skipping merlin_oracle.do"
    exit
}

tempname BV3
local checks 0

// ======================================================================== //
// A: additional log-scale timescale with offset() on the reference
// ======================================================================== //
stexcess (age, df(2) noorthog time2(df(2) offset(ageoff) noorthog)) ///
         (age, df(2) noorthog), indicator(excess)
local llA = e(ll)
local kA_r   "`e(knotsref)'"
local kA_rt2 "`e(knotsref_t2)'"
local kA_e   "`e(knotsexc)'"
matrix `BV3' = e(b)

qui merlin (_t age rcs(_t, knots(`kA_rt2') offset(ageoff) log)        ///
                   rcs(_t, knots(`kA_r') log)                         ///
                , family(user, llf(merlin_stexcess_logl) failure(_d)) ///
                  timevar(_t))                                        ///
           (    age rcs(_t, knots(`kA_e') log)                        ///
                , family(null, reffailure(1)) timevar(_t))            ///
           , indicator(excess) chintpoints(30) nogen                  ///
             from(`BV3') iterate(0) evaltype(gf1) search(off)
di as txt "[A] time2(offset, log):    v2 ll = " %14.8f `llA' ///
    "   merlin ll = " %14.8f e(ll)
assert reldif(e(ll), `llA') < 1e-8
local ++checks

// full merlin maximisation from a perturbed start reaches the same optimum
matrix B0 = `BV3'
forvalues j = 1/`=colsof(B0)' {
    matrix B0[1, `j'] = B0[1, `j'] * 0.9 - 0.01
}
qui merlin (_t age rcs(_t, knots(`kA_rt2') offset(ageoff) log)        ///
                   rcs(_t, knots(`kA_r') log)                         ///
                , family(user, llf(merlin_stexcess_logl) failure(_d)) ///
                  timevar(_t))                                        ///
           (    age rcs(_t, knots(`kA_e') log)                        ///
                , family(null, reffailure(1)) timevar(_t))            ///
           , indicator(excess) chintpoints(30) nogen                  ///
             from(B0) evaltype(gf1) search(off)
di as txt "[A] merlin re-maximised:   ll = " %14.8f e(ll)
assert reldif(e(ll), `llA') < 1e-8
matrix BM = e(b)
forvalues j = 1/`=colsof(BM)' {
    assert reldif(BM[1, `j'], `BV3'[1, `j']) < 2e-4
}
local ++checks

// prediction semantics: merlin's stexcess hazard/chazard (per-row indicator).
// merlin wires merlin_p_stexcess_h/_ch itself when e(cmd2) == "stexcess",
// which needs an e-class helper (the v1 wrapper sets it the same way).
capture program drop _stx_set_cmd2
program _stx_set_cmd2, eclass
    ereturn local cmd2 stexcess
    ereturn local hfunction1  merlin_p_stexcess_h
    ereturn local chfunction1 merlin_p_stexcess_ch
end
_stx_set_cmd2
tempvar mh mch
qui predict double `mh',  hazard  outcome(1)
qui predict double `mch', chazard outcome(1) chintpoints(50)  // = Mata core
// merlin re-maximised: move the Mata core to merlin's optimum for comparison
stexcess (age, knots(`kA_r') noorthog                                 ///
               time2(knots(`kA_rt2') offset(ageoff) noorthog))        ///
         (age, knots(`kA_e') noorthog), indicator(excess) from(BM)
predict vh,  hazard
predict vch, chazard
assert reldif(vh,  `mh')  < 1e-6 if !missing(vh)
assert reldif(vch, `mch') < 1e-6 if !missing(vch)
local ++checks

// ======================================================================== //
// B: tvctime (identity-scale tvc) + additional timescale with its own tvc
// ======================================================================== //
stexcess (age, df(2) noorthog)                                        ///
         (age, df(2) noorthog tvc(age) dftvc(2) tvctime               ///
               time2(df(1) offset(ageoff) tvc(age) dftvc(1) noorthog)) ///
         , indicator(excess)
local llB = e(ll)
local kB_r   "`e(knotsref)'"
local kB_e   "`e(knotsexc)'"
local kB_et2 "`e(knotsexc_t2)'"
matrix `BV3' = e(b)

// the Mata core's internally sited tvc knots (patient events; identity scale
// for the tvctime spline, log(t + ageoff) for the time2 tvc)
mata {
    t   = st_data(., "_t");  d = st_data(., "_d")
    ind = st_data(., "excess");  off = st_data(., "ageoff")
    ev  = select(t,       (d :== 1) :& (ind :== 1))
    evo = select(t + off, (d :== 1) :& (ind :== 1))
    st_local("kB_tvc1", invtokens(strtrim(strofreal(
        _stx_knots(ev, 2)', "%21.0g"))))
    st_local("kB_t2tvc", invtokens(strtrim(strofreal(
        _stx_knots(ln(evo), 1)', "%21.0g"))))
}

qui merlin (_t age rcs(_t, knots(`kB_r') log)                         ///
                , family(user, llf(merlin_stexcess_logl) failure(_d)) ///
                  timevar(_t))                                        ///
           (    age age#rcs(_t, knots(`kB_tvc1'))                     ///
                    rcs(_t, knots(`kB_et2') offset(ageoff) log)       ///
                    age#rcs(_t, knots(`kB_t2tvc') offset(ageoff) log) ///
                    rcs(_t, knots(`kB_e') log)                        ///
                , family(null, reffailure(1)) timevar(_t))            ///
           , indicator(excess) chintpoints(30) nogen                  ///
             from(`BV3') iterate(0) evaltype(gf1) search(off)
di as txt "[B] tvctime + time2 tvc:   v2 ll = " %14.8f `llB' ///
    "   merlin ll = " %14.8f e(ll)
assert reldif(e(ll), `llB') < 1e-8
local ++checks

// ======================================================================== //
// C: moffset() + identity-scale (time) additional timescale
// ======================================================================== //
stexcess (age, df(2) noorthog time2(df(1) moffset(mvar) time noorthog)) ///
         (age, df(2) noorthog), indicator(excess)
local llC = e(ll)
local kC_r   "`e(knotsref)'"
local kC_rt2 "`e(knotsref_t2)'"
local kC_e   "`e(knotsexc)'"
matrix `BV3' = e(b)

qui merlin (_t age rcs(_t, knots(`kC_rt2') moffset(mvar))             ///
                   rcs(_t, knots(`kC_r') log)                         ///
                , family(user, llf(merlin_stexcess_logl) failure(_d)) ///
                  timevar(_t))                                        ///
           (    age rcs(_t, knots(`kC_e') log)                        ///
                , family(null, reffailure(1)) timevar(_t))            ///
           , indicator(excess) chintpoints(30) nogen                  ///
             from(`BV3') iterate(0) evaltype(gf1) search(off)
di as txt "[C] time2(moffset, time):  v2 ll = " %14.8f `llC' ///
    "   merlin ll = " %14.8f e(ll)
assert reldif(e(ll), `llC') < 1e-8
local ++checks

// ======================================================================== //
// D: offset() on the main (baseline) timescale
// ======================================================================== //
stexcess (age, df(2) offset(ageoff) noorthog)                         ///
         (age, df(2) noorthog), indicator(excess)
local llD = e(ll)
local kD_r "`e(knotsref)'"
local kD_e "`e(knotsexc)'"
matrix `BV3' = e(b)

qui merlin (_t age rcs(_t, knots(`kD_r') offset(ageoff) log)          ///
                , family(user, llf(merlin_stexcess_logl) failure(_d)) ///
                  timevar(_t))                                        ///
           (    age rcs(_t, knots(`kD_e') log)                        ///
                , family(null, reffailure(1)) timevar(_t))            ///
           , indicator(excess) chintpoints(30) nogen                  ///
             from(`BV3') iterate(0) evaltype(gf1) search(off)
di as txt "[D] baseline offset():     v2 ll = " %14.8f `llD' ///
    "   merlin ll = " %14.8f e(ll)
assert reldif(e(ll), `llD') < 1e-8
local ++checks

// ======================================================================== //
// E: noconstant -- internal consistency (no merlin analogue needed)
// ======================================================================== //
stexcess (age, df(2) noorthog)(age, df(2) noorthog noconstant), ///
    indicator(excess)
local names : colnames e(b)
assert !`: list posof "exc:_cons" in names'
predict en, netsurvival ci
assert !missing(en[1])
local ++checks

// ---- standardisation with an offset timescale: per-individual fallback
// path must equal the mean of row-wise predictions ----
stexcess (age, df(2) noorthog time2(df(2) offset(ageoff) noorthog)) ///
         (age, df(2) noorthog), indicator(excess)
tempvar tfix
gen double `tfix' = 2.5
predict s_row, survival timevar(`tfix') at(excess 1)
qui su s_row if e(sample)
local smean = r(mean)
predict s_std, survival standardise timevar(`tfix') at(excess 1) ci
assert reldif(s_std[1], `smean') < 1e-8
local ++checks

di as txt _n "merlin_oracle.do completed: `checks' checks passed."
