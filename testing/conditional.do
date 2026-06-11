// conditional.do -- certify ltruncated() (conditional predictions) through
// exact identities with unconditional predictions:
//   S(t|t0) = S(t)/S(t0);  H(t|t0) = H(t)-H(t0);
//   RMST(t|t0) = (RMST(t)-RMST(t0))/S(t0);
//   standardised: mean S(t)/mean S(t0); rmst integrates that ratio.
// Identities that mix quadrature panels -- one integral over (t0,t] vs the
// difference of two over (0,t] -- hold to quadrature accuracy (~1e-6
// relative), so those tolerances are 1e-4; same-panel identities are tight.

clear all
set seed 8128

local N = 2500
set obs `=2*`N''
gen byte excess = (_n > `N')
gen double age = rnormal()
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

range tt 1 4.8 20
gen double t0 = 1                                 // conditioning time
gen double t0g = tt                               // for evaluating S(t0)

// ---- covariate-pattern identities ----
predict sc, survival ltruncated(t0) timevar(tt) ci at(age 1 excess 1)
predict su, survival timevar(tt) at(age 1 excess 1)
predict su0, survival timevar(t0) at(age 1 excess 1)
assert reldif(sc, su/su0) < 1e-4 if !missing(sc)
assert sc[1] == 1                                 // t == t0 limit
assert !missing(sc_lci[2]) & sc_lci[2] < sc[2] & sc[2] < sc_uci[2]

predict nc, netsurvival ltruncated(t0) timevar(tt) at(age 1)
predict nu, netsurvival timevar(tt) at(age 1)
predict nu0, netsurvival timevar(t0) at(age 1)
assert reldif(nc, nu/nu0) < 1e-4 if !missing(nc)

predict cc, chazard ltruncated(t0) timevar(tt) ci at(age 1 excess 1)
predict cu, chazard timevar(tt) at(age 1 excess 1)
predict cu0, chazard timevar(t0) at(age 1 excess 1)
assert reldif(cc, cu - cu0) < 1e-4 if !missing(cc) & tt > 1
predict lcc, logchazard ltruncated(t0) timevar(tt) at(age 1 excess 1)
assert reldif(lcc, ln(cc)) < 1e-10 if !missing(lcc) & tt > 1

predict cifc, cif ltruncated(t0) timevar(tt) ci at(age 1 excess 1)
assert reldif(cifc, 1 - sc) < 1e-10 if !missing(cifc)

predict rc, rmst ltruncated(t0) timevar(tt) ci at(age 1 excess 1)
predict ru, rmst timevar(tt) at(age 1 excess 1)
predict ru0, rmst timevar(t0) at(age 1 excess 1)
assert reldif(rc, (ru - ru0)/su0) < 1e-4 if !missing(rc) & tt > 1
predict tlc, timelost ltruncated(t0) timevar(tt) at(age 1 excess 1)
assert reldif(tlc, (tt - t0) - rc) < 1e-8 if !missing(tlc)

// conditional contrast
predict srk, sratio at1(age 1 excess 1) at2(age 0 excess 1) ///
    ltruncated(t0) timevar(tt) ci
predict sc0, survival ltruncated(t0) timevar(tt) at(age 0 excess 1)
assert reldif(srk, sc/sc0) < 1e-8 if !missing(srk) & tt > 1

// rows with t < t0 are missing
gen double tbad = 0.5
predict sb, survival ltruncated(t0) timevar(tbad) at(age 1 excess 1)
assert missing(sb[1])

// ---- standardised marginal-conditional identities ----
predict msc, survival standardise ltruncated(t0) timevar(tt) ci at(excess 1)
predict msu, survival standardise timevar(tt) at(excess 1)
predict msu0, survival standardise timevar(t0g) at(excess 1)
assert reldif(msc, msu/msu0[1]) < 1e-4 if !missing(msc)
assert msc[1] == 1
predict mcifc, cif standardise ltruncated(t0) timevar(tt) ci at(excess 1)
assert reldif(mcifc, 1 - msc) < 1e-10 if !missing(mcifc)

// standardised conditional rmst vs a dense trapezoid of the conditional
// standardised survival curve
range ttd 1 3 400
predict mscd, survival standardise ltruncated(t0) timevar(ttd) at(excess 1)
qui gen double trap = ///
    (mscd + mscd[_n-1])/2*(ttd - ttd[_n-1]) if _n > 1 & !missing(mscd)
qui su trap if ttd <= 3
local rmst_trap = r(sum)
gen double t3 = 3
predict mrc, rmst standardise ltruncated(t0) timevar(t3) ci at(excess 1)
assert reldif(mrc[1], `rmst_trap') < 2e-4
predict mtlc, timelost standardise ltruncated(t0) timevar(t3) at(excess 1)
assert reldif(mtlc[1], (3 - 1) - mrc[1]) < 1e-8

// ---- errors ----
rcof "predict zz, hazard ltruncated(t0)" == 198
rcof "predict zz, excesshazard ltruncated(t0)" == 198
rcof "predict zz, chazard standardise ltruncated(t0)" == 198

di as txt _n "conditional.do completed."
