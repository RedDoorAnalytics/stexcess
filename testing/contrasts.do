// contrasts.do -- the hazard/CIF/RMST difference and ratio prediction
// statistics (hdifference, hratio, cifdifference, cifratio, rmstdifference,
// rmstratio) plus the excess-hazard ratio (excesshratio). basic_run.do covers
// sdifference/sratio; this file certifies the remaining contrasts via the
// defining identity:
//   Xdifference at1(A) at2(B) == [X at(A)] - [X at(B)]
//   Xratio      at1(A) at2(B) == [X at(A)] / [X at(B)]
// computed from the same fit, so agreement is to machine precision.
// Self-contained; needs only the package on the adopath.

clear all
set seed 9001

local N = 3000
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

range tt 0 5 50

// the two covariate patterns contrasted (totals, differing in age)
local A "age 1 excess 1"
local B "age 0 excess 1"

// ---- direct quantities at each pattern ----
predict hA, hazard timevar(tt) at(`A')
predict hB, hazard timevar(tt) at(`B')
predict cA, cif    timevar(tt) at(`A')
predict cB, cif    timevar(tt) at(`B')
predict rA, rmst   timevar(tt) at(`A')
predict rB, rmst   timevar(tt) at(`B')

// ---- hazard difference / ratio ----
predict hd, hdifference at1(`A') at2(`B') timevar(tt) ci
predict hr, hratio      at1(`A') at2(`B') timevar(tt) ci
assert reldif(hd, hA - hB) < 1e-9 if !missing(hd)
assert reldif(hr, hA / hB) < 1e-9 if !missing(hr)
assert !missing(hd_lci[2]) & !missing(hr_lci[2])      // CIs produced

// ---- CIF difference / ratio (ratio is 0/0 at t=0 -> missing) ----
predict cd, cifdifference at1(`A') at2(`B') timevar(tt) ci
predict cr, cifratio      at1(`A') at2(`B') timevar(tt) ci
assert reldif(cd, cA - cB) < 1e-9 if !missing(cd)
assert reldif(cr, cA / cB) < 1e-9 if !missing(cr) & tt > 0
assert missing(cr) if tt == 0                          // 0/0 limit left missing

// ---- RMST difference / ratio ----
predict rd, rmstdifference at1(`A') at2(`B') timevar(tt) ci
predict rr, rmstratio      at1(`A') at2(`B') timevar(tt) ci
assert reldif(rd, rA - rB) < 1e-9 if !missing(rd)
assert reldif(rr, rA / rB) < 1e-9 if !missing(rr) & tt > 0

// ---- excess-hazard ratio / difference (excess-only; ignore the indicator) --
predict eA, excesshazard timevar(tt) at(`A')
predict eB, excesshazard timevar(tt) at(`B')
predict er, excesshratio      at1(`A') at2(`B') timevar(tt) ci
predict ed, excesshdifference at1(`A') at2(`B') timevar(tt) ci
assert reldif(er, eA / eB) < 1e-9 if !missing(er)
assert reldif(ed, eA - eB) < 1e-9 if !missing(ed)
assert !missing(er_lci[2]) & !missing(er_uci[2])      // CIs produced
assert !missing(ed_lci[2]) & !missing(ed_uci[2])

// ---- self-contrast: difference == 0, ratio == 1 ----
predict hd0, hdifference at1(`A') at2(`A') timevar(tt)
predict hr1, hratio      at1(`A') at2(`A') timevar(tt)
predict er1, excesshratio      at1(`A') at2(`A') timevar(tt)
predict ed0, excesshdifference at1(`A') at2(`A') timevar(tt)
assert abs(hd0) < 1e-10 if !missing(hd0)
assert reldif(hr1, 1) < 1e-10 if !missing(hr1)
assert reldif(er1, 1) < 1e-10 if !missing(er1)
assert abs(ed0) < 1e-10 if !missing(ed0)

di as txt _n "contrasts.do completed."
