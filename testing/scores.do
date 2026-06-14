// scores.do -- certify predict, scores: one variable per e(b) coefficient,
// the per-record contribution to the (weighted) gradient over the
// estimation sample. The decisive checks are exact: the score cross-product
// reproduces the robust and cluster sandwiches (V_oim B V_oim), base/omitted
// factor-variable columns are exactly 0, and the reference/excess block
// structure matches the estimating equations. A loose colsum check confirms
// the scores sum to ~0 (the gradient; its raw scale is parameter-dependent
// because the spline coefficients differ by orders of magnitude).

clear all
set seed 90210

local N = 3000
set obs `=2*`N''
gen long id  = ceil(_n/3)                      // clusters spanning records
gen byte excess = (_n > `N')
gen double age  = rnormal()
gen byte grp = 1 + floor(3*runiform())
gen double w = 0.5 + 2*runiform()
gen double _E = -ln(runiform())
gen double _lr = 0.30*age + 0.2*(grp==2) + 0.4*(grp==3)
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
gen byte esamp = 0

// ======================================================================== //
// A: joint -- robust and cluster sandwiches reproduced exactly; the excess
// block of a control record's score is exactly 0
// ======================================================================== //
qui stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
matrix Adef = e(V)                             // OIM bread (default e(V))
local kfull = colsof(e(b))
local pr = e(dfref) + 2                         // age + rcs(df) + cons = df+2
predict sc*, scores
unab scv : sc*
assert `: word count `scv'' == `kfull'
qui replace esamp = e(sample)

// control records contribute nothing to the excess block (cols pr+1..kfull)
forvalues j = `=`pr'+1'/`kfull' {
    local v : word `j' of `scv'
    qui su `v' if esamp & excess == 0
    assert r(max) == 0 & r(min) == 0
}
// loose gradient check (scale is parameter-dependent under the spline basis)
foreach v of local scv {
    qui su `v' if esamp
    assert abs(r(mean)*r(N)) < 1e-2
}

qui count if esamp
local n = r(N)
mata {
    S  = st_data(., tokens(st_local("scv")), "esamp")
    A  = st_matrix("Adef")
    st_matrix("Vrep", A * ((`n'/(`n'-1)) * cross(S, S)) * A)
    cl = st_data(., "id", "esamp")
    o  = order(cl, 1)
    SC = panelsum(S[o, .], panelsetup(cl[o], 1))
    Gc = rows(SC)
    st_matrix("Vrepc", A * ((Gc/(Gc-1)) * cross(SC, SC)) * A)
}
qui stexcess (age, df(3))(age, df(3)), indicator(excess) vce(robust) nolog
mata: st_local("d", strofreal(mreldif(st_matrix("e(V)"), st_matrix("Vrep"))))
assert `d' < 1e-7
qui stexcess (age, df(3))(age, df(3)), indicator(excess) vce(cluster id) nolog
mata: st_local("dc", strofreal(mreldif(st_matrix("e(V)"), st_matrix("Vrepc"))))
assert `dc' < 1e-7
di as txt "[A] joint: robust reldif `d', cluster reldif `dc'"

// ======================================================================== //
// B: factor variables -- base/omitted columns exactly 0; robust reproduced
// ======================================================================== //
qui stexcess (age i.grp, df(3))(age, df(3)), indicator(excess) nolog
matrix Adef = e(V)
predict fsc*, scores
qui replace esamp = e(sample)
local cn : colfullnames e(b)
local j = 1
foreach v of varlist fsc* {
    local nm : word `j' of `cn'
    if strpos("`nm'", "b.grp") | strpos("`nm'", "o.") {
        qui su `v' if esamp
        assert r(max) == 0 & r(min) == 0       // base/omitted score == 0
    }
    local ++j
}
unab fscv : fsc*
qui count if esamp
local n = r(N)
mata {
    S = st_data(., tokens(st_local("fscv")), "esamp")
    st_matrix("Vrep", st_matrix("Adef") * ///
        ((`n'/(`n'-1)) * cross(S, S)) * st_matrix("Adef"))
}
qui stexcess (age i.grp, df(3))(age, df(3)), indicator(excess) vce(robust) nolog
mata: st_local("df", strofreal(mreldif(st_matrix("e(V)"), st_matrix("Vrep"))))
assert `df' < 1e-7
di as txt "[B] factor variables: base cols zero, robust reldif `df'"

// ======================================================================== //
// C: two-stage -- control rows score only the reference block, patient rows
// only the excess block (the stacked estimating equations)
// ======================================================================== //
qui stexcess (age, df(3))(age, df(3)), indicator(excess) twostage nolog
local kfull = colsof(e(b))
local pr = e(dfref) + 2
capture drop sc*
predict sc*, scores
unab scv : sc*
qui replace esamp = e(sample)
forvalues j = `=`pr'+1'/`kfull' {              // exc block, control rows == 0
    local v : word `j' of `scv'
    qui su `v' if esamp & excess == 0
    assert r(max) == 0 & r(min) == 0
}
forvalues j = 1/`pr' {                          // ref block, patient rows == 0
    local v : word `j' of `scv'
    qui su `v' if esamp & excess == 1
    assert r(max) == 0 & r(min) == 0
}
di as txt "[C] twostage: block-diagonal score structure confirmed"

// ======================================================================== //
// D: pweights -- weighted scores (w*s) reproduce the robust sandwich. The
// oim bread is the iweight e(V) (iw and pw share the pseudo-likelihood and
// Hessian; only the variance estimator differs), so the pweight sandwich is
// bread * (n/(n-1)) cross(scores) * bread, exactly as for the unweighted case.
// ======================================================================== //
stset survtime [iw=w], failure(died)
qui stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
matrix Aiw = e(V)                              // oim bread at the weighted fit
stset survtime [pw=w], failure(died)
qui stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
capture drop sc*
predict sc*, scores
qui replace esamp = e(sample)
unab scv : sc*
foreach v of local scv {                       // still the (weighted) gradient
    qui su `v' if esamp
    assert abs(r(mean)*r(N)) < 1e-2
}
qui count if esamp
local n = r(N)
mata {
    S = st_data(., tokens(st_local("scv")), "esamp")
    st_matrix("Vrep", st_matrix("Aiw") * ///
        ((`n'/(`n'-1)) * cross(S, S)) * st_matrix("Aiw"))
}
mata: st_local("dpw", strofreal(mreldif(st_matrix("e(V)"), st_matrix("Vrep"))))
assert `dpw' < 1e-7
di as txt "[D] pweights: scores reproduce robust sandwich, reldif `dpw'"

// fweights: the scores remain the per-record gradient contributions (sum to
// the gradient), but their cross-product is the w^2-sum, NOT the fweight
// (expanded-data) sandwich -- so it is the unweighted/pweight/cluster
// sandwich that scores reproduce, not the fweight one (see the help).
gen int fw = 1 + floor(3*runiform())
stset survtime [fw=fw], failure(died)
qui stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
capture drop sc*
predict sc*, scores
qui replace esamp = e(sample)
unab scv : sc*
foreach v of local scv {
    qui su `v' if esamp
    assert abs(r(mean)*r(N)) < 1e-2
}
di as txt "[D] fweights: scores are the gradient contributions"
stset survtime, failure(died)

// ======================================================================== //
// E: errors
// ======================================================================== //
qui stexcess (age, df(3))(age, df(3)), indicator(excess) nolog
rcof "predict s1 s2, scores" == 198            // wrong number of names
rcof "margins, predict(hazard)" == 322         // margins not supported

di as txt _n "scores.do completed."
