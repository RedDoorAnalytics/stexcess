// vce.do -- certify vce(robust) and vce(cluster) through exact identities:
// robust == pweights with unit weights; cluster with singleton clusters ==
// robust; duplicate-record clusters == single-copy robust (the cluster
// sandwich collapses the copies); fweights + cluster == the expanded data.

clear all
set seed 60902

local N = 2000
set obs `=2*`N''
gen long id = _n
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

// reference knots for the algebraic identities (raw basis, fixed knots, so
// reparametrisations and knot-siting conventions cannot interfere)
qui stexcess (age, df(2))(age, df(2)), indicator(excess) nolog
local kr "`e(knotsref)'"
local ke "`e(knotsexc)'"
local spec (age, knots(`kr') noorthog)(age, knots(`ke') noorthog)

// ---- vce(robust) == pweights with unit weights (same formula) ----
stexcess `spec', indicator(excess) vce(robust) nolog
assert "`e(vce)'" == "robust" & "`e(vcetype)'" == "Robust"
matrix BR = e(b)
matrix VR = e(V)
gen double one = 1
stset survtime [pw=one], failure(died)
stexcess `spec', indicator(excess) nolog
assert mreldif(e(b), BR) < 1e-10 & mreldif(e(V), VR) < 1e-8

// ---- cluster with singleton clusters == robust ----
stset survtime, failure(died)
stexcess `spec', indicator(excess) vce(cluster id) nolog
assert "`e(vce)'" == "cluster" & "`e(vcetype)'" == "Robust"
assert "`e(clustvar)'" == "id" & e(N_clust) == 2*`N'
assert mreldif(e(b), BR) < 1e-10 & mreldif(e(V), VR) < 1e-8

// ---- duplicate-record clusters == single-copy robust ----
preserve
expand 3
stset survtime, failure(died)
stexcess `spec', indicator(excess) vce(cluster id) nolog
assert e(N_clust) == 2*`N'
assert mreldif(e(b), BR) < 1e-8 & mreldif(e(V), VR) < 1e-6
restore

// ---- fweights + cluster == expanded data + cluster ----
gen byte fw = 1 + floor(3*runiform())
gen long cl = ceil(id/5)                      // clusters spanning records
stset survtime [fw=fw], failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) vce(cluster cl) nolog
matrix BFW = e(b)
matrix VFW = e(V)
local ncl = e(N_clust)
preserve
expand fw
stset survtime, failure(died)
stexcess (age, df(2))(age, df(2)), indicator(excess) vce(cluster cl) nolog
assert e(N_clust) == `ncl'
assert mreldif(e(b), BFW) < 1e-8 & mreldif(e(V), VFW) < 1e-6
restore

// ---- string cluster variables; twostage + cluster ----
stset survtime, failure(died)
gen str12 sid = "c" + string(cl)
stexcess `spec', indicator(excess) vce(cluster sid) nolog
assert e(N_clust) == `ncl' & "`e(clustvar)'" == "sid"
stexcess `spec', indicator(excess) twostage vce(cluster cl) nolog
assert "`e(vce)'" == "cluster" & e(N_clust) == `ncl'
predict s2, netsurvival ci at(age 0)
assert !missing(s2[1])

// ---- errors ----
rcof "stexcess `spec', indicator(excess) vce(bootstrap)" == 198
stset survtime [iw=one], failure(died)
rcof "stexcess `spec', indicator(excess) vce(robust)" == 101
stset survtime, failure(died)

di as txt _n "vce.do completed."
