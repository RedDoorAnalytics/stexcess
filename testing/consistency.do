// consistency.do -- standard postestimation tooling must work off e():
// test, testparm, lincom, estat ic, estimates store/table/replay.

clear all
capture sysuse stexcess_example, clear
if _rc {
    di as err "stexcess_example.dta not found along the adopath"
    exit 601
}
stset stime, failure(died)
stexcess (age female, df(3))(age female i.stage, df(3)), ///
    indicator(patient) nolog

test [ref]age
assert r(p) < . & r(df) == 1
testparm i.stage, equation(exc)
assert r(p) < . & r(df) == 2
lincom [exc]3.stage - [exc]2.stage
assert r(estimate) < . & r(se) < .
lincom [exc]age, eform
estat ic
matrix S = r(S)
assert reldif(S[1, 3], e(ll)) < 1e-12 & S[1, 4] == e(rank)

estimates store m1
stexcess (age female, df(2))(age female, df(2)), indicator(patient) nolog
estimates store m2
estimates table m1 m2, keep(ref:age exc:age) se
estimates restore m1
stexcess                                      // replay from restored e()
lrtest m2 m1, force                           // nested spline df, same terms

di as txt _n "consistency.do completed."
