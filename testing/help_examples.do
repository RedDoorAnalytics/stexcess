// help_examples.do -- run every example from the help files verbatim, so
// the documented examples are guaranteed to work. Uses the packaged example
// dataset (sysuse when installed; falls back to the repo copy when the
// driver has put data/ on the adopath).

clear all
capture sysuse stexcess_example, clear
if _rc {
    di as err "stexcess_example.dta not found along the adopath"
    exit 601
}

// ---- help stexcess: Examples ----
stset stime, failure(died)
stexcess (age female, df(3))(age female, df(3)), indicator(patient)
stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)
stexcess (age female, df(3))(age female, df(3) tvc(age) dftvc(2)), ///
    indicator(patient)
stexcess (female, df(3) time2(df(3) offset(age)))(age female, df(3)), ///
    indicator(patient)
stexcess (age female, df(3))(age female, df(3)), indicator(patient) twostage

// ---- help stexcess postestimation: Examples ----
sysuse stexcess_example, clear
stset stime, failure(died)
stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)
predict h, hazard ci
predict sn, netsurvival ci
range tt 0 5 100
predict sn3, netsurvival at(age 70 female 1 stage 3) timevar(tt) ci
line sn3 sn3_lci sn3_uci tt, sort
predict ms, survival standardise timevar(tt) ci at(patient 1)
predict rsr, sratio at1(age 70 female 0 patient 1) ///
    at2(age 70 female 0 patient 0) timevar(tt) ci
predict eh, excesshazard standardise timevar(tt) ci

assert !missing(h[1]) & !missing(sn[1]) & !missing(ms[2]) & !missing(rsr[2])

// stronger than non-missing: documented identities hold on the fitted model
predict s_id,  survival    at(age 70 female 1 stage 3) timevar(tt)
predict c_id,  cif         at(age 70 female 1 stage 3) timevar(tt)
predict ch_id, chazard     at(age 70 female 1 stage 3) timevar(tt)
predict lc_id, logchazard  at(age 70 female 1 stage 3) timevar(tt)
assert reldif(c_id, 1 - s_id) < 1e-9 if !missing(c_id)
assert reldif(lc_id, ln(ch_id)) < 1e-9 if !missing(lc_id) & tt > 0
qui su ms if tt > 0                             // tt==0 is the filled limit S=1
assert r(min) > 0 & r(max) < 1                  // standardised survival in (0,1)

di as txt _n "help_examples.do completed."
