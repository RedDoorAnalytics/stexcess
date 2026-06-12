// run_all.do -- run the full test battery. From the repository root:
//   stata-mp -q -b do testing/run_all.do
//
// Needs the merlin package for merlin_oracle.do (it skips itself if absent).
//
// On-demand validations not run here:
//   testing/weights_mc.do      pweight sandwich Monte Carlo (~3 minutes)
//   stexcess-dev repo          Python/JAX-oracle parity (testing/
//                              mata_parity.do, needs that repo's venv) and
//                              the v1-command cross-checks (parity1/3 fit
//                              v1 WITHOUT this repo's ado/ on the adopath,
//                              then parity2/4 compare)

clear all
adopath ++ "`c(pwd)'/ado"
adopath ++ "`c(pwd)'/data"          // so sysuse finds the example dataset
mata: mata mlib index

do "testing/basic_run.do"
do "testing/errors.do"
do "testing/weights.do"
do "testing/vce.do"
do "testing/conditional.do"
do "testing/multiple_timescales.do"
do "testing/factor_vars.do"
do "testing/help_examples.do"
do "testing/consistency.do"
do "testing/merlin_oracle.do"

di as txt _n "run_all.do: all test files completed."
