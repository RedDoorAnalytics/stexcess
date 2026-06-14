{smcl}
{* *! version 2.0.0 stexcess (Mata core)}{...}
{vieweralsosee "stexcess postestimation" "help stexcess_postestimation"}{...}
{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:stexcess} {hline 2}}Modelled excess hazard models{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{cmd:stexcess} {cmd:(}{it:reference_model}{cmd:)} {cmd:(}{it:excess_model}{cmd:)}
{ifin}{cmd:,} {opth ind:icator(varname)} [{it:options}]

{phang2}where {it:reference_model} and {it:excess_model} are each specified as{p_end}

{phang3}[{varlist}] [{cmd:,} {it:model_options}]{p_end}

{phang2}The data must be {help stset}. Delayed entry (left truncation,
{cmd:stset, enter()}) is supported: each record contributes the hazard
integral over its own ({it:t0}, {it:t}]. Weights are taken from the
{cmd:stset} specification: {cmd:fweight}s reproduce the expanded dataset
exactly (including default knot placement), {cmd:pweight}s give
design-consistent estimates with a robust (sandwich) variance, and
{cmd:iweight}s weight the log likelihood with a model-based variance.
{opt vce(robust)} and {opt vce(cluster clustvar)} request (cluster-)robust
sandwich variances for the joint fit, with the usual finite-sample factors;
{opt vce(cluster clustvar)} also clusters the {opt twostage} sandwich.
{it:varlist} may contain factor variables and interactions
(e.g. {cmd:i.stage} or {cmd:c.age##i.sex}); base and omitted levels are
shown in the coefficient table with zero coefficients. Factor variables are
not allowed inside {opt tvc()} or as offset variables.{p_end}

{synoptset 27 tabbed}{...}
{synopthdr:options}
{synoptline}
{synopt :{opth ind:icator(varname)}}required; 1 = excess (patient) record, 0 = reference (control){p_end}
{synopt :{opt two:stage}}two-stage estimation with sandwich (robust) standard errors{p_end}
{synopt :{opt vce(robust)}, {opt vce(cluster clustvar)}}robust or cluster-robust standard errors{p_end}
{synopt :{opt chint:points(#)}}Gauss-Legendre nodes for cumulative hazards (default 30){p_end}
{synopt :{opt from(matname)}}starting values (one column per model parameter){p_end}
{synopt :{opt iter:ate(#)}}maximum number of Newton-Raphson iterations{p_end}
{synopt :{opt tol:erance(#)}}coefficient-vector convergence tolerance{p_end}
{synopt :{opt ltol:erance(#)}}log-likelihood convergence tolerance{p_end}
{synopt :{opt nrtol:erance(#)}}scaled-gradient convergence tolerance{p_end}
{synopt :{opt eform}}display exponentiated coefficients{p_end}
{synopt :{opt nolog}}suppress the log-likelihood iteration log{p_end}
{synopt :{opt level(#)}}confidence level{p_end}
{synoptline}

{synoptset 27 tabbed}{...}
{synopthdr:model_options}
{synoptline}
{synopt :{opt df(#)}}degrees of freedom for the baseline spline (default 3){p_end}
{synopt :{opth knots(numlist)}}baseline knot locations (transformed scale, including boundaries){p_end}
{synopt :{opt time}}baseline spline in time rather than log time{p_end}
{synopt :{opth off:set(varname)}}added to {cmd:_t} before the baseline spline is formed{p_end}
{synopt :{opth moff:set(varname)}}subtracted from {cmd:_t} before the baseline spline is formed{p_end}
{synopt :{opth tvc(varlist)}}covariates with time-varying effects{p_end}
{synopt :{opth dftvc(numlist)}}df for each time-varying effect; required with {opt tvc()}{p_end}
{synopt :{opt tvctime}}time-varying effect splines in time rather than log time{p_end}
{synopt :{cmd:time2(}{it:mt_opts}{cmd:)}}additional timescale; up to {cmd:time5()}{p_end}
{synopt :{opt noorth:og}}do not orthogonalise this component's spline bases{p_end}
{synopt :{opt nocons:tant}}omit the constant term{p_end}
{synoptline}

{synoptset 27 tabbed}{...}
{synopthdr:mt_opts}
{synoptline}
{synopt :{opt df(#)}}degrees of freedom for the timescale spline{p_end}
{synopt :{opth knots(numlist)}}knot locations (transformed scale); one of {opt df()}/{opt knots()} is required{p_end}
{synopt :{opth offset(varname)}}added to {cmd:_t} to define the timescale (e.g. age at diagnosis for attained age){p_end}
{synopt :{opth moffset(varname)}}subtracted from {cmd:_t} to define the timescale (resets the clock){p_end}
{synopt :{opt time}}timescale spline in time rather than log time{p_end}
{synopt :{opth tvc(varlist)}}covariates with time-varying effects on this timescale{p_end}
{synopt :{opth dftvc(numlist)}}df for each such effect; required with {opt tvc()}{p_end}
{synopt :{opt tvctime}}those effect splines in time rather than log time{p_end}
{synopt :{opt noorthog}}do not orthogonalise this timescale's spline bases{p_end}
{synoptline}


{title:Description}

{pstd}
{cmd:stexcess} fits a relative-survival / excess-hazard model in which the
expected (reference) rate is itself a fitted hazard model estimated from a
control cohort, rather than taken from external population life tables. The
reference and excess hazards are flexible parametric models -- restricted
cubic splines on the log-hazard scale -- estimated jointly; the reference
parameters are shared across the control records and the excess records'
total hazard. The syntax follows the original (merlin-based) stexcess v1.

{pstd}
Estimation runs entirely inside Stata: the designs are built once in Mata and
maximised with Mata's {helpb mf_optimize:optimize()} -- the Newton-Raphson
engine underneath {helpb ml} -- using an exact analytic log likelihood,
gradient and Hessian. Predictions and standardisation use analytic
delta-method Jacobians. There are no dependencies outside Stata. Starting
values for the joint fit come from first fitting the reference model to the
control records alone, with the excess parameters started at zero apart from
a crude event-rate intercept; {opt from()} overrides this.

{pstd}
By default both components are estimated jointly by maximum likelihood, so
patient (excess) records also contribute information about the reference
hazard. With {opt twostage} the reference model is instead fitted to the
control records only (stage 1) and the excess model is then fitted to the
patient records with the reference parameters held fixed (stage 2), analogous
to supplying a known background hazard. This insulates the reference fit from
any misspecification of the excess model. Standard errors come from a stacked
M-estimation sandwich variance that propagates stage-1 uncertainty into the
excess parameters; all postestimation (CIs, standardisation, contrasts) uses
this joint variance automatically.


{title:Timescales and knots}

{pstd}
Each component's baseline is a restricted cubic spline in log time (or time,
with {opt time}), optionally shifted by {opt offset()}/{opt moffset()} before
the transform. {cmd:time2()}-{cmd:time5()} add further timescales as splines
in {cmd:_t} {it:+ offset - moffset}, each on its own log or natural scale --
e.g. {cmd:time2(df(3) offset(agediag))} adds an attained-age timescale when
age at diagnosis is in {cmd:agediag}. Time-varying effects are covariate x
spline interactions on the relevant timescale; main-timescale {opt tvc()}
splines are functions of {cmd:_t} alone (no offset), matching v1.

{pstd}
Knots in any {opt knots()} option are on the {bf:transformed} scale of that
spline (log scale unless {opt time} is given), in ascending order, and
include the two boundary knots, so a spline with {it:df} degrees of freedom
has {it:df}+1 knots. Explicit knots override df-based placement. This is the
same scale and format in which the original merlin-based {cmd:stexcess}
stores its knots (e.g. {cmd:e(knots_1_2_1)}), so those can be passed straight
through to reproduce a v1 fit. Default knot siting differs from v1 in two
small ways: knots are placed by interpolated (rather than rounded
order-statistic) centiles, and the excess component's knots are sited from
patient ({it:indicator} = 1) events only, whereas v1 uses all events --
explicit knots are how v1 fits are reproduced exactly. The knots actually
used are stored in {cmd:e(knotsref)}, {cmd:e(knotsexc)} and, for additional
timescales, {cmd:e(knotsref_t}{it:#}{cmd:)} / {cmd:e(knotsexc_t}{it:#}{cmd:)}.


{title:Examples}

{pstd}The examples use the package's simulated dataset: a cancer cohort
with matched population controls and 5-year follow-up. It loads from
GitHub below; {cmd:net get stexcess} downloads a local copy.{p_end}

{phang2}{stata `"use https://raw.githubusercontent.com/RedDoorAnalytics/stexcess/main/data/stexcess_example.dta, clear"':. use stexcess_example.dta, clear  (from GitHub)}{p_end}
{phang2}{stata "stset stime, failure(died)":. stset stime, failure(died)}{p_end}

{pstd}Joint flexible parametric excess hazard model{p_end}

{phang2}{stata "stexcess (age female, df(3))(age female, df(3)), indicator(patient)":. stexcess (age female, df(3))(age female, df(3)), indicator(patient)}{p_end}

{pstd}Stage effects on the excess hazard, via factor variables{p_end}

{phang2}{stata "stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)":. stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)}{p_end}

{pstd}A time-varying excess effect of age{p_end}

{phang2}{stata "stexcess (age female, df(3))(age female, df(3) tvc(age) dftvc(2)), indicator(patient)":. stexcess (age female, df(3))(age female, df(3) tvc(age) dftvc(2)), indicator(patient)}{p_end}

{pstd}Attained age as an additional reference timescale (age at diagnosis
enters as the offset, so the timescale is age + follow-up time){p_end}

{phang2}{stata "stexcess (female, df(3) time2(df(3) offset(age)))(age female, df(3)), indicator(patient)":. stexcess (female, df(3) time2(df(3) offset(age)))(age female, df(3)), indicator(patient)}{p_end}

{pstd}Two-stage estimation: reference fitted to the controls only, with a
stacked sandwich variance{p_end}

{phang2}{stata "stexcess (age female, df(3))(age female, df(3)), indicator(patient) twostage":. stexcess (age female, df(3))(age female, df(3)), indicator(patient) twostage}{p_end}


{title:Postestimation}

{pstd}See {help stexcess_postestimation:stexcess postestimation}.


{title:Stored results}

{synoptset 18 tabbed}{...}
{p2col 5 18 22 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(N_ref)}}number of reference (control) records{p_end}
{synopt:{cmd:e(N_exc)}}number of excess (patient) records{p_end}
{synopt:{cmd:e(N_fail)}}number of failures{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom (covariate coefficients){p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(dfref)}, {cmd:e(dfexc)}}baseline spline df{p_end}
{synopt:{cmd:e(chintpoints)}}quadrature nodes{p_end}
{synopt:{cmd:e(converged)}}1 if the optimiser converged{p_end}
{p2col 5 18 22 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:stexcess}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(method)}}{cmd:joint} or {cmd:twostage}{p_end}
{synopt:{cmd:e(vce)}}{cmd:oim}, {cmd:robust} or {cmd:cluster}{p_end}
{synopt:{cmd:e(clustvar)}, {cmd:e(N_clust)}}cluster variable and number of clusters{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust} with {opt twostage} or {cmd:pweight}s{p_end}
{synopt:{cmd:e(wtype)}, {cmd:e(wexp)}}weight type and variable, when stset with weights{p_end}
{synopt:{cmd:e(indicator)}}excess indicator variable{p_end}
{synopt:{cmd:e(refvars)}, {cmd:e(excvars)}}covariates in each equation{p_end}
{synopt:{cmd:e(knotsref)}, {cmd:e(knotsexc)}}baseline knots used (transformed scale){p_end}
{synopt:{cmd:e(knotsref_t}{it:#}{cmd:)}, {cmd:e(knotsexc_t}{it:#}{cmd:)}}additional-timescale knots used{p_end}
{synopt:{cmd:e(atvars)}}variables that may appear in predict's {cmd:at()} options{p_end}
{p2col 5 18 22 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector ({cmd:ref:} and {cmd:exc:} equations){p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}

{pstd}
Following {helpb streg}, the record counts {cmd:e(N_ref)}, {cmd:e(N_exc)} and
{cmd:e(N_fail)} are sums of the weights when the data are {cmd:stset} with
weights (so {cmd:e(N_ref)} + {cmd:e(N_exc)} need not equal {cmd:e(N)} under
{cmd:iweight}s or {cmd:pweight}s). {cmd:e(N)} is the sum of the weights under
{cmd:fweight}s and the physical observation count otherwise.


{title:Author}

{p 5 12 2}{bf:Michael J. Crowther}{p_end}
{p 5 12 2}Red Door Analytics AB{p_end}
{p 5 12 2}Stockholm, Sweden{p_end}
{p 5 12 2}michael.crowther@reddooranalytics.se{p_end}
