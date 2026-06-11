{smcl}
{* *! version 2.0.0 stexcess postestimation (Mata core)}{...}
{vieweralsosee "stexcess" "help stexcess"}{...}
{title:Title}

{p2colset 5 33 35 2}{...}
{p2col :{cmd:stexcess postestimation} {hline 2}}Predictions after stexcess{p_end}
{p2colreset}{...}


{title:Syntax for predict}

{p 8 16 2}
{cmd:predict} {newvar} {ifin} [{cmd:,} {it:statistic} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt :{opt h:azard}}hazard function (the default){p_end}
{synopt :{opt ch:azard}}cumulative hazard function{p_end}
{synopt :{opt logch:azard}}log cumulative hazard function{p_end}
{synopt :{opt surv:ival}}survivor function{p_end}
{synopt :{opt cif}}cumulative incidence function (1 - survival){p_end}
{synopt :{opt rmst}}restricted mean survival time within (0,{it:t}]{p_end}
{synopt :{opt timel:ost}}time lost due to an event within (0,{it:t}]{p_end}
{synopt :{opt nets:urvival}}net (excess-only) survival exp(-H_exc){p_end}
{synopt :{opt exc:esshazard}}excess hazard function{p_end}
{synopt :{opt rmstn:et}}net restricted mean survival time{p_end}
{synopt :{opt hdiff:erence}}difference in hazard functions, {opt at1()} vs {opt at2()}{p_end}
{synopt :{opt sdiff:erence}}difference in survival functions{p_end}
{synopt :{opt cifdiff:erence}}difference in cumulative incidence functions{p_end}
{synopt :{opt rmstdiff:erence}}difference in restricted mean survival times{p_end}
{synopt :{opt hr:atio}}ratio of hazard functions{p_end}
{synopt :{opt sr:atio}}ratio of survival functions{p_end}
{synopt :{opt cifr:atio}}ratio of cumulative incidence functions{p_end}
{synopt :{opt rmstr:atio}}ratio of restricted mean survival times{p_end}
{synoptline}

{synoptset 22 tabbed}{...}
{synopthdr:options}
{synoptline}
{synopt :{opt at(at_spec)}}set covariates, the indicator or offset variables to fixed values{p_end}
{synopt :{opt at1(at_spec)}, {opt at2(at_spec)}}as {opt at()}, for difference/ratio statistics{p_end}
{synopt :{opt zero:s}}set all covariates and the indicator to zero{p_end}
{synopt :{opt stand:ardise}}regression-standardised (g-formula) prediction over the estimation sample{p_end}
{synopt :{opt ci}}delta-method confidence interval in {it:newvar}{cmd:_lci}/{it:newvar}{cmd:_uci}{p_end}
{synopt :{opth ti:mevar(varname)}}evaluate predictions at these time points (default {cmd:_t}){p_end}
{synopt :{opt level(#)}}confidence level{p_end}
{synoptline}

{phang}{it:at_spec} is {it:varname} {it:#} [{it:varname} {it:#} ...]{p_end}


{title:Description}

{pstd}
Predictions follow the original (merlin-based) stexcess: covariates {bf:and
the excess indicator} default to their observed, row-by-row values, so for a
patient record ({it:indicator} = 1) {cmd:hazard} is the total hazard
{it:h_ref + h_exc}, while for a control record it is the reference hazard.
{opt at()} overrides any of them with fixed values -- in particular the
indicator itself: {cmd:at(}{it:indicator}{cmd: 1)} forces total quantities
and {cmd:at(}{it:indicator}{cmd: 0)} reference-only quantities, regardless of
each record's arm. Variables usable inside {opt at()} are the {bf:underlying}
model covariates (for factor variables and interactions, the variable itself,
e.g. {cmd:at(stage 3)} -- the indicator/interaction terms are recomputed from
the overridden values), the indicator and any timescale offset variables
(stored in {cmd:e(atvars)}). {opt zeros} sets the underlying covariates and
the indicator (but not offset variables) to zero; {opt at()} wins where both
apply.

{pstd}
The net statistics ({opt netsurvival}, {opt excesshazard}, {opt rmstnet})
involve only the excess component and ignore the indicator.

{pstd}
All confidence intervals are analytic delta-method intervals on a suitable
transformed scale (log for hazards, cumulative hazards and RMST-type
quantities; complementary log-log for survival-type quantities; identity for
differences).


{title:Standardisation}

{pstd}
{opt standardise} averages the requested statistic over the estimation-sample
covariate distribution (the g-formula), with delta-method CIs; when the data
are stset with weights, the average is weighted accordingly. It may be
combined with {opt at()}/{opt zeros}, which then pin selected variables
across the whole population -- counterfactual standardisation. For example,
after {cmd:stexcess (age)(age), indicator(patient)}:

{phang2}{cmd:. range tt 0 5 50}{p_end}
{phang2}{cmd:. predict ms, survival standardise timevar(tt) ci at(patient 1)}{p_end}

{pstd}
gives the marginal all-cause survival had everyone been a patient, averaged
over the observed age distribution. Difference/ratio statistics and
{opt logchazard} are not available with {opt standardise}.


{title:Remarks}

{pstd}
Predictions are unconditional from time 0, also after fitting with delayed
entry. Out-of-sample prediction is allowed: predictions are computed at
{opt timevar()} for every observation in the {it:if/in} sample using that
row's (possibly overridden) covariate, indicator and offset values.
Quantities that are undefined at {it:t} = 0 on the log-time scale are filled
in with their limits (survival-type statistics 1; cumulative-hazard,
CIF and RMST-type statistics 0).


{title:Examples}

{pstd}Setup, using the simulated dataset shipped with the package{p_end}

{phang2}{stata `"use https://raw.githubusercontent.com/RedDoorAnalytics/stexcess/main/data/stexcess_example.dta, clear"':. use stexcess_example.dta, clear  (from GitHub)}{p_end}
{phang2}{stata "stset stime, failure(died)":. stset stime, failure(died)}{p_end}
{phang2}{stata "stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)":. stexcess (age female, df(3))(age female i.stage, df(3)), indicator(patient)}{p_end}

{pstd}Hazard and net survival, observed covariates and arm{p_end}

{phang2}{stata "predict h, hazard ci":. predict h, hazard ci}{p_end}
{phang2}{stata "predict sn, netsurvival ci":. predict sn, netsurvival ci}{p_end}

{pstd}Net survival for a 70-year-old woman with stage III disease, on a
time grid{p_end}

{phang2}{stata "range tt 0 5 100":. range tt 0 5 100}{p_end}
{phang2}{stata "predict sn3, netsurvival at(age 70 female 1 stage 3) timevar(tt) ci":. predict sn3, netsurvival at(age 70 female 1 stage 3) timevar(tt) ci}{p_end}
{phang2}{stata "line sn3 sn3_lci sn3_uci tt, sort":. line sn3 sn3_lci sn3_uci tt, sort}{p_end}

{pstd}Marginal (standardised) all-cause survival had everyone been a
patient, averaged over the observed covariate distribution{p_end}

{phang2}{stata "predict ms, survival standardise timevar(tt) ci at(patient 1)":. predict ms, survival standardise timevar(tt) ci at(patient 1)}{p_end}

{pstd}Relative survival ratio (total vs reference) for 70-year-old men{p_end}

{phang2}{stata "predict rsr, sratio at1(age 70 female 0 patient 1) at2(age 70 female 0 patient 0) timevar(tt) ci":. predict rsr, sratio at1(age 70 female 0 patient 1) at2(age 70 female 0 patient 0) timevar(tt) ci}{p_end}

{pstd}Standardised excess hazard{p_end}

{phang2}{stata "predict eh, excesshazard standardise timevar(tt) ci":. predict eh, excesshazard standardise timevar(tt) ci}{p_end}
