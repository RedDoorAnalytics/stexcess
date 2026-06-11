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
each record's arm. Variables usable inside {opt at()} are the model
covariates, the indicator and any timescale offset variables (stored in
{cmd:e(atvars)}). {opt zeros} sets covariates and the indicator (but not
offset variables) to zero; {opt at()} wins where both apply.

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
covariate distribution (the g-formula), with delta-method CIs. It may be
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
Out-of-sample prediction is allowed: predictions are computed at
{opt timevar()} for every observation in the {it:if/in} sample using that
row's (possibly overridden) covariate, indicator and offset values.
Quantities that are undefined at {it:t} = 0 on the log-time scale are filled
in with their limits (survival-type statistics 1; cumulative-hazard,
CIF and RMST-type statistics 0).


{title:Examples}

{phang}Total and net survival for a patient aged 60:{p_end}
{phang2}{cmd:. predict s_all, survival at(age 60 patient 1) ci}{p_end}
{phang2}{cmd:. predict s_net, netsurvival at(age 60) ci}{p_end}

{phang}Relative survival ratio (total vs reference) at age 60:{p_end}
{phang2}{cmd:. predict rsr, sratio at1(age 60 patient 1) at2(age 60 patient 0) ci}{p_end}

{phang}Standardised excess hazard over the patient covariate distribution:{p_end}
{phang2}{cmd:. predict eh, excesshazard standardise timevar(tt) ci}{p_end}
