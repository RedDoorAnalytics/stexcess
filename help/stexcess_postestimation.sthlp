{smcl}
{* *! version 1.0.0  06feb2023}{...}
{vieweralsosee "stexcess" "help stexcess"}{...}
{vieweralsosee "merlin" "help merlin"}{...}
{viewerjumpto "Syntax" "stexcess_postestimation##syntax"}{...}
{viewerjumpto "Description" "stexcess_postestimation##description"}{...}
{viewerjumpto "Options" "stexcess_postestimation##options"}{...}
{viewerjumpto "Remarks" "stexcess_postestimation##remarks"}{...}
{viewerjumpto "Examples" "stexcess_postestimation##examples"}{...}

{marker syntax}{...}
{title:Syntax for predict}

{pstd}
Syntax for predictions following a {helpb stexcess} model

{p 8 16 2}
{cmd:predict}
{it:newvarname}
{ifin} [{cmd:,}
{it:{help stexcess_postestimation##statistic:statistic}}
{it:{help stexcess_postestimation##opts_table:options}}]


{marker statistic}{...}
{synoptset 22 tabbed}{...}
{synopthdr:statistic}
{synoptline}
{synopt :{opt surv:ival}}survivor function{p_end}
{synopt :{opt cif}}cumulative incidence function{p_end}
{synopt :{opt h:azard}}hazard function{p_end}
{synopt :{opt ch:azard}}cumulative hazard function{p_end}
{synopt :{opt logch:azard}}log cumulative hazard function{p_end}
{synopt :{opt rmst}}restricted mean survival time, within (0,{it:t}]{p_end}
{synopt :{opt timel:ost}}time lost due to an event, within (0,{it:t}]{p_end}
{synopt :{opt hdiff:erence}}difference in hazard functions{p_end}
{synopt :{opt sdiff:erence}}difference in survival functions{p_end}
{synopt :{opt cifdiff:erence}}difference in cumulative incidence functions{p_end}
{synopt :{opt rmstdiff:erence}}difference in restricted mean survival functions{p_end}
{synopt :{opt hr:atio}}ratio of hazard functions{p_end}
{synopt :{opt sr:atio}}ratio of survival functions{p_end}
{synopt :{opt cifr:atio}}ratio of cumulative incidence functions{p_end}
{synopt :{opt rmstr:atio}}ratio of restricted mean survival functions{p_end}
{synoptline}

{marker opts_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr:options}
{synoptline}
{syntab:Main}
{synopt :{opt at(at_spec)}}specify covariate values for prediction{p_end}
{synopt :{opt zero:s}}set all covariates to zero{p_end}
{synopt :{opt at1(at_spec)}}specify covariate values for prediction; for use with difference and ratio predictions{p_end}
{synopt :{opt at2(at_spec)}}specify covariate values for prediction; for use with difference and ratio predictions{p_end}
{synopt :{opt ci}}calculate confidence intervals{p_end}
{synopt :{cmd:timevar(}{varname}{cmd:)}}calculate predictions at specified time-points{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:predict} is a standard postestimation command of Stata.
This entry concerns use of {cmd:predict} after {helpb stexcess}.

{pstd}
{cmd:predict} after {cmd:stexcess} creates new variables containing
observation-by-observation values of estimated observed response variables,
linear predictions of observed response variables, or other such functions.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang} 
{cmd:survival} calculates the survival function. 

{phang} 
{cmd:cif} calculates the cumulative incidence function at time {it:t}, where {it:t} is the time at which predictions are made. 
This is 1 - survival. 

{phang} 
{cmd:hazard} calculates the hazard function at time {it:t}, where {it:t} is the time at which predictions are made. 

{phang}
{cmd:chazard} calculates the cumulative hazard function at time {it:t}, where {it:t} is the time at which predictions are made. 

{phang}
{cmd:logchazard} calculates the log of the cumulative hazard function at time {it:t}, where {it:t} is the time at which predictions are made. 

{phang} 
{cmd:rmst} calculates the restricted mean survival time, which is the integral of the survival function within the interval 
(0,{it:t}], where {it:t} is the time at which predictions are made. 

{phang} 
{cmd:timelost} calculates the time lost due to the event occuring, within the interval (0,{it:t}]. 
This is the integral of the {cmd:cif} between (0,{it:t}].

{phang} 
{cmd:hdifference} calculates the difference in hazard function at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:sdifference} calculates the difference in survival function at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:cifdifference} calculates the difference in cumulative incidence function at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:rmstdifference} calculates the difference in restricted mean survival time at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:hratio} calculates the ratio of hazard functions at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:sratio} calculates the ratio of survival functions at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:cifratio} calculates the ratio of cumulative incidence functions at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang} 
{cmd:rmstratio} calculates the ratio of restricted mean survival times at time {it:t}, where {it:t} is the time at which predictions are made, 
across the covariate patterns specified in {cmd:at1()} and {cmd:at2()}. 

{phang}
{opt at(varname # [ varname # ...])} requests that the covariates specified by 
the listed {it:varname}(s) be set to the listed {it:#} values. For example,
{cmd:at(trt 1 age 50)} would evaluate predictions at {cmd:trt} = 1 and
{cmd:age} = 50. This is a useful way to obtain out of sample predictions. Other covariates in your model, but 
{bf:not} included in {cmd:at()} will be set to their observed values, i.e. the values in your dataset. Note that if {cmd:at()} 
is used together with {cmd:zeros}, all covariates not listed in {cmd:at()} are set to zero. See also {cmd:zeros}.

{phang}
{opt zeros} sets all covariates to zero. See also {cmd:at()}. Note, any response variables will be skipped, i.e. not set 
to zero, so if a response variable for one model is included as a covariate in another - it will {it:not} be set to zero. Also note 
that it {cmd:at1()} and {cmd:at2()} are specified, then {cmd:zeros} applies to both.

{phang}
{opt at1(varname # [ varname # ...])} does the same as {cmd:at()} but for use in conjunction with {cmd:?difference} or 
{cmd:?ratio} predictions.

{phang}
{opt at2(varname # [ varname # ...])} does the same as {cmd:at()} but for use in conjunction with {cmd:?difference} or 
{cmd:?ratio} predictions.

{phang}
{cmd:ci} specifies that confidence intervals are calculated for the predicted {it:statistic}. The multivariate delta 
method (i.e. {cmd:predictnl}) is used for all calculations, except when a {cmd:family(cox)} model has been fitted, 
in which case bootstrapping is used. The calculated confidence intervals are generated in {it:newvarname_lci} 
and {it:newvarname_uci}.

{phang}
{cmd:timevar(}{varname}{cmd:)} calculate predictions at specified time-points. 
For survival models, the default is to calculate predictions at the response times. 
For a {cmd:merlin} model where a {cmd:timevar()} was specified, then the default will use the original 
{cmd:timevar()}. This option overides it.{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
Out-of-sample prediction is allowed for all {cmd:predict} options.


{marker examples}{...}
{title:Examples}

{phang}Fit a flexible parametric excess hazard model:{p_end}
{cmd:    . use simdata, clear}
{cmd:    . stset stime, failure(died)}
{cmd:    . stmerlin (age sex, df(3))(age sex, df(2)), indicator(cancer)}

{phang}Predict the survival function for the control group:{p_end}
{cmd:    . predict s1, survival at(cancer 0)}

{phang}Predict the survival function for the cancer group:{p_end}
{cmd:    . predict s2, survival at(cancer 1)}
