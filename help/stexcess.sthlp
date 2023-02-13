{smcl}
{* *! version 1.0.0  ?????2023}{...}
{vieweralsosee "stexcess postestimation" "help stexcess_postestimation"}{...}
{vieweralsosee "merlin" "help merlin"}{...}
{vieweralsosee "stmerlin" "help stmerlin"}{...}
{title:Title}

{p2colset 5 17 19 2}{...}
{p2col:{helpb stexcess} {hline 2}}estimation of a flexible parametric excess 
hazard model with a modelled expected rate{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 12 2}
{cmd:stexcess} {bf:(}{it:reference_model}{bf:)}{bf:(}{it:excess_model}{bf:)} 
{ifin} , 
{bf:indicator}({varname}) [{bf:,} {help stexcess##options:{it:options}}]

{phang2}
where {it:control_model} and {it:excess_model} are specified with
{p_end}

{phang3}        
[{it:{help merlin_models:indepsyntax}}] , 
{help stexcess##model_options:{it:model_options}}
{p_end}
        
{phang3}and {it:{help merlin_models:indepsyntax}} is a {helpb merlin} 
linear predictor, which can be anything from a simple {varlist}, to 
directly specifying spline or fractional polynomial functions of continuous 
covariates.{p_end}
	
{phang2}You must {cmd:stset} your data before using {cmd:stexcess}; 
see {manhelp stset ST}.{p_end}


{synoptset 27}{...}
{marker options}{...}
{synopthdr:options}
{synoptline}
{synopt:{opth ind:icator(varname)}}an indicator variable that must be coded 
as 0 for control observations and 1 for excess observations{p_end}
{synopt:{opt chintp:oints(#)}}defines the number of Gauss-Legendre 
quadrature points used when calculating the cumulative hazard function{p_end}
{synopt :{opt eform}}display exponentiated coefficients{p_end}
{synopt:{it:maximize_opts}}control the maximization process{p_end}
{synopt:{it:display_opts}}display options{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 27}{...}
{marker model_options}{...}
{synopthdr:model options}
{synoptline}
{synopt :{opt nocons:tant}}omit the constant term{p_end}
{synopt :{opth df(#)}}degrees of freedom for the baseline function with 
{cmd:rp} or {cmd:rcs} models; see details{p_end}
{synopt :{opt knots(knots_list)}}knot locations for the baseline function 
with {cmd:rp} or {cmd:rcs} models; see details{p_end}
{synopt :{opth tvc(varlist)}}time-dependent effects{p_end}
{synopt :{opth dftvc(numlist)}}degrees of freedom for each time-dependent 
effect{p_end}
{synopt :{opt tvctime}}use splines of time rather than log time for 
time-dependent effects{p_end}
{synopt :{opt noorth:og}}turns off the default orthogonalisation of any 
spline terms{p_end}
{synopt:{bf:time#(}{help stmerlin##mt_opts:{it:mt_opts}})}define two to five 
additional timescales modelled with restricted cubic 
splines, specified with {cmd:time2({help stmerlin##mt_opts:{it:mt_opts}})}, 
with a maximum of {cmd:time5({help stmerlin##mt_opts:{it:mt_opts}})}{p_end}
{synoptline}
{p2colreset}{...}

{synoptset 27}{...}
{marker mt_opts}{...}
{synopthdr:multiple timescale options}
{synoptline}
{synopt :{opth offset(varname)}}defines the offset to be added to {cmd:_t} 
to define the additional timescale{p_end}
{synopt :{opth moffset(varname)}}defines the offset to be taken away 
("minused") from {cmd:_t} to define the additional timescale{p_end}
{synopt :{opth df(#)}}degrees of freedom for timescale spline function{p_end}
{synopt :{opt knots(knots_list)}}knot locations for the timescale spline 
function{p_end}
{synopt :{opt time}}use splines of time rather than log time for 
the timescale{p_end}
{synopt :{opth tvc(varlist)}}time-dependent effects on the additional 
timescale{p_end}
{synopt :{opth dftvc(numlist)}}degrees of freedom for each time-dependent 
effect{p_end}
{synopt :{opt tvctime}}use splines of time rather than log time for 
time-dependent effects{p_end}
{synopt :{opt noorth:og}}turns off the default orthogonalisation of any 
spline terms{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:stexcess} fits modelled excess survival models, using restricted cubic 
splines on the log hazard scale. Multiple timescales can also be included. 
Time-dependent effects can be specified using restricted cubic splines 
through options, or in alternative ways using {helpb merlin}'s linear 
predictor syntax. For predictions available post-estimation, see 
{helpb stexcess postestimation}.
{p_end}

{pstd}
The {helpb merlin} command fits an extremely broad class of mixed effects 
regression models for linear, non-linear and user-defined outcomes. For 
full details and many tutorials, take a look at the accompanying 
website: {browse "https://reddooranalytics.se/products/merlin":{bf:reddooranalytics.se/products/merlin}}.
{p_end}


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt indicator(varname)} an indicator variable that must be coded 
as 0 for control observations and 1 for excess observations.
{p_end}

{phang}
{opt chintpoints(varname)} defines the number of Gauss-Legendre 
quadrature points used when calculating the total cumulative hazard function. 
A higher number of quadrature points increases the accuracy, but is 
computationally slower.
{p_end}

{phang}
{opt eform} displays exponentiated coefficients for the main linear 
predictors.
{p_end}

{phang}
{it:maximize_opts} control the maxi`ation process.
{p_end}

{phang}
{it:display_opts} display options.
{p_end}

{marker model_options}{...}
{dlgtab:Model options}

{phang}{opt noconstant} suppresses the constant (intercept) term in the linear predictor.

{phang}{opt df(#)} degrees of freedom for the baseline log [cumulative] hazard function, i.e. number of restricted cubic 
spline terms when using {cmd:distribution(rp)} or {cmd:distribution(rcs)}. Internal knots are placed at centiles of the 
event times. Boundary knots are placed at the minimum and maximum event times.

{phang}{opt knots(knots_list)} either:

{phang2}defines the knot locations for the spline functions used to model the baseline 
log [cumulative] hazard function when using {cmd:distribution(rp)} or {cmd:distribution(rcs)}. Must include boundary 
knots. Knots should be specified in increasing order.

{phang2}defines the knot locations (cut-points) of the baseline function for {cmd:distribution(pwexponential)}. Knots should 
be specified in increasing order.

{phang}{opt tvc(varlist)} specifies the variables that have time-dependent effects. Time-dependent effects are fitted 
using restricted cubic splines of time or log time (the default). The degrees of freedom are specified using the 
{cmd:dftvc()} option. Note, {cmd:tvc()}s are not supported with generalised gamma, log normal or log logistic models.

{phang}{opt dftvc(numlist)} degrees of freedom for the time-dependent effects specified in {cmd:tvc()}. If only one number is 
specified, then the same degrees of freedom are applied to all {cmd:tvc()}s, otherwise, a number must be specified for each.

{phang}{opt tvctime} specified that restricted cubic splines of time are used to model time-dependent effects, rather than 
the default of log time.

{phang}{opt noorthog} suppresses orthogonal transformation of spline variables.

{marker multitime_details}{...}
{dlgtab:Multiple timescales}

{phang}
{opt offset(varname)} defines the offset to be added to {cmd:_t} to define the additional timescale. If time since diagnosis 
was the main timescale, and you wish to add attained age as a second timescale, the {cmd:offset()} would contain age at diagnosis.

{phang}
{opt moffset(varname)} defines the offset to be taken away ("minused") from {cmd:_t} to define the additional timescale, i.e. to 
reset the clock.

{phang}{opt df(#)} degrees of freedom for the additional timescale function, i.e. number of restricted cubic 
spline terms. Internal knots are placed at centiles of the event times. Boundary knots are placed at the minimum and maximum event times.

{phang}{opt knots(knots_list)} defines the knot locations for the spline functions used to model the additional timescale. Must include 
boundary knots. Knots should be specified in increasing order.

{phang}{opt time} specifies that restricted cubic splines of time are used to model the additional timescale, rather than the default of log time.

{phang}{opt tvc(varlist)} specifies the variables that have time-dependent effects on the additional timescale. Time-dependent effects 
are fitted using restricted cubic splines of time or log time (the default). The degrees of freedom are specified using the 
{cmd:dftvc()} option. Note, {cmd:tvc()}s are not supported with generalised gamma, log normal or log logistic models.

{phang}{opt dftvc(numlist)} degrees of freedom for the time-dependent effects specified in {cmd:tvc()}. If only one number is 
specified, then the same degrees of freedom are applied to all {cmd:tvc()}s, otherwise, a number must be specified for each.

{phang}{opt tvctime} specified that restricted cubic splines of time are used to model time-dependent effects on the additional timescale, 
rather than the default of log time.

{phang}{opt noorthog} suppresses orthogonal transformation of spline variables (additional timescale and any time-dependent effects).


{title:Examples}

{phang}Fit a flexible parametric excess hazard model:{p_end}
{cmd:    . use simdata, clear}
{cmd:    . stset stime, failure(died)}
{cmd:    . stmerlin (age sex, df(3))(age sex, df(2)), indicator(cancer)}

{phang}Model time-dependent effects:{p_end}
{cmd:    . stmerlin (age sex, df(3) tvc(age) dftvc(1))}
{cmd:               (age sex, df(2)), indicator(cancer)}


{title:Author}

{p 5 12 2}
{bf:Michael J. Crowther}{p_end}
{p 5 12 2}
Red Door Analytics AB{p_end}
{p 5 12 2}
Stockholm, Sweden{p_end}
{p 5 12 2}
michael@reddooranalytics.se{p_end}
{p 5 12 2}
{browse "https://reddooranalytics.se":{bf:reddooranalytics.se}}
