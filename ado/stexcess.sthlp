{smcl}
{* *! version 1.0.0 03feb2022}{...}
{title:Title}

{p2colset 5 15 18 2}{...}
{p2col :{hi:nnplot} {hline 2}}plot an artificial neural network{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{phang2}
{cmd: nnplot} {cmd:,} {it:options}

{synoptset 29 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opth outputs(#)}}number of outputs{p_end}
{synopt:{opth inputs(#)}}number of inputs{p_end}
{synopt:{opth hlayers(#)}}number of hidden layers{p_end}
{synopt:{opth hnodes(numlist)}}number of nodes per hidden layer{p_end}
{synopt:{opt colors(list)}}list of colors to use in each layer{p_end}
{synopt:{opt inlabels(list)}}labels for inputs{p_end}
{synopt:{opt outlabels(list)}}labels for outputs{p_end}
{synopt:{opt hlabel(string)}}labels for hidden layers{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:nnplot} plots a schematic diagram of an example neural network. {cmd:nnplot} allows any number of outputs (responses), 
inputs (covariates), hidden layers, and nodes in each layer. All connections are assumed between layer 
{it:i} and layer {it:i+1}. 
{p_end}
	

{marker options}{...}
{title:Options}

{phang}
{opt outputs(#)} specifies the number of outputs in the network.

{phang}
{opt inputs(#)} specifies the number of inputs in the network.

{phang}
{opt hlayers(#)} number of hidden layers in the network; default is {cmd:hlayers(1)}

{phang}
{opt hnodes(numlist)} number of nodes per hidden layer; default is {cmd:hnodes(1)}. If you only specify one number, 
and {cmd:hlayers()} is > 1, then {cmd:neuralnet} assumes that number of nodes in all hidden layers. 

{phang}
{opt colors(list)} override the default color list, which is {cmd:colors(navy forest_green maroon purple khaki)}

{phang}
{opt inlabels(list)} override the default input labels, which are {cmd:x#}

{phang}
{opt outlabels(list)} override the default output labels, which are {cmd:y#}

{phang}
{opt hlabel(string)} override the default hidden layer label, which is {cmd:Hidden layer}, with the number added on the end.


{marker examples}{...}
{title:Example: Three-layer neural network}

{pstd}{cmd:nnplot, outputs(1) inputs(5) hlayers(1) hnodes(3)}{p_end}


{title:Author}

{pstd}Michael J. Crowther{p_end}
{pstd}Red Door Analytics AB{p_end}
{pstd}Stockholm, Sweden{p_end}
{pstd}E-mail: {browse "mailto:michael@reddooranalytics.se":michael@reddooranalytics.se}{p_end}

{phang}
Please report any errors you may find.{p_end}
