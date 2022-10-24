
local drive /Users/Michael/Documents/reddooranalytics/products/stexcess

cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/ado"

pr drop _all
mata mata clear

clear 
set seed 725
set obs 1000
gen id1 = _n
gen cancer = runiform()>0.5
gen age = rnormal(0,5)

survsim stime died , 	hazard(0.01 :* exp(0.001 :* age :+ log({t})) :+ exp(-1 :* cancer :+ 0.002 :* cancer :* age))  ///
                        maxt(10)

// mata:
// real matrix modexcess(transmorphic gml, real matrix t)
// {
// 	haz_expect = exp(merlin_util_xzb(gml,t))
//         haz_excess = exp(merlin_util_xzb_mod(gml,2))
// 	return(log(haz_expect :+ haz_excess))
// }
// end

stset stime, f(died)

stexcess        (age, df(3))            ///
                (age, df(3))            ///
                , indicator(cancer)


// merlin  (stime age fp(stime, pow(0)), timevar(stime) family(user, failure(died) loghf(modexcess)))  ///
//         (cancer cancer#age, nocons family(null))

