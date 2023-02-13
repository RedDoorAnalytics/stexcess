//local drive Z:/
local drive /Users/Michael/Documents
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

local drive /Users/Michael/Documents/reddooranalytics/products/stexcess

cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/ado"
adopath ++ "`drive'/help"

pr drop _all
mata mata clear

clear 
set seed 725
set obs 5000
gen id1 = _n
gen cancer = runiform()>0.5
gen age = rnormal(0,5)

survsim stime died , 	hazard(exp(log(0.01) :+0.001 :* age :+          ///
                                log({t})) :+ cancer :* exp(log(0.02) :+ ///
                                0.002 :* age :+                         ///
                                0.1 :* log({t})))                       ///
                        maxt(10) // nodes(150)

stset stime, f(died)

// - sync indicator in merlin struct                                    DONE
// - replay sync                                                        DONE
// - e(predictnotok) ...                                                DONE
// - equation names in output                                           DONE
// - predict help file                                                  DONE
// - multiple timescale sync & example                                  DONE
// - 
      
stexcess (age, df(1) noorthog)          ///
         (age, df(1) noorthog)          ///
         ,                              ///
         indicator(cancer)
         
stexcess
         
// predict test1, mu
         
predict s1, survival
predict s2, survival at(cancer 0) //standardise
predict s3, survival at(cancer 1 age 45)
