local drive /Users/michael/My Drive/software
cd "`drive'/merlin"
adopath ++ "`drive'/merlin"
adopath ++ "`drive'/merlin/merlin"
adopath ++ "`drive'/merlin/stmerlin"
clear all

tr:do ./build/buildmlib.do
mata mata clear

local drive /Users/michael/My Drive/software/stexcess

cd "`drive'"
adopath ++ "`drive'"
adopath ++ "`drive'/ado"
adopath ++ "`drive'/help"

pr drop _all
mata mata clear

clear 
set seed 725
set obs 10000
gen id1 = _n
gen cancer = runiform()>0.5
gen age = rnormal(0,10)
gen tes0 = 0

survsim stime died , 	hazard(exp(log(0.01) :+0.01 :* age :+          ///
                                log({t})) :+ cancer :* exp(log(0.02) :+ ///
                                0.002 :* age :+                         ///
                                0.1 :* log({t})))                       ///
                        maxt(10) // nodes(150)

stset stime, f(died)

timer clear
timer on 1
merlin (_t age rcs(_t, df(1) log event) , 	///
	family(user, loghf(merlin_stexcess_logh) ///
		failure(_d)) timevar(_t)) 		///
       ( rcs(_t, df(1)   log   event) , family(null, reffailure(1)) ///
	timevar(_t)) if _st==1 & !missing(cancer) , 	///
		mordred nogen indicator(cancer)
timer off 1
timer on 2
merlin (_t age rcs(_t, df(1) log event) , 	///
	family(user, llf(merlin_stexcess_logl) ///
		failure(_d)) timevar(_t)) 		///
       ( rcs(_t, df(1)   log   event) , family(null, reffailure(1)) ///
	timevar(_t)) if _st==1 & !missing(cancer) , 	///
		mordred nogen indicator(cancer) evaltype(gf0) 
timer off 2	
timer on 3
merlin (_t age rcs(_t, df(1) log event) , 	///
	family(user, llf(merlin_stexcess_logl) ///
		failure(_d)) timevar(_t)) 		///
       ( rcs(_t, df(1)   log   event) , family(null, reffailure(1)) ///
	timevar(_t)) if _st==1 & !missing(cancer) , 	///
		mordred nogen indicator(cancer) evaltype(gf1) 
timer off 3
timer list

/*
. timer list
   1:      9.74 /        1 =       9.7380
   2:      8.62 /        1 =       8.6250
   3:      7.13 /        1 =       7.1330
*/

exit 
// - sync indicator in merlin struct                                    DONE
// - replay sync                                                        DONE
// - e(predictnotok) ...                                                DONE
// - equation names in output                                           DONE
// - predict help file                                                  DONE
// - multiple timescale sync & example                                  DONE
// - 
      
stexcess (age, df(1) noorthog)       ///
         (, df(1) noorthog)          ///
         ,                              ///
         indicator(cancer)
         
stexcess
         
// predict test1, mu
cap drop time
range time 0 5 100
predict s1, survival at(cancer 0) timevar(time)
predict s2, survival at(cancer 1) timevar(time)


scatter s1 s2 time
