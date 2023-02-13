
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

survsim stime died , 	hazard(exp(log(0.01) :+0.01 :* (age :+ {t}) :+   ///
                                log({t})) :+ cancer :* exp(log(0.02) :+  ///
                                0.02 :* age :+                           ///
                                0.1 :* log({t})))                        ///
                        maxt(10) //nodes(150)

stset stime, f(died)
  
stexcess (, df(1) noorthog time2(df(1) offset(age) time noorthog))      ///
         (age, df(1) noorthog)    ///
         ,                        ///
         indicator(cancer) 
         
predict s1, survival
predict s2, relsurv at(cancer 0 age 45)

predict s3, sratio at1(cancer 1 age 45) at2(cancer 0 age 45)
