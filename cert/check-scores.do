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
                        maxt(10)

stset stime, f(died)

stexcess (age, df(1) noorthog)  ///
         (age, df(1) noorthog)  ///
         ,                      ///
         indicator(cancer) evaltype(gf0)
mat b1 = e(b)
mat v1 = e(V)

stexcess (age, df(1) noorthog)  ///
         (age, df(1) noorthog)  ///
         ,                      ///
         indicator(cancer) evaltype(gf1)
mat b2 = e(b)
mat v2 = e(V)

assert mreldif( b1 , b2 ) < 1E-6
assert mreldif( v1 , v2 ) < 1E-7
