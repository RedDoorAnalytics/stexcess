
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
                        maxt(10) 

stset stime, f(died)

stexcess (age, df(1) noorthog)  ///
         (age, df(1) noorthog)  ///
         ,                      ///
         indicator(cancer) 

assert `"`e(cmd2)'"'            == `"stexcess"'
assert `"`e(predictnotok)'"'    == `"mu eta mudifference etadifference etaratio muratio"'
assert `"`e(chintpoints)'"'     == `"30"'
assert `"`e(indicator)'"'       == `"cancer"'
assert `"`e(nap2)'"'            == `"0"'
assert `"`e(ndistap2)'"'        == `"0"'
assert `"`e(constant2)'"'       == `"1"'
assert `"`e(knots_2_2_1)'"'     == `"-3.8339042 2.302535"'
assert `"`e(cmplabels2)'"'      == `"age rcs()"'
assert `"`e(Nvars_2)'"'         == `"1 1"'
assert `"`e(timevar2)'"'        == `"_t"'
assert `"`e(family2)'"'         == `"null"'
assert `"`e(nap1)'"'            == `"0"'
assert `"`e(ndistap1)'"'        == `"0"'
assert `"`e(constant1)'"'       == `"1"'
assert `"`e(knots_1_2_1)'"'     == `"-3.8339042 2.302535"'
assert `"`e(cmplabels1)'"'      == `"age rcs()"'
assert `"`e(Nvars_1)'"'         == `"1 1"'
assert `"`e(llfunction1)'"'     == `"merlin_stexcess_logl"'
assert `"`e(timevar1)'"'        == `"_t"'
assert `"`e(failure1)'"'        == `"_d"'
assert `"`e(response1)'"'       == `"_t _d"'
assert `"`e(family1)'"'         == `"user"'
assert `"`e(allvars)'"'         == `"_d _t age"'
assert `"`e(title)'"'           == `"Excess hazard model"'
assert `"`e(modellabels)'"'     == `"reference excess"'
assert `"`e(cmd)'"'             == `"merlin"'
assert `"`e(hasopts)'"'         == `"1"'
assert `"`e(from)'"'            == `"1"'
assert `"`e(predict)'"'         == `"merlin_p"'
assert `"`e(opt)'"'             == `"moptimize"'
assert `"`e(vce)'"'             == `"oim"'
assert `"`e(user)'"'            == `"merlin_gf()"'
assert `"`e(crittype)'"'        == `"log likelihood"'
assert `"`e(ml_method)'"'       == `"gf1"'
assert `"`e(singularHmethod)'"' == `"m-marquardt"'
assert `"`e(technique)'"'       == `"nr"'
assert `"`e(which)'"'           == `"max"'
assert `"`e(properties)'"'      == `"b V"'

assert         e(rank)       == 6
assert         e(N)          == 5000
assert         e(k)          == 6
assert         e(k_eq)       == 6
assert         e(noconstant) == 0
assert         e(consonly)   == 1
assert         e(k_dv)       == 0
assert         e(converged)  == 1
assert         e(rc)         == 0
assert         e(k_autoCns)  == 0
assert reldif( e(ll)          , -8556.741072162085) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 2
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] =  .0016721647096895
mat T_b[1,2] =  .9683304471765269
mat T_b[1,3] = -4.506772750104814
mat T_b[1,4] = -.0142980256539864
mat T_b[1,5] =  .1331646204558134
mat T_b[1,6] = -4.035418277551483
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1E-8
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =    .00003228768158
mat T_V[1,2] =  1.31036919125e-06
mat T_V[1,3] = -2.14197116057e-07
mat T_V[1,4] = -.0000560703418006
mat T_V[1,5] =  -.000011340917284
mat T_V[1,6] =  -.000020357319092
mat T_V[2,1] =  1.31036919125e-06
mat T_V[2,2] =  .0026400515852382
mat T_V[2,3] = -.0047112647052271
mat T_V[2,4] =  6.34000059944e-06
mat T_V[2,5] = -.0009792651105189
mat T_V[2,6] =  .0021291565855443
mat T_V[3,1] = -2.14197116057e-07
mat T_V[3,2] = -.0047112647052271
mat T_V[3,3] =   .009350492320766
mat T_V[3,4] = -.0000540567668836
mat T_V[3,5] =  .0009335554664661
mat T_V[3,6] = -.0047713714798007
mat T_V[4,1] = -.0000560703418006
mat T_V[4,2] =  6.34000059944e-06
mat T_V[4,3] = -.0000540567668836
mat T_V[4,4] =  .0003685090294888
mat T_V[4,5] =  .0000994185574262
mat T_V[4,6] =  .0001969422350609
mat T_V[5,1] =  -.000011340917284
mat T_V[5,2] = -.0009792651105189
mat T_V[5,3] =  .0009335554664661
mat T_V[5,4] =  .0000994185574262
mat T_V[5,5] =  .0062160341531639
mat T_V[5,6] = -.0043156863486805
mat T_V[6,1] =  -.000020357319092
mat T_V[6,2] =  .0021291565855443
mat T_V[6,3] = -.0047713714798007
mat T_V[6,4] =  .0001969422350609
mat T_V[6,5] = -.0043156863486805
mat T_V[6,6] =  .0133357990737088
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1E-8
mat drop C_V T_V

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] =  -.000133916036714
mat T_gradient[1,2] = -.0000809319083315
mat T_gradient[1,3] = -.0000399607672723
mat T_gradient[1,4] =  .0000422478175737
mat T_gradient[1,5] =  .0000179053762835
mat T_gradient[1,6] = -.0000127702030061
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1E-8
mat drop C_gradient T_gradient
