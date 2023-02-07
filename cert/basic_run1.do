
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
assert `"`e(predictnotok)'"'    == `"mu"'
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
assert `"`e(loghfunction1)'"'   == `"merlin_stexcess_logh"'
assert `"`e(timevar1)'"'        == `"_t"'
assert `"`e(failure1)'"'        == `"_d"'
assert `"`e(response1)'"'       == `"_t _d"'
assert `"`e(family1)'"'         == `"user"'
assert `"`e(allvars)'"'         == `"_d _t age"'
assert `"`e(title)'"'           == `"Excess hazard model"'
assert `"`e(cmd)'"'             == `"merlin"'
assert `"`e(hasopts)'"'         == `"1"'
assert `"`e(from)'"'            == `"1"'
assert `"`e(predict)'"'         == `"merlin_p"'
assert `"`e(opt)'"'             == `"moptimize"'
assert `"`e(vce)'"'             == `"oim"'
assert `"`e(user)'"'            == `"merlin_gf()"'
assert `"`e(crittype)'"'        == `"log likelihood"'
assert `"`e(ml_method)'"'       == `"gf0"'
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
assert reldif( e(ll)          , -8556.741072163895) <  1E-8
assert         e(k_eq_model) == 0
assert         e(Nmodels)    == 2
assert         e(Nlevels)    == 1

qui {
mat T_b = J(1,6,0)
mat T_b[1,1] =  .0016722318190681
mat T_b[1,2] =   .968330123611433
mat T_b[1,3] = -4.506771190683557
mat T_b[1,4] = -.0142985163622496
mat T_b[1,5] =  .1331617357998256
mat T_b[1,6] =  -4.03542056116722
}
matrix C_b = e(b)
assert mreldif( C_b , T_b ) < 1e-06
mat drop C_b T_b

qui {
mat T_V = J(6,6,0)
mat T_V[1,1] =   .000032287578259
mat T_V[1,2] =  1.30996921976e-06
mat T_V[1,3] = -2.13883269231e-07
mat T_V[1,4] = -.0000560702661379
mat T_V[1,5] = -.0000113391474606
mat T_V[1,6] = -.0000203585945067
mat T_V[2,1] =  1.30996921976e-06
mat T_V[2,2] =  .0026400637408468
mat T_V[2,3] = -.0047112889405032
mat T_V[2,4] =  6.34270286376e-06
mat T_V[2,5] = -.0009792889873146
mat T_V[2,6] =  .0021292105326087
mat T_V[3,1] = -2.13883269231e-07
mat T_V[3,2] = -.0047112889405032
mat T_V[3,3] =  .0093505363525647
mat T_V[3,4] = -.0000540604351974
mat T_V[3,5] =  .0009336040504441
mat T_V[3,6] = -.0047714714863026
mat T_V[4,1] = -.0000560702661379
mat T_V[4,2] =  6.34270286376e-06
mat T_V[4,3] = -.0000540604351974
mat T_V[4,4] =  .0003685107009205
mat T_V[4,5] =  .0000994094381325
mat T_V[4,6] =  .0001969542018765
mat T_V[5,1] = -.0000113391474606
mat T_V[5,2] = -.0009792889873146
mat T_V[5,3] =  .0009336040504441
mat T_V[5,4] =  .0000994094381325
mat T_V[5,5] =  .0062161267956487
mat T_V[5,6] = -.0043158360015108
mat T_V[6,1] = -.0000203585945067
mat T_V[6,2] =  .0021292105326087
mat T_V[6,3] = -.0047714714863026
mat T_V[6,4] =  .0001969542018765
mat T_V[6,5] = -.0043158360015108
mat T_V[6,6] =  .0133360896374532
}
matrix C_V = e(V)
assert mreldif( C_V , T_V ) < 1e-06
mat drop C_V T_V

qui {
mat T_gradient = J(1,6,0)
mat T_gradient[1,1] = -.0001385775937531
mat T_gradient[1,2] = -.0001519475629455
mat T_gradient[1,3] = -.0000606568581677
mat T_gradient[1,4] =  .0012616871374724
mat T_gradient[1,5] =  .0001386189127942
mat T_gradient[1,6] =  .0002201901997686
}
matrix C_gradient = e(gradient)
assert mreldif( C_gradient , T_gradient ) < 1e-06
mat drop C_gradient T_gradient
