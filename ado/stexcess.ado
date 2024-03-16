*! version 1.1.0  16mar2024

/*
History
v1.1.0
- new evaluator calculatign full log likelihood
- analytics scores now used
*/


program define stexcess, eclass sortpreserve properties(st)
        version 17
        local copy0 `0'
	syntax [anything] [if] [in] , [*]
	if replay() {
                if (`"`e(cmd2)'"'!="stexcess") error 301
                Replay `copy0'
                exit
	}
	else Estimate `copy0'
	ereturn local cmd2 stexcess
	ereturn local cmdline2 stexcess `copy0'
end

program Estimate, eclass
        version 17
        
        st_is 2 analysis
        
        cap which merlin
        if _rc>0 {
                di as txt "The {bf:merlin} package is required; type {bf:ssc install merlin}"
                exit 198
        }
        
        _parse expand cmd glob : 0
        if `cmd_n'>2 | `cmd_n'<2 {
                di as error "control and excess models required"
                exit 198
        }
        local if `"`glob_if'"'
        local in `"`glob_in'"'
        local opts `"`glob_op'"'
        
        // parse global opts
        
        local 0 , `glob_op'
        syntax [if] [in] ,      INDicator(varname)              ///
                                [                               ///
                                        CHINTPoints(passthru)   ///
                                        FROM(passthru)          ///
                                        DEBUG                   ///
					EVALtype(passthru)	///
                                        *                       ///
                                ]

        if "`debug'"!="" {
                local noisily noisily
        }
	if "`evaltype'"=="" {
		local evaltype evaltype(gf1)
	}
                                
        global ind `indicator'          //!! fix - add to merlin struct
        
        //parse control
        _parsemodel `cmd_1'
        local clp1 `"`s(clp)'"'
        local nocons1 "`s(noconstant)'"
        local sv_clp1 `"`s(sv_clp)'"'
        local sv_opts1   `"`s(sv_opts)'"'
        
        //parse excess
        _parsemodel `cmd_2'
        local clp2 `"`s(clp)'"'
        local nocons2 "`s(noconstant)'"
        
        // sample
        if "`if'"!="" {
                local if if `if' & _st==1 & !missing(`indicator')
        }
        else {
                local if if _st==1 & !missing(`indicator')
        }
        local samp `if' `in'
        
        //check indicator
        qui count if `indicator'==0 
        if `r(N)'==0 {
                di as error "No control observations"
                exit 198
        }
        qui count if `indicator'==1 
        if `r(N)'==0 {
                di as error "No excess observations"
                exit 198
        }
        
        //left truncation
        qui su _t0 `samp'
        local hasltrunc = `r(max)' > 0
        if `hasltrunc' {
                local ltrunc ltruncated(_t0)
        }
        
        if "`from'"=="" {
                di ""
                di as txt "Obtaining starting values:"
                //starting values
                tempname init init1 init2
                cap `noisily' {
                        stmerlin `clp1' `samp' & `indicator'==0      	///
                                , dist(exp)                             ///
                                  nogen                                 ///
                                  `nocons1'                             //
                }
                if _rc>0 {
                        di as error "Error in obtaining starting values"
                        exit `_rc'
                }
                mat `init1' = e(b)
                          
                tempvar exprate
                qui predict `exprate' , hazard outcome(1)
                       
                cap `noisily' {
                        stmerlin `clp2' `samp' & `indicator'==1         ///
                                , dist(exp)                             ///
                                  nogen                                 ///
                                  bhazard(`exprate')                    ///
                                  showmerlin                            ///
                                  `nocons2'                             //
                }
                if _rc>0 {
                        di as error "Error in obtaining starting values"
                        exit `_rc'
                }
                mat `init2' = e(b)
                mat `init' = `init1',`init2'
                local from from(`init')
        }
        
        // Fit
        merlin 	(_t     `clp1'                                  ///
                        ,                                       ///
                        family(user,                            ///
                                llf(merlin_stexcess_logl)     	///
                                failure(_d)                     ///
                                `ltrunc')                       ///
                        timevar(_t)                             ///
                        `nocons1'                               ///
                )                                               ///
                (       `clp2'                                  ///
                        ,                                       ///
                        family(null, reffailure(1))             ///
                        timevar(_t)                             ///
                        `nocons2'                               ///
                )                                               ///
                `samp'                                          ///
                ,                                               ///
                mordred nogen                                   ///
                indicator(`indicator')                          ///
                modellabels("reference excess")                 ///
                `from'                                          ///
                `chintpoints'                                   ///
		`evaltype'					///
                `options'					//
                
        ereturn local predictnotok mu eta mudifference ///
                        etadifference etaratio muratio
        

end

program _parsemodel, sclass
        syntax [anything] ,                             ///
                        [                               ///
                                DF(passthru)	        ///
                                KNOTS(passthru)	        ///
                                NOORTHog	        ///
                                TIME                    ///
							///
				OFFset(passthru)	///
				MOFFset(passthru)	///
                                                        ///
                                TVC(varlist)	        ///
                                DFTvc(numlist)	        ///
                                TVCTIME		        ///
                                                        ///
                                TIME2(string)	        ///
                                TIME3(string)	        ///
                                TIME4(string)	        ///
                                TIME5(string)	        ///
                                NOCONStant              ///
                        ]                               //
        
        local varlist `"`anything'"'
        local nocons `noconstant'
        
        sreturn local sv_clp `"`varlist'"'
        local sv_opts `"`df' `knots' `noorthog' `time' `offset' `moffset'"'
        local sv_opts `"`sv_opts' tvc(`tvc') dftvc(`dftvc') `tvctime'"'
        forvalues i=2/5 {
                local sv_opts `"`sv_opts' time`i'(`time`i'')"'
        }
        sreturn local sv_opts `"`sv_opts'"'
        
        // baseline
        if "`noorthog'"==""     local orth orthog
        if "`time'"==""         local log log
        
        local rcsbase rcs(_t, `df' `knots' `orth' `log' ///
			`offset' `moffset' event)
        
        // tvcs
        if "`tvc'"!="" {
                if "`dftvc'"=="" {
                        di as error "dftvc() required"
                        exit 198
                }
                
                local tvclog 
                if "`tvctime'"=="" {
                        local tvclog log
                }
                
                local tvcorthog
                if "`noorthog'"=="" {
                        local tvcorthog orthog
                }
                
                local Ntvc : word count `tvc'
                local Ntvcdf : word count `dftvc'
                if `Ntvcdf''>1 & `Ntvcdf'!=`Ntvc' {
                        di as error "Number of dftvc() elements do not match tvc()"
                        exit 198
                }
                
                if `Ntvcdf'==1 {
                        foreach var in `tvc' {
                                local tvcs `tvcs' `var'#rcs(_t, df(`dftvc') event `tvclog' `tvcorthog')
                        }
                }
                else {
                        local ind = 1
                        foreach var in `tvc' {
                                local dftvcv : word `ind' of `dftvc'
                                local tvcs `tvcs' `var'#rcs(_t, df(`dftvcv') event `tvclog' `tvcorthog')
                                local ind = `ind' + 1
                        }
                }
        }
        
        // additional timescales
        forvalues j=2/5 {
                if "`time`j''"!="" {
                        _merlin_parse_timescale , `time`j''
                        local multitimes `multitimes' `s(multitimes)'
                }
        }
        
        sreturn local clp `varlist' `tvcs' `multitimes' `rcsbase'
        sreturn local noconstant `noconstant'
end

program _merlin_parse_timescale, sclass
        version 17
        syntax ,[				        ///
                        OFFset(passthru)		///
                        MOFFset(passthru)		///
                        DF(numlist max=1 int >0)	///
                        KNOTS(numlist min=2)		///
                        NOORTHog			///
                        TIME				///
                                                        ///
                        TVC(varlist)			///
                        DFTvc(numlist)			///
                        TVCTime				///
                ]                                       //
        
        
        if "`noorthog'"==""     local orthog orthog
        if "`time'"==""         local log log
        if "`tvctime'"==""      local tvctime log
        
        if "`df'"=="" {
                local knotsopt knots(`knots')
                local Ndf = `: word count `knots'' - 1 
        }
        else {
                local knotsopt df(`df')
                local Ndf = `df'
                local event event
        }
        
        local multitimes rcs(_t, `offset' `moffset' `event' `knotsopt' `orthog' `log')

        if "`tvc'"!="" {
                
                local Ntvc : word count `tvc'
                local Ntvcdf : word count `dftvc'
                
                if `Ntvcdf'>1 & `Ntvcdf'!=`Ntvc' {
                        di as error "Number of dftvc() elements do not match tvc()"
                        exit 198
                }
                
                if `Ntvcdf'==1 {
                        foreach var in `tvc' {
                                local multitimes `multitimes' `var'#rcs(_t, `offset' `moffset' `event' df(`dftvc') `orthog' `tvctime')
                        }
                }
                else {
                        local ind = 1
                        foreach var in `tvc' {
                                local dftvci : word `ind' of `dftvc'
                                local multitimes `multitimes' `var'#rcs(_t, `offset' `moffset' `event' df(`dftvci') `orthog' `tvctime')
                                local ind = `ind' + 1
                        }
                }

        }
        
        sreturn local multitimes `multitimes'
end

program Replay
	merlin `0'
end
