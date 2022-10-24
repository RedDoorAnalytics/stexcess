*! version 1.0.0  ?????2022

program define stexcess, eclass sortpreserve properties(st)
        version 17
        local copy0 `0'
	syntax [anything] [if] [in] , [*]
	if replay() & "`options'"=="" {
		if (`"`e(cmd2)'"'!="stmerlin") error 301
		Replay `copy0'
		exit
	}
	else Estimate `copy0'
	ereturn local cmd2 stexcess
	ereturn local cmdline2 stexcess `copy0'


end

program Estimate, eclass
        version 17
        
        local GLOB      INDicator(varname)      ///
                                                //
        
        
        _parse expand cmd glob : 0
        if `cmd_n'>2 | `cmd_n'<2 {
                di "control and excess models required"
                exit 1986
        }
        local if `"`glob_if'"'
        local in `"`glob_in'"'
        local opts `"`glob_op'"'
        
        forvalues i=1/2 {
                
                local 0 `"`cmd_`i''"'
                syntax anything ,                       ///
                                [                       ///
                                        DF(passthru)	///
                                        KNOTS(passthru)	///
                                        NOORTHog	///
                                        TIME            ///
                                                        ///
                                        TVC(passthru)	///
                                        DFTvc(passthru)	///
                                        TVCTIME		///
                                                        ///
                                        TIME2(string)	///
                                        TIME3(string)	///
                                        TIME4(string)	///
                                        TIME5(string)	///
                                        NOCONStant      ///
                                ]                       //
                
                local varlist_`i' `anything'
                local nocons`i' `noconstant'
                
                // baseline
                if "`noorthog'"==""     local orth orthog
                if "`time'"==""         local log log
                
		local rcsbase`i' rcs(_t, `df' `knots' `orth' `log' event)
                
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
			
			local Ntvc`i' : word count `tvc'
			local Ntvcdf`i' : word count `dftvc'
			if `Ntvcdf`i'''>1 & `Ntvcdf`i''!=`Ntvc`i'' {
				di as error "Number of dftvc() elements do not match tvc()"
				exit 198
			}
			
			if `Ntvcdf`i''==1 {
				foreach var in `tvc' {
					local tvcs`i' `tvcs`i'' `var'#rcs(_t, df(`dftvc') event `tvclog' `tvcorthog')
				}
			}
			else {
				local ind = 1
				foreach var in `tvc' {
					local dftvcv : word `ind' of `dftvc'
					local tvcs`i' `tvcs`i'' `var'#rcs(_t, df(`dftvcv') event `tvclog' `tvcorthog')
					local ind = `ind' + 1
				}
			}
                }
                
                // additional timescales
                forvalues j=2/5 {
                        _merlin_parse_timescale , `time`j''
                        local multitimes`i' `multitimes`i'' `s(multitimes)'
                }
                
        }
        
        // parse global opts
        
        local 0 , `glob_op'
        syntax , `GLOB'
        
        // sample
        tempvar touse
        mark `touse' `if' `in'
        
        //left truncation
        qui su _t0 if `touse'==1 & _st==1
        local hasltrunc = `r(max)' > 0
        if `hasltrunc' {
                local ltrunc ltruncated(_t0)
        }
        
        // Fit
        
        merlin 	(_t                                             ///
                        `varlist1'                              ///
                        `tvcs1'                                 ///
                        `multitimes1'                           ///
                        `rcsbase1'                              ///
                        ,                                       ///
                        family(user,    loghf(_stexcess_logh)   ///
                                        failure(_d)             ///
                                        `ltrunc'                ///
                                        )                       ///
                        timevar(_t)                             ///
                        `nocons1'                               ///
                )                                               ///
                (       `varlist2'                              ///
                        `tvcs2'                                 ///
                        `multitimes2'                           ///
                        `rcsbase2'                              ///
                        ,                                       ///
                        family(null)                            ///
                        timevar(_t)                             ///
                        noconstant                              ///
                )                                               ///
                if _st==1                                       ///
                `from' 						///
                `debug'						///
                `level' 					///
                `nolog'						///
                `evaltype'					///
                `mlopts'					//

        
        

end

program Replay
	merlin `0'
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
                local event event
        }
        else {
                local knotsopt df(`df')
                local Ndf = `df'
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
        
        sreturn local timescale
end

version 17

mata:
real matrix _stexcess_logh(transmorphic gml, real matrix t)
{
	haz_expect = exp(merlin_util_xzb(gml,t))
        haz_excess = exp(merlin_util_xzb_mod(gml,2))
	return(log(haz_expect :+ haz_excess))
}
end
