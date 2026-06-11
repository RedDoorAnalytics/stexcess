*! stexcess v2.0.0  (Mata core)  -- modelled excess hazard models
*! Syntax follows the original (merlin-based) stexcess v1.1.1:
*!   stexcess (refspec)(excspec) [if][in], indicator(varname) [options]
*! where each model spec is  [varlist][, model_options]  with per-component
*! df()/knots()/noorthog/time/offset()/moffset()/tvc()/dftvc()/tvctime/
*! time2()-time5()/noconstant. Estimation and postestimation run in Mata
*! (library lstexcess.mlib, source ado/stexcess_mata.do): designs built once,
*! maximised with optimize() -- the Newton-Raphson engine beneath ml -- using
*! an exact analytic d2 evaluator. No dependencies outside Stata.

program stexcess, eclass properties(st)
    version 19.5

    if replay() {
        if "`e(cmd)'" != "stexcess" error 301
        Display `0'
        exit
    }
    Estimate `0'
    ereturn local cmdline `"stexcess `0'"'
end

program Estimate, eclass
    st_is 2 analysis                       // require stset

    // ---- split (refspec)(excspec) [if][in], options ----
    _parse expand stxcmd stxglob : 0
    if `stxcmd_n' != 2 {
        di as err "reference and excess model specifications required, as in"
        di as err "    stexcess (refspec)(excspec), indicator(varname) ..."
        exit 198
    }

    local 0 `"`stxglob_if' `stxglob_in', `stxglob_op'"'
    syntax [if] [in], INDicator(varname numeric) ///
        [ CHINTpoints(integer 30) TWOstage Level(cilevel) FROM(name) ///
          ITERate(numlist max=1 integer >=0) TOLerance(numlist max=1 >0) ///
          LTOLerance(numlist max=1 >0) NRTOLerance(numlist max=1 >0) ///
          EFORM noLOG DEBUG EVALtype(string) ]
    // debug and evaltype() are accepted for v1 compatibility; the Mata core
    // always uses its exact analytic d2 evaluator, so they are no-ops

    if `chintpoints' < 1 {
        di as err "chintpoints() must be >= 1"
        exit 198
    }
    if "`from'" != "" {
        confirm matrix `from'
    }

    // ---- per-component model specs ----
    foreach c in ref exc {
        if "`c'" == "ref" {
            local spec `"`stxcmd_1'"'
            local eqn reference
        }
        else {
            local spec `"`stxcmd_2'"'
            local eqn excess
        }
        _parsemodel `eqn' `spec'
        local `c'vars `s(vars)'
        local `c'cons `s(cons)'
        local `c'tslist `s(tslist)'
        foreach j of local `c'tslist {
            foreach o in df knots log orthog off moff tvc dftvc tvclog {
                local `c'_ts`j'_`o' `"`s(ts`j'_`o')'"'
            }
        }
    }

    // ---- sample / requirements ----
    if `"`: char _dta[st_w]'"' != "" {
        di as err "stexcess does not support weights " ///
            `"(data are stset with `: char _dta[st_w]')"'
        exit 101
    }
    marksample touse
    markout `touse' _t _d `indicator' `refvars' `excvars'
    foreach c in ref exc {
        foreach j of local `c'tslist {
            markout `touse' ``c'_ts`j'_off' ``c'_ts`j'_moff' ``c'_ts`j'_tvc'
        }
    }
    qui count if `touse'
    if r(N) == 0 error 2000
    local nobs = r(N)

    qui count if !inlist(`indicator', 0, 1) & `touse'
    if r(N) {
        di as err "indicator(`indicator') must be coded 0 (reference) " ///
            "or 1 (excess) on the estimation sample"
        exit 450
    }
    qui count if `indicator' == 0 & `touse'
    local nref = r(N)
    local nexc = `nobs' - `nref'
    qui count if _d == 1 & `touse'
    local nfail = r(N)
    qui count if `indicator' == 1 & _d == 1 & `touse'
    if r(N) == 0 {
        di as err "no events among excess (indicator = 1) records"
        exit 2000
    }
    if "`twostage'" != "" & `nref' == 0 {
        di as err "twostage requires reference (indicator = 0) records"
        exit 2000
    }

    // drop collinear covariates (constant covariates are collinear with the
    // spline intercept); tvc variables enter the design through their spline
    // interactions, so a constant tvc variable cannot just be dropped
    foreach c in ref exc {
        local eqn = cond("`c'" == "ref", "reference", "excess")
        if "``c'vars'" != "" & "``c'cons'" == "1" {
            qui _rmcoll ``c'vars' if `touse', forcedrop
            local keep `r(varlist)'
            if "`keep'" != "``c'vars'" {
                local dropped : list `c'vars - keep
                di as txt "note: `dropped' omitted from the `eqn' " ///
                    "equation because of collinearity"
                local `c'vars `keep'
            }
        }
        foreach j of local `c'tslist {
            foreach v of local `c'_ts`j'_tvc {
                qui su `v' if `touse'
                if r(sd) == 0 | r(sd) >= . {
                    di as err "tvc variable `v' is constant " ///
                        "in the `eqn' equation"
                    exit 459
                }
            }
        }
    }

    // ---- hand parameters to the Mata fit driver via macros/scalars ----
    local _stx_t       _t
    local _stx_t0      _t0
    local _stx_d       _d
    local _stx_ind     `indicator'
    local _stx_touse   `touse'
    local _stx_from    `from'
    local _stx_nolog   `log'
    local _stx_iterate `iterate'
    local _stx_ptol    `tolerance'
    local _stx_vtol    `ltolerance'
    local _stx_nrtol   `nrtolerance'
    foreach c in ref exc {
        local _stx_`c'vars   ``c'vars'
        local _stx_`c'_cons  ``c'cons'
        local _stx_`c'_tslist ``c'tslist'
        foreach j of local `c'tslist {
            foreach o in df knots log orthog off moff tvc dftvc tvclog {
                local _stx_`c'_ts`j'_`o' ``c'_ts`j'_`o''
            }
        }
    }
    local _stx_nnodes  `chintpoints'
    local _stx_twostage = ("`twostage'" != "")
    tempname b V
    local _stx_bmat `b'
    local _stx_Vmat `V'

    st_show
    mata: _stx_fit()

    // ---- assemble e() and post ----
    local names `_stx_names'
    matrix colnames `b' = `names'
    matrix colnames `V' = `names'
    matrix rownames `V' = `names'

    ereturn post `b' `V', esample(`touse') obs(`nobs') depname(_t)
    ereturn scalar ll        = `_stx_ll'
    ereturn scalar k         = `_stx_k'
    ereturn scalar df_m     = `: word count `refvars'' + ///
                              `: word count `excvars''
    ereturn scalar converged = `_stx_conv'
    ereturn scalar ic        = `_stx_iter'
    ereturn scalar iterations = `_stx_iter'
    ereturn scalar chintpoints = `chintpoints'
    mata: st_local("vrank", strofreal(rank(st_matrix("e(V)"))))
    ereturn scalar rank      = `vrank'
    foreach c in ref exc {
        local dfb ``c'_ts1_df'
        if "`dfb'" == "" {
            local dfb = `: word count ``c'_ts1_knots'' - 1
        }
        ereturn scalar df`c' = `dfb'
    }
    ereturn scalar N_ref     = `nref'
    ereturn scalar N_exc     = `nexc'
    ereturn scalar N_fail    = `nfail'
    ereturn local  knotsref  "`_stx_kref'"
    ereturn local  knotsexc  "`_stx_kexc'"
    foreach c in ref exc {
        foreach j of local `c'tslist {
            if `j' > 1 {
                ereturn local knots`c'_t`j' "`_stx_k`c'_t`j''"
            }
        }
    }
    ereturn local  atvars    "`_stx_atvars'"
    ereturn local  method    = cond("`twostage'" != "", "twostage", "joint")
    if "`twostage'" != "" {
        ereturn local vce     "robust"
        ereturn local vcetype "Robust"
    }
    else ereturn local vce "oim"
    ereturn local  title     "Modelled excess hazard model"
    local mnotok Hazard CHazard LOGCHazard SURVival CIF RMST TIMELost ///
        NETSurvival EXCesshazard RMSTNet HDIFFerence SDIFFerence ///
        CIFDIFFerence RMSTDIFFerence HRatio SRatio CIFRatio RMSTRatio ///
        STANDardise
    ereturn local marginsnotok "`mnotok'"
    ereturn local  indicator "`indicator'"
    ereturn local  refvars   "`refvars'"
    ereturn local  excvars   "`excvars'"
    ereturn local  tvc       "`exc_ts1_tvc'"
    ereturn local  tvcref    "`ref_ts1_tvc'"
    ereturn local  predict   "stexcess_p"
    ereturn local  cmd       "stexcess"


    Display, level(`level') `eform'
end

program Display
    syntax [, Level(cilevel) EFORM ]
    local eopt = cond("`eform'" != "", `"eform("exp(b)")"', "")
    di ""
    di as txt "`e(title)'" ///
        _col(49) as txt "Number of obs     =" as res %10.0fc e(N)
    di _col(49) as txt "No. of failures   =" as res %10.0fc e(N_fail)
    di _col(49) as txt "Reference records =" as res %10.0fc e(N_ref)
    di _col(49) as txt "Excess records    =" as res %10.0fc e(N_exc)
    di as txt "Log likelihood = " as res %10.4f e(ll)
    if "`e(method)'" == "twostage" {
        di as txt "Two-stage estimation: reference fitted to controls " ///
            "only; stacked sandwich variance."
    }
    if e(converged) == 0 di as err "convergence not achieved"
    ereturn display, level(`level') `eopt'
    di as txt "Equations: {bf:ref} = reference (control) hazard, " ///
        "{bf:exc} = excess hazard"
end

// parse one model spec: [varlist][, model_options]; results in s().
// Timescale 1 is the baseline; time2()-time5() add timescales as splines in
// (t + offset - moffset) on the log (default) or natural (time) scale.
program _parsemodel, sclass
    sreturn clear
    gettoken eqn 0 : 0                     // "reference" | "excess"

    syntax [varlist(numeric default=none)] , ///
        [ DF(numlist max=1 integer >0) KNOTS(numlist ascending min=2) ///
          NOORTHog TIME OFFset(varname numeric) MOFFset(varname numeric) ///
          TVC(varlist numeric) DFTvc(numlist integer >0) TVCTIME ///
          TIME2(string) TIME3(string) TIME4(string) TIME5(string) ///
          NOCONStant ]

    if "`df'" == "" & "`knots'" == "" local df 3
    if "`df'" != "" & "`knots'" != "" {
        di as err "`eqn' model: only one of df() and knots() is allowed"
        exit 198
    }
    _expand_dftvc `eqn' "`tvc'" "`dftvc'"
    local dftvc `s(dftvc)'

    sreturn local vars       "`varlist'"
    sreturn local cons       = cond("`noconstant'" == "", "1", "0")
    sreturn local ts1_df     "`df'"
    sreturn local ts1_knots  "`knots'"
    sreturn local ts1_log    = cond("`time'" == "", "1", "0")
    sreturn local ts1_orthog = cond("`noorthog'" == "", "1", "0")
    sreturn local ts1_off    "`offset'"
    sreturn local ts1_moff   "`moffset'"
    sreturn local ts1_tvc    "`tvc'"
    sreturn local ts1_dftvc  "`dftvc'"
    sreturn local ts1_tvclog = cond("`tvctime'" == "", "1", "0")

    local tslist 1
    forvalues j = 2/5 {
        if `"`time`j''"' != "" {
            local tslist `tslist' `j'
            local 0 `", `time`j''"'
            syntax , [ OFFset(varname numeric) MOFFset(varname numeric) ///
                DF(numlist max=1 integer >0) KNOTS(numlist ascending min=2) ///
                NOORTHog TIME TVC(varlist numeric) ///
                DFTvc(numlist integer >0) TVCTime ]
            if "`df'" == "" & "`knots'" == "" {
                di as err "`eqn' model: time`j'() requires df() or knots()"
                exit 198
            }
            if "`df'" != "" & "`knots'" != "" {
                di as err "`eqn' model: time`j'(): only one of df() and " ///
                    "knots() is allowed"
                exit 198
            }
            _expand_dftvc "`eqn' time`j'()" "`tvc'" "`dftvc'"
            sreturn local ts`j'_df     "`df'"
            sreturn local ts`j'_knots  "`knots'"
            sreturn local ts`j'_log    = cond("`time'" == "", "1", "0")
            sreturn local ts`j'_orthog = cond("`noorthog'" == "", "1", "0")
            sreturn local ts`j'_off    "`offset'"
            sreturn local ts`j'_moff   "`moffset'"
            sreturn local ts`j'_tvc    "`tvc'"
            sreturn local ts`j'_dftvc  "`s(dftvc)'"
            sreturn local ts`j'_tvclog = cond("`tvctime'" == "", "1", "0")
        }
    }
    sreturn local tslist "`tslist'"
end

// dftvc() is required with tvc(); one value applies to all tvc variables,
// otherwise one value per variable (v1 behaviour)
program _expand_dftvc, sclass
    args what tvc dftvc
    if "`tvc'" == "" {
        sreturn local dftvc ""
        exit
    }
    if "`dftvc'" == "" {
        di as err "`what' model: dftvc() required when tvc() is specified"
        exit 198
    }
    local ntvc : word count `tvc'
    local ndf  : word count `dftvc'
    if `ndf' == 1 {
        local out
        forvalues i = 1/`ntvc' {
            local out `out' `dftvc'
        }
        sreturn local dftvc "`out'"
        exit
    }
    if `ndf' != `ntvc' {
        di as err "`what' model: number of dftvc() elements does not " ///
            "match tvc()"
        exit 198
    }
    sreturn local dftvc "`dftvc'"
end
