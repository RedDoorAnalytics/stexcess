*! stexcess_p v2.0.0 -- postestimation predictions for stexcess (Mata core)
*! v1-style statistics and semantics: covariates and the excess indicator are
*! taken from the data row by row (so hazard = h_ref + ind*h_exc), overridable
*! with at()/zeros -- e.g. at(<indicator> 0) gives reference-only quantities.
*! Net (excess-only) statistics, standardisation and delta-method CIs are
*! computed in Mata (lstexcess.mlib) with analytic Jacobians.

program stexcess_p
    version 19.5

    syntax newvarname [if] [in], [ ///
        SURVival CIF Hazard CHazard LOGCHazard RMST TIMELost ///
        NETSurvival EXCesshazard RMSTNet ///
        HDIFFerence SDIFFerence CIFDIFFerence RMSTDIFFerence ///
        HRatio SRatio CIFRatio RMSTRatio ///
        STANDardise ///
        AT(string) AT1(string) AT2(string) ZEROs ///
        TImevar(varname numeric) ///
        CI Level(cilevel) ]

    if "`e(cmd)'" != "stexcess" error 301
    // NB: at() must be declared BEFORE at1()/at2() -- syntax binds a typed
    // option to the first declared match, and "at" abbreviates "at1".

    // ---- resolve the requested statistic (exactly one) ----
    local stats survival cif hazard chazard logchazard rmst timelost ///
        netsurvival excesshazard rmstnet hdifference sdifference ///
        cifdifference rmstdifference hratio sratio cifratio rmstratio
    local q ""
    foreach opt of local stats {
        if "``opt''" != "" local q `q' `opt'
    }
    if `: word count `q'' > 1 {
        di as err "specify only one prediction statistic"
        exit 198
    }
    if "`q'" == "" local q hazard
    local stat `q'

    // statistic -> Mata quantity (+ contrast kind)
    local kind ""
    if inlist("`stat'", "hdifference", "sdifference", "cifdifference", ///
        "rmstdifference") local kind difference
    if inlist("`stat'", "hratio", "sratio", "cifratio", "rmstratio") ///
        local kind ratio
    local quantity "`stat'"
    if "`kind'" != "" {
        local quantity = cond(substr("`stat'", 1, 1) == "h", "hazard", ///
            cond(substr("`stat'", 1, 1) == "s", "survival", ///
            cond(substr("`stat'", 1, 3) == "cif", "cif", "rmst")))
    }
    else {
        if "`stat'" == "netsurvival" local quantity netsurv
    }
    local docontrast = ("`kind'" != "")

    if "`timevar'" == "" local timevar _t
    // output sample = requested if/in restricted to obs with a timevar value
    marksample touse, novarlist
    markout `touse' `timevar'

    // ---- option combinations ----
    if `docontrast' {
        if `"`at'"' != "" {
            di as err "at() may not be combined with difference/ratio " ///
                "statistics; use at1() and at2()"
            exit 198
        }
        if "`standardise'" != "" {
            di as err "difference/ratio statistics are not available " ///
                "with standardise"
            exit 198
        }
    }
    else if `"`at1'"' != "" | `"`at2'"' != "" {
        di as err "at1()/at2() require a difference or ratio statistic"
        exit 198
    }
    if "`standardise'" != "" & "`stat'" == "logchazard" {
        di as err "logchazard is not available with standardise"
        exit 198
    }

    // ---- validate at-style specs: known names, name-value pairs ----
    local atok `e(atvars)'
    foreach a in at at1 at2 {
        if `"``a''"' != "" {
            local newspec ""
            tokenize `"``a''"'
            while "`1'" != "" {
                capture unab 1 : `1'
                if !`: list 1 in atok' {
                    di as err "`a'(): `1' is not a covariate, offset or " ///
                        "indicator variable of the fitted model"
                    exit 198
                }
                if "`2'" == "" {
                    di as err "`a'() requires variable-value pairs"
                    exit 198
                }
                capture confirm number `2'
                if _rc {
                    di as err "`a'(): invalid value for `1'"
                    exit 198
                }
                local _stx_`a'_`1' = `2'
                local newspec `newspec' `1' `2'
                macro shift 2
            }
            local `a' `newspec'
        }
    }

    // ---- output variables: built as tempvars so a failed prediction
    // leaves nothing behind, renamed into place on success ----
    if "`ci'" != "" {
        confirm new variable `varlist'_lci `varlist'_uci
    }
    tempvar out lci uci
    local _stx_out `out'
    qui gen double `out' = .
    if "`ci'" != "" {
        local _stx_lci `lci'
        local _stx_uci `uci'
        qui gen double `lci' = .
        qui gen double `uci' = .
    }

    local _stx_timevar  `timevar'
    local _stx_touse    `touse'
    local _stx_quantity `quantity'
    local _stx_kind     `kind'
    local _stx_zeros    `zeros'
    local _stx_ci       `ci'

    if `docontrast' {
        mata: _stx_contrast(`level')
        local qlabel "stexcess `stat'"
    }
    else if "`standardise'" != "" {
        // standardisation population = estimation sample (with at()/zeros
        // overrides applied across the population)
        tempvar pop
        qui gen byte `pop' = e(sample)
        local _stx_poptouse `pop'
        mata: _stx_standsurv(`level')
        local qlabel "standardised `stat'"
    }
    else {
        mata: _stx_predict(`level')
        local qlabel "stexcess `stat'"
    }
    rename `out' `varlist'
    if "`ci'" != "" {
        rename `lci' `varlist'_lci
        rename `uci' `varlist'_uci
        local _stx_lci `varlist'_lci
        local _stx_uci `varlist'_uci
    }

    // limits at t = 0 (log-time splines cannot be evaluated there)
    local one  survival netsurvival sratio
    local zero chazard cif rmst timelost rmstnet sdifference ///
        cifdifference rmstdifference
    if `: list stat in one' | `: list stat in zero' {
        local v = cond(`: list stat in one', 1, 0)
        qui replace `varlist' = `v' if `timevar' == 0 & `touse'
        if "`ci'" != "" {
            qui replace `_stx_lci' = `v' if `timevar' == 0 & `touse'
            qui replace `_stx_uci' = `v' if `timevar' == 0 & `touse'
        }
    }

    label var `varlist' "`qlabel'"
    if "`ci'" != "" {
        label var `_stx_lci' "lower `=100*`level''% CI"
        label var `_stx_uci' "upper `=100*`level''% CI"
    }
end
