// stexcess_mata.do -- Mata library source for stexcess (pure-Stata core).
// Compile with ado/build_lstexcess.do -> ado/lstexcess.mlib (shipped artifact).
//
// Numerics mirror the certified Python reference implementation
// (python/stexcess), which is retained for dev-time certification only, and
// extend it with v1 (merlin-based stexcess) features: multiple timescales
// (time2()-time5(), each a spline in t + offset - moffset on the log or
// natural scale), baseline offsets, identity-scale splines (time/tvctime),
// per-variable tvc dfs and noconstant.
//
// Estimation: Mata optimize() (the Newton-Raphson engine beneath ml) with an
// exact analytic d2 evaluator. Predictions: analytic delta-method Jacobians.
// Prediction semantics follow v1/merlin: covariates (and the excess
// indicator) are taken from the data row by row, overridable with at()/zeros,
// so e.g. hazard = h_ref(x_i) + ind_i * h_exc(x_i).

version 19.5

mata:
mata set matastrict on

// ========================================================================= //
// structures
// ========================================================================= //

struct _stx_rcs {                 // one restricted-cubic-spline specification
    real colvector knots          // transformed scale, ascending, df+1 entries
    real matrix R                 // orthogonalisation factor; 0x0 if raw basis
    real scalar logs              // 1: basis in ln(v); 0: basis in v
}

struct _stx_ts {                  // one timescale of a component
    real scalar num               // user-facing number: 1 = main, 2..5
    string scalar offvar, moffvar // offset variables ("" if none)
    struct _stx_rcs scalar base
    string rowvector tvcnames     // variables with tv effects on this scale
    struct _stx_rcs rowvector tvcspecs   // aligned to tvcnames
}

struct _stx_comp {                // one hazard component (reference or excess)
    string rowvector covnames     // included covariate terms, design order
    string rowvector covfull      // all terms incl. base/omitted (display)
    real rowvector covincl       // 1 if covfull[j] is estimated
    string rowvector datnames     // covnames + tvc-only variables (data cols)
    string rowvector datvars      // variables holding datnames' values (fv
                                  // terms resolve to fvrevar tempvars)
    struct _stx_ts rowvector ts   // ts[1] is the main (baseline) timescale
    real scalar cons              // 1 if _cons included
}

struct _stx_model {               // fitted model, stored for postestimation
    struct _stx_comp scalar ref, exc
    string scalar indvar          // excess indicator variable name
    string scalar wvar, wtype     // stset weights ("" if none)
    real rowvector b              // included parameters only
    real matrix V                 //   "
    real rowvector bsel           // their column positions within e(b),
                                  // whose layout includes base/omitted terms
    real scalar pr                // number of reference-equation parameters
}

struct _stx_edata {               // estimation blocks for the d2 evaluators
    real matrix CDev, CWq         // control records: reference design
    real colvector Ccw, Cd
    real matrix RDev, RWq         // patient records: reference design
    real colvector Rcw
    real matrix EDev, EWq         // patient records: excess design
    real colvector Ecw, Pd
    real scalar pr
    real colvector hrfix          // stage 2: patient reference hazard (fixed)
    real scalar ceR               //   "    : patient reference cumhaz total
}

// ========================================================================= //
// utilities
// ========================================================================= //

// user-facing error: clean one-line message + return code, no Mata traceback
void _stx_error(real scalar rc, string scalar msg)
{
    errprintf("%s\n", msg)
    exit(rc)
}

// ========================================================================= //
// quadrature, knots, spline bases
// ========================================================================= //

// Gauss-Legendre nodes/weights on [-1,1] via Golub-Welsch
real matrix _stx_gl(real scalar G)
{
    real matrix Jm, X
    real rowvector L
    real colvector k, b
    real scalar i

    if (G == 1) return((0, 2))
    k = 1::(G - 1)
    b = k :/ sqrt(4 :* k:^2 :- 1)
    Jm = J(G, G, 0)
    for (i = 1; i < G; i++) {
        Jm[i, i + 1] = b[i]
        Jm[i + 1, i] = b[i]
    }
    X = .
    L = .
    symeigensystem(Jm, X, L)
    return((L', 2 :* (X[1, .]:^2)'))   // nodes; weights = 2*phi1^2
}

// knots at equally spaced centiles of the (transformed) event values
// (numpy.percentile with linear interpolation, as in splines.default_knots).
// With weights, positions are taken in the expanded order-statistic space
// (cumulative weights), so integer fweights match the expanded data exactly
real colvector _stx_knots(real colvector v, real scalar df,
    | real colvector wt)
{
    real colvector s, w, cw, out
    real matrix sv
    real scalar n, j, pos, lo, f, W, vlo, vhi

    if (rows(v) == 0) _stx_error(2000, "no events to site spline knots")
    if (args() < 3 | rows(wt) == 0) {
        s = sort(v, 1)
        n = rows(s)
        out = J(df + 1, 1, .)
        for (j = 0; j <= df; j++) {
            pos = (n - 1) * (j / df)
            lo  = floor(pos)
            f   = pos - lo
            out[j + 1] = (lo + 1 >= n ? s[n]
                : s[lo + 1] + f * (s[lo + 2] - s[lo + 1]))
        }
        return(out)
    }
    sv = sort((v, wt), 1)
    s  = sv[., 1]
    w  = sv[., 2]
    cw = runningsum(w)
    W  = cw[rows(cw)]
    out = J(df + 1, 1, .)
    for (j = 0; j <= df; j++) {
        pos = (W - 1) * (j / df)
        lo  = floor(pos)
        f   = pos - lo
        // expanded element at 0-based index i is s[min{j : cw[j] > i}]
        vlo = s[sum(cw :<= lo) + 1]
        vhi = (lo + 1 >= W ? s[rows(s)] : s[sum(cw :<= lo + 1) + 1])
        out[j + 1] = vlo + f * (vhi - vlo)
    }
    return(out)
}

// strictly-increasing knot check (catches duplicate sited knots when event
// times are heavily tied, and bad user-supplied knot lists)
void _stx_checkknots(real colvector kn, string scalar what)
{
    if (rows(kn) < 2) _stx_error(459, what + ": at least 2 knots (df >= 1) required")
    if (missing(kn)) _stx_error(459, what + ": knots contain missing values")
    if (min(kn[|2 \ rows(kn)|] - kn[|1 \ rows(kn) - 1|]) <= 0) {
        _stx_error(459, what + ": knots not strictly increasing " +
            "(too few distinct event times?)")
    }
}

real colvector _stx_cube(real colvector u)
{
    return((u :* (u :> 0)):^3)
}

// raw RCS basis, no intercept: col 1 linear, cols 2..d restricted cubics
real matrix _stx_rcsbasis(real colvector x, real colvector k)
{
    real scalar d, j, lam, kmin, kmax
    real matrix B

    d = rows(k) - 1
    kmin = k[1]
    kmax = k[d + 1]
    B = J(rows(x), d, .)
    B[., 1] = x
    for (j = 2; j <= d; j++) {
        lam = (kmax - k[j]) / (kmax - kmin)
        B[., j] = _stx_cube(x :- k[j]) - lam :* _stx_cube(x :- kmin) -
            (1 - lam) :* _stx_cube(x :- kmax)
    }
    return(B)
}

// economy-QR orthogonalisation factor of [1, rcs], positive diagonal;
// x is already on the spec's transformed scale. With weights the rows are
// sqrt(w)-scaled, which reproduces the QR of fweight-expanded data
real matrix _stx_orthogR(real colvector x, real colvector kn,
    | real colvector wt)
{
    real matrix B1, H, R1
    real rowvector tau
    real colvector s

    B1 = (J(rows(x), 1, 1), _stx_rcsbasis(x, kn))
    if (args() == 3 & rows(wt)) B1 = sqrt(wt) :* B1
    H = .
    tau = .
    R1 = .
    hqrd(B1, H, tau, R1)
    s = sign(diagonal(R1))
    s = s + (s :== 0)
    return(R1 :* s)
}

// (optionally orthogonalised) basis at natural-scale input v; the spec's
// log/identity transform is applied here. icpt=1 keeps the leading constant
// column. ln(v) for v <= 0 is missing and propagates by design.
real matrix _stx_basis(real colvector v, struct _stx_rcs scalar spec,
    real scalar icpt)
{
    real matrix B
    real colvector x

    x = (spec.logs ? ln(v) : v)
    B = (J(rows(x), 1, 1), _stx_rcsbasis(x, spec.knots))
    if (rows(spec.R)) B = B * luinv(spec.R)   // B @ inv(R)
    return(icpt ? B : B[|1, 2 \ rows(B), cols(B)|])
}

// ========================================================================= //
// design construction
// ========================================================================= //

real scalar _stx_datidx(struct _stx_comp scalar C, string scalar name)
{
    real rowvector idx

    idx = select(1..cols(C.datnames), C.datnames :== name)
    if (!cols(idx)) _error(3498, "tvc variable not among the data columns")
    return(idx[1])
}

// component design at natural-scale times t with data X (rows aligned,
// columns = C.datnames) and per-timescale offsets OFF (rows aligned, one
// column per C.ts entry). Column order mirrors v1's linear predictor:
// [covariates] [main-ts tvc blocks] [ts2 base, ts2 tvcs] ... [baseline rcs]
// [_cons if cons]. Main-timescale tvc splines take no offset (as in v1).
real matrix _stx_cols(struct _stx_comp scalar C, real colvector t,
    real matrix X, real matrix OFF)
{
    real matrix D, B
    real colvector v
    real scalar n, s, k

    n = rows(t)
    D = J(n, 0, 0)
    if (cols(C.covnames)) D = X[|1, 1 \ n, cols(C.covnames)|]
    for (k = 1; k <= cols(C.ts[1].tvcnames); k++) {
        B = _stx_basis(t, C.ts[1].tvcspecs[k], 0)
        D = (D, X[., _stx_datidx(C, C.ts[1].tvcnames[k])] :* B)
    }
    for (s = 2; s <= cols(C.ts); s++) {
        v = t + OFF[., s]
        D = (D, _stx_basis(v, C.ts[s].base, 0))
        for (k = 1; k <= cols(C.ts[s].tvcnames); k++) {
            B = _stx_basis(v, C.ts[s].tvcspecs[k], 0)
            D = (D, X[., _stx_datidx(C, C.ts[s].tvcnames[k])] :* B)
        }
    }
    B = _stx_basis(t + OFF[., 1], C.ts[1].base, 1)
    D = (D, B[|1, 2 \ n, cols(B)|])
    if (C.cons) D = (D, B[., 1])
    return(D)
}

// "eq:term" labels in design order
string scalar _stx_names(struct _stx_comp scalar C, string scalar eq)
{
    string scalar out, tslab
    real scalar s, j, i

    out = ""
    for (j = 1; j <= cols(C.covfull); j++) {
        out = out + " " + eq + ":" + C.covfull[j]
    }
    for (j = 1; j <= cols(C.ts[1].tvcnames); j++) {
        for (i = 1; i <= rows(C.ts[1].tvcspecs[j].knots) - 1; i++) {
            out = out + " " + eq + ":" + C.ts[1].tvcnames[j] +
                "_tvc" + strofreal(i)
        }
    }
    for (s = 2; s <= cols(C.ts); s++) {
        tslab = "_t" + strofreal(C.ts[s].num)
        for (i = 1; i <= rows(C.ts[s].base.knots) - 1; i++) {
            out = out + " " + eq + ":" + tslab + "_rcs" + strofreal(i)
        }
        for (j = 1; j <= cols(C.ts[s].tvcnames); j++) {
            for (i = 1; i <= rows(C.ts[s].tvcspecs[j].knots) - 1; i++) {
                out = out + " " + eq + ":" + C.ts[s].tvcnames[j] +
                    tslab + "_tvc" + strofreal(i)
            }
        }
    }
    for (j = 1; j <= rows(C.ts[1].base.knots) - 1; j++) {
        out = out + " " + eq + ":_rcs" + strofreal(j)
    }
    if (C.cons) out = out + " " + eq + ":_cons"
    return(out)
}

// number of parameters of one component
real scalar _stx_npar(struct _stx_comp scalar C)
{
    real scalar p, s, j

    p = cols(C.covnames) + rows(C.ts[1].base.knots) - 1 + C.cons
    for (s = 1; s <= cols(C.ts); s++) {
        if (s > 1) p = p + rows(C.ts[s].base.knots) - 1
        for (j = 1; j <= cols(C.ts[s].tvcnames); j++) {
            p = p + rows(C.ts[s].tvcspecs[j].knots) - 1
        }
    }
    return(p)
}

// number of parameters of one component in the FULL e(b) layout, which
// also carries base/omitted factor-variable terms with zero coefficients
real scalar _stx_kfull(struct _stx_comp scalar C)
{
    return(_stx_npar(C) - cols(C.covnames) + cols(C.covfull))
}

// full-layout column positions of the component's included parameters
real rowvector _stx_bsel(struct _stx_comp scalar C, real scalar off)
{
    real rowvector sel
    real scalar j, prest

    sel = J(1, 0, .)
    for (j = 1; j <= cols(C.covfull); j++) {
        if (C.covincl[j]) sel = (sel, off + j)
    }
    prest = _stx_npar(C) - cols(C.covnames)
    if (prest) sel = (sel, off :+ cols(C.covfull) :+ (1..prest))
    return(sel)
}

// event design + record-major stacked quadrature design over (t0, t];
// wrec are record weights (stset weights at fit time, 1 otherwise), folded
// into the per-record quadrature weights
void _stx_qdesign(real colvector t, real colvector t0, real matrix X,
    real matrix OFF, struct _stx_comp scalar C, real colvector nd,
    real colvector w, real colvector wrec, real matrix Dev, real matrix Wq,
    real colvector cw)
{
    real matrix U
    real colvector half
    real scalar n, G

    n = rows(t)
    G = rows(nd)
    Dev = _stx_cols(C, t, X, OFF)
    half = 0.5 :* (t - t0)
    U = (0.5 :* (t + t0)) :+ half * nd'
    Wq = _stx_cols(C, vec(U'),
        (cols(X) ? X # J(G, 1, 1) : J(n * G, 0, 0)), OFF # J(G, 1, 1))
    cw = (wrec :* half) # w
}

// per-record sums over each record's G quadrature rows (missing propagates,
// so predictions at t <= 0 come back missing rather than silently 0)
real colvector _stx_rowsumG(real colvector v, real scalar G)
{
    return(rowsum(colshape(v, G), 1))
}

// nodewise design sums: out[i,j] = sum_g e[(i-1)G+g] * Wq[(i-1)G+g, j]
real matrix _stx_qsum(real colvector e, real matrix Wq, real scalar G)
{
    real matrix out
    real scalar j

    out = J(rows(e) / G, cols(Wq), .)
    for (j = 1; j <= cols(Wq); j++) {
        out[., j] = rowsum(colshape(e :* Wq[., j], G), 1)
    }
    return(out)
}

// ========================================================================= //
// offsets and component specs (driven by the .ado macro contract)
// ========================================================================= //

// combined offset (off - moff) over touse rows; prediction-time at()
// overrides are applied by the wrapper as data replacements, so values are
// always read straight from the variables
real colvector _stx_offdata(string scalar offvar, string scalar moffvar,
    string scalar touse, real scalar n)
{
    real colvector off

    off = J(n, 1, 0)
    if (offvar != "")  off = off + st_data(., offvar, touse)
    if (moffvar != "") off = off - st_data(., moffvar, touse)
    return(off)
}

// per-timescale offset matrix for one component over touse rows
real matrix _stx_offmat(struct _stx_comp scalar C, string scalar touse,
    real scalar n)
{
    real matrix OFF
    real scalar s

    OFF = J(n, cols(C.ts), 0)
    for (s = 1; s <= cols(C.ts); s++) {
        if (C.ts[s].offvar != "" | C.ts[s].moffvar != "") {
            OFF[., s] = _stx_offdata(C.ts[s].offvar, C.ts[s].moffvar,
                touse, n)
        }
    }
    return(OFF)
}

// build one component spec from the wrapper's macro contract:
//   _stx_<tag>vars, _stx_<tag>_cons, _stx_<tag>_tslist and, per timescale j,
//   _stx_<tag>_ts<j>_{df,knots,log,orthog,off,moff,tvc,dftvc,tvclog}.
// Knots are sited from the transformed timescale values at this component's
// knot-siting event records (evmask over touse rows); the orthogonalisation
// factor uses the same sample.
void _stx_compspec(struct _stx_comp scalar C, string scalar tag,
    string scalar touse, real colvector t, real colvector evmask,
    real colvector wt, string scalar what)
{
    string scalar p, kn_s
    string rowvector tsnums
    real colvector t_ev, off, v_ev, w_ev
    real rowvector dftvc
    real scalar n, s, j, df

    n = rows(t)
    w_ev = select(wt, evmask)
    C.covnames = tokens(st_local("_stx_" + tag + "vars"))
    C.covfull  = tokens(st_local("_stx_" + tag + "covfull"))
    C.covincl  = strtoreal(tokens(st_local("_stx_" + tag + "covincl")))
    if (!cols(C.covfull)) {                   // plain-variable fallback
        C.covfull = C.covnames
        C.covincl = J(1, cols(C.covnames), 1)
    }
    C.datnames = C.covnames
    // fv terms resolve to fvrevar tempvars supplied by the wrapper;
    // tvc-only variables (appended below) hold their own values
    C.datvars = tokens(st_local("_stx_" + tag + "covmap"))
    if (!cols(C.datvars)) C.datvars = C.covnames
    C.cons = st_local("_stx_" + tag + "_cons") != "0"
    tsnums = tokens(st_local("_stx_" + tag + "_tslist"))
    C.ts = J(1, cols(tsnums), _stx_ts())

    // NB: fields are assigned directly into C.ts[s] -- copying structs
    // through a scalar temp would alias and corrupt earlier entries
    for (s = 1; s <= cols(tsnums); s++) {
        C.ts[s].num = strtoreal(tsnums[s])
        p = "_stx_" + tag + "_ts" + tsnums[s] + "_"
        C.ts[s].offvar  = st_local(p + "off")
        C.ts[s].moffvar = st_local(p + "moff")
        C.ts[s].base.logs = st_local(p + "log") != "0"

        off = _stx_offdata(C.ts[s].offvar, C.ts[s].moffvar, touse, n)
        t_ev = select(t + off, evmask)
        kn_s = st_local(p + "knots")
        df = strtoreal(st_local(p + "df"))
        v_ev = (C.ts[s].base.logs ? ln(t_ev) : t_ev)
        C.ts[s].base.knots = (kn_s == "" ? _stx_knots(v_ev, df, w_ev)
                                         : strtoreal(tokens(kn_s))')
        _stx_checkknots(C.ts[s].base.knots, what +
            (s > 1 ? " time" + tsnums[s] + "()" : "") + " baseline spline")
        C.ts[s].base.R = (st_local(p + "orthog") != "0"
            ? _stx_orthogR(v_ev, C.ts[s].base.knots, w_ev) : J(0, 0, .))

        C.ts[s].tvcnames = tokens(st_local(p + "tvc"))
        if (cols(C.ts[s].tvcnames)) {
            dftvc = strtoreal(tokens(st_local(p + "dftvc")))
            C.ts[s].tvcspecs = J(1, cols(C.ts[s].tvcnames), _stx_rcs())
            for (j = 1; j <= cols(C.ts[s].tvcnames); j++) {
                C.ts[s].tvcspecs[j].logs = st_local(p + "tvclog") != "0"
                // main-timescale tvc splines take no offset (v1 behaviour)
                v_ev = (s == 1 ? select(t, evmask) : t_ev)
                if (C.ts[s].tvcspecs[j].logs) v_ev = ln(v_ev)
                C.ts[s].tvcspecs[j].knots = _stx_knots(v_ev, dftvc[j], w_ev)
                _stx_checkknots(C.ts[s].tvcspecs[j].knots, what + " tvc(" +
                    C.ts[s].tvcnames[j] + ") spline")
                C.ts[s].tvcspecs[j].R = (st_local(p + "orthog") != "0"
                    ? _stx_orthogR(v_ev, C.ts[s].tvcspecs[j].knots, w_ev)
                    : J(0, 0, .))
                if (!anyof(C.datnames, C.ts[s].tvcnames[j])) {
                    C.datnames = (C.datnames, C.ts[s].tvcnames[j])
                    C.datvars  = (C.datvars,  C.ts[s].tvcnames[j])
                }
            }
        }
    }
}

// validate timescale positivity for log-scale splines over (t0, t]
void _stx_checkts(struct _stx_comp scalar C, string scalar touse,
    real colvector t, real colvector t0, string scalar what)
{
    real colvector off
    real scalar s, needpos, j

    for (s = 1; s <= cols(C.ts); s++) {
        needpos = C.ts[s].base.logs
        for (j = 1; j <= cols(C.ts[s].tvcnames); j++) {
            if (s > 1) needpos = needpos | C.ts[s].tvcspecs[j].logs
        }
        if (!needpos) continue
        if (C.ts[s].offvar == "" & C.ts[s].moffvar == "") continue
        off = _stx_offdata(C.ts[s].offvar, C.ts[s].moffvar, touse, rows(t))
        if (missing(off)) {
            _stx_error(459, what + ": offset variables contain missing values")
        }
        if (min(t + off) <= 0 | min(t0 + off) < 0) {
            _stx_error(459, what + " time" + strofreal(C.ts[s].num) +
                "(): t + offset must be > 0 (log scale); use the time option?")
        }
    }
}

// ========================================================================= //
// d2 evaluators (exact ll, gradient, Hessian)
// ========================================================================= //

// joint likelihood over both parameter blocks
void _stx_eval_joint(real scalar todo, real rowvector b,
    struct _stx_edata scalar D, real scalar lnf, real rowvector g,
    real matrix H)
{
    real colvector thr, the, eC, eR, eE, hr, he, pi, w2
    real matrix Hrr, Hre, Hee
    real scalar pr, k

    pr = D.pr
    k  = cols(b)
    thr = b[(1..pr)]'
    the = b[((pr + 1)..k)]'

    eC = D.Ccw :* exp(D.CWq * thr)
    eR = D.Rcw :* exp(D.RWq * thr)
    eE = D.Ecw :* exp(D.EWq * the)
    hr = exp(D.RDev * thr)
    he = exp(D.EDev * the)

    lnf = cross(D.Cd, D.CDev * thr) - sum(eC) +
          cross(D.Pd, ln(hr + he)) - sum(eR) - sum(eE)
    if (missing(lnf) | todo < 1) return

    pi = hr :/ (hr + he)
    g = (cross(D.CDev, D.Cd) - cross(D.CWq, eC) +
         cross(D.RDev, D.Pd :* pi) - cross(D.RWq, eR) \
         cross(D.EDev, D.Pd :* (1 :- pi)) - cross(D.EWq, eE))'
    if (todo < 2) return

    w2  = D.Pd :* pi :* (1 :- pi)
    Hrr = cross(D.RDev, w2, D.RDev) -
          cross(D.CWq, eC, D.CWq) - cross(D.RWq, eR, D.RWq)
    Hre = -cross(D.RDev, w2, D.EDev)
    Hee = cross(D.EDev, w2, D.EDev) - cross(D.EWq, eE, D.EWq)
    H = makesymmetric((Hrr, Hre \ Hre', Hee))
}

// two-stage, stage 1: reference model over control records only
void _stx_eval_ref(real scalar todo, real rowvector b,
    struct _stx_edata scalar D, real scalar lnf, real rowvector g,
    real matrix H)
{
    real colvector thr, eC

    thr = b'
    eC = D.Ccw :* exp(D.CWq * thr)
    lnf = cross(D.Cd, D.CDev * thr) - sum(eC)
    if (missing(lnf) | todo < 1) return
    g = (cross(D.CDev, D.Cd) - cross(D.CWq, eC))'
    if (todo < 2) return
    H = makesymmetric(-cross(D.CWq, eC, D.CWq))
}

// two-stage, stage 2: excess model over patients, reference hazard fixed
void _stx_eval_exc(real scalar todo, real rowvector b,
    struct _stx_edata scalar D, real scalar lnf, real rowvector g,
    real matrix H)
{
    real colvector the, eE, he, pi, w2

    the = b'
    eE = D.Ecw :* exp(D.EWq * the)
    he = exp(D.EDev * the)
    lnf = cross(D.Pd, ln(D.hrfix + he)) - D.ceR - sum(eE)
    if (missing(lnf) | todo < 1) return
    pi = D.hrfix :/ (D.hrfix + he)
    g = (cross(D.EDev, D.Pd :* (1 :- pi)) - cross(D.EWq, eE))'
    if (todo < 2) return
    w2 = D.Pd :* pi :* (1 :- pi)
    H = makesymmetric(cross(D.EDev, w2, D.EDev) - cross(D.EWq, eE, D.EWq))
}

// one optimize() run; returns the handle so results can be queried
transmorphic _stx_optimize(pointer(real function) scalar fn,
    struct _stx_edata scalar D, real rowvector b0, real scalar trace,
    | real scalar maxiter)
{
    transmorphic S
    real scalar v

    S = optimize_init()
    optimize_init_evaluator(S, fn)
    optimize_init_evaluatortype(S, "d2")
    optimize_init_argument(S, 1, D)
    optimize_init_params(S, b0)
    optimize_init_technique(S, "nr")
    optimize_init_tracelevel(S, (trace ? "value" : "none"))
    optimize_init_valueid(S, "log likelihood")
    // user maximize options from the wrapper (empty locals -> defaults);
    // an explicit maxiter argument (the starting-values pre-fit) wins
    v = strtoreal(st_local("_stx_ptol"))
    if (v < .) optimize_init_conv_ptol(S, v)
    v = strtoreal(st_local("_stx_vtol"))
    if (v < .) optimize_init_conv_vtol(S, v)
    v = strtoreal(st_local("_stx_nrtol"))
    if (v < .) optimize_init_conv_nrtol(S, v)
    v = strtoreal(st_local("_stx_iterate"))
    if (args() == 5)  optimize_init_conv_maxiter(S, maxiter)
    else if (v < .)   optimize_init_conv_maxiter(S, v)
    (void) _optimize(S)
    return(S)
}

// ========================================================================= //
// model store (process-level, like e() for postestimation)
// ========================================================================= //

void _stx_putmodel(struct _stx_model scalar M)
{
    external pointer(struct _stx_model scalar) scalar STX_FIT

    STX_FIT = &M
}

struct _stx_model scalar _stx_getmodel()
{
    external pointer(struct _stx_model scalar) scalar STX_FIT

    if (STX_FIT == NULL) {
        _stx_error(301, "no stexcess fit in memory; run (or rerun) stexcess first")
    }
    return(*STX_FIT)
}

// fetch the stored model and verify it matches the active e() results --
// guards against predicting from a stale store after, e.g., estimates restore
struct _stx_model scalar _stx_usemodel()
{
    struct _stx_model scalar M
    real matrix eb

    M = _stx_getmodel()
    eb = st_matrix("e(b)")
    // two steps: Mata | does not short-circuit, and subscripting requires
    // the columns to exist
    if (cols(eb) < max(M.bsel)) {
        _stx_error(301, "fit in memory does not match e(b) " +
            "(estimates restore?); rerun stexcess")
    }
    if (mreldif(eb[M.bsel], M.b) > 1e-12) {
        _stx_error(301, "fit in memory does not match e(b) " +
            "(estimates restore?); rerun stexcess")
    }
    return(M)
}

// ========================================================================= //
// fit driver -- consumes the .ado macro/scalar contract
// ========================================================================= //

void _stx_fit()
{
    struct _stx_edata scalar D
    struct _stx_model scalar M
    transmorphic S, S2
    string scalar touse, fromname
    real colvector t, t0, d, ind, nd, w, cm, pm, eC, eR, eE, hr, he, pi, w2
    real colvector wgt, wd, thr, the
    real matrix Spr, Spe, Brob
    real matrix Xr, Xe, OFFr, OFFe, glm, Dev, Wq, V, A, B, Sc, Sp, Ainv
    real matrix A11, A21, A22
    real rowvector b0, b
    real colvector cw
    real scalar n, G, k, kfull, rate, twostage, ll, conv, iter, pe, s, j
    real scalar trace
    real rowvector bf
    real matrix Vf

    touse = st_local("_stx_touse")
    t   = st_data(., st_local("_stx_t"), touse)
    t0  = st_data(., st_local("_stx_t0"), touse)
    d   = st_data(., st_local("_stx_d"), touse)
    ind = st_data(., st_local("_stx_ind"), touse)
    n   = rows(t)
    M.wvar  = st_local("_stx_wvar")
    M.wtype = st_local("_stx_wtype")
    wgt = (M.wvar == "" ? J(n, 1, 1) : st_data(., M.wvar, touse))
    wd  = wgt :* d
    twostage = st_local("_stx_twostage") == "1"
    trace = st_local("_stx_nolog") == ""

    // input validation
    if (min(t) <= 0) _stx_error(459, "all exit times must be > 0")
    if (min(t0) < 0 | min(t - t0) <= 0) _stx_error(459, "require 0 <= t0 < t for every record")
    cm = ind :== 0
    pm = ind :== 1
    if (sum(cm) + sum(pm) < n) _stx_error(450, "indicator must be coded 0 (reference) / 1 (excess)")
    if (!sum(pm :& (d :== 1))) _stx_error(2000, "no excess events: cannot site the excess baseline knots")
    if (twostage & !sum(cm)) _stx_error(2000, "twostage requires reference (control) records")

    // component specs: ref knots from ALL events, exc from patient events
    _stx_compspec(M.ref, "ref", touse, t, d :== 1, wgt, "reference")
    _stx_compspec(M.exc, "exc", touse, t, (d :== 1) :& pm, wgt, "excess")
    _stx_checkts(M.ref, touse, t, t0, "reference")
    _stx_checkts(M.exc, touse, t, t0, "excess")
    M.indvar = st_local("_stx_ind")

    // store the knots used (transformed scale) for e()
    for (s = 1; s <= cols(M.ref.ts); s++) {
        st_local("_stx_kref" + (s > 1 ? "_t" + strofreal(M.ref.ts[s].num) : ""),
            invtokens(strtrim(strofreal(M.ref.ts[s].base.knots', "%21.0g"))))
    }
    for (s = 1; s <= cols(M.exc.ts); s++) {
        st_local("_stx_kexc" + (s > 1 ? "_t" + strofreal(M.exc.ts[s].num) : ""),
            invtokens(strtrim(strofreal(M.exc.ts[s].base.knots', "%21.0g"))))
    }


    G = strtoreal(st_local("_stx_nnodes"))
    glm = _stx_gl(G)
    nd = glm[., 1]
    w  = glm[., 2]

    Xr = (cols(M.ref.datvars) ? st_data(., M.ref.datvars, touse) : J(n, 0, 0))
    Xe = (cols(M.exc.datvars) ? st_data(., M.exc.datvars, touse) : J(n, 0, 0))
    OFFr = _stx_offmat(M.ref, touse, n)
    OFFe = _stx_offmat(M.exc, touse, n)

    // estimation blocks
    Dev = .
    Wq = .
    cw = .
    // weights enter the likelihood only through the per-record event
    // multipliers (Cd/Pd = w*d) and quadrature weights (cw), so the d2
    // evaluators and the two-stage sandwich are weight-correct as they stand
    _stx_qdesign(select(t, cm), select(t0, cm),
        (cols(Xr) ? select(Xr, cm) : J(sum(cm), 0, 0)), select(OFFr, cm),
        M.ref, nd, w, select(wgt, cm), Dev, Wq, cw)
    D.CDev = Dev
    D.CWq  = Wq
    D.Ccw  = cw
    D.Cd   = select(wd, cm)
    _stx_qdesign(select(t, pm), select(t0, pm),
        (cols(Xr) ? select(Xr, pm) : J(sum(pm), 0, 0)), select(OFFr, pm),
        M.ref, nd, w, select(wgt, pm), Dev, Wq, cw)
    D.RDev = Dev
    D.RWq  = Wq
    D.Rcw  = cw
    _stx_qdesign(select(t, pm), select(t0, pm),
        (cols(Xe) ? select(Xe, pm) : J(sum(pm), 0, 0)), select(OFFe, pm),
        M.exc, nd, w, select(wgt, pm), Dev, Wq, cw)
    D.EDev = Dev
    D.EWq  = Wq
    D.Ecw  = cw
    D.Pd   = select(wd, pm)
    D.pr   = cols(D.CDev)
    pe = cols(D.EDev)
    k  = D.pr + pe
    if (missing(D.CDev) | missing(D.CWq) | missing(D.RDev) | missing(D.RWq) |
        missing(D.EDev) | missing(D.EWq)) {
        _stx_error(459, "design contains missing values " +
            "(offset variables missing, or log of a non-positive timescale?)")
    }

    // crude exponential start values, or user-supplied from()
    b0 = J(1, k, 0)
    if (twostage) {
        rate = max((sum(select(wd, cm)) /
                    sum(select(wgt :* (t - t0), cm)), 1e-4))
    }
    else {
        rate = max((sum(wd) / sum(wgt :* (t - t0)), 1e-4))
    }
    if (M.ref.cons) b0[D.pr] = ln(rate)
    rate = max((sum(select(wd, pm)) /
                sum(select(wgt :* (t - t0), pm)), 1e-4))
    if (M.exc.cons) b0[k] = ln(0.5 * rate)
    // bsel maps the included parameters into the full e(b) layout (which
    // also carries base/omitted factor-variable terms as zero coefficients)
    M.bsel = (_stx_bsel(M.ref, 0), _stx_bsel(M.exc, _stx_kfull(M.ref)))
    kfull = _stx_kfull(M.ref) + _stx_kfull(M.exc)
    fromname = st_local("_stx_from")
    if (fromname != "") {
        if (cols(st_matrix(fromname)) != kfull) {
            _stx_error(198, "from(): matrix must have " + strofreal(kfull) +
                " columns (one per e(b) column)")
        }
        b0 = st_matrix(fromname)[M.bsel]
    }

    if (!twostage) {
        // starting values: fit the reference model to the control records
        // alone (cheap: half the data, half the parameters, exact d2) and
        // start the excess block at zero apart from its crude event-rate
        // intercept; falls back to the crude values if the pre-fit fails
        if (fromname == "" & sum(D.Cd) > 0) {
            if (trace) printf("\n{txt}Obtaining starting values:\n")
            S = _stx_optimize(&_stx_eval_ref(), D, b0[(1..D.pr)], 0, 25)
            if (optimize_result_errorcode(S) == 0 &
                !missing(optimize_result_params(S))) {
                b0[(1..D.pr)] = optimize_result_params(S)
            }
        }
        if (trace) printf("\n{txt}Fitting full model:\n")
        S = _stx_optimize(&_stx_eval_joint(), D, b0, trace)
        b = optimize_result_params(S)
        V = optimize_result_V_oim(S)
        ll = optimize_result_value(S)
        conv = optimize_result_converged(S)
        iter = optimize_result_iterations(S)

        // pweights: robust sandwich A^-1 B A^-1 with B from the per-record
        // weighted scores (the weights are already inside Cd/Pd/cw)
        if (M.wtype == "pweight") {
            thr = b[(1..D.pr)]'
            the = b[((D.pr + 1)..k)]'
            eC = D.Ccw :* exp(D.CWq * thr)
            eR = D.Rcw :* exp(D.RWq * thr)
            eE = D.Ecw :* exp(D.EWq * the)
            pi = exp(D.RDev * thr) :/
                 (exp(D.RDev * thr) + exp(D.EDev * the))
            Sc  = D.Cd :* D.CDev - _stx_qsum(eC, D.CWq, G)
            Spr = (D.Pd :* pi) :* D.RDev - _stx_qsum(eR, D.RWq, G)
            Spe = (D.Pd :* (1 :- pi)) :* D.EDev - _stx_qsum(eE, D.EWq, G)
            Brob = cross((Sc, J(rows(Sc), pe, 0) \ Spr, Spe),
                         (Sc, J(rows(Sc), pe, 0) \ Spr, Spe))
            V = V * Brob * V
            V = 0.5 :* (V + V')
        }
    }
    else {
        // stage 1: reference on controls; stage 2: excess with ref fixed
        if (trace) printf("\n{txt}Fitting reference model (controls only):\n")
        S = _stx_optimize(&_stx_eval_ref(), D, b0[(1..D.pr)], trace)
        eR = D.Rcw :* exp(D.RWq * optimize_result_params(S)')
        D.hrfix = exp(D.RDev * optimize_result_params(S)')
        D.ceR = sum(eR)
        if (trace) printf("\n{txt}Fitting excess model (reference fixed):\n")
        S2 = _stx_optimize(&_stx_eval_exc(), D, b0[((D.pr + 1)..k)], trace)
        b = (optimize_result_params(S), optimize_result_params(S2))
        ll = optimize_result_value(S) + optimize_result_value(S2)
        conv = optimize_result_converged(S) & optimize_result_converged(S2)
        iter = optimize_result_iterations(S) + optimize_result_iterations(S2)

        // stacked M-estimation sandwich V = A^-1 B A^-T
        eC = D.Ccw :* exp(D.CWq * b[(1..D.pr)]')
        eE = D.Ecw :* exp(D.EWq * b[((D.pr + 1)..k)]')
        hr = D.hrfix
        he = exp(D.EDev * b[((D.pr + 1)..k)]')
        pi = hr :/ (hr + he)
        w2 = D.Pd :* pi :* (1 :- pi)
        A11 = cross(D.CWq, eC, D.CWq)
        A21 = cross(D.EDev, w2, D.RDev)
        A22 = cross(D.EWq, eE, D.EWq) - cross(D.EDev, w2, D.EDev)
        A = (A11, J(D.pr, pe, 0) \ A21, A22)
        Sc = D.Cd :* D.CDev - _stx_qsum(eC, D.CWq, G)
        Sp = (D.Pd :* (1 :- pi)) :* D.EDev - _stx_qsum(eE, D.EWq, G)
        // the score rows carry one factor of the weights; pweights keep it
        // squared (sum of (w s)(w s)'), fweights/iweights need the
        // expanded-data sum of w s s'
        if (M.wtype != "pweight" & M.wvar != "") {
            Sc = Sc :/ sqrt(select(wgt, cm))
            Sp = Sp :/ sqrt(select(wgt, pm))
        }
        B = (cross(Sc, Sc), J(D.pr, pe, 0) \ J(pe, D.pr, 0), cross(Sp, Sp))
        Ainv = luinv(A)
        V = Ainv * B * Ainv'
        V = 0.5 :* (V + V')
    }

    M.b  = b
    M.pr = D.pr
    M.V  = V
    _stx_putmodel(M)

    bf = J(1, kfull, 0)
    bf[M.bsel] = b
    Vf = J(kfull, kfull, 0)
    Vf[M.bsel, M.bsel] = V
    st_matrix(st_local("_stx_bmat"), bf)
    st_matrix(st_local("_stx_Vmat"), Vf)
    st_local("_stx_ll",   strofreal(ll, "%21.0g"))
    st_local("_stx_k",    strofreal(k))
    st_local("_stx_conv", strofreal(conv))
    st_local("_stx_iter", strofreal(iter))
    st_local("_stx_names",
        strtrim(_stx_names(M.ref, "ref") + _stx_names(M.exc, "exc")))
}

// ========================================================================= //
// delta-method CIs
// ========================================================================= //

void _stx_delta(real colvector est, real matrix Jc, real matrix V,
    string scalar transform, real scalar level,
    real colvector lci, real colvector uci)
{
    real colvector se, S, zt, se_t, l, u
    real scalar z

    se = rowsum((Jc * V) :* Jc)
    se = sqrt(se :* (se :> 0))         // clip negative variances from rounding
    z = invnormal(0.5 + level / 200)
    if (transform == "log") {
        se_t = se :/ est
        lci = est :* exp(-z :* se_t)
        uci = est :* exp( z :* se_t)
        return
    }
    if (transform == "cloglog" | transform == "cif") {
        // cif: est = 1 - S with se(est) = se(S); bound S then complement
        S = (transform == "cif" ? 1 :- est : est)
        S = rowmin((rowmax((S, J(rows(S), 1, 1e-12))),
                    J(rows(S), 1, 1 - 1e-12)))
        zt = ln(-ln(S))
        se_t = abs(1 :/ (S :* ln(S))) :* se
        l = exp(-exp(zt :+ z :* se_t))     // lower S bound
        u = exp(-exp(zt :- z :* se_t))     // upper S bound
        if (transform == "cif") {
            lci = 1 :- u
            uci = 1 :- l
        }
        else {
            lci = l
            uci = u
        }
        return
    }
    lci = est - z :* se
    uci = est + z :* se
}

string scalar _stx_transform(string scalar q)
{
    if (q == "survival" | q == "netsurv") return("cloglog")
    if (q == "cif") return("cif")
    if (q == "logchazard") return("identity")
    return("log")   // hazards, cumulative hazards, rmst, timelost
}

// quantities that involve the per-record excess indicator
real scalar _stx_needind(string scalar q)
{
    return(q == "hazard" | q == "chazard" | q == "logchazard" |
           q == "survival" | q == "cif" | q == "rmst" | q == "timelost")
}

// ========================================================================= //
// row-based prediction data
// ========================================================================= //

// the wrapper rebuilds factor-variable terms with fvrevar at predict time
// (so at()/zeros replacements of the underlying variables propagate) and
// passes the resulting variable maps in _stx_refmap/_stx_excmap
void _stx_mapinfo()
{
    struct _stx_model scalar M

    M = _stx_usemodel()
    st_local("_stx_refdat", invtokens(M.ref.datnames))
    st_local("_stx_excdat", invtokens(M.exc.datnames))
}

// covariate/offset/indicator blocks for predictions over touse rows
void _stx_preddata(struct _stx_model scalar M, string scalar touse,
    real scalar m, real scalar needind,
    real matrix Xr, real matrix Xe, real matrix OFFr, real matrix OFFe,
    real colvector ind)
{
    string rowvector mr, me

    mr = tokens(st_local("_stx_refmap"))
    me = tokens(st_local("_stx_excmap"))
    if (cols(mr) != cols(M.ref.datnames) |
        cols(me) != cols(M.exc.datnames)) {
        _error(3498, "internal: prediction variable map misaligned")
    }
    Xr = (cols(mr) ? st_data(., mr, touse) : J(m, 0, 0))
    Xe = (cols(me) ? st_data(., me, touse) : J(m, 0, 0))
    OFFr = _stx_offmat(M.ref, touse, m)
    OFFe = _stx_offmat(M.exc, touse, m)
    if (needind) {
        if (_st_varindex(M.indvar) >= .) {
            _stx_error(111, "the indicator variable " + M.indvar +
                " was not found; supply it with at()")
        }
        ind = st_data(., M.indvar, touse)
    }
    else ind = J(m, 1, 1)
}

// ========================================================================= //
// covariate-pattern quantities with analytic Jacobians (per-row data)
// ========================================================================= //

// hazard/cumhaz/survival families row by row; est (m x 1) and, when doJ,
// Jc (m x k) on the natural scale. ind weights the excess component.
void _stx_rowquant(struct _stx_model scalar M, real colvector times,
    real matrix Xr, real matrix Xe, real matrix OFFr, real matrix OFFe,
    real colvector ind, string scalar q, real scalar G, real scalar doJ,
    real colvector est, real matrix Jc)
{
    real matrix rDev, rWq, eDev, eWq, Jr, Je, glm
    real colvector rcw, ecw, er, ee, hr, he, Hr, He, ch, zero
    real colvector nd, w
    real scalar m, pr, k, needH, needref

    m = rows(times)
    pr = M.pr
    k = cols(M.b)
    glm = _stx_gl(G)
    nd = glm[., 1]
    w  = glm[., 2]
    zero = J(m, 1, 0)
    needref = !(q == "netsurv" | q == "excesshazard")
    needH = !(q == "hazard" | q == "excesshazard")

    rDev = rWq = eDev = eWq = .
    rcw = ecw = .
    hr = Hr = J(m, 1, 0)
    Jr = J(m, pr, 0)
    if (needref) {
        _stx_qdesign(times, zero, Xr, OFFr, M.ref, nd, w,
            J(m, 1, 1), rDev, rWq, rcw)
        hr = exp(rDev * M.b[(1..pr)]')
        if (needH) {
            er = rcw :* exp(rWq * M.b[(1..pr)]')
            Hr = _stx_rowsumG(er, G)
            if (doJ) Jr = _stx_qsum(er, rWq, G)
        }
    }
    _stx_qdesign(times, zero, Xe, OFFe, M.exc, nd, w,
        J(m, 1, 1), eDev, eWq, ecw)
    he = exp(eDev * M.b[((pr + 1)..k)]')
    if (needH) {
        ee = ecw :* exp(eWq * M.b[((pr + 1)..k)]')
        He = _stx_rowsumG(ee, G)
        if (doJ) Je = _stx_qsum(ee, eWq, G)
    }

    Jc = J(0, 0, .)
    if (q == "hazard" | q == "excesshazard") {
        est = (q == "hazard" ? hr + ind :* he : he)
        if (!doJ) return
        Jc = J(m, k, 0)
        if (q == "hazard") {
            Jc[|1, 1 \ m, pr|] = hr :* rDev
            Jc[|1, pr + 1 \ m, k|] = (ind :* he) :* eDev
        }
        else Jc[|1, pr + 1 \ m, k|] = he :* eDev
        return
    }
    if (q == "chazard" | q == "logchazard") {
        ch = Hr + ind :* He
        est = (q == "logchazard" ? ln(ch) : ch)
        if (!doJ) return
        Jc = (Jr, ind :* Je)
        if (q == "logchazard") Jc = Jc :/ ch
        return
    }
    if (q == "survival" | q == "cif") {
        ch = exp(-(Hr + ind :* He))            // S(t)
        est = (q == "cif" ? 1 :- ch : ch)
        if (!doJ) return
        Jc = (-ch :* Jr, -ch :* (ind :* Je))
        if (q == "cif") Jc = -Jc
        return
    }
    if (q == "netsurv") {
        est = exp(-He)
        if (!doJ) return
        Jc = (J(m, pr, 0), -est :* Je)
        return
    }
    _error(3498, "unknown quantity: " + q)
}

// RMST(tau) = int_0^tau S(u) du by outer Gauss-Legendre over the survival
// curve (inner nodes = G, outer = 40); timelost = tau - RMST;
// rmstnet integrates net survival
void _stx_rmstrows(struct _stx_model scalar M, real colvector taus,
    real matrix Xr, real matrix Xe, real matrix OFFr, real matrix OFFe,
    real colvector ind, string scalar q, real scalar G, real scalar doJ,
    real colvector est, real matrix Jc)
{
    real matrix glm, Jo, U, Xrx, Xex
    real colvector no, wo, Sflat, half_out
    real scalar m, Mo, j

    m = rows(taus)
    glm = _stx_gl(40)
    no = glm[., 1]
    wo = glm[., 2]
    Mo = 40
    half_out = 0.5 :* taus
    U = half_out * (no :+ 1)'              // (m x Mo) outer evaluation times

    Xrx = (cols(Xr) ? Xr # J(Mo, 1, 1) : J(m * Mo, 0, 0))
    Xex = (cols(Xe) ? Xe # J(Mo, 1, 1) : J(m * Mo, 0, 0))
    Sflat = .
    Jo = .
    _stx_rowquant(M, vec(U'), Xrx, Xex, OFFr # J(Mo, 1, 1),
        OFFe # J(Mo, 1, 1), ind # J(Mo, 1, 1),
        (q == "rmstnet" ? "netsurv" : "survival"), G, doJ, Sflat, Jo)
    est = half_out :* (colshape(Sflat, Mo) * wo)
    if (q == "timelost") est = taus - est
    if (!doJ) return
    Jc = J(m, cols(M.b), .)
    for (j = 1; j <= cols(M.b); j++) {
        Jc[., j] = half_out :* (colshape(Jo[., j], Mo) * wo)
    }
    if (q == "timelost") Jc = -Jc
}

// dispatch: any quantity over per-row data (rmst inner nodes = G)
void _stx_quantity(struct _stx_model scalar M, real colvector times,
    real matrix Xr, real matrix Xe, real matrix OFFr, real matrix OFFe,
    real colvector ind, string scalar q, real scalar G, real scalar doJ,
    real colvector est, real matrix Jc)
{
    if (q == "rmst" | q == "rmstnet" | q == "timelost") {
        _stx_rmstrows(M, times, Xr, Xe, OFFr, OFFe, ind, q, G, doJ, est, Jc)
    }
    else {
        _stx_rowquant(M, times, Xr, Xe, OFFr, OFFe, ind, q, G, doJ, est, Jc)
    }
}

// ========================================================================= //
// standardisation (g-formula) with analytic Jacobians
// ========================================================================= //

// per-individual hazard/cumhaz (and Jacobian blocks) for ONE component at
// ONE time over a population chunk. Elements whose timescale carries no
// offset have node bases shared across individuals (evaluated G times, not
// c*G) -- the fast path that makes standardisation cheap at scale; offset
// timescales fall back to per-individual bases.
void _stx_std_comp(struct _stx_comp scalar C, real colvector th,
    real scalar tj, real matrix X, real matrix OFF, real colvector nd,
    real colvector w, real scalar needh, real scalar needH, real scalar doJ,
    real colvector h, real colvector H, real matrix Jh, real matrix JH)
{
    pointer(struct _stx_rcs scalar) scalar psp
    real matrix E, W, EW, B, JhM, JHM, LP
    real colvector u, lph, gk, z, Ew, voff
    real scalar G, c, half, ncov, nel, e, s, kk, p, pos, l, shared, isb1
    real scalar zidx, sx, nk, kk0, kend
    pointer(real matrix) rowvector cEV, cQ
    real rowvector czi, cp, cshared, cb1

    G = rows(nd)
    c = rows(OFF)
    half = 0.5 * tj
    u = half :+ half :* nd                        // nodes on (0, tj]
    ncov = cols(C.covnames)

    // enumerate spline elements in design order (ts1 tvcs, ts2.., ts1 base)
    nel = cols(C.ts[1].tvcnames)
    for (s = 2; s <= cols(C.ts); s++) nel = nel + 1 + cols(C.ts[s].tvcnames)
    nel = nel + 1                                  // ts1 base (+ cons)
    cEV = cQ = J(1, nel, NULL)
    czi = cp = cshared = cb1 = J(1, nel, 0)

    lph = (ncov ? X[|1, 1 \ c, ncov|] * th[|1 \ ncov|] : J(c, 1, 0))
    LP  = (needH ? J(c, G, 0) :+ lph : J(0, 0, .))

    pos = ncov
    e = 0
    for (s = 1; s <= cols(C.ts) + 1; s++) {
        // pass order: s=1 -> ts1 tvcs; 2..nts -> additional ts; last -> ts1 base
        isb1 = s == cols(C.ts) + 1
        sx = (isb1 ? 1 : s)
        voff = (sx > 1 | isb1 ? OFF[., sx] : J(c, 1, 0))
        shared = (!isb1 & sx == 1 ? 1 : !any(voff :!= 0))
        nk = cols(C.ts[sx].tvcnames)
        kk0  = (isb1 ? 0 : (sx == 1 ? 1 : 0))      // 0 slot = the ts base
        kend = (isb1 ? 0 : nk)
        for (kk = kk0; kk <= kend; kk++) {
            e = e + 1
            if (kk == 0) {                         // a baseline element
                // NB: a pointer, because struct-scalar assignment copies into
                // the variable's current referent and would corrupt the model
                psp = &(C.ts[sx].base)
                zidx = 0
            }
            else {
                psp = &(C.ts[sx].tvcspecs[kk])
                zidx = _stx_datidx(C, C.ts[sx].tvcnames[kk])
            }
            p = rows(psp->knots) - 1 + (isb1 ? 1 : 0) // ts1 base keeps cons col
            // basis order for ts1 base: [cons, rcs..]; params: [rcs.., cons];
            // gk is padded so gk[p] (the cons coefficient) is 0 under nocons
            if (isb1) {
                gk = th[|pos + 1 \ pos + p - 1|] \ 0
                if (C.cons) gk[p] = th[pos + p]
            }
            else gk = th[|pos + 1 \ pos + p|]
            if (needh) {
                if (shared) {
                    cEV[e] = &(_stx_basis(J(1, 1, tj), *psp, isb1))
                    lph = lph + (isb1
                        ? J(c, 1, (*cEV[e])[|1, 2 \ 1, p|] * gk[|1 \ p - 1|] +
                            (*cEV[e])[1] * gk[p])
                        : (zidx ? X[., zidx] :* ((*cEV[e]) * gk)[1]
                                : J(c, 1, ((*cEV[e]) * gk)[1])))
                }
                else {
                    cEV[e] = &(_stx_basis(J(c, 1, tj) + voff, *psp, isb1))
                    lph = lph + (isb1
                        ? (*cEV[e])[|1, 2 \ c, p|] * gk[|1 \ p - 1|] +
                            (*cEV[e])[., 1] * gk[p]
                        : (zidx ? X[., zidx] :* ((*cEV[e]) * gk)
                                : (*cEV[e]) * gk))
                }
            }
            if (needH) {
                if (shared) {
                    cQ[e] = &(_stx_basis(u, *psp, isb1))
                    z = (isb1
                        ? (*cQ[e])[|1, 2 \ G, p|] * gk[|1 \ p - 1|] +
                            (*cQ[e])[., 1] * gk[p]
                        : (*cQ[e]) * gk)               // G x 1
                    LP = LP + (zidx ? X[., zidx] * z' : J(c, 1, 1) * z')
                }
                else {
                    // node times per individual: row i = u' + voff_i
                    cQ[e] = &(_stx_basis(vec((J(c, 1, 1) * u' :+ voff)'),
                        *psp, isb1))
                    z = (isb1
                        ? (*cQ[e])[|1, 2 \ c * G, p|] * gk[|1 \ p - 1|] +
                            (*cQ[e])[., 1] * gk[p]
                        : (*cQ[e]) * gk)               // cG x 1
                    LP = LP + (zidx ? X[., zidx] :* colshape(z, G)
                                    : colshape(z, G))
                }
            }
            czi[e] = zidx
            cp[e] = p
            cshared[e] = shared
            cb1[e] = isb1
            pos = pos + (isb1 ? p - 1 + C.cons : p)
        }
    }

    if (needh) h = exp(lph)
    if (needH) {
        E = exp(LP)
        Ew = E * w
        H = half :* Ew
    }
    if (!doJ) return

    // Jacobian blocks, walked in the same design order
    if (needh) JhM = (ncov ? h :* X[|1, 1 \ c, ncov|] : J(c, 0, 0))
    if (needH) {
        W = J(c, 1, w')
        EW = E :* W
        JHM = (ncov ? H :* X[|1, 1 \ c, ncov|] : J(c, 0, 0))
    }
    for (e = 1; e <= nel; e++) {
        p = cp[e]
        zidx = czi[e]
        isb1 = cb1[e]
        if (needh) {
            if (cshared[e]) {
                B = (zidx ? (h :* X[., zidx]) * (*cEV[e]) : h * (*cEV[e]))
            }
            else {
                B = (zidx ? (h :* X[., zidx]) :* (*cEV[e]) : h :* (*cEV[e]))
            }
            JhM = (JhM, (isb1
                ? (C.cons ? (B[|1, 2 \ c, p|], B[., 1]) : B[|1, 2 \ c, p|])
                : B))
        }
        if (needH) {
            if (cshared[e]) {
                B = E * (w :* (*cQ[e]))                // c x p
                B = (zidx ? (half :* X[., zidx]) :* B : half :* B)
            }
            else {
                B = J(c, p, .)
                for (l = 1; l <= p; l++) {
                    B[., l] = half :*
                        rowsum(EW :* colshape((*cQ[e])[., l], G), 1)
                }
                if (zidx) B = X[., zidx] :* B
            }
            JHM = (JHM, (isb1
                ? (C.cons ? (B[|1, 2 \ c, p|], B[., 1]) : B[|1, 2 \ c, p|])
                : B))
        }
    }
    if (needh) Jh = JhM
    if (needH) JH = JHM
}

// sum of the quantity (and its Jacobian) over a population chunk, per time
void _stx_std_chunk(struct _stx_model scalar M, real colvector times,
    real matrix Xr, real matrix Xe, real matrix OFFr, real matrix OFFe,
    real colvector ind, real colvector wp, string scalar q,
    real colvector nd, real colvector w, real scalar doJ,
    real colvector gsum, real matrix Jsum)
{
    real colvector thr, the, hr, he, Hr, He, vals, S
    real matrix Jhr, Jhe, JHr, JHe
    real scalar m, j, pr, k, needh, needref

    m = rows(times)
    pr = M.pr
    k = cols(M.b)
    thr = M.b[(1..pr)]'
    the = M.b[((pr + 1)..k)]'
    needh   = q == "hazard" | q == "excesshazard"
    needref = !(q == "netsurv" | q == "excesshazard")

    hr = he = Hr = He = .
    Jhr = Jhe = JHr = JHe = .
    for (j = 1; j <= m; j++) {
        if (needref) {
            _stx_std_comp(M.ref, thr, times[j], Xr, OFFr, nd, w,
                needh, !needh, doJ, hr, Hr, Jhr, JHr)
        }
        _stx_std_comp(M.exc, the, times[j], Xe, OFFe, nd, w,
            needh, !needh, doJ, he, He, Jhe, JHe)

        if (q == "hazard")             vals = hr + ind :* he
        else if (q == "excesshazard")  vals = he
        else if (q == "chazard")       vals = Hr + ind :* He
        else if (q == "survival" | q == "cif" | q == "netsurv") {
            S = (q == "netsurv" ? exp(-He) : exp(-(Hr + ind :* He)))
            vals = (q == "cif" ? 1 :- S : S)
        }
        else _error(3498, "unknown quantity: " + q)

        if (hasmissing(vals)) {                     // e.g. grid time <= 0
            gsum[j] = .
            if (doJ) Jsum[j, .] = J(1, k, .)
            continue
        }
        gsum[j] = gsum[j] + sum(wp :* vals)
        if (!doJ) continue

        if (needh) {
            if (needref) {
                Jsum[|j, 1 \ j, pr|] = Jsum[|j, 1 \ j, pr|] +
                    colsum(wp :* Jhr)
                Jsum[|j, pr + 1 \ j, k|] = Jsum[|j, pr + 1 \ j, k|] +
                    colsum((wp :* ind) :* Jhe)
            }
            else {
                Jsum[|j, pr + 1 \ j, k|] = Jsum[|j, pr + 1 \ j, k|] +
                    colsum(wp :* Jhe)
            }
        }
        else if (q == "chazard") {
            Jsum[|j, 1 \ j, pr|] = Jsum[|j, 1 \ j, pr|] + colsum(wp :* JHr)
            Jsum[|j, pr + 1 \ j, k|] = Jsum[|j, pr + 1 \ j, k|] +
                colsum((wp :* ind) :* JHe)
        }
        else {                       // survival family: d exp(-H) = -S dH
            if (q == "cif") S = -S   // d cif = +S dH
            if (needref) {
                Jsum[|j, 1 \ j, pr|] = Jsum[|j, 1 \ j, pr|] +
                    colsum((-wp :* S) :* JHr)
                Jsum[|j, pr + 1 \ j, k|] = Jsum[|j, pr + 1 \ j, k|] +
                    colsum((-wp :* S) :* (ind :* JHe))
            }
            else {
                Jsum[|j, pr + 1 \ j, k|] = Jsum[|j, pr + 1 \ j, k|] +
                    colsum((-wp :* S) :* JHe)
            }
        }
    }
}

// standardised estimate over a population: averaged value + Jacobian
void _stx_standest(struct _stx_model scalar M, real colvector times,
    real matrix Xrpop, real matrix Xepop, real matrix OFFrpop,
    real matrix OFFepop, real colvector indpop, real colvector wpop,
    string scalar q, real scalar G, real scalar chunk, real scalar doJ,
    real colvector est, real matrix Jc)
{
    real matrix glm, Jsum
    real colvector nd, w, gsum
    real scalar N, W, lo, hi, m, k

    m = rows(times)
    k = cols(M.b)
    N = rows(OFFrpop)
    W = sum(wpop)
    glm = _stx_gl(G)
    nd = glm[., 1]
    w  = glm[., 2]
    gsum = J(m, 1, 0)
    Jsum = J(m, k, 0)
    for (lo = 1; lo <= N; lo = lo + chunk) {
        hi = min((lo + chunk - 1, N))
        _stx_std_chunk(M, times,
            (cols(Xrpop) ? Xrpop[|lo, 1 \ hi, .|] : J(hi - lo + 1, 0, 0)),
            (cols(Xepop) ? Xepop[|lo, 1 \ hi, .|] : J(hi - lo + 1, 0, 0)),
            OFFrpop[|lo, 1 \ hi, .|], OFFepop[|lo, 1 \ hi, .|],
            indpop[|lo \ hi|], wpop[|lo \ hi|], q, nd, w, doJ, gsum, Jsum)
    }
    est = gsum :/ W
    Jc = Jsum :/ W
}

// ========================================================================= //
// postestimation drivers -- macro contract from stexcess_p.ado
// ========================================================================= //

// results are stashed in Mata and written out by _stx_flush() once the
// wrapper has restored the data (at()/zeros run under preserve, so anything
// st_store'd before the restore would be rolled back)
void _stx_stash(real colvector est, real matrix Jc,
    struct _stx_model scalar M, string scalar transform, real scalar level)
{
    external real colvector STX_O_est, STX_O_lci, STX_O_uci
    real colvector lci, uci

    STX_O_est = est
    STX_O_lci = STX_O_uci = J(0, 1, .)
    if (st_local("_stx_ci") == "") return
    lci = .
    uci = .
    _stx_delta(est, Jc, M.V, transform, level, lci, uci)
    STX_O_lci = lci
    STX_O_uci = uci
}

void _stx_flush()
{
    external real colvector STX_O_est, STX_O_lci, STX_O_uci
    string scalar touse

    touse = st_local("_stx_touse")
    st_store(., st_local("_stx_out"), touse, STX_O_est)
    if (st_local("_stx_ci") != "") {
        st_store(., st_local("_stx_lci"), touse, STX_O_lci)
        st_store(., st_local("_stx_uci"), touse, STX_O_uci)
    }
    STX_O_est = STX_O_lci = STX_O_uci = J(0, 1, .)
}

// plain predict: observed-row covariates/indicator (the wrapper has already
// applied any at()/zeros overrides as data replacements)
void _stx_predict(real scalar level)
{
    struct _stx_model scalar M
    string scalar touse, q
    real colvector times, est, ind
    real matrix Jc, Xr, Xe, OFFr, OFFe
    real scalar doJ

    M = _stx_usemodel()
    touse = st_local("_stx_touse")
    times = st_data(., st_local("_stx_timevar"), touse)
    q = st_local("_stx_quantity")
    doJ = st_local("_stx_ci") != ""
    Xr = Xe = OFFr = OFFe = .
    ind = .
    _stx_preddata(M, touse, rows(times), _stx_needind(q),
        Xr, Xe, OFFr, OFFe, ind)
    est = .
    Jc = .
    _stx_quantity(M, times, Xr, Xe, OFFr, OFFe, ind, q, 50, doJ, est, Jc)
    _stx_stash(est, Jc, M, _stx_transform(q), level)
}

// contrast: the wrapper applies at1(), calls _stx_cside1(), restores,
// applies at2(), then _stx_cside2() combines and stashes
void _stx_cside1()
{
    external real colvector STX_C_e1
    external real matrix STX_C_J1
    struct _stx_model scalar M
    string scalar touse, q
    real colvector times, e1, ind
    real matrix J1, Xr, Xe, OFFr, OFFe
    real scalar doJ

    M = _stx_usemodel()
    touse = st_local("_stx_touse")
    times = st_data(., st_local("_stx_timevar"), touse)
    q = st_local("_stx_quantity")
    doJ = st_local("_stx_ci") != ""
    Xr = Xe = OFFr = OFFe = .
    ind = .
    e1 = .
    J1 = .
    _stx_preddata(M, touse, rows(times), _stx_needind(q),
        Xr, Xe, OFFr, OFFe, ind)
    _stx_quantity(M, times, Xr, Xe, OFFr, OFFe, ind, q, 50, doJ, e1, J1)
    STX_C_e1 = e1
    STX_C_J1 = (doJ ? J1 : J(0, 0, .))
}

void _stx_cside2(real scalar level)
{
    external real colvector STX_C_e1
    external real matrix STX_C_J1
    struct _stx_model scalar M
    string scalar touse, q, kind
    real colvector times, e2, est, ind
    real matrix J2, Jc, Xr, Xe, OFFr, OFFe
    real scalar doJ

    M = _stx_usemodel()
    touse = st_local("_stx_touse")
    times = st_data(., st_local("_stx_timevar"), touse)
    q = st_local("_stx_quantity")
    kind = st_local("_stx_kind")
    doJ = st_local("_stx_ci") != ""
    Xr = Xe = OFFr = OFFe = .
    ind = .
    e2 = .
    J2 = .
    Jc = .
    _stx_preddata(M, touse, rows(times), _stx_needind(q),
        Xr, Xe, OFFr, OFFe, ind)
    _stx_quantity(M, times, Xr, Xe, OFFr, OFFe, ind, q, 50, doJ, e2, J2)
    if (kind == "ratio") {
        est = STX_C_e1 :/ e2
        if (doJ) Jc = (STX_C_J1 :* e2 - STX_C_e1 :* J2) :/ (e2 :^ 2)
        _stx_stash(est, Jc, M, "log", level)
    }
    else {
        est = STX_C_e1 - e2
        if (doJ) Jc = STX_C_J1 - J2
        _stx_stash(est, Jc, M, "identity", level)
    }
    STX_C_e1 = J(0, 1, .)
    STX_C_J1 = J(0, 0, .)
}

// standardised predict: average over e(sample); at()/zeros overrides are
// applied across the population by the wrapper (counterfactual
// standardisation); n_nodes = 50 (40 for rmst, outer 40), chunked at 4000
void _stx_standsurv(real scalar level)
{
    struct _stx_model scalar M
    string scalar touse, pop, q
    real colvector times, est, taus, half_out, Sflat, no, wo, ind, wpop
    real matrix Jc, Xr, Xe, OFFr, OFFe, glm, U, Jo
    real scalar doJ, m, Mo, j, npop

    M = _stx_usemodel()
    touse = st_local("_stx_touse")
    pop = st_local("_stx_poptouse")
    times = st_data(., st_local("_stx_timevar"), touse)
    q = st_local("_stx_quantity")
    doJ = st_local("_stx_ci") != ""
    npop = rows(st_data(., pop, pop))
    Xr = Xe = OFFr = OFFe = .
    ind = .
    _stx_preddata(M, pop, npop,
        _stx_needind(q == "rmst" | q == "timelost" ? "survival" : q),
        Xr, Xe, OFFr, OFFe, ind)
    wpop = (M.wvar == "" ? J(npop, 1, 1) : st_data(., M.wvar, pop))
    est = .
    Jc = .

    if (q == "rmst" | q == "rmstnet" | q == "timelost") {
        // outer GL over the standardised survival curve (inner = outer = 40)
        taus = times
        m = rows(taus)
        glm = _stx_gl(40)
        no = glm[., 1]
        wo = glm[., 2]
        Mo = 40
        half_out = 0.5 :* taus
        U = half_out * (no :+ 1)'
        Sflat = .
        Jo = .
        _stx_standest(M, vec(U'), Xr, Xe, OFFr, OFFe, ind, wpop,
            (q == "rmstnet" ? "netsurv" : "survival"), 40, 4000, doJ,
            Sflat, Jo)
        est = half_out :* (colshape(Sflat, Mo) * wo)
        if (q == "timelost") est = taus - est
        if (doJ) {
            Jc = J(m, cols(M.b), .)
            for (j = 1; j <= cols(M.b); j++) {
                Jc[., j] = half_out :* (colshape(Jo[., j], Mo) * wo)
            }
            if (q == "timelost") Jc = -Jc
        }
        _stx_stash(est, Jc, M, "log", level)
        return
    }
    _stx_standest(M, times, Xr, Xe, OFFr, OFFe, ind, wpop, q, 50, 4000, doJ,
        est, Jc)
    _stx_stash(est, Jc, M, _stx_transform(q), level)
}

end
