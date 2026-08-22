*===================================================================
* auc_samplesize.do
*
* Sample size for estimating a C-statistic (AUC of a ROC curve) to a
* specified precision (i.e. so that the two-sided (1-alpha)% CI for
* the AUC has a target half-width).
*
* METHOD
* ------
* Uses the variance of the empirical AUC derived by Hanley & McNeil
* (1982, Radiology 143:29-36), which is algebraically identical to
* the variance of the Mann-Whitney U statistic:
*
*     Var(AUC) = [ Q0 + (n1-1)*Q1 + (n2-1)*Q2 ] / (n1*n2)
*
*     Q0 = AUC*(1-AUC)
*     Q1 = AUC/(2-AUC) - AUC^2
*     Q2 = 2*AUC^2/(1+AUC) - AUC^2
*
*     n1 = number of "case"/diseased subjects
*     n2 = number of "control"/non-diseased subjects
*
* Sample size for a target half-width d (i.e. 95% CI = AUC +/- d) is
* obtained by setting Var(AUC) = (d / z)^2 and solving for n1, n2 for
* a chosen case:control allocation. This is the same variance formula
* used in Hanley & McNeil's original tables and in other software
* that does closed-form AUC precision calculations (e.g. the R
* package "presize", function prec_auc()).
*
* IMPORTANT / LIMITATIONS (please read)
* --------------------------------------
* - This is a large-sample Wald-type approximation, not an exact
*   result. It becomes unreliable for very small n or AUC close to
*   1 (or 0.5 with small n), because the normal approximation to the
*   sampling distribution of AUC breaks down and the true CI can be
*   asymmetric/truncated at 1.
* - It assumes the AUC value you supply is correct (a priori
*   assumption) - as with any sample-size calculation, results are
*   only as good as this assumption.
* - "halfwidth" here means the *margin of error*: CI = AUC +/- d. If
*   you have a target for the *full* CI width w, supply halfwidth =
*   w/2.
* - I am not aware of a built-in Stata command (as of Stata 18) that
*   performs this specific calculation (Stata's `power` and `ciwidth`
*   suites do not list an AUC/ROC/C-statistic method in their
*   documented method lists). This program implements the formula
*   directly and transparently so every step can be checked.
*
* USAGE
* -----
*   . do auc_samplesize.do        // defines the program
*   . aucsize, auc(0.75) halfwidth(0.05)
*   . aucsize, auc(0.75) halfwidth(0.05) prev(0.3) level(90)
*
* OPTIONS
*   auc(#)        hypothesised/assumed AUC (C-statistic), 0 < auc < 1
*   halfwidth(#)  desired margin of error for the CI, > 0
*                 (full CI width = 2*halfwidth)
*   prev(#)       case fraction n1/(n1+n2), 0 < prev < 1
*                 default 0.5 (equal numbers of cases and controls)
*   level(#)      confidence level in percent, default 95
*
* Returns (via `return list' after running):
*   r(n1), r(n2), r(n), r(halfwidth_achieved)
*===================================================================

capture program drop aucsize
program define aucsize, rclass
    version 14
    syntax , Auc(real) HALFwidth(real) [ Prev(real 0.5) Level(real 95) ]

    * ---- input checks ----
    if `auc' <= 0 | `auc' >= 1 {
        display as error "auc() must be strictly between 0 and 1"
        exit 198
    }
    if `halfwidth' <= 0 {
        display as error "halfwidth() must be > 0"
        exit 198
    }
    if `prev' <= 0 | `prev' >= 1 {
        display as error "prev() must be strictly between 0 and 1"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        display as error "level() must be strictly between 0 and 100"
        exit 198
    }

    local alpha = 1 - `level'/100
    local z = invnormal(1 - `alpha'/2)

    * ---- Hanley-McNeil (1982) variance components ----
    local q0 = `auc'*(1-`auc')
    local q1 = `auc'/(2-`auc') - `auc'^2
    local q2 = 2*`auc'^2/(1+`auc') - `auc'^2

    * ---- target variance implied by the desired half-width ----
    local vtarget = (`halfwidth'/`z')^2

    * ---- closed-form (continuous-n) starting value ----
    * With n1 = prev*n, n2 = (1-prev)*n, Var(AUC)=vtarget reduces to
    * the quadratic  A*n^2 + B*n + C = 0 :
    local A = `vtarget'*`prev'*(1-`prev')
    local B = -(`prev'*`q1' + (1-`prev')*`q2')
    local C = `q1' + `q2' - `q0'
    local disc = `B'^2 - 4*`A'*`C'
    if `disc' >= 0 & `A' != 0 {
        local n0 = (-`B' + sqrt(`disc')) / (2*`A')
    }
    else {
        local n0 = 4   // fallback; the search loop below will find the real answer
    }

    * ---- integer search: increase n until the *achieved* half-width  ----
    * ---- (computed with rounded, integer n1/n2) is <= the target     ----
    * ---- half-width. This is the objective success criterion.       ----
    local n = max(2, floor(`n0') - 5)
    local ok = 0
    local iter = 0
    while `ok' == 0 {
        local iter = `iter' + 1
        if `iter' > 2000000 {
            display as error "search did not converge -- check inputs"
            exit 498
        }
        local n1 = ceil(`prev'*`n')
        local n2 = `n' - `n1'
        if `n2' < 1 {
            local n2 = 1
            local n1 = `n' - 1
        }
        local var = (`q0' + (`n1'-1)*`q1' + (`n2'-1)*`q2') / (`n1'*`n2')
        local achieved = `z'*sqrt(`var')
        if `achieved' <= `halfwidth' {
            local ok = 1
        }
        else {
            local n = `n' + 1
        }
    }

    local ntotal = `n1' + `n2'

    display as text "{hline 62}"
    display as text "Sample size for a target precision of the AUC (C-statistic)"
    display as text "Method: Hanley & McNeil (1982) variance formula"
    display as text "{hline 62}"
    display as text "Assumed AUC                    = " as result %6.3f `auc'
    display as text "Confidence level                = " as result %5.1f `level' "%"
    display as text "Target half-width (+/-)         = " as result %6.4f `halfwidth'
    display as text "Case fraction n1/(n1+n2)         = " as result %5.3f `prev'
    display as text "{hline 62}"
    display as text "Required n1 (cases)             = " as result `n1'
    display as text "Required n2 (controls)          = " as result `n2'
    display as text "Required total n                = " as result `ntotal'
    display as text "Achieved half-width at this n    = " as result %7.5f `achieved'
    display as text "{hline 62}"

    return scalar n1 = `n1'
    return scalar n2 = `n2'
    return scalar n  = `ntotal'
    return scalar halfwidth_achieved = `achieved'
end

*-------------------------------------------------------------------
* Example (comment out or edit as needed):
* aucsize, auc(0.75) halfwidth(0.05)
* return list
*-------------------------------------------------------------------
