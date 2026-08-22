*===================================================================
* auc_samplesize.do
*
* Sample size for estimating a C-statistic (AUC of a ROC curve) to a
* specified precision (i.e. so that the two-sided (1-alpha)% CI for
* the AUC has a target half-width), computed across a RANGE of case
* fractions prev = n1/(n1+n2) rather than a single value -- e.g. to
* see how total n trades off against the cases:controls split.
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
* - For each case fraction, n is found by exact integer search
*   (increase n until the discrete achieved half-width is <= the
*   target; n-1 always fails) -- this loop logic is unchanged from
*   the single-prev version, just run once per requested prev value.
*
* USAGE
* -----
*   . do auc_samplesize.do        // defines the program
*   . aucsize, auc(0.75) halfwidth(0.05)
*   . aucsize, auc(0.75) halfwidth(0.05) prev(0.1(0.1)0.9) level(90)
*   . aucsize, auc(0.75) halfwidth(0.05) prev(0.2 0.3 0.5)
*
* OPTIONS
*   auc(#)        hypothesised/assumed AUC (C-statistic), 0 < auc < 1
*   halfwidth(#)  desired margin of error for the CI, > 0
*                 (full CI width = 2*halfwidth)
*   prev(numlist) one or more case fractions n1/(n1+n2), each in (0,1);
*                 accepts Stata numlist syntax including ranges, e.g.
*                 prev(0.1(0.1)0.9) or prev(0.2 0.3 0.5); default "0.5"
*   level(#)      confidence level in percent, default 95
*
* Returns: r(results) -- a (#prev x 5) matrix with columns
*          prev, n1, n2, n, halfwidth_achieved
*===================================================================

capture program drop aucsize
program define aucsize, rclass
    version 14
    syntax , Auc(real) HALFwidth(real) [ PREV(string) Level(real 95) ]

    * ---- input checks ----
    if `auc' <= 0 | `auc' >= 1 {
        display as error "auc() must be strictly between 0 and 1"
        exit 198
    }
    if `halfwidth' <= 0 {
        display as error "halfwidth() must be > 0"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        display as error "level() must be strictly between 0 and 100"
        exit 198
    }
    if `"`prev'"' == "" {
        local prev "0.5"
    }

    local alpha = 1 - `level'/100
    local z = invnormal(1 - `alpha'/2)

    * ---- Hanley-McNeil (1982) variance components ----
    * (depend only on auc, so computed once, outside the loop over
    *  case fractions)
    local q0 = `auc'*(1-`auc')
    local q1 = `auc'/(2-`auc') - `auc'^2
    local q2 = 2*`auc'^2/(1+`auc') - `auc'^2

    * ---- target variance implied by the desired half-width ----
    local vtarget = (`halfwidth'/`z')^2

    * ---- expand the requested case-fraction list and validate ----
    local plist ""
    foreach p of numlist `prev' {
        if `p' <= 0 | `p' >= 1 {
            display as error "each value in prev() must be strictly between 0 and 1 (got `p')"
            exit 198
        }
        local plist "`plist' `p'"
    }
    local nrows : word count `plist'
    matrix results = J(`nrows', 5, .)

    display as text "{hline 66}"
    display as text "Sample size across a range of case fractions (n1/(n1+n2))"
    display as text "Target precision for the AUC (Hanley & McNeil 1982 variance)"
    display as text "{hline 66}"
    display as text "auc = " as result %5.3f `auc' as text "   halfwidth = " as result %6.4f `halfwidth' ///
        as text "   level = " as result %5.1f `level' "%"
    display as text "{hline 66}"
    display as text %8s "prev" _col(14) "n1" _col(24) "n2" _col(34) "n" _col(46) "achieved halfwidth"
    display as text "{hline 66}"

    local row = 0
    foreach p of local plist {
        local row = `row' + 1

        * closed-form (continuous-n) starting value: with n1=p*n,
        * n2=(1-p)*n, Var(AUC)=vtarget reduces to A*n^2+B*n+C=0
        local A = `vtarget'*`p'*(1-`p')
        local B = -(`p'*`q1' + (1-`p')*`q2')
        local C = `q1' + `q2' - `q0'
        local disc = `B'^2 - 4*`A'*`C'
        if `disc' >= 0 & `A' != 0 {
            local n0 = (-`B' + sqrt(`disc')) / (2*`A')
        }
        else {
            local n0 = 4   // fallback; the search loop below will find the real answer
        }

        * exact integer search: this is the objective success
        * criterion -- smallest n at this case fraction such that the
        * achieved half-width (integer n1,n2) is <= the target
        local n = max(2, floor(`n0') - 5)
        local ok = 0
        local iter = 0
        while `ok' == 0 {
            local iter = `iter' + 1
            if `iter' > 2000000 {
                display as error "search did not converge for prev=`p' -- check inputs"
                exit 498
            }
            local n1 = ceil(`p'*`n')
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

        display as result %8.3f `p' _col(14) %4.0f `n1' _col(24) %4.0f `n2' _col(34) %5.0f `ntotal' _col(48) %7.5f `achieved'

        matrix results[`row',1] = `p'
        matrix results[`row',2] = `n1'
        matrix results[`row',3] = `n2'
        matrix results[`row',4] = `ntotal'
        matrix results[`row',5] = `achieved'
    }
    display as text "{hline 66}"

    matrix colnames results = prev n1 n2 n halfwidth_achieved
    return matrix results = results
end

*-------------------------------------------------------------------
* Example (comment out or edit as needed):
* aucsize, auc(0.75) halfwidth(0.05) prev(0.1(0.1)0.9)
* return list
* matrix list r(results)
*-------------------------------------------------------------------
