*===================================================================
* auc_test_samplesize.do
*
* Sample size for a ONE-SIDED test of a single AUC (C-statistic)
* against a fixed margin:
*
*     H0: AUC <= auc0        vs        Ha: AUC > auc0
*
* at one-sided type I error alpha, with power (1-beta) to detect a
* true AUC of auc1 (auc1 > auc0), computed across a RANGE of case
* fractions prev = n1/(n1+n2) rather than a single value -- e.g. to
* see how total n trades off against the cases:controls split.
*
* THIS IS NOT the same calculation as auc_samplesize.do (a separate
* file), which sizes a study for the *precision* of a CI (an
* estimation problem, no power/beta term).
*
* THIS ALSO IS NOT a comparison of two different diagnostic tests'
* AUCs (e.g. new test vs reference test) -- that needs the covariance
* between two AUCs (paired, DeLong-type) or the variance of a
* difference of two independent AUCs, which is a different file. Ask
* if that's what you actually need.
*
* METHOD
* ------
* Hanley & McNeil (1982, Radiology 143:29-36) large-sample variance
* of the empirical AUC:
*
*     Var(AUC) = [ Q0 + (n1-1)*Q1 + (n2-1)*Q2 ] / (n1*n2)
*     Q0 = AUC*(1-AUC),  Q1 = AUC/(2-AUC) - AUC^2,  Q2 = 2*AUC^2/(1+AUC) - AUC^2
*
* Wald test Z = (AUChat - auc0)/SE, rejected when Z > z_(1-alpha).
* The rejection boundary uses the variance evaluated AT auc0 (Var0);
* power against the true value auc1 uses the variance evaluated AT
* auc1 (Var1) -- these differ because the Hanley-McNeil variance
* depends on the AUC value itself. This "unpooled null/alt variance"
* structure (same as used for, e.g., unpooled two-proportion
* non-inferiority tests) gives the sample-size equation, solved here
* separately for each case fraction in the requested range:
*
*     z_alpha * SE0(n) + z_beta * SE1(n) = auc1 - auc0
*
* For each prev value, n is found by exact integer search (increase
* n until the discrete criterion above holds; n-1 always fails it --
* this was checked in a Python mirror of this exact loop logic before
* writing the Stata version). I have NOT checked these numbers
* against a published worked example/table -- I don't have a
* citable one -- so treat this as internally consistent with the
* stated formula, not as validated against an external source.
*
* LIMITATIONS
* -----------
* - Large-sample Wald approximation; unreliable for small n or AUC
*   near 1 (or near 0.5 with very small n).
* - Assumes auc0 and auc1 are both correctly specified a priori.
* - One-sided upper test only (Ha: AUC > auc0).
*
* USAGE
* -----
*   . do auc_test_samplesize.do
*   . auc1test, auc0(0.5) auc1(0.75) alpha(0.025) power(0.80) prev(0.1(0.1)0.9)
*   . auc1test, auc0(0.7) auc1(0.8) prev(0.2 0.3 0.5)
*   . auc1test, auc0(0.5) auc1(0.75)                       // prev defaults to 0.5 only
*
* OPTIONS
*   auc0(#)    null / margin AUC value being tested against, 0<auc0<1
*   auc1(#)    anticipated true AUC used for power, auc0 < auc1 < 1
*   alpha(#)   one-sided type I error, default 0.025
*   power(#)   desired power, default 0.80
*   prev(numlist) one or more case fractions n1/(n1+n2), each in (0,1);
*              accepts Stata numlist syntax including ranges, e.g.
*              prev(0.1(0.1)0.9) or prev(0.2 0.3 0.5); default "0.5"
*
* Returns: r(results) -- a (#prev x 5) matrix with columns
*          prev, n1, n2, n, lhs_achieved (see display header for lhs)
*===================================================================

capture program drop auc1test
program define auc1test, rclass
    version 14
    syntax , Auc0(real) Auc1(real) [ Alpha(real 0.025) Power(real 0.80) PREV(string) ]

    * ---- input checks ----
    if `auc0' <= 0 | `auc0' >= 1 {
        display as error "auc0() must be strictly between 0 and 1"
        exit 198
    }
    if `auc1' <= 0 | `auc1' >= 1 {
        display as error "auc1() must be strictly between 0 and 1"
        exit 198
    }
    if `auc1' <= `auc0' {
        display as error "auc1() must be greater than auc0() (this file tests Ha: AUC > auc0)"
        exit 198
    }
    if `alpha' <= 0 | `alpha' >= 1 {
        display as error "alpha() must be strictly between 0 and 1"
        exit 198
    }
    if `power' <= 0 | `power' >= 1 {
        display as error "power() must be strictly between 0 and 1"
        exit 198
    }
    if `"`prev'"' == "" {
        local prev "0.5"
    }

    local z_a = invnormal(1 - `alpha')
    local z_b = invnormal(`power')
    local delta = `auc1' - `auc0'

    * ---- Hanley-McNeil variance components at auc0 and at auc1 ----
    * (these depend only on auc0/auc1, so computed once, outside the
    *  loop over case fractions)
    local q0_0 = `auc0'*(1-`auc0')
    local q1_0 = `auc0'/(2-`auc0') - `auc0'^2
    local q2_0 = 2*`auc0'^2/(1+`auc0') - `auc0'^2

    local q0_1 = `auc1'*(1-`auc1')
    local q1_1 = `auc1'/(2-`auc1') - `auc1'^2
    local q2_1 = 2*`auc1'^2/(1+`auc1') - `auc1'^2

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

    display as text "{hline 74}"
    display as text "Sample size across a range of case fractions (n1/(n1+n2))"
    display as text "H0: AUC <= auc0  vs  Ha: AUC > auc0   (Hanley-McNeil variance)"
    display as text "{hline 74}"
    display as text "auc0 = " as result %5.3f `auc0' as text "   auc1 = " as result %5.3f `auc1' ///
        as text "   one-sided alpha = " as result %6.4f `alpha' as text "   power = " as result %5.3f `power'
    display as text "{hline 74}"
    display as text %8s "prev" _col(14) "n1" _col(24) "n2" _col(34) "n" _col(48) "z_a*SE0+z_b*SE1"
    display as text "{hline 74}"

    local row = 0
    foreach p of local plist {
        local row = `row' + 1

        * closed-form (continuous-n, leading-order) starting value
        local W0 = `q1_0'/(1-`p') + `q2_0'/`p'
        local W1 = `q1_1'/(1-`p') + `q2_1'/`p'
        local n0 = ((`z_a'*sqrt(`W0') + `z_b'*sqrt(`W1'))/`delta')^2

        * exact integer search: this is the objective success
        * criterion -- smallest n at this case fraction such that
        * z_a*SE0(n) + z_b*SE1(n) <= delta
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
            local var0 = (`q0_0' + (`n1'-1)*`q1_0' + (`n2'-1)*`q2_0') / (`n1'*`n2')
            local var1 = (`q0_1' + (`n1'-1)*`q1_1' + (`n2'-1)*`q2_1') / (`n1'*`n2')
            local lhs = `z_a'*sqrt(`var0') + `z_b'*sqrt(`var1')
            if `lhs' <= `delta' {
                local ok = 1
            }
            else {
                local n = `n' + 1
            }
        }
        local ntotal = `n1' + `n2'

        display as result %8.3f `p' _col(14) %4.0f `n1' _col(24) %4.0f `n2' _col(34) %5.0f `ntotal' _col(46) %9.5f `lhs'

        matrix results[`row',1] = `p'
        matrix results[`row',2] = `n1'
        matrix results[`row',3] = `n2'
        matrix results[`row',4] = `ntotal'
        matrix results[`row',5] = `lhs'
    }
    display as text "{hline 74}"

    matrix colnames results = prev n1 n2 n lhs_achieved
    return matrix results = results
    return scalar delta = `delta'
end

*-------------------------------------------------------------------
* Example (comment out or edit as needed):
* auc1test, auc0(0.5) auc1(0.75) alpha(0.025) power(0.80) prev(0.1(0.1)0.9)
* return list
* matrix list r(results)
*-------------------------------------------------------------------
