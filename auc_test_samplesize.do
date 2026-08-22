*===================================================================
* auc_test_samplesize.do
*
* Sample size for a ONE-SIDED test of a single AUC (C-statistic)
* against a fixed margin:
*
*     H0: AUC <= auc0        vs        Ha: AUC > auc0
*
* at one-sided type I error alpha, with power (1-beta) to detect a
* true AUC of auc1 (auc1 > auc0). auc0 is whatever fixed margin you
* are testing against (e.g. 0.5 for "better than chance", or a
* pre-specified minimum-acceptable/non-inferiority value).
*
* THIS IS NOT the same calculation as auc_samplesize.do (the earlier
* file), which sized a study for the *precision* of a CI (an
* estimation problem, no power/beta term). This file sizes a study
* for a *hypothesis test with specified power*, which is a different
* and generally larger sample size for the same alpha.
*
* THIS ALSO IS NOT a comparison of two different diagnostic tests'
* AUCs (e.g. new test vs reference test). That is a different, more
* complex problem: if the two AUCs are estimated on the SAME subjects
* (paired/correlated ROC curves) you need the covariance between the
* two AUCs (e.g. via DeLong's method); if on independent samples you
* need the variance of the difference of two independent AUCs. Ask
* if that is what you actually need and I will write that separately
* rather than adapt this file to it.
*
* METHOD
* ------
* Uses the Hanley & McNeil (1982, Radiology 143:29-36) large-sample
* variance of the empirical AUC (same formula as auc_samplesize.do):
*
*     Var(AUC) = [ Q0 + (n1-1)*Q1 + (n2-1)*Q2 ] / (n1*n2)
*     Q0 = AUC*(1-AUC)
*     Q1 = AUC/(2-AUC) - AUC^2
*     Q2 = 2*AUC^2/(1+AUC) - AUC^2
*
* The test statistic is Wald-type, Z = (AUChat - auc0)/SE, rejected
* when Z > z_(1-alpha). The rejection boundary uses the variance
* evaluated AT auc0 (Var0); power against the true value auc1 uses
* the variance evaluated AT auc1 (Var1) -- these differ because the
* Hanley-McNeil variance depends on the AUC value itself. This
* "unpooled/unequal variance under H0 vs H1" structure is the
* standard form for one-sided tests of this kind (the same structure
* used for, e.g., unpooled two-proportion non-inferiority tests), and
* leads to the sample-size equation:
*
*     z_alpha * SE0(n) + z_beta * SE1(n) = auc1 - auc0
*
* where SE0(n) = sqrt(Var(auc0; n1,n2)), SE1(n) = sqrt(Var(auc1; n1,n2)).
*
* I solved this numerically (search over integer n) rather than in
* closed form, and checked the implementation logic (in Python,
* mirroring this exact code) against several worked cases to confirm
* it (a) meets the alpha/power criterion and (b) is the minimum n
* that does so (n-1 fails). I have NOT checked it against a published
* worked example/table specific to AUC sample size -- I don't have a
* citable one on hand, so treat the numeric agreement as "internally
* consistent with the stated formula," not as "matches a known
* published number."
*
* LIMITATIONS
* -----------
* - Large-sample Wald approximation; unreliable for small n or AUC
*   near 1 (or near 0.5 with very small n), for the same reasons as
*   auc_samplesize.do.
* - Assumes auc0 and auc1 are both correctly specified a priori.
* - One-sided upper test only (Ha: AUC > auc0). If your alternative
*   runs the other way, this file as written is not set up for that.
*
* USAGE
* -----
*   . do auc_test_samplesize.do
*   . auc1test, auc0(0.5) auc1(0.75) alpha(0.025) power(0.80)
*   . auc1test, auc0(0.7) auc1(0.8) prev(0.3)
*
* OPTIONS
*   auc0(#)   null / margin AUC value being tested against, 0<auc0<1
*   auc1(#)   anticipated true AUC used for power, auc0 < auc1 < 1
*   alpha(#)  one-sided type I error, default 0.025
*   power(#)  desired power, default 0.80
*   prev(#)   case fraction n1/(n1+n2), default 0.5
*
* Returns: r(n1) r(n2) r(n) r(lhs_achieved) r(delta)
*===================================================================

capture program drop auc1test
program define auc1test, rclass
    version 14
    syntax , Auc0(real) Auc1(real) [ Alpha(real 0.025) Power(real 0.80) Prev(real 0.5) ]

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
    if `prev' <= 0 | `prev' >= 1 {
        display as error "prev() must be strictly between 0 and 1"
        exit 198
    }

    local z_a = invnormal(1 - `alpha')
    local z_b = invnormal(`power')
    local delta = `auc1' - `auc0'

    * ---- Hanley-McNeil variance components at auc0 and at auc1 ----
    local q0_0 = `auc0'*(1-`auc0')
    local q1_0 = `auc0'/(2-`auc0') - `auc0'^2
    local q2_0 = 2*`auc0'^2/(1+`auc0') - `auc0'^2

    local q0_1 = `auc1'*(1-`auc1')
    local q1_1 = `auc1'/(2-`auc1') - `auc1'^2
    local q2_1 = 2*`auc1'^2/(1+`auc1') - `auc1'^2

    * ---- closed-form (continuous-n, leading-order) starting value ----
    * approximation: Var_k(n) ~= [Q1k/(1-prev) + Q2k/prev] / n for large n
    local W0 = `q1_0'/(1-`prev') + `q2_0'/`prev'
    local W1 = `q1_1'/(1-`prev') + `q2_1'/`prev'
    local n0 = ((`z_a'*sqrt(`W0') + `z_b'*sqrt(`W1'))/`delta')^2

    * ---- integer search: increase n until the EXACT discrete        ----
    * ---- criterion  z_a*SE0(n) + z_b*SE1(n) <= delta  holds. This   ----
    * ---- is the objective success criterion for this calculation.  ----
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

    display as text "{hline 62}"
    display as text "Sample size: one-sided test of a single AUC against a margin"
    display as text "H0: AUC <= auc0   vs   Ha: AUC > auc0"
    display as text "Method: Hanley & McNeil (1982) variance, unpooled null/alt SE"
    display as text "{hline 62}"
    display as text "auc0 (margin / null)            = " as result %6.3f `auc0'
    display as text "auc1 (anticipated true AUC)      = " as result %6.3f `auc1'
    display as text "One-sided alpha                  = " as result %6.4f `alpha'
    display as text "Power                            = " as result %5.3f `power'
    display as text "Case fraction n1/(n1+n2)          = " as result %5.3f `prev'
    display as text "{hline 62}"
    display as text "Required n1 (cases)              = " as result `n1'
    display as text "Required n2 (controls)           = " as result `n2'
    display as text "Required total n                 = " as result `ntotal'
    display as text "z_a*SE0 + z_b*SE1 at this n       = " as result %7.5f `lhs'
    display as text "  (must be <= auc1-auc0 = " as result %6.4f `delta' as text ")"
    display as text "{hline 62}"

    return scalar n1 = `n1'
    return scalar n2 = `n2'
    return scalar n  = `ntotal'
    return scalar lhs_achieved = `lhs'
    return scalar delta = `delta'
end

*-------------------------------------------------------------------
* Example (comment out or edit as needed):
auc1test, auc0(0.75) auc1(0.85) alpha(0.025) power(0.80)
* return list
*-------------------------------------------------------------------
