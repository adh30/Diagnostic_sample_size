*===================================================================
* dta_samplesize.do
*
* Sample size for a diagnostic accuracy study, to estimate sensitivity
* and specificity to specified precision, across a RANGE of disease
* prevalence values.
*
* THIS IS A DIFFERENT DESIGN SCENARIO from auc_samplesize.do and
* auc_test_samplesize.do. In those files, prev() was an allocation
* ratio (n1/(n1+n2)) that YOU, the researcher, choose directly (a
* case-control-style design: you decide how many cases and how many
* controls to enrol). HERE, prev() is the disease PREVALENCE in the
* population you are sampling consecutively/cross-sectionally from --
* you do not control it, you can only anticipate a plausible value
* (or range of values, hence looping over a numlist) and plan total
* enrolment N so that, at that prevalence, you expect enough diseased
* AND enough non-diseased subjects to hit your precision targets for
* both sensitivity and specificity simultaneously.
*
* METHOD
* ------
* Buderer NMF. "Statistical methodology: I. Incorporating the
* prevalence of disease into the sample size calculation for
* sensitivity and specificity." Acad Emerg Med. 1996;3(9):895-900.
*
*     n_se = z^2 * Se*(1-Se) / me_se^2      (diseased subjects needed
*                                             to estimate Se to +/- me_se)
*     n_sp = z^2 * Sp*(1-Sp) / me_sp^2      (non-diseased subjects needed
*                                             to estimate Sp to +/- me_sp)
*     N_sensitivity = n_se / prev
*     N_specificity = n_sp / (1-prev)
*     N_total       = max(N_sensitivity, N_specificity)
*
* ROUNDING -- CORRECTED after checking against Akoglu's calculator
* -------------------------------------------------------------------
* An earlier version of this file rounded UP (ceil) at n_se/n_sp, then
* rounded UP again at N_sensitivity/N_specificity -- a double-ceiling
* that made totals too large. I verified the correct convention
* directly against Akoglu's online calculator (the tool built to
* accompany the user's-guide paper this method comes from --
* https://turkjemergmed.com/calculator, spreadsheet at
* https://development.b4bynd.com/sse4das/, "Single-Test Design, new
* diagnostic tests" tab), by entering values into it and reading its
* output:
*     Se=50%,   prev=30%, me=10%  -> N(se)=320,  N(sp)=137
*     Se=74.6%, prev=30%, me=5%   -> N(se)=971,  N(sp)=416
*     Se=85.7%, prev=30%, me=5%   -> N(se)=628
*     Sp=62.4%, prev=30%, me=5%   -> N(sp)=515
* (the calculator's single "Sens or Spec" field is used once per run;
* since n_se depends only on Se and n_sp only on Sp, running it once
* with Se and once with Sp and reading the respective output is
* equivalent to giving it two different Se/Sp values.)
*
* All six values match EXACTLY (and only) when: (1) n_se and n_sp are
* kept as continuous, unrounded numbers all the way through, with NO
* intermediate rounding, and (2) the single rounding step, at the very
* end, is round-to-nearest -- NOT ceiling. E.g. for Se=74.6%: raw n_se
* = 291.158 (never rounded); N_sensitivity = 291.158/0.30 = 970.53,
* which rounds to 971 (matches). My original ceil()-based version gave
* n_se=ceil(291.158)=292, then N_sensitivity=ceil(292/0.30)=974 --
* three too many, and NOT what the calculator you checked against
* reports. This file now reproduces the calculator's numbers exactly
* (I re-derived and re-checked all four cases above in Python before
* rewriting the Stata code).
*
* CONSEQUENCE OF ROUND-TO-NEAREST (be aware of this trade-off): unlike
* a ceiling-based n, round-to-nearest does not guarantee the achieved
* margin of error is <= your target -- about half the time it will be
* fractionally worse than requested (e.g. 320 subjects here gives an
* achieved margin very slightly ABOVE 10%, not below). The shortfall
* is always under 1 subject's worth of precision, so it is unlikely to
* matter in practice, but it is a real difference from a "round up to
* guarantee precision" convention, and I am matching Akoglu's
* calculator's convention here because that is what you are
* cross-checking against, not because round-to-nearest is more
* conservative -- it is slightly less conservative.
*
* LIMITATIONS
* -----------
* - Large-sample normal-approximation CI for a proportion (Wald-type),
*   same as the standard formula n=z^2 p(1-p)/d^2 applied twice. This
*   degrades for Se or Sp very close to 0 or 1, or for small n.
* - Assumes Se, Sp, and prevalence are all correctly anticipated a
*   priori; as with any sample-size calculation the output is only as
*   good as those assumptions.
* - Assumes consecutive/cross-sectional sampling (prevalence in your
*   sample matches the population prevalence you specify). If instead
*   you are free to choose your own cases:controls ratio (a
*   case-control design), this is not the right file -- ask if you
*   want that version instead.
*
* USAGE
* -----
*   . do dta_samplesize.do
*   . dtasize, se(0.746) sp(0.746) me(0.05) level(95) prev(0.1(0.05)0.4)
*   . dtasize, se(0.857) sp(0.624) me(0.05) me_sp(0.08) prev(0.2 0.3)
*
* OPTIONS
*   se(#)        anticipated sensitivity, 0<se<1
*   sp(#)        anticipated specificity, 0<sp<1
*   me(#)        marginal error (+/-) for BOTH Se and Sp CIs, unless
*                me_sp() is also given, > 0
*   me_sp(#)     optional separate marginal error for Sp only;
*                defaults to me() if not specified
*   level(#)     confidence level in percent, default 95
*                (equivalently, type I error alpha = 1-level/100)
*   prev(numlist) one or more anticipated disease-prevalence values,
*                each in (0,1); accepts numlist ranges, e.g.
*                prev(0.1(0.05)0.4) or prev(0.2 0.3 0.5)
*
* Returns: r(results) -- a (#prev x 6) matrix with columns
*          prev, N_sensitivity, N_specificity, N_total,
*          n_diseased_raw, n_nondiseased_raw
*===================================================================

capture program drop dtasize
program define dtasize, rclass
    version 14
    syntax , Se(real) Sp(real) ME(real) [ ME_sp(real 0) Level(real 95) PREV(string) ]

    * ---- input checks ----
    if `se' <= 0 | `se' >= 1 {
        display as error "se() must be strictly between 0 and 1"
        exit 198
    }
    if `sp' <= 0 | `sp' >= 1 {
        display as error "sp() must be strictly between 0 and 1"
        exit 198
    }
    if `me' <= 0 {
        display as error "me() must be > 0"
        exit 198
    }
    if `me_sp' == 0 {
        local me_sp = `me'
    }
    else if `me_sp' <= 0 {
        display as error "me_sp() must be > 0"
        exit 198
    }
    if `level' <= 0 | `level' >= 100 {
        display as error "level() must be strictly between 0 and 100"
        exit 198
    }
    if `"`prev'"' == "" {
        display as error "prev() is required -- specify one or more anticipated prevalence values"
        exit 198
    }

    local alpha = 1 - `level'/100
    local z = invnormal(1 - `alpha'/2)

    * ---- Buderer (1996): raw (UNROUNDED) group sizes for Se and Sp ----
    * kept as continuous numbers -- do not round here (see header note)
    local n_se_raw = (`z'^2) * `se'*(1-`se') / `me'^2
    local n_sp_raw = (`z'^2) * `sp'*(1-`sp') / `me_sp'^2

    * ---- expand the requested prevalence list and validate ----
    local plist ""
    foreach p of numlist `prev' {
        if `p' <= 0 | `p' >= 1 {
            display as error "each value in prev() must be strictly between 0 and 1 (got `p')"
            exit 198
        }
        local plist "`plist' `p'"
    }
    local nrows : word count `plist'
    matrix results = J(`nrows', 6, .)

    display as text "{hline 84}"
    display as text "Sample size for a diagnostic accuracy study (Buderer 1996)"
    display as text "Precision targets for sensitivity and specificity, across a range of prevalence"
    display as text "(rounding matches Akoglu's calculator: round-to-nearest at the final step only)"
    display as text "{hline 84}"
    display as text "Se = " as result %5.3f `se' as text "   Sp = " as result %5.3f `sp' ///
        as text "   me(Se) = " as result %5.3f `me' as text "   me(Sp) = " as result %5.3f `me_sp' ///
        as text "   level = " as result %5.1f `level' "%"
    display as text "raw n(diseased), unrounded    = " as result %9.3f `n_se_raw'
    display as text "raw n(nondiseased), unrounded = " as result %9.3f `n_sp_raw'
    display as text "{hline 84}"
    display as text %8s "prev" _col(14) "N_sens" _col(24) "N_spec" _col(34) "N_total" _col(46) "binding"
    display as text "{hline 84}"

    local row = 0
    foreach p of local plist {
        local row = `row' + 1
        local Nd = round(`n_se_raw'/`p')
        local Nn = round(`n_sp_raw'/(1-`p'))
        if `Nd' >= `Nn' {
            local N = `Nd'
            local binding "sensitivity"
        }
        else {
            local N = `Nn'
            local binding "specificity"
        }

        display as result %8.3f `p' _col(14) %5.0f `Nd' _col(24) %5.0f `Nn' _col(34) %6.0f `N' _col(46) "`binding'"

        matrix results[`row',1] = `p'
        matrix results[`row',2] = `Nd'
        matrix results[`row',3] = `Nn'
        matrix results[`row',4] = `N'
        matrix results[`row',5] = `n_se_raw'
        matrix results[`row',6] = `n_sp_raw'
    }
    display as text "{hline 84}"

    matrix colnames results = prev N_sensitivity N_specificity N_total n_diseased_raw n_nondiseased_raw
    return matrix results = results
    return scalar n_diseased_raw = `n_se_raw'
    return scalar n_nondiseased_raw = `n_sp_raw'
end

*-------------------------------------------------------------------
* Verification example (matches Akoglu's calculator exactly):
* dtasize, se(0.746) sp(0.746) me(0.05) level(95) prev(0.3)
*   -> should give N_sensitivity=971, N_specificity=416
dtasize, se(0.75) sp(0.75) me(0.1) level(95) prev(0.25(0.05)0.75)

* return list
* matrix list r(results)
*-------------------------------------------------------------------
