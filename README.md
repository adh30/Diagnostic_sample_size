# Diagnostic_sample_size

Stata `.do` files implementing sample-size calculations for diagnostic-accuracy and AUC (C-statistic) study designs. Each file is self-contained: running it with `do` defines a Stata program (command) that can then be called with the options described in that file's own header comments. Currently in development - various comments to self in readme.md

Repository: https://github.com/adh30/Diagnostic_sample_size

## Contents

| File | Command defined | Purpose |
|---|---|---|
| `auc_samplesize.do` | `aucsize` | Sample size to estimate a single AUC to a target confidence-interval half-width, for one specified case fraction (`prev`). |
| `auc_samplesize_1.do` | `aucsize` | Same as above, but computed across a range/list of case fractions (`prev` accepts a numlist) instead of a single value; returns a matrix of results instead of scalars. |
| `auc_test_samplesize.do` | `auc1test` | Sample size for a one-sided test of a single AUC against a fixed margin (H0: AUC ≤ `auc0` vs Ha: AUC > `auc0`), for one specified case fraction. |
| `auc_test_samplesize_1.do` | `auc1test` | Same test as above, computed across a range/list of case fractions. |
| `dta_samplesize_2.do` | `dtasize` | Sample size for a diagnostic accuracy study (sensitivity and specificity) across a range of anticipated disease-prevalence values. Per-group sample size uses the Wald (normal-approximation) formula, Buderer (1996) prevalence adjustment, with a single round-to-nearest at the final step. |
| `dta_samplesize_4.do` | `dtasize` (plus helper `_cp_minn`) | Same design scenario as `dta_samplesize_2.do`, but the per-group sample size for Se/Sp is found by exact integer search against the Clopper-Pearson (1934) exact binomial confidence interval, then combined via Buderer (1996) and rounded up (ceiling) throughout. Deliberately more conservative than `dta_samplesize_2.do`. |

`aucsize` and the two `dtasize` files are two different designs, not two versions of the same calculation: in the `auc_*.do` files, `prev` is an allocation ratio you choose (a case-control-style design — you decide how many cases and controls to enrol). In the `dta_*.do` files, `prev` is the disease prevalence in the population you sample from — you can only anticipate it, not choose it.

## Requirements

- `auc_samplesize.do`, `auc_samplesize_1.do`, `auc_test_samplesize.do`, `auc_test_samplesize_1.do`, and `dta_samplesize_2.do` each declare `version 14`.
- `dta_samplesize_4.do` declares `version 10` for both `_cp_minn` and `dtasize`, and uses the `invibeta()`/`ibeta()` scalar functions. I did not check which Stata version first introduced `invibeta()`, so if you are running Stata older than 14, confirm that function is available in your version before relying on this file.

## Usage

Each file must be run with `do` to define its command before that command can be called.

```stata
. do auc_samplesize.do
. aucsize, auc(0.75) halfwidth(0.05)
. aucsize, auc(0.75) halfwidth(0.05) prev(0.3) level(90)

. do auc_samplesize_1.do
. aucsize, auc(0.75) halfwidth(0.05) prev(0.1(0.1)0.9)
. return list
. matrix list r(results)

. do auc_test_samplesize.do
. auc1test, auc0(0.5) auc1(0.75) alpha(0.025) power(0.80)

. do auc_test_samplesize_1.do
. auc1test, auc0(0.5) auc1(0.75) alpha(0.025) power(0.80) prev(0.1(0.1)0.9)

. do dta_samplesize_2.do
. dtasize, se(0.746) sp(0.746) me(0.05) level(95) prev(0.1(0.05)0.4)

. do dta_samplesize_4.do
. dtasize, se(0.746) sp(0.746) me(0.05) level(95) prev(0.1(0.05)0.4)
```

Full option lists are documented in the header comment of each file (`OPTIONS` section) and are not repeated here to avoid the two copies drifting out of sync.

## Methods and references

- Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology. 1982;143(1):29-36. — variance of the empirical AUC, used in `auc_samplesize.do`, `auc_samplesize_1.do`, `auc_test_samplesize.do`, `auc_test_samplesize_1.do`.
- Buderer NMF. Statistical methodology: I. Incorporating the prevalence of disease into the sample size calculation for sensitivity and specificity. Acad Emerg Med. 1996;3(9):895-900. — used in both `dta_samplesize_2.do` and `dta_samplesize_4.do`.
- Clopper CJ, Pearson ES. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika. 1934;26(4):404-413. — exact binomial CI used for the per-group sample size in `dta_samplesize_4.do`.

## Notes on the files (as found)

- **Same command name reused across paired files.** `aucsize` is defined in both `auc_samplesize.do` and `auc_samplesize_1.do`; `auc1test` in both `auc_test_samplesize.do` and `auc_test_samplesize_1.do`; `dtasize` in both `dta_samplesize_2.do` and `dta_samplesize_4.do`. Each program starts with `capture program drop <name>`, so sourcing both files in the same session leaves only the version loaded last in effect.
- **Active (uncommented) example calls.** `auc_test_samplesize.do` (line 196), `dta_samplesize_2.do` (line 219), and `dta_samplesize_4.do` (line 292) each end with an example call to their command that is not commented out, so it executes automatically when the file is run with `do`. The equivalent example blocks in `auc_samplesize.do`, `auc_samplesize_1.do`, and `auc_test_samplesize_1.do` are fully commented out.
- **Header comment in `dta_samplesize_2.do` and `dta_samplesize_4.do` says `dta_samplesize.do`.** The comment block at the top of each file begins `* dta_samplesize.do`, which does not match either file's actual name on disk.
- **Trailing lines in `dta_samplesize_4.do` (lines 290-292).** The last few lines before the closing rule read:
  ```
   * Example:
  . * dtasize, se(0.746) sp(0.746) me(0.05) level(95) prev(0.1(0.05)0.4)
  . dtasize, se(0.85) sp(0.85) me(0.10) level(95) prev(0.25(0.05)0.75)
  ```
  Two of these lines carry a leading `. ` before the content (consistent with text copied from a Stata results log, where each command is echoed with a `.` prompt), and one comment line has a leading space before its `*`. I have not verified how Stata's do-file parser handles a line beginning with `.` in this position, so I am flagging it rather than asserting it does or doesn't cause a problem — worth checking if you run this file as-is.
- Each file's own header comment documents its method, options, return values, and stated limitations (large-sample/Wald approximation validity, a priori assumptions, one-sided vs two-sided scope, etc.) in more detail than is repeated here.
