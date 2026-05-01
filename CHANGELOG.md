# Changelog

## v1.3.0

### Bug fix: z-score sign opposite to beta

The internal z-score used for weighted combination had its sign flipped relative
to the corresponding beta. The function `z_transform_fess()` was computing
`Phi^{-1}(p/2) * sign(beta)`, which returns the lower-tail quantile — always
negative — and therefore the opposite sign from `beta`. The correct formula is
`Phi^{-1}(1 - p/2) * sign(beta)` (upper-tail quantile). The displayed `z`
column in the output was therefore sign-inconsistent with `beta` for all
meta-analysed variants.

P-values were unaffected (computed from `|z|`), but the `z` output column and
any downstream use of it were wrong.

Credit: Brady Ryan (University of Michigan) identified this bug.

### Bug fix: negative correlation estimates clamped to zero

Tetrachoric correlation estimates for near-zero overlap could come out slightly
negative due to sampling noise. A negative off-diagonal entry in the correlation
matrix makes `zse < 1`, producing anti-conservative corrected p-values, and can
also make the variance-covariance matrix for beta non-positive-definite.

Negative off-diagonal correlations are now clamped to 0 by default. The original
behaviour can be restored with `--no-cap-correlations`.

### Build: upgraded to C++14

The Makefile now uses `-std=c++14`. Boost ≥ 1.82 emits a deprecation warning
when compiled under C++11; C++14 silences it. No source changes were required.

### Output column names clarified

The output column formerly named `p` is now `p_corrected`, and the column
formerly named `p_fess` is now `p_stouffer`. The README descriptions of these
columns have also been corrected: `p_stouffer` is a Stouffer sample-size
weighted Z-score (not inverse-variance), and `beta`/`se` are derived from the
IVW method (paired with `p_wald`). Numeric output is unchanged.

## v1.2.0

### Breaking change: `--use-beta-sign` removed

The sign-of-beta dichotomization method (Southam et al. 2017) is now always
used for correlation matrix estimation. The old p-value transform method
(`b_transform(p) = Phi^{-1}(1-p) <= 0`) severely underestimated the
inter-study correlation and has been removed along with the `--use-beta-sign`
flag that was introduced in v1.1.0.

### Bug fix: missing factor of 2 in variance formula

`produit_reciproque_asymetrique()` computes `sum_{k<l} w_k w_l rho_{kl}` by
iterating over the upper triangle of the symmetric correlation matrix. However,
the full off-diagonal sum is `2 * sum_{k<l} w_k w_l rho_{kl}`. The missing
factor of 2 caused the corrected z-score standard error and the beta standard
error to be underestimated, inflating type 1 error rates under sample overlap.

Credit: Brady Ryan (University of Michigan) identified this bug.

### Bug fix: Digby tetrachoric exponent

The Digby approximation for tetrachoric correlation used an exponent of 0.75,
but the correct value is pi/4 (0.7854). This caused systematic negative bias
in correlation estimates (e.g. bias of -0.02 at rho=0.5 with the old exponent
vs -0.003 with pi/4).

The default is now pi/4. A new `--digby-exponent` CLI option allows overriding
this value (e.g. `--digby-exponent 0.75` to restore the old behaviour).

### Code cleanup

Removed ~300 lines of dead code: commented-out legacy meta-analysis body,
2-study prototype with undefined `tetrachoric()` function call, and the
now-unused `b_transform()` and `add_minimum()` functions.

## v1.1.0

### Bug fix: dichotomization method for correlation estimation

The original implementation used `b_transform(p) = Phi^{-1}(1-p) <= 0` to dichotomize
p-values before computing tetrachoric correlation between studies. This ignores the
direction of effect and severely underestimates the inter-study correlation,
leading to undercorrection for sample overlap.

The paper (Southam et al. 2017) describes a different procedure:
`z_k = Phi^{-1}(p_k/2) * sgn(beta_k)`, then dichotomize on the sign of the
resulting z-score. This method correctly recovers the inter-study correlation
under the null hypothesis.

**Null simulation results** (100,000 variants under the null):

| True rho | Old method rho-hat | New method rho-hat | Old lambda_GC | New lambda_GC |
|----------|--------------------|--------------------|---------------|---------------|
| 0.0      | -0.005             | 0.006              | 1.017         | 1.011         |
| 0.1      | 0.008              | 0.097              | 1.090         | 1.042         |
| 0.3      | 0.054              | 0.284              | 1.258         | 1.129         |
| 0.5      | 0.159              | 0.475              | 1.375         | 1.181         |
| 0.7      | 0.360              | 0.670              | 1.401         | 1.185         |

The new method is enabled with `--use-beta-sign` and will become the default in a
future version. The old method remains the default for backwards compatibility.

Credit: Brady Ryan (University of Michigan) identified the discrepancy between the
paper and the implementation.

### Bug fix: allele harmonization in first pass

The first pass (correlation matrix estimation) did not harmonize alleles between
input files. When effect alleles were swapped between studies, the sign of beta
was incorrect, producing wrong tetrachoric correlation estimates. The first pass
now detects allele flips and negates beta accordingly, mirroring the allele
harmonization already present in the second pass (meta-analysis).

### Bug fix: missing return in `initialise()`

The `initialise()` function (`int inline initialise(int, char**)`) had no `return`
statement. This is undefined behavior in C++ and caused crashes (double-free) on
some compilers with `-O3`.

### Improvement: tiered arbitrary precision arithmetic

The previous implementation had a single fallback from `long double` to
`cpp_dec_float<200>`, which overflowed for p-values below ~1e-200.

All four transform functions (`b_transform`, `z_transform` x2, `z_transform_fess`)
now use a 4-tier precision cascade:

```
long double -> cpp_dec_float<200> -> cpp_dec_float<500> -> cpp_dec_float<1000> -> safe fallback
```

The final fallback caps z-scores at +/-38 (equivalent to p ~ 1e-315), which is
reached only for the most extreme p-values in GWAS (e.g. GPC3 height locus at
p = 2.8e-1135).

### Build changes

- Makefile now uses system Boost (`/usr/include`, `lib/`) instead of hardcoded
  Sanger cluster paths. Original paths preserved as comments.
- Removed `-static` flag (dynamic linking). For static builds, pass
  `CFLAGS="-O3 -std=c++11 -static"` to make.
