# Changelog

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
