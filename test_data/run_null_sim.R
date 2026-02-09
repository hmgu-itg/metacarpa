library(data.table)

# ============================================================
# Null simulation: compare current vs beta-sign method
# ============================================================

METACARPA_BIN <- "../src/metacarpa"
LIB_PATH <- "../src/lib"

lambda_gc <- function(p) {
  p <- as.numeric(p)
  p <- p[!is.na(p) & p > 0 & p < 1]
  chi2 <- qnorm(p / 2)^2
  median(chi2) / qchisq(0.5, 1)
}

simulate_gwas <- function(n_variants, rho, n1 = 10000, n2 = 10000, seed = 42) {
  set.seed(seed)
  Sigma <- matrix(c(1, rho, rho, 1), 2, 2)
  z <- MASS::mvrnorm(n_variants, mu = c(0, 0), Sigma = Sigma)

  se1 <- rep(1 / sqrt(n1), n_variants)
  se2 <- rep(1 / sqrt(n2), n_variants)
  beta1 <- z[, 1] * se1
  beta2 <- z[, 2] * se2
  p1 <- 2 * pnorm(-abs(z[, 1]))
  p2 <- 2 * pnorm(-abs(z[, 2]))

  chrs <- rep(1:22, length.out = n_variants)
  pos  <- ave(seq_len(n_variants), chrs, FUN = function(x) seq_along(x) * 1000)

  common <- data.table(
    SNPID = paste0(chrs, ":", pos, ":A:G"),
    RSID  = paste0("rs", seq_len(n_variants)),
    CHR   = chrs,
    POS   = as.integer(pos)
  )

  study1 <- cbind(common, data.table(
    EFFECT_ALLELE = "A", OTHER_ALLELE = "G",
    EFFECT_ALLELE_FREQ = 0.3,
    BETA = beta1, SE = se1, P = p1, N = n1
  ))

  study2 <- cbind(common, data.table(
    EFFECT_ALLELE = "A", OTHER_ALLELE = "G",
    EFFECT_ALLELE_FREQ = 0.3,
    BETA = beta2, SE = se2, P = p2, N = n2
  ))

  setorder(study1, CHR, POS)
  setorder(study2, CHR, POS)

  list(study1 = study1, study2 = study2, z = z, rho_true = rho)
}

run_metacarpa <- function(sim, method = "current") {
  tmpdir <- tempdir()
  f1 <- file.path(tmpdir, "sim_study1.tsv")
  f2 <- file.path(tmpdir, "sim_study2.tsv")
  fout <- file.path(tmpdir, paste0("sim_output_", method, ".tsv"))

  fwrite(sim$study1, f1, sep = "\t")
  fwrite(sim$study2, f2, sep = "\t")

  flag <- if (method == "betasign") "--use-beta-sign" else ""

  cmd <- sprintf(
    "LD_LIBRARY_PATH=%s %s -I %s -I %s -O %s -t '\t' -c 3 -q 4 -u 5 -v 6 -a 7 -b 8 -s 9 -p 10 -n 11 -i 2 %s 2>&1",
    LIB_PATH, METACARPA_BIN, f1, f2, fout, flag
  )

  output <- system(cmd, intern = TRUE)

  mat_files <- list.files(tmpdir, pattern = paste0("sim_output_", method, ".*matrix\\.txt$"),
                          full.names = TRUE)
  rho_est <- NA
  if (length(mat_files) > 0) {
    mat_lines <- readLines(mat_files[1])
    rho_est <- as.numeric(mat_lines[3])
    file.remove(mat_files)
  }

  result <- NULL
  if (file.exists(fout)) {
    result <- fread(fout)
    pcols <- c("p_wald", "p_corrected", "p_stouffer")
    for (col in pcols) {
      if (col %in% names(result)) {
        result[, (col) := as.numeric(get(col))]
      }
    }
    file.remove(fout)
  }

  list(rho_est = rho_est, result = result, log = output)
}

# ============================================================
# Run simulations
# ============================================================
rho_values <- c(0.0, 0.1, 0.2, 0.3, 0.5, 0.7)
n_variants <- 100000

sim_results <- data.table(
  rho_true = numeric(),
  rho_current = numeric(),
  rho_betasign = numeric(),
  lambda_wald_current = numeric(),
  lambda_wald_betasign = numeric(),
  lambda_stouffer = numeric(),
  t1e_wald_current = numeric(),
  t1e_wald_betasign = numeric(),
  t1e_stouffer = numeric()
)

for (rho in rho_values) {
  cat(sprintf("\n=== rho = %.1f ===\n", rho))
  sim <- simulate_gwas(n_variants, rho, seed = 123)

  res_cur <- run_metacarpa(sim, "current")
  res_bs  <- run_metacarpa(sim, "betasign")

  cat(sprintf("  Estimated rho: current = %.4f, beta-sign = %.4f\n",
              res_cur$rho_est, res_bs$rho_est))

  lam_cur <- lambda_gc(res_cur$result$p_wald)
  lam_bs  <- lambda_gc(res_bs$result$p_wald)
  lam_unc <- lambda_gc(res_cur$result$p_stouffer)

  t1e_cur <- mean(res_cur$result$p_wald < 0.05, na.rm = TRUE)
  t1e_bs  <- mean(res_bs$result$p_wald < 0.05, na.rm = TRUE)
  t1e_unc <- mean(as.numeric(res_cur$result$p_stouffer) < 0.05, na.rm = TRUE)

  cat(sprintf("  Lambda GC (p_wald): current = %.3f, beta-sign = %.3f, uncorrected = %.3f\n",
              lam_cur, lam_bs, lam_unc))
  cat(sprintf("  Type 1 error (p_wald < 0.05): current = %.4f, beta-sign = %.4f, uncorrected = %.4f\n",
              t1e_cur, t1e_bs, t1e_unc))

  sim_results <- rbind(sim_results, data.table(
    rho_true = rho,
    rho_current = res_cur$rho_est,
    rho_betasign = res_bs$rho_est,
    lambda_wald_current = lam_cur,
    lambda_wald_betasign = lam_bs,
    lambda_stouffer = lam_unc,
    t1e_wald_current = t1e_cur,
    t1e_wald_betasign = t1e_bs,
    t1e_stouffer = t1e_unc
  ))
}

cat("\n\n==========================================\n")
cat("        SIMULATION RESULTS SUMMARY        \n")
cat("==========================================\n\n")
print(sim_results)

# Save results
fwrite(sim_results, "null_sim_results.csv")
cat("\nResults saved to null_sim_results.csv\n")

# ============================================================
# Generate plots
# ============================================================
pdf("null_sim_plots.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))

# 1. Estimated vs true rho
plot(sim_results$rho_true, sim_results$rho_current, type = "b",
     pch = 19, col = "#D55E00", lwd = 2,
     xlab = expression("True " * rho), ylab = expression("Estimated " * rho),
     main = "Correlation recovery",
     xlim = c(0, 0.8), ylim = c(0, 0.8))
lines(sim_results$rho_true, sim_results$rho_betasign, type = "b",
      pch = 17, col = "#0072B2", lwd = 2)
abline(0, 1, col = "grey40", lty = 2)
legend("topleft", legend = c("Current", "Beta-sign"),
       col = c("#D55E00", "#0072B2"), pch = c(19, 17), lwd = 2, bty = "n")

# 2. Lambda GC
yr <- range(c(sim_results$lambda_wald_current, sim_results$lambda_wald_betasign,
              sim_results$lambda_stouffer), na.rm = TRUE)
plot(sim_results$rho_true, sim_results$lambda_wald_current, type = "b",
     pch = 19, col = "#D55E00", lwd = 2,
     xlab = expression("True " * rho), ylab = expression("Lambda"[GC]),
     main = "Genomic inflation",
     xlim = c(0, 0.8), ylim = yr)
lines(sim_results$rho_true, sim_results$lambda_wald_betasign, type = "b",
      pch = 17, col = "#0072B2", lwd = 2)
lines(sim_results$rho_true, sim_results$lambda_stouffer, type = "b",
      pch = 15, col = "grey60", lwd = 2)
abline(h = 1, col = "grey40", lty = 2)
legend("topleft", legend = c("Current", "Beta-sign", "Uncorrected"),
       col = c("#D55E00", "#0072B2", "grey60"), pch = c(19, 17, 15), lwd = 2, bty = "n")

# 3. Type 1 error at alpha = 0.05
yr <- range(c(sim_results$t1e_wald_current, sim_results$t1e_wald_betasign,
              sim_results$t1e_stouffer), na.rm = TRUE)
plot(sim_results$rho_true, sim_results$t1e_wald_current, type = "b",
     pch = 19, col = "#D55E00", lwd = 2,
     xlab = expression("True " * rho), ylab = "Type 1 error rate",
     main = expression("Type 1 error (" * alpha * " = 0.05)"),
     xlim = c(0, 0.8), ylim = yr)
lines(sim_results$rho_true, sim_results$t1e_wald_betasign, type = "b",
      pch = 17, col = "#0072B2", lwd = 2)
lines(sim_results$rho_true, sim_results$t1e_stouffer, type = "b",
      pch = 15, col = "grey60", lwd = 2)
abline(h = 0.05, col = "grey40", lty = 2)
legend("topleft", legend = c("Current", "Beta-sign", "Uncorrected"),
       col = c("#D55E00", "#0072B2", "grey60"), pch = c(19, 17, 15), lwd = 2, bty = "n")

dev.off()
cat("Plots saved to null_sim_plots.pdf\n")
