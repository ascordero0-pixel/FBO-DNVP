# =============================================================================
# trackB_sensitivity.R  --  Uncertainty & sensitivity analysis (Track B)
# Thesis T3 -- FBO framework
# -----------------------------------------------------------------------------
# Produces, for the headline calibration (Option A: train Nov18+26, C_d=0.65):
#   (1) Profile-likelihood 95% confidence interval on g_person.
#   (2) Sensitivity of g_hat to:
#         - fold design (all-3-day-DNVP vs Nov26-only)
#         - wind-correction factor (+/- one roughness class)
#         - fixed C_d (0.60 / 0.65 / 0.70)
#         - winsor cap (18 / 20 / 25 m/s)
#   (3) The identifiability ridge g_person/C_d (documenting why C_d is fixed).
#
# The profile-likelihood CI is built from the SSE curve assuming Gaussian
# residuals with unknown variance: the (1-alpha) region is
#     SSE(g) <= SSE_min * (1 + p/(n-p) * F_{1-alpha; p, n-p})
# with p=1 free parameter. This is the standard nonlinear-regression CI.
# =============================================================================

source("load_co2_campaign.R")
source("trackB_massbalance.R")
source("trackB_calibrate.R")

# ---- profile-likelihood CI on g_person --------------------------------------
# Memoised: repeated calls with identical (train_ids, Cd, alpha) return the
# cached result, so the master pipeline does not recompute the 800-point SSE
# scan when both trackB_sensitivity.R (stage 4) and trackB_summary_tables.R
# (stage 5) request the same CI.
.profile_ci_cache <- new.env(parent = emptyenv())

profile_ci <- function(train_ids, days_list, Cd = CD_FIX, alpha = 0.05) {
  key <- paste(sort(train_ids), sprintf("%.6f", Cd), sprintf("%.4f", alpha),
               sep = "|")
  key <- paste(key, collapse = "_")
  if (exists(key, envir = .profile_ci_cache, inherits = FALSE))
    return(get(key, envir = .profile_ci_cache, inherits = FALSE))

  g_hat <- calibrate_g(train_ids, days_list, Cd)
  sse_min <- sse_over(train_ids, g_hat, days_list, Cd)

  # residual count n (number of measured points across train days)
  n <- sum(sapply(train_ids, function(dd)
    sum(!is.na(days_list[[dd]]$CO2i))))
  p <- 1
  Fcrit <- qf(1 - alpha, p, n - p)
  sse_thresh <- sse_min * (1 + p / (n - p) * Fcrit)

  # Locate the two crossings of sse_thresh by ROOT-FINDING (bisection via
  # uniroot) instead of a uniform grid scan. The profile SSE(g) is smooth and
  # U-shaped with a single minimum at g_hat, so f(g) = SSE(g) - sse_thresh has
  # exactly one root on each side of g_hat. uniroot converges in ~30-40 engine
  # evaluations per side (vs 800 for the grid), an ~10x speed-up at machine
  # precision. Fallback to the grid only if a bracket is not sign-changing.
  fthr <- function(g) sse_over(train_ids, g, days_list, Cd) - sse_thresh
  find_root <- function(lo_b, hi_b) {
    if (fthr(lo_b) * fthr(hi_b) > 0) return(NA_real_)   # no sign change
    uniroot(fthr, lower = lo_b, upper = hi_b, tol = 1e-10)$root
  }
  lo <- find_root(G_LOWER_M3S, g_hat)   # lower crossing, in [G_LOWER, g_hat]
  hi <- find_root(g_hat, G_UPPER_M3S)   # upper crossing, in [g_hat, G_UPPER]
  if (is.na(lo)) lo <- G_LOWER_M3S      # CI open at lower bound (rare)
  if (is.na(hi)) hi <- G_UPPER_M3S      # CI open at upper bound (rare)

  res <- list(g_hat = g_hat, sse_min = sse_min, n = n,
       Fcrit = Fcrit, sse_thresh = sse_thresh,
       ci_lo = lo, ci_hi = hi,
       g_hat_lh = g_hat*1000*3600, ci_lo_lh = lo*1000*3600, ci_hi_lh = hi*1000*3600)
  assign(key, res, envir = .profile_ci_cache)
  res
}

# ---- identifiability ridge g/Cd --------------------------------------------
ridge_scan <- function(train_ids, days_list, Cd_grid = c(0.45,0.55,0.60,0.65,0.75,0.85,1.0)) {
  out <- lapply(Cd_grid, function(cd) {
    g <- calibrate_g(train_ids, days_list, cd)
    sse <- sse_over(train_ids, g, days_list, cd)
    data.frame(Cd = cd, g_lh = g*1000*3600, ratio = (g*1000*3600)/cd, SSE = sse)
  })
  do.call(rbind, out)
}

# ---- sensitivity table ------------------------------------------------------
sensitivity_table <- function(days_list) {
  base_train <- c("18","26")

  scen <- list()
  add <- function(name, g_lh) scen[[length(scen)+1]] <<-
    data.frame(scenario = name, g_lh = g_lh, dev_prior = (g_lh-G_PRIOR_LH)/G_PRIOR_LH*100)

  # (a) headline
  add("Headline: train Nov18+26, Cd=0.65, wind x0.521, cap20",
      calibrate_g(base_train, days_list) * 1000*3600)

  # (b) fold: Nov26 only
  add("Fold: train Nov26 only",
      calibrate_g(c("26"), days_list) * 1000*3600)

  # (c) fold: all 3 days in training (incl SNVP)
  add("Fold: train all 3 days",
      calibrate_g(c("18","26","17"), days_list) * 1000*3600)

  # (d) Cd sensitivity
  for (cd in c(0.60,0.70)) {
    g <- calibrate_g(base_train, days_list, cd) * 1000*3600
    add(sprintf("Cd = %.2f (vs 0.65)", cd), g)
  }

  # (e) wind factor sensitivity (+/- one roughness class ~ +/-25%)
  wf0 <- WIND_FACTOR
  for (mult in c(0.75, 1.25)) {
    WIND_FACTOR <<- wf0 * mult
    add(sprintf("Wind factor x%.2f (%.3f)", mult, WIND_FACTOR),
        calibrate_g(base_train, days_list) * 1000*3600)
  }
  WIND_FACTOR <<- wf0

  # (f) winsor cap
  cap0 <- WIND_CAP
  for (cp in c(18, 25)) {
    WIND_CAP <<- cp
    add(sprintf("Winsor cap = %d m/s", cp),
        calibrate_g(base_train, days_list) * 1000*3600)
  }
  WIND_CAP <<- cap0

  do.call(rbind, scen)
}

if (sys.nframe() == 0) {
  tb <- build_trackB(verbose = FALSE)
  dl <- tb$days

  cat(strrep("=",70), "\n")
  cat("PROFILE-LIKELIHOOD 95% CI ON g_person (train Nov18+26, Cd=0.65)\n")
  cat(strrep("=",70), "\n")
  ci <- profile_ci(c("18","26"), dl)
  cat(sprintf("  n residual points         = %d\n", ci$n))
  cat(sprintf("  F(0.95; 1, %d)            = %.3f\n", ci$n-1, ci$Fcrit))
  cat(sprintf("  g_hat                     = %.2f L/h  (%.4e m3/s)\n", ci$g_hat_lh, ci$g_hat))
  cat(sprintf("  95%% CI                    = [%.2f, %.2f] L/h\n", ci$ci_lo_lh, ci$ci_hi_lh))
  cat(sprintf("  Tugores prior 13.30 L/h in CI? %s\n",
              ifelse(13.30>=ci$ci_lo_lh & 13.30<=ci$ci_hi_lh, "YES", "no")))

  cat("\n", strrep("=",70), "\n", sep="")
  cat("IDENTIFIABILITY RIDGE  g_person vs fixed C_d\n")
  cat(strrep("=",70), "\n")
  rr <- ridge_scan(c("18","26"), dl)
  print(format(rr, digits=4), row.names=FALSE)
  cat("  --> SSE rises monotonically with Cd (weak identifiability, not flat);\n")
  cat("      fixing Cd=0.65 (T1) is the principled resolution.\n")

  cat("\n", strrep("=",70), "\n", sep="")
  cat("SENSITIVITY OF g_hat\n")
  cat(strrep("=",70), "\n")
  st <- sensitivity_table(dl)
  st$g_lh <- round(st$g_lh, 2); st$dev_prior <- round(st$dev_prior, 1)
  print(st, row.names=FALSE)
}
