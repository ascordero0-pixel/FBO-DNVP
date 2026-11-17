# =============================================================================
# trackB_calibrate.R  --  Bayesian-style single-parameter calibration (Track B)
# Thesis T3 -- FBO framework
# -----------------------------------------------------------------------------
# IDENTIFIABILITY STRATEGY (documented, per project methodology):
#   The pair {g_person, C_d} is only weakly identifiable: in the ventilation-
#   dominated regime the steady-state balance gives
#       C_ss - C_ext  ~  N * g_person / (C_d * A_eq * sqrt(...))  ~  g_person / C_d
#   i.e. a likelihood ridge along g_person/C_d = const. This is the same class
#   of degeneracy resolved elsewhere in the project (e.g. Q_inf<->C_d in Track A).
#   Resolution: FIX C_d = 0.65 (building aerodynamic property, from T1) and
#   calibrate ONLY g_person (the occupational quantity that Track B exists to
#   estimate). The Tugores et al. (2024) prior (13.3 L/h at age 11) then serves
#   as an INDEPENDENT cross-validation of the calibrated value, not as a constraint.
#
# LOO-CV design (T1 monitoring campaign):
#   TRAIN = Nov 18 + Nov 26  (DNVP)
#   TEST  = Dec 17           (SNVP, scenario S2, WWeqR const 2.16%)
#   Per-day metrics also reported to expose the model's validity domain
#   (Nov 18 = storm day, most demanding case).
# =============================================================================

source("load_co2_campaign.R")
source("trackB_massbalance.R")

# ---- Tugores prior ----------------------------------------------------------
G_PRIOR_LH   <- 13.3                      # L/h per person, 11-yr interpolation
G_PRIOR_M3S  <- G_PRIOR_LH / 1000 / 3600  # m3/s
G_LOWER_M3S  <- 1e-6                       # bounds for optimiser
G_UPPER_M3S  <- 1.5e-5

# ---- objective: SSE over a set of day ids ----------------------------------
sse_over <- function(day_ids, g_person, days_list, Cd = CD_FIX) {
  tot <- 0
  for (dd in day_ids) {
    sim <- simulate_day(days_list[[dd]], g_person, Cd)
    r <- sim$Csim - sim$Cmeas
    tot <- tot + sum(r^2, na.rm = TRUE)
  }
  tot
}

# ---- 1-D optimisation of g_person (C_d fixed) ------------------------------
calibrate_g <- function(train_ids, days_list, Cd = CD_FIX) {
  opt <- optimize(function(g) sse_over(train_ids, g, days_list, Cd),
                  interval = c(G_LOWER_M3S, G_UPPER_M3S), tol = 1e-12)
  opt$minimum
}

# ---- driver ----------------------------------------------------------------
run_calibration <- function() {
  tb <- build_trackB(verbose = TRUE)
  days_list <- tb$days

  TRAIN <- c("18", "26")   # DNVP
  TEST  <- c("17")         # SNVP / S2

  cat("\n", strrep("=", 68), "\n", sep = "")
  cat("TRACK B CALIBRATION  --  C_d fixed at ", CD_FIX,
      " (T1); calibrate g_person\n", sep = "")
  cat(strrep("=", 68), "\n", sep = "")

  g_hat <- calibrate_g(TRAIN, days_list)
  g_hat_lh <- g_hat * 1000 * 3600

  cat(sprintf("\n  g_person (calibrated) = %.4e m3/s = %.2f L/h per person\n",
              g_hat, g_hat_lh))
  cat(sprintf("  Tugores prior         = %.4e m3/s = %.2f L/h (independent)\n",
              G_PRIOR_M3S, G_PRIOR_LH))
  cat(sprintf("  deviation from prior  = %+.1f %%\n",
              (g_hat_lh - G_PRIOR_LH) / G_PRIOR_LH * 100))

  # metrics
  report <- function(ids, label) {
    all_sim <- do.call(rbind, lapply(ids, function(dd)
      simulate_day(days_list[[dd]], g_hat)))
    m <- gof(all_sim)
    cat(sprintf("  %-16s NMBE=%+6.2f%%  CV-RMSE=%6.2f%%  RMSE=%4.0f ppm  n=%d\n",
                label, m$NMBE, m$CVRMSE, m$RMSE, m$n))
    m
  }

  cat("\n  --- Goodness of fit (ASHRAE Guideline 14) ---\n")
  m_tr <- report(TRAIN, "TRAIN Nov18+26:")
  m_te <- report(TEST,  "TEST  Dec17/S2:")
  cat(sprintf("  Generalisation gap (CV-RMSE test-train) = %+.2f pp\n",
              m_te$CVRMSE - m_tr$CVRMSE))

  cat("\n  --- Per-day breakdown (validity domain) ---\n")
  for (dd in c("18","26","17")) {
    md <- gof(simulate_day(days_list[[dd]], g_hat))
    tag <- c("18"="DNVP storm","26"="DNVP moderate","17"="SNVP moderate")[[dd]]
    cat(sprintf("    day %s (%-14s): NMBE=%+6.2f%%  CV-RMSE=%6.2f%%  n=%d\n",
                dd, tag, md$NMBE, md$CVRMSE, md$n))
  }

  invisible(list(g_hat = g_hat, g_hat_lh = g_hat_lh,
                 train = m_tr, test = m_te))
}

if (sys.nframe() == 0) {
  res <- run_calibration()
}
