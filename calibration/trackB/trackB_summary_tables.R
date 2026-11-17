# =============================================================================
# trackB_summary_tables.R  --  Calibration summary tables (Track B)
# Thesis T3 -- FBO framework
# -----------------------------------------------------------------------------
# Produces the most decision-relevant tables to CONFIRM the calibration, each
# exported to CSV (for the manuscript) and printed to console. Tables:
#
#   T1  Calibrated parameter + provenance + 95% CI
#   T2  Goodness-of-fit by fold (ASHRAE Guideline 14: NMBE, CV-RMSE, RMSE)
#   T3  Per-day breakdown (validity domain) with mean wind and regime
#   T4  MSE decomposition (bias / variance / phase) per day  [Theil / Kobayashi]
#   T5  Sensitivity of g_person to methodological choices
#
# The MSE decomposition (T4) is included because aggregate metrics (NMBE,
# CV-RMSE) can mask temporal phase error. For simulated s_i and measured m_i:
#
#   MSE = (1/n) Σ (s_i - m_i)^2
#       = (mean_s - mean_m)^2                      [SB  : squared bias]
#       + (sd_s - sd_m)^2                          [SDSD: variance mismatch]
#       + 2 * sd_s * sd_m * (1 - r)                [LCS : lack of phase corr.]
#
# where sd_ are population standard deviations (÷n) and r is the Pearson
# correlation. SB+SDSD+LCS = MSE exactly. LCS isolates the error attributable
# to imperfect temporal synchronisation (e.g. NDIR sensor lag, mixing latency),
# which is physically distinct from a g_person bias (SB).
# =============================================================================

source("load_co2_campaign.R")
source("trackB_massbalance.R")
source("trackB_calibrate.R")
source("trackB_sensitivity.R")

G_HAT      <- 4.0943e-06      # m3/s  (headline, Option A)
G_HAT_LH   <- G_HAT*1000*3600 # 14.74 L/h
OUTDIR     <- "runs"
if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---- helper: full GoF incl. R2, slope, and MSE decomposition ---------------
gof_full <- function(sim_df) {
  ok <- !is.na(sim_df$Cmeas)
  s <- sim_df$Csim[ok]; m <- sim_df$Cmeas[ok]; n <- length(m)
  mbar <- mean(m); sbar <- mean(s)
  # ASHRAE metrics
  nmbe   <- sum(s - m) / (n * mbar) * 100
  rmse   <- sqrt(mean((s - m)^2))
  cvrmse <- rmse / mbar * 100
  # regression stats
  r   <- suppressWarnings(cor(s, m))
  R2  <- r^2
  slope <- cov(m, s) / var(m)           # simulated regressed on measured
  # population sd (÷n, consistent with MSE decomposition)
  sd_s <- sqrt(mean((s - sbar)^2))
  sd_m <- sqrt(mean((m - mbar)^2))
  MSE  <- mean((s - m)^2)
  SB   <- (sbar - mbar)^2
  SDSD <- (sd_s - sd_m)^2
  LCS  <- 2 * sd_s * sd_m * (1 - r)
  list(n=n, mean_meas=mbar, mean_sim=sbar, NMBE=nmbe, CVRMSE=cvrmse, RMSE=rmse,
       R2=R2, slope=slope, r=r, MSE=MSE, SB=SB, SDSD=SDSD, LCS=LCS,
       SB_pct=SB/MSE*100, SDSD_pct=SDSD/MSE*100, LCS_pct=LCS/MSE*100)
}

build_tables <- function() {
  tb <- build_trackB(verbose = FALSE)
  dl <- tb$days
  TRAIN <- c("18","26"); TEST <- c("17")
  regime <- c("18"="DNVP storm","26"="DNVP moderate","17"="SNVP moderate")

  sim_day <- function(dd) simulate_day(dl[[dd]], G_HAT)
  sim_set <- function(ids) do.call(rbind, lapply(ids, sim_day))

  mean_wind <- function(dd) {
    w <- dl[[dd]]$WindS; round(mean(w[!is.na(w)]), 1)
  }

  # ---- T1: parameter + CI ----
  ci <- profile_ci(TRAIN, dl)
  T1 <- data.frame(
    Parameter = c("g_person (calibrated)","g_person prior (Tugores 2024)",
                  "C_d (fixed)","Q_inf (fixed)","C_ext (fixed)",
                  "Delta_Cp (fixed)","H_stack (fixed)","A_wall (fixed)",
                  "wind factor u_site/u10"),
    Value = c(sprintf("%.2f L/h",G_HAT_LH), "13.30 L/h","0.65","0.64 ren/h",
              "420 ppm","0.75","0.40 m","15.67 m2","0.521"),
    CI95 = c(sprintf("[%.2f, %.2f] L/h",ci$ci_lo_lh,ci$ci_hi_lh),
             "range 9.63-13.29 (age 10)","-","-","-","-","-","-","z=11 m urban"),
    Source = c("Track B calibration","Tugores et al. 2024","T1","Track A",
               "default (not measured)","EN 16798-7 B.3.3.4","window elevations",
               "T1","EN 16798-7 6.4.2.2"),
    stringsAsFactors = FALSE)

  # ---- T2: GoF by fold ----
  g_tr <- gof_full(sim_set(TRAIN)); g_te <- gof_full(sim_set(TEST))
  T2 <- data.frame(
    Fold = c("TRAIN (Nov18+26, DNVP)","TEST (Dec17, SNVP/S2)"),
    n = c(g_tr$n, g_te$n),
    NMBE_pct = round(c(g_tr$NMBE, g_te$NMBE), 2),
    CVRMSE_pct = round(c(g_tr$CVRMSE, g_te$CVRMSE), 2),
    RMSE_ppm = round(c(g_tr$RMSE, g_te$RMSE), 0),
    R2 = round(c(g_tr$R2, g_te$R2), 3),
    slope = round(c(g_tr$slope, g_te$slope), 3),
    stringsAsFactors = FALSE)
  # ASHRAE G14 hourly pass flags
  T2$ASHRAE_NMBE_ok  <- abs(T2$NMBE_pct) <= 10
  T2$ASHRAE_CVRMSE_ok <- T2$CVRMSE_pct <= 30

  # ---- T3: per-day validity domain ----
  T3 <- do.call(rbind, lapply(c("26","17","18"), function(dd) {
    g <- gof_full(sim_day(dd))
    data.frame(Day=dd, Regime=regime[[dd]], mean_u10=mean_wind(dd),
               n=g$n, NMBE_pct=round(g$NMBE,2), CVRMSE_pct=round(g$CVRMSE,2),
               R2=round(g$R2,3), stringsAsFactors=FALSE)
  }))

  # ---- T4: MSE decomposition per day ----
  T4 <- do.call(rbind, lapply(c("26","17","18"), function(dd) {
    g <- gof_full(sim_day(dd))
    data.frame(Day=dd, Regime=regime[[dd]],
               MSE_ppm2=round(g$MSE,0),
               RMSE_ppm=round(sqrt(g$MSE),0),
               bias_SB_pct=round(g$SB_pct,1),
               var_SDSD_pct=round(g$SDSD_pct,1),
               phase_LCS_pct=round(g$LCS_pct,1),
               r=round(g$r,3), stringsAsFactors=FALSE)
  }))

  # ---- T5: sensitivity ----
  T5 <- sensitivity_table(dl)
  T5$g_lh <- round(T5$g_lh, 2); T5$dev_prior <- round(T5$dev_prior, 1)
  names(T5) <- c("Scenario","g_person_Lh","dev_from_prior_pct")

  # ---- write CSVs ----
  write.csv(T1, file.path(OUTDIR,"TableB1_parameters.csv"), row.names=FALSE)
  write.csv(T2, file.path(OUTDIR,"TableB2_gof_by_fold.csv"), row.names=FALSE)
  write.csv(T3, file.path(OUTDIR,"TableB3_validity_domain.csv"), row.names=FALSE)
  write.csv(T4, file.path(OUTDIR,"TableB4_mse_decomposition.csv"), row.names=FALSE)
  write.csv(T5, file.path(OUTDIR,"TableB5_sensitivity.csv"), row.names=FALSE)

  # ---- print ----
  line <- function() cat(strrep("-",72),"\n")
  cat("\n"); line(); cat("TABLE B1 - Calibrated parameters & provenance\n"); line()
  print(T1, row.names=FALSE)
  cat("\n"); line(); cat("TABLE B2 - Goodness of fit by fold (ASHRAE G14)\n"); line()
  print(T2, row.names=FALSE)
  cat("\n"); line(); cat("TABLE B3 - Per-day validity domain\n"); line()
  print(T3, row.names=FALSE)
  cat("\n"); line(); cat("TABLE B4 - MSE decomposition (bias/variance/phase)\n"); line()
  print(T4, row.names=FALSE)
  cat("  [SB=squared bias, SDSD=variance mismatch, LCS=lack of phase correlation;\n")
  cat("   SB+SDSD+LCS=100%. High LCS => temporal phase error dominates.]\n")
  cat("\n"); line(); cat("TABLE B5 - Sensitivity of g_person\n"); line()
  print(T5, row.names=FALSE)

  cat("\nCSV files written: TableB1..TableB5_*.csv\n")
  invisible(list(T1=T1,T2=T2,T3=T3,T4=T4,T5=T5))
}

if (sys.nframe() == 0) tabs <- build_tables()
