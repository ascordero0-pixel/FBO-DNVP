# ============================================================================
# MASTER_FIS_F1.R  --  Orchestrator for T3 FBO PHASE 1 (Bayesian Optimization)
# ----------------------------------------------------------------------------
# Runs, in dependency order, the EnergyPlus-coupled Phase-1 pipeline that turns
# the validated Phase-0 engine (see F0/MASTER_FIS.R) into the optimized DNVP
# set-points for the 9 scenarios. FAIL-FAST: aborts on the first error rather
# than propagating wrong values downstream.
#
# Pipeline (this master, EnergyPlus REQUIRED for the heavy stages):
#   VALIDATE    re-run the Phase-0 gates (FIS + g1..g5 + control law)   [fast]
#   AQI GATE    make_aqi_ext.R   outdoor worst-ICA series -> Schedule:File (C:/EP_work)
#   SMOKE       one evaluate_U() call -> validate the IDF->E+->U coupling
#   SEED        LHS 50x9 -> F1/seed_9scen_results.csv        [EXPENSIVE ~450 E+ runs]
#   REPRICE     reprice.R  re-score U with Cost-active terciles (in place, backups)
#   BO          bo_driver.R run_bo()  GP(Matern5/2+ARD)+EI loop  [EXPENSIVE ~150 E+ runs]
#   PRODUCTS    build_best_x()  argmax-U per scenario -> F1/best_x_9scen_v2.csv
#   MONITOR     bo_monitor.R  monitor_bo()   (read-only, run in a SECOND session)
#
# Project layout: ONE FIS.Rproj at the root. Reusable engine in core/ ; Phase-1
# drivers + IDF in F1/ ; input CSVs in data/ ; EnergyPlus I/O redirected to the
# ASCII path C:/EP_work (non-ASCII OneDrive paths break E+). Paths via `here`.
#
# Canonical constants (locked; do NOT edit casually):
#   V = 126.7031 m3 | station 41091009 (Torneo) | outdoor year 2024 |
#   x = (tau_CO2[800,1200], delta_hys[50,200], d_slot[5,30], f_open[0.10,f_max(c)],
#        tau_AQ[1,6]) ; North Axis {NE=0, N=329, E=59} ; f_max {w154 .16, w216 .22, w410 .42}.
#   FIS anchoring (Table B7): g1[550/800] g2 ramp 5-10% g3 terciles 5.40/12.02
#   g4 terciles 6.44/8.44 g5[2.5/3.5] Cost[1/3,2/3] U{0.20/0.50/0.85}.
#   FIS correction (2026-08): Cost sub-FIS rule SC7 (g3=High, g4=Low) M->H -- energy is no longer
#   masked by low switching (surfaced by the Phase-2 NSGA-II benchmark). ZERO impact on best_x /
#   seed / BO optima (the masked corner was never sampled) -> no re-run needed. Validation 105/105.
#
# STATE: Phase 1 is ALREADY RUN (seed + BO closed, best_x exists). This master is
# the REPRODUCIBILITY artifact: the heavy toggles default to FALSE so a plain
# source() only re-validates the engine and (optionally) rebuilds best_x. Flip a
# toggle to re-run that stage from scratch.
#
# Usage:  open FIS.Rproj, then  source(here::here("F1", "MASTER_FIS_F1.R"))
#         (pause OneDrive before any stage that launches EnergyPlus)
# ============================================================================

library(here)

# ---- 0. Configuration + toggles --------------------------------------------
CORE_DIR <- here::here("core")
F1_DIR   <- here::here("F1")
DATA_DIR <- here::here("data")
V_ZONE   <- 126.7031
SEED_CSV <- file.path(F1_DIR, "seed_9scen_results.csv")
BO_DIR   <- file.path(F1_DIR, "bo")
BEST_X   <- file.path(BO_DIR, "best_x_9scen_v2.csv")   # canonical location [repo patch P4]

# Feature toggles. Heavy (EnergyPlus) stages are OFF by default: a plain source()
# only validates the engine and rebuilds best_x from the existing seed+BO CSVs.
RUN_VALIDATE  <- TRUE     # re-run Phase-0 unit gates (fast, no E+)
RUN_AQI       <- FALSE    # rebuild AQI_ext series -> C:/EP_work (needs data/, no E+)
RUN_SMOKE     <- FALSE    # one evaluate_U() coupling test (1 annual E+ run)
RUN_SEED      <- FALSE    # LHS 50x9  -- EXPENSIVE (~450 E+ runs)
RUN_SEED_SIGNCHECK <- TRUE # (only within SEED) 1 E+ run to validate the North Axis sign first
RUN_REPRICE   <- FALSE    # re-score U with Cost-active terciles (modifies CSVs in place)
RUN_BO        <- FALSE    # Bayesian optimization loop -- EXPENSIVE (~150 E+ runs)
RUN_BEST_X    <- TRUE     # assemble best_x_9scen_v2.csv from seed + BO (fast, no E+)

# scenario space (mirrors bo_driver.R)
ORIENT <- c(NE = 0, N = 329, E = 59); FMAX <- c(w154 = 0.16, w216 = 0.22, w410 = 0.42)
ALL9   <- as.vector(t(outer(names(ORIENT), names(FMAX), paste, sep = "_")))
VARS   <- c("tau_CO2", "delta_hys", "d_slot", "f_open", "tau_AQ")

# ---- FAIL-FAST guard #1: required scripts + inputs present ------------------
.need_core <- file.path(CORE_DIR, c("fis_scoring.R", "indicators.R", "g5.R",
                                    "control_law.R", "load_outdoor.R"))
.need_f1   <- file.path(F1_DIR, c("fbo_eval.R", "make_aqi_ext.R", "bo_driver.R", "reprice.R"))
.miss <- c(.need_core, .need_f1)[!file.exists(c(.need_core, .need_f1))]
if (length(.miss)) stop("Missing script(s): ", paste(basename(.miss), collapse = ", "),
                        "\n  (open the project via FIS.Rproj so `here` resolves the root)")

cat("== MASTER_FIS_F1: starting ==\n\n")

# ---- 1. Source Phase-1 modules (define only; heavy loops are gated below) ---
# fbo_eval.R sources the Phase-0 core (indicators/g5/fis_scoring/load_outdoor)
# and defines fbo_setup() + evaluate_U(). make_aqi_ext.R and bo_driver.R only
# DEFINE functions when sourced. reprice.R is NOT sourced here (it re-scores on
# source) -- it is handled under RUN_REPRICE below.
cat("[1] sourcing Phase-1 modules ...\n")
source(file.path(CORE_DIR, "control_law.R"))  # run_tests_control() for the gate [repo patch P1]
source(file.path(F1_DIR, "fbo_eval.R"))        # fbo_setup(), evaluate_U()
source(file.path(F1_DIR, "make_aqi_ext.R"))    # make_aqi_ext()
source(file.path(F1_DIR, "bo_driver.R"))       # run_bo(), SCN, LO, ...
source(file.path(F1_DIR, "bo_monitor.R"))      # monitor_bo()  (read-only helper)

# ---- 2. Validation gates (FAIL-FAST, Phase-0 engine) -----------------------
if (RUN_VALIDATE) {
  cat("\n[2] validation gates (abort on first failure)\n")
  cat("--- scoring FIS (expect 105) ---\n");      run_tests()
  cat("--- indicators g1..g4 (expect 16) ---\n"); run_tests_indicators()
  cat("--- g5 (expect 15) ---\n");                run_tests_g5()
  cat("--- control law (expect 12) ---\n");       run_tests_control()
} else cat("\n[2] validation skipped (RUN_VALIDATE = FALSE)\n")

# ---- 3. AQI_ext control-gate series ----------------------------------------
if (RUN_AQI) {
  cat("\n[3] building AQI_ext series -> C:/EP_work ...\n")
  make_aqi_ext()                                 # writes AQI_ext_2024_5min.csv (8760 h)
} else cat("\n[3] AQI_ext skipped (RUN_AQI = FALSE; already built)\n")

# ---- 4. Smoke test: IDF -> EnergyPlus -> U coupling ------------------------
if (RUN_SMOKE) {
  cat("\n[4] smoke test: one evaluate_U() call (1 annual E+ run) ...\n")
  P_IDF <- c(file.path(F1_DIR, "v4_FBO_annual_prod.idf"), "C:/EP_work/v4_FBO_annual_prod.idf")
  P_IDF <- P_IDF[file.exists(P_IDF)][1]
  P_EPW <- list.files(DATA_DIR, pattern = "\\.epw$", full.names = TRUE)[1]
  fbo_setup(base_idf = P_IDF, epw = P_EPW, work_dir = "C:/EP_work/FBO", V = V_ZONE)
  r <- evaluate_U(list(tau_CO2 = 900, delta_hys = 100, d_slot = 10, f_open = 0.28, tau_AQ = 5),
                  north = ORIENT[["N"]], tag = "smoke_master")
  cat(sprintf("    U=%.4f | g1=%.0f | p_III=%.2f | g3=%.2f g4=%.2f g5=%.2f\n",
              r$U, r$g1, r$g2[2], r$g3, r$g4, r$g5))
} else cat("\n[4] smoke test skipped (RUN_SMOKE = FALSE)\n")

# ---- 5. LHS seed (50 x 9) --  EXPENSIVE (~450 annual E+ runs) ---------------
# seed_driver.R defines check_north_sign() and run_seed(). run_seed() is RESUMABLE
# (skips done tags, incremental save) and auto-builds AQI_ext if missing. It writes
# F1/seed_9scen_results.csv. The seed is ALREADY BUILT, so this stays OFF by default.
if (RUN_SEED) {
  .sd <- file.path(F1_DIR, "seed_driver.R")
  if (!file.exists(.sd)) stop("[5] seed_driver.R not found in F1/")
  cat("\n[5] LHS seed 50x9 -> F1/seed_9scen_results.csv (EXPENSIVE) ...\n")
  source(.sd)                                    # defines check_north_sign(), run_seed()
  if (RUN_SEED_SIGNCHECK) { cat("    North Axis sign check (1 E+ run) ...\n"); check_north_sign() }
  run_seed()                                     # full sweep; resumable (re-source master to continue)
} else cat("\n[5] LHS seed skipped (RUN_SEED = FALSE; seed already built)\n")

# ---- 6. Reprice: re-score U with Cost-active terciles ----------------------
# reprice.R re-scores the U column of the seed and every F1/bo/bo_<scn>.csv using
# the corrected g3/g4 terciles (Cost activated), keeping the physical g1..g5.
# It CREATES backups *_placeholderU.csv and REWRITES the CSVs IN PLACE.
if (RUN_REPRICE) {
  cat("\n[6] reprice: re-scoring U (Cost-active terciles; backups + in-place) ...\n")
  source(file.path(F1_DIR, "reprice.R"))
} else cat("\n[6] reprice skipped (RUN_REPRICE = FALSE; U already Cost-active)\n")

# ---- 7. Bayesian optimization loop --  EXPENSIVE ---------------------------
if (RUN_BO) {
  cat("\n[7] Bayesian optimization (GP Matern5/2+ARD + EI) -- EXPENSIVE ...\n")
  stopifnot("falta seed CSV" = file.exists(SEED_CSV))
  run_bo()                                        # los 9 escenarios; ver bo_driver.R para n_max/STALL_MAX
} else cat("\n[7] BO skipped (RUN_BO = FALSE; BO already run -> F1/bo/)\n")

# ---- 8. Product: assemble best_x (argmax U per scenario) -------------------
# best_x = the max-U evaluated point per scenario over (seed + BO). Rebuilds
# F1/best_x_9scen_v2.csv deterministically from the current CSVs. Uses the FIS
# to RE-SCORE U (so the file is always Cost-active/consistent with fis_scoring).
build_best_x <- function(out = BEST_X) {
  sf <- build_scoring_fis(default_variables())
  keep <- c("scn", VARS, "g1", "g2_pII", "g2_pIII", "g3", "g4", "g5")
  s <- read.csv(SEED_CSV, stringsAsFactors = FALSE)
  rows <- list(s[, keep])
  for (f in list.files(BO_DIR, pattern = "^bo_.*\\.csv$", full.names = TRUE)) {
    b <- read.csv(f, stringsAsFactors = FALSE); if (all(keep %in% names(b))) rows[[length(rows)+1]] <- b[, keep]
  }
  d <- do.call(rbind, rows); d <- d[stats::complete.cases(d[, c(VARS, "g1")]), ]
  d$U <- vapply(seq_len(nrow(d)), function(i) scoring_utility(
    sf, g1 = d$g1[i], g2 = c(d$g2_pII[i], d$g2_pIII[i]),
    g3 = d$g3[i], g4 = d$g4[i], g5 = d$g5[i]), numeric(1))
  bx <- do.call(rbind, lapply(ALL9, function(scn) {
    ds <- d[d$scn == scn, ]; if (!nrow(ds)) return(NULL)
    r <- ds[which.max(ds$U), ]
    data.frame(scn = scn, r[VARS], U = round(r$U, 4), g1 = round(r$g1, 1),
               p_II = round(r$g2_pII, 2), p_III = round(r$g2_pIII, 2),
               g3 = round(r$g3, 1), g4 = round(r$g4, 1), stringsAsFactors = FALSE)
  }))
  write.csv(bx, out, row.names = FALSE)
  cat(sprintf("    best_x -> %s  (%d escenarios; mejor U = %.4f en %s)\n",
              out, nrow(bx), max(bx$U), bx$scn[which.max(bx$U)]))
  invisible(bx)
}
if (RUN_BEST_X) {
  cat("\n[8] assembling best_x (argmax U per scenario) ...\n")
  BEST <- build_best_x()
} else cat("\n[8] best_x skipped (RUN_BEST_X = FALSE)\n")

# ---- 9. Done ---------------------------------------------------------------
cat("\n== MASTER_FIS_F1: pipeline complete ==\n")
cat("   MONITOR (2nd R session, read-only): source(here::here('F1','bo_monitor.R')); monitor_bo()\n")
cat("   next: Phase 2 benchmark -> source(here::here('F2','MASTER_FIS_F2.R'))\n")
