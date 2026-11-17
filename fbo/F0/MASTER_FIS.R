# ============================================================================
# MASTER_FIS.R  --  Orchestrator for the T3 FBO "FIS" project (Phase-0, pure R)
# ----------------------------------------------------------------------------
# Runs, in dependency order, every script built for the scoring / indicator
# layer of the Fuzzy-Bayesian-Optimization framework, gating each stage behind
# its own validation suite (FAIL-FAST: the pipeline aborts on the first error
# rather than propagating wrong values downstream).
#
# Pipeline (this master, EnergyPlus NOT required):
#   ENGINE      fis_scoring.R   hierarchical Mamdani scoring FIS   g1..g5 -> U
#               indicators.R    direct indicators                 series -> g1..g4
#               g5.R            outdoor-import indicator           series + Cout -> g5
#               control_law.R   crisp DNVP opening law (Eq. 9)     state -> s(t)
#               load_outdoor.R  Torneo-2024 pollutant-series loader
#   VALIDATE    105 + 13 + 15 + 12 unit checks + loader self-checks
#   DATA        build the continuous 5-min outdoor series (Cout_PM25/NO2/O3)
#   PRODUCTS    Section 2.7 membership-function figure ; g5 sensitivity analysis
#
# Project layout (Option 1): ONE FIS.Rproj at the root. Reusable modules in
# core/ ; input CSVs in data/ ; this driver in F0/ ; outputs in runs/. Paths are
# resolved with `here` from the .Rproj root, so this runs from any subfolder.
#
# Canonical constants locked earlier this project (do not edit casually):
#   zone volume V = 126.7031 m3 (IDF v4) ; station 41091009 (Torneo, Seville) ;
#   outdoor year = 2024 (last validated) ; ICA sub-index = Orden TEC/351/2019 ;
#   theta_rm clamped to the EN 16798-1 [10,30] applicability range.
#
# 2026-07 update (evaluator closure, see claude/T3_Fase1_*):
#   * fis_scoring.R  : g2 fail-ramp RE-ANCHORED [5%,10%] (5% = EN Cat III limit).
#   * indicators.R   : Palanca 2 -- air-speed comfort extension (+dtheta*s on g2 upper).
#   These change the g2 membership figure (Section 2.7) and the g2 row of Table 6.
#
# 2026-08 update (benchmark closure, see claude/T3_Fase1_evaluador_U_correccion_SC7.md):
#   * fis_scoring.R  : Cost sub-FIS rule SC7 (g3=High, g4=Low)  M -> H. The Phase-2 NSGA-II
#     benchmark surfaced that LOW switching (g4) masked HIGH overcooling energy (g3); this
#     restores energy in the Cost. ZERO blast radius: best_x, seed and BO optima unchanged
#     (the masked corner was never sampled). Validation still 105/105 (T1 (H,L)->H, (M,L)->M).
#
# Outside this master (EnergyPlus-coupled, run in the local environment):
#   EMS-embedded control law (Eq. 12, native) + monitored-day smoke test ;
#   LHS (50x9) closing the g3/g4 terciles ; the FBO loop (Phase 1).
#
# Usage:  open FIS.Rproj, then  source(here::here("F0", "MASTER_FIS.R"))
# ============================================================================

library(here)   # resolves all paths from the FIS.Rproj root


# ---- 0. Configuration ------------------------------------------------------
CORE_DIR <- here::here("core")  # folder holding the reusable engine modules
DATA_DIR <- here::here("data")  # folder holding the three Torneo 2024 CSVs
STATION  <- "41091009"          # Torneo (Seville); the loader filters to it
YEAR     <- 2024                # last consolidated/validated year of the network
V_ZONE   <- 126.7031            # canonical zone volume [m3]

# Feature toggles: set to FALSE to skip the slower "product" stages.
RUN_FIGURES     <- TRUE         # membership-function plots (fast)
RUN_SENSITIVITY <- TRUE         # g5 sensitivity sweep (~15 s, reads the CSVs)

# Outdoor CSV files. If your PM2.5 is the cleaned copy, change its name here;
# either way load_outdoor() keeps only station 41091009 internally.
CSV <- list(
  PM25 = file.path(DATA_DIR, "PM25_HH_2024_SEV_Torneo.csv"),
  NO2  = file.path(DATA_DIR, "NO2_HH_2024_SEV_Torneo.csv"),
  O3   = file.path(DATA_DIR, "O3_HH_2024_SEV_Torneo.csv"))

OUT_DIR <- here::here("runs"); if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# FAIL-FAST guard #1: all engine scripts and input files must be present.
.engine  <- c("fis_scoring.R", "indicators.R", "g5.R", "control_law.R", "load_outdoor.R")
.eng_paths <- file.path(CORE_DIR, .engine)
.missing_scripts <- .engine[!file.exists(.eng_paths)]
if (length(.missing_scripts))
  stop("Missing engine script(s) in core/: ", paste(.missing_scripts, collapse = ", "),
       "\n  (open the project via FIS.Rproj so `here` resolves the root)")
.missing_csv <- unlist(CSV)[!file.exists(unlist(CSV))]
if (length(.missing_csv))
  stop("Missing input CSV(s): ", paste(.missing_csv, collapse = ", "),
       "\n  (place the Torneo CSVs in ", normalizePath(DATA_DIR, mustWork = FALSE), " or edit DATA_DIR)")

cat("== MASTER_FIS: starting  (station", STATION, "| year", YEAR, ") ==\n\n")


# ---- 1. Engine modules -----------------------------------------------------
# Sourcing only DEFINES functions; the sys.nframe() guards inside each file keep
# their test harnesses from auto-running here, so we control validation below.
cat("[1] sourcing engine modules ...\n")
source(file.path(CORE_DIR, "fis_scoring.R"))   # scoring FIS: MF primitives, rule bases, ScoringFIS
source(file.path(CORE_DIR, "indicators.R"))    # g1 (CO2 P95), g2 (comfort ladder), g3 (EDM), g4 (effort)
source(file.path(CORE_DIR, "g5.R"))            # per-species dynamic balance + worst ICA sub-index
source(file.path(CORE_DIR, "control_law.R"))   # crisp opening law (Eq. 9): state -> window state s(t)
source(file.path(CORE_DIR, "load_outdoor.R"))  # wide-CSV parsing, gap-filling, 5-min expansion


# ---- 2. Validation gates (FAIL-FAST) --------------------------------------
# Each suite calls stop() on any failure, which aborts the whole master. A clean
# run prints "RESULT: N passed, 0 failed" for each layer.
cat("\n[2] validation gates (abort on first failure)\n")
cat("--- scoring FIS (expect 105 passed) ---\n");        run_tests()
cat("--- indicators g1..g4 (expect 16 passed) ---\n");   run_tests_indicators()  # +3 Palanca-2 checks
cat("--- g5 (expect 15 passed) ---\n");                  run_tests_g5()
cat("--- control law (expect 12 passed) ---\n");         run_tests_control()


# ---- 3. Outdoor data pipeline ---------------------------------------------
# Build the continuous full-year 5-min outdoor series that g5.R consumes.
cat("\n[3] building outdoor series (Torneo", YEAR, ") ...\n")
outdoor <- load_outdoor(CSV, station = STATION, year = YEAR,
                        five_min = TRUE, verbose = TRUE)

# FAIL-FAST guard #2: the series must be complete (leap year, no residual NA).
stopifnot(nrow(outdoor$hourly)  == 8784,
          nrow(outdoor$fivemin) == 8784 * 12,
          !anyNA(outdoor$fivemin[c("Cout_PM25", "Cout_NO2", "Cout_O3")]))
cat("    outdoor 5-min series OK:", nrow(outdoor$fivemin), "rows, no NA\n")
# `outdoor$fivemin` (timestamp + Cout_*) is the object to cbind with the
# EnergyPlus Q(t) and the occupancy mask when the coupling is wired up.


# ---- 4. Product: Section 2.7 membership-function figure -------------------
if (RUN_FIGURES) {
  cat("\n[4] rendering Section 2.7 membership-function figure -> runs/ ...\n")
  .mf <- file.path(CORE_DIR, "plot_mf_panel.R")
  if (file.exists(.mf)) source(.mf)  # T3_Fig_MF_primaries_3panel (png+pdf), plotmath labels
  else cat("    plot_mf_panel.R not found in core/; figure skipped\n")
} else cat("\n[4] figure skipped (RUN_FIGURES = FALSE)\n")


# ---- 5. Product: g5 sensitivity analysis ----------------------------------
if (RUN_SENSITIVITY) {
  cat("\n[5] g5 sensitivity to assumed loss rates k_s -> runs/ ...\n")
  Sys.setenv(TORNEO_DIR = DATA_DIR)          # sensitivity_g5.R reads the CSVs from here
  source(file.path(CORE_DIR, "sensitivity_g5.R"))  # T3_g5_sensitivity (png + pdf) + console report
} else cat("\n[5] sensitivity skipped (RUN_SENSITIVITY = FALSE)\n")


# ---- 6. Done ---------------------------------------------------------------
cat("\n== MASTER_FIS: pipeline complete ==\n")
cat("   validated engine (FIS + g1..g5) + outdoor data + products in ",
    normalizePath(OUT_DIR), "\n", sep = "")
cat("   next (EnergyPlus, local env): EMS-embedded control law (Eq. 12)",
    "+ smoke test -> LHS (g3/g4 terciles) -> FBO loop (Phase 1)\n")
