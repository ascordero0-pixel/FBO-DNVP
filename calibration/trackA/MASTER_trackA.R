# =============================================================================
# MASTER_trackA.R  --  Master runner for the Track A (passive thermal) closure
# =============================================================================
# Reproduces, in a single top-to-bottom pass, every calibrated metric and figure
# reported in manuscript Section 3.1 and Appendix A.1 for the case-study
# classroom. Fail-fast: geometry is audited before any metric is trusted.
#
# Locked environment : DesignBuilder 7.3.1.003 / EnergyPlus 9.4.0.002
#                      Python 3.13 + eppy 0.5.69, RStudio + reticulate
# Canonical geometry : floor area = 45.2511 m2, volume = 126.7031 m3
# Run period         : 1 Nov - 20 Dec 2021, timestep 12 (5 min), 50 free-running nights
#
# ---- CANONICAL (ANCHORED) TRACK A VALUES ------------------------------------
#   v3 (single-zone, RETAINED model):
#       NMBE_T / CV-RMSE_T = -5.41 % / 8.04 %  (5-min,  n = 5399)
#                          = -5.42 % / 8.04 %  (hourly, n = 451)   -> ASHRAE G14 PASS
#       full-campaign NMBE = +0.12 %
#   v6 (two-zone, REJECTED - annex negative result, RUN_V6 optional):
#       NMBE_T / CV-RMSE_T = +6.05 % / 9.21 %  (hourly)            -> G14 PASS
#       full-campaign NMBE = +12.43 % ; rigid offset v6-v3 = +2.1 C (11.5 pp)
#   tau : per-night exponential fit, median = 38.1 h (R2 median 0.974, 43 nights)
#         (pooled lumped fit = 57.4 h, R2 = 0.04 -> outlier estimator, annex contrast)
#   UA  = 99.3 W/K  (direct-exterior: 3 facades incl. corrected north Wall_3
#                    + 9 windows; ceiling & floor are interior CTF boundaries.
#                    T1 survey value 200.7 W/K counts the roof as direct loss
#                    and is a NON-comparable upper bound.)
#   Infiltration = 0.64 ren/h global / 0.41 ren/h nocturnal
#   Calibration point : crack_scale = 1, cap_mult = 1, q_heat = 20 W/m2 (STEP)
#   Occupancy : 0.44 p/m2 (~20 occ), metabolic 105 W/person, radiant frac 0.30,
#               weekday schedule 09:10-14:00
#
# ---- PUBLISHABLE IDF --------------------------------------------------------
#   T1_aula_base_94_v4.idf  is the self-contained IDF: the inherited T1 export
#   tagged the north facade (Wall_3) with a spurious ExteriorNaturalVentedCavity;
#   v4 persists the Wall_3 -> Outdoors correction in the file, so a reviewer
#   reproduces -5.41 % by opening the IDF directly (no runtime normalization
#   needed). Proven equivalent to v3 + runtime normalize (both give -5.41/8.04).
#
# ---- REPRODUCIBILITY NOTES (read before running) ----------------------------
#   * Handles 'ce'/'cfg' live in .STATE and DIE on RStudio restart (<pointer:0x0>).
#     Always start from a CLEAN session and run this file top to bottom.
#   * The loader flag '.CALIB_SERIES_LOADED' can survive a restart while 'D'
#     degrades to an atomic vector -> "$ operator is invalid for atomic vectors".
#     This master clears it defensively before the figure phase.
#   * Pause OneDrive sync before executing (file locks corrupt EnergyPlus runs).
#   * cfg$work_dir may read empty via reticulate on a stale handle; PATHS below
#     pins runs_dir = "runs" explicitly.
# =============================================================================

# ---- run options ------------------------------------------------------------
RUN_V6 <- TRUE   # regenerate the two-zone model (needed for Fig A.2). Set FALSE
                 # to reproduce only the retained single-zone (v3) numbers.

# =============================================================================
# PHASE 0 - environment setup
# =============================================================================
source("run_calibration.R")     # defines .STATE and the phase functions
phase1_setup()                  # builds ce / cfg inside .STATE
phase2_check_inputs()           # verifies EPW / IDF / IDD / EPLUS / XLSX present

# pin paths explicitly (do NOT rely on cfg$work_dir, which can read empty)
PATHS <- list(runs_dir = "runs",
              xlsx_track_a = "T3_Calibration_TRACK_A_50_nights.xlsx")
stopifnot(dir.exists(PATHS$runs_dir), file.exists(PATHS$xlsx_track_a))

# =============================================================================
# PHASE 1 - publishable IDF + geometry/envelope audit (fail-fast)
# =============================================================================
# Generate the corrected, self-contained IDF (v4) from v3, then audit geometry.
# verify_trackA_geometry() stops on the first FAIL, so no downstream metric is
# ever computed on an unverified envelope.
source("make_v4.R")                  # writes T1_aula_base_94_v4.idf (verified)
source("verify_trackA_geometry.R")
verify_trackA_geometry(idf_v3 = "T1_aula_base_94_v3.idf",
                       idf_v4 = "T1_aula_base_94_v4.idf",
                       htm    = file.path(PATHS$runs_dir, "v3_trackA", "eplustbl.htm"),
                       strict = TRUE)   # -> expect 12/12 checks passed, UA = 99.3

# =============================================================================
# PHASE 2 - definitive single-zone (v3) metrics  [ANCHOR]
# =============================================================================
# trackA_5min_metrics.R normalizes Wall_3 at runtime and MUST reproduce the
# -5.41 % / 8.04 % anchor. If it returns -4.70 %, the normalization was skipped
# -> STOP and investigate before proceeding.
source("trackA_5min_metrics.R")      # NMBE/CV-RMSE (5-min + hourly), tau pooled, UA
source("tau_per_night.R")            # per-night tau (median 38.1 h) + exports
                                     #   runs/tau_per_night.csv (43 nights)

# =============================================================================
# PHASE 3 - envelope conductance UA (corrected, from the htm)
# =============================================================================
# UA_aula_v4_fixed.R REPLACES the deprecated UA_aula_v3.R (which had two bugs:
# it missed Wall_3 and inflated the interior ceiling U). Reads the authoritative
# EnergyPlus Envelope Summary and returns the direct-exterior UA = 99.3 W/K.
source("UA_aula_v4_fixed.R")
ua_aula(file.path(PATHS$runs_dir, "v3_trackA", "eplustbl.htm"))

# =============================================================================
# PHASE 4 - two-zone model (v6)  [OPTIONAL - annex negative result]
# =============================================================================
# The two-zone variant is a REJECTED candidate (it overshoots the nocturnal bias
# by a rigid +2.1 C). It is not part of the retained model, but its anchored
# series is required for Appendix A Fig A.2. close_trackA_v6.R normalizes Wall_3
# in BOTH zones and its v3 control arm must reproduce -5.41 %; it self-verifies
# (stops if any NatVentCavity/OSCM remains). Backup the previous v6 series first.
if (isTRUE(RUN_V6)) {
  if (file.exists(file.path(PATHS$runs_dir, "v6_trackA", "eplusout.csv")) &&
      !file.exists(file.path(PATHS$runs_dir, "v6_trackA", "eplusout_prev_BACKUP.csv")))
    file.copy(file.path(PATHS$runs_dir, "v6_trackA", "eplusout.csv"),
              file.path(PATHS$runs_dir, "v6_trackA", "eplusout_prev_BACKUP.csv"))
  source("close_trackA_v6.R")        # -> v6 nocturnal +6.05 % / 9.21 % (hourly)
} else {
  message("[MASTER] RUN_V6 = FALSE -> skipping two-zone model (Fig A.2 will be stale).")
}

# =============================================================================
# PHASE 5 - manuscript figures (anchored via load_calibrated_series.R)
# =============================================================================
# Defensive reset of the loader flag: after a restart it may be TRUE while D is
# a stale atomic vector. Force a clean reload so D is a proper data.frame.
rm(list = intersect(c("D", "df", "s6", ".CALIB_SERIES_LOADED"), ls()))
source("load_calibrated_series.R")   # anchors v3 (-5.41) and v6 (+6.05); [OK] if both pass

source("fig5_representative_night.R")     # Fig 5  (Section 3.1.2)
source("figA_calibration_v3v6.R")         # Fig A.2 (Appendix A) - needs v6 (Phase 4)
source("figA3_tau.R")                     # Fig A.3 (Appendix A) - needs Phase 2 CSV

message("\n[MASTER] Track A closure complete. ",
        "Metrics + tables + figures reproduced from a clean session.")

# =============================================================================
# DEPRECATED / RETIRED (do NOT re-add to the canonical flow):
#   * closure_trackA_v3_v6.R : returns -4.70 % (does NOT normalize Wall_3);
#                              superseded by figA_calibration_v3v6.R (anchored).
#   * UA_aula_v3.R           : two counting bugs -> 156 W/K; use UA_aula_v4_fixed.R.
#   * Config "crack 1.5-3 / cap 15" and "UA 200.7" and metrics "-4.70/7.65":
#                              all retired; anchored values are in the header above.
# =============================================================================
