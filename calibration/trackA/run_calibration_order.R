# =====================================================================
# run_calibration_order.R — action log (run line by line)
# Do NOT source() this whole file: results are inspected between phases
# and decisions are made along the way. This is the reproducible
# bookkeeping of the calibration.
# =====================================================================

# --- Startup: engine already appended (v2 + Wall_3 patch + saveas UTF-8) ---
source("run_calibration.R")     # defines everything, runs nothing
phase1_setup()                  # link R<->Python, load engine, verify v2 ext.
phase2_check_inputs()           # EPW / IDF v3 / workbook / paths
phase3_test_run()               # ONE isolated run (circuit already validated)

# --- PILOT — 8 runs (quick validation, already passed) ----------------
resA_pilot <- phase4_track_a(crack_scale  = c(0.5, 1.0),
                             cap_mult     = c(1.0, 15.0),
                             heating_w_m2 = c(0.0, 30.0))
print(resA_pilot)

# --- FULL SWEEP — 120 runs (6 x 5 x 4) --------------------------------
# NOTE: ~50 min. OneDrive PAUSED. Do not close RStudio or let the PC sleep.
resA <- phase4_track_a()
table(resA$accept_triple_gate)                          # 99 TRUE / 21 FALSE
print(resA[resA$accept_triple_gate %in% c(TRUE, "True"), ])   # the passers

# --- PHYSICAL CHECK of infiltration (final calibration point) ---------
# Chosen point: crack=1.0 (template nominal), cap=1 (default), heat=20.
source("check_infiltration.R")
res_inf <- check_infiltration_ach(crack = 1.0, cap = 1, heat = 20)
# -> mean 0.638 / nocturnal 0.409 ach (consistent with a 1946 building)

# =====================================================================
# Track A calibration ends here. Track B (CO2 mass-balance) is a
# separate, self-contained pipeline: MASTER_trackB.R in the TRACK B
# project folder (formerly phases 5-7 of run_calibration.R, now removed).
# =====================================================================