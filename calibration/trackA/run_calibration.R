# =====================================================================
# run_calibration.R  —  T3 Digital-Twin Calibration · MASTER ORCHESTRATOR
# ---------------------------------------------------------------------
# Seville classroom (T1 model). EnergyPlus 9.4.0 / eppy / reticulate.
# IDF base = v3 (corrected geometry). Track A uses the v2 engine
# extension (TrackAGridV2 / calibrate_track_a_v2) appended to
# calibration_engine.py.
#
# USAGE (manual, one phase at a time):
#   source("run_calibration.R")
#   phase1_setup()
#   phase2_check_inputs()
#   phase3_test_run()
#   # pilot (8 runs):
#   resA <- phase4_track_a(crack_scale=c(0.5,1.0), cap_mult=c(1.0,15.0),
#                          heating_w_m2=c(0.0,30.0))
#   # full grid (120 runs):
#   resA <- phase4_track_a()
#
# Track B (CO2) is calibrated separately: see MASTER_trackB.R in the
# TRACK B project folder. It is not orchestrated from this file.
#
# To move machines: edit ONLY the PATHS block below.
# =====================================================================


# ---------------------------------------------------------------------
# CONFIG — all paths in one place. Edit here and nowhere else.
# ---------------------------------------------------------------------
PATHS <- list(
  python_exe = "C:/Users/usuario/AppData/Local/Programs/Python/Python313/python.exe",
  engine_dir = ".",                                   # holds calibration_engine.py
  idd        = "C:/EnergyPlusV9-4-0/Energy+.idd",
  eplus_exe  = "C:/EnergyPlusV9-4-0/energyplus.exe",
  epw        = "weather/Seville_measured_NovDec2021.epw",
  idf_base   = "T1_aula_base_94_v3.idf",              # <-- v3 (corrected geometry)
  # OUTPUTS: runs/ inside the project. Works because OneDrive sync is PAUSED.
  # If OneDrive resumes, E+ temp files may get locked mid-write — re-pause
  # sync, or switch this to a non-OneDrive path like "C:/EP_runs/T3".
  runs_dir   = "runs",
  xlsx_track_a = "T3_Calibration_TRACK_A_50_nights.xlsx"
)

# Internal handles populated by phase1 (do not edit).
.STATE <- new.env()


# ---------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------
.ok   <- function(msg) cat("  [OK] ", msg, "\n", sep = "")
.info <- function(msg) cat("  [..] ", msg, "\n", sep = "")
.die  <- function(msg) stop("\n  [STOP] ", msg, "\n", call. = FALSE)
.head <- function(msg) cat("\n=== ", msg, " ===\n", sep = "")

.is_onedrive <- function(p) grepl("onedrive", p, ignore.case = TRUE)

# Robust pandas -> R data.frame conversion.
# reticulate sometimes hands back a live pandas object instead of an R
# data.frame; as.data.frame() then fails. This converts via to_dict("list")
# (a plain dict of column->list, which py_to_r always converts cleanly).
.py_df_to_r <- function(obj) {
  if (inherits(obj, "data.frame")) return(obj)
  d <- tryCatch(reticulate::py_to_r(obj$to_dict("list")), error = function(e) NULL)
  if (!is.null(d)) {
    d <- lapply(d, function(col) {
      v <- unlist(col, use.names = FALSE)
      if (is.null(v)) NA else v
    })
    return(as.data.frame(d, stringsAsFactors = FALSE, optional = TRUE))
  }
  as.data.frame(obj, stringsAsFactors = FALSE)   # last resort
}


# ---------------------------------------------------------------------
# PHASE 1 — environment: link, engine, metrics self-test
# ---------------------------------------------------------------------
phase1_setup <- function() {
  .head("PHASE 1 — environment setup")

  if (!requireNamespace("reticulate", quietly = TRUE))
    .die("Package 'reticulate' not installed. Run install.packages('reticulate').")
  library(reticulate)

  if (Sys.getenv("RETICULATE_PYTHON") == "") {
    if (!file.exists(PATHS$python_exe))
      .die(paste0("Python not found at: ", PATHS$python_exe))
    use_python(PATHS$python_exe, required = TRUE)
  }

  cfg_py <- py_config()
  .info(paste0("Python: ", cfg_py$python))
  if (!grepl("^3\\.13", cfg_py$version)) .die("Python is not 3.13 — wrong interpreter.")
  for (m in c("eppy", "numpy", "pandas"))
    if (!py_module_available(m)) .die(paste0("Python module '", m, "' not visible to R."))
  .ok("Python 3.13 + eppy + numpy + pandas all visible.")

  if (!file.exists(file.path(PATHS$engine_dir, "calibration_engine.py")))
    .die(paste0("calibration_engine.py not found in: ", PATHS$engine_dir))
  ce <- import_from_path("calibration_engine", path = PATHS$engine_dir)

  # Sanity: confirm the Track A v2 extension is present in the engine.
  if (!py_has_attr(ce, "calibrate_track_a_v2") || !py_has_attr(ce, "TrackAGridV2"))
    .die(paste0("Engine loaded but the Track A v2 extension is missing ",
                "(TrackAGridV2 / calibrate_track_a_v2). Append ",
                "calibration_engine_v2_trackA_append.py to calibration_engine.py."))

  cfg <- ce$EPConfig(idd_path = PATHS$idd, eplus_exe = PATHS$eplus_exe,
                     epw_path = PATHS$epw, work_dir = PATHS$runs_dir)
  .ok(paste0("Engine loaded (v3 + Track A v2 ext). E+ ", cfg$engine_version,
             " / DesignBuilder ", cfg$designbuilder_version))

  ce$`_selftest`()
  .ok("Metrics self-test passed (sign convention + zero-error).")

  assign("ce", ce, envir = .STATE)
  assign("cfg", cfg, envir = .STATE)
  invisible(TRUE)
}


# ---------------------------------------------------------------------
# PHASE 2 — verify all inputs exist and paths are sane
# ---------------------------------------------------------------------
phase2_check_inputs <- function() {
  .head("PHASE 2 — input verification")
  if (!exists("ce", envir = .STATE)) .die("Run phase1_setup() first.")

  req <- c(EPW = PATHS$epw, IDF = PATHS$idf_base,
           IDD = PATHS$idd, EPLUS = PATHS$eplus_exe,
           XLSX_A = PATHS$xlsx_track_a)
  for (nm in names(req)) {
    if (!file.exists(req[[nm]])) .die(paste0(nm, " not found at: ", req[[nm]]))
    .ok(paste0(nm, " present: ", req[[nm]]))
  }

  if (!dir.exists(PATHS$runs_dir)) {
    .info(paste0("Creating runs dir: ", PATHS$runs_dir))
    dir.create(PATHS$runs_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (.is_onedrive(normalizePath(PATHS$runs_dir, mustWork = FALSE)))
    .info(paste0("runs_dir is inside OneDrive (", PATHS$runs_dir,
                 ") — OK ONLY while sync is PAUSED."))
  else
    .ok(paste0("runs_dir outside OneDrive: ", PATHS$runs_dir))

  invisible(TRUE)
}


# ---------------------------------------------------------------------
# PHASE 3 — one isolated EnergyPlus run (smoke test of the whole circuit)
# NOTE: this does NOT apply the Track A patches (no heating / no EMS-off /
# no Wall_3 normalisation) — it only verifies the IDF->E+->read circuit.
# ---------------------------------------------------------------------
phase3_test_run <- function(crack_scale = 1.0, cap_mult = 1.0) {
  .head("PHASE 3 — isolated test run")
  if (!exists("cfg", envir = .STATE)) .die("Run phase1_setup() first.")
  ce  <- get("ce",  envir = .STATE)
  cfg <- get("cfg", envir = .STATE)

  .info(paste0("Editing IDF (crack_scale=", crack_scale,
               ", cap_mult=", cap_mult, ") and launching EnergyPlus..."))
  idf <- ce$load_idf(cfg, PATHS$idf_base)
  ce$set_timestep(idf, 12L)
  ce$scale_afn_crack_coeffs(idf, crack_scale)
  ce$set_capacitance_multiplier(idf, cap_mult)
  ce$add_output_variables(idf)
  run_idf <- ce$save_idf(idf, file.path(PATHS$runs_dir, "TEST_run.idf"))

  out_dir <- ce$run_energyplus(cfg, run_idf, "TEST_run")
  .ok(paste0("EnergyPlus run completed. Output dir: ", out_dir))

  sim <- ce$read_eplus_variable(out_dir, "Zone Mean Air Temperature")
  # robust read: avoid forcing the whole pandas frame
  if (inherits(sim, "data.frame")) {
    sim_r <- sim
  } else {
    sim_r <- data.frame(
      timestamp = as.character(py_to_r(sim$timestamp$astype("str")$tolist())),
      value     = as.numeric(py_to_r(sim$value$tolist())),
      stringsAsFactors = FALSE)
  }
  .ok(paste0("Read Zone Mean Air Temperature: ", nrow(sim_r), " timesteps."))
  cat("  First rows:\n"); print(utils::head(sim_r, 3))
  .ok("Full circuit verified: IDF edit -> E+ run -> read. Ready for sweeps.")
  invisible(sim_r)
}


# ---------------------------------------------------------------------
# PHASE 4 — Track A (thermal) grid sweep, v2 engine (+heating, triple gate)
# Pass crack_scale / cap_mult / heating_w_m2 to override the default grid
# (e.g. for a small pilot). Leave them NULL for the full 6x5x4 = 120 sweep.
# ---------------------------------------------------------------------
phase4_track_a <- function(p = 1L, eps_slope = 0.06, radiant_fraction = 0.30,
                           crack_scale = NULL, cap_mult = NULL, heating_w_m2 = NULL) {
  .head("PHASE 4 — Track A thermal calibration (v2: +heating, triple gate)")
  if (!exists("cfg", envir = .STATE)) .die("Run phase1_setup() first.")
  ce  <- get("ce",  envir = .STATE); cfg <- get("cfg", envir = .STATE)
  if (!requireNamespace("readxl", quietly = TRUE)) .die("Install 'readxl'.")

  # --- measured indoor T (unchanged loader) ---
  .info("Loading measured indoor T (Track A workbook)...")
  d <- readxl::read_excel(PATHS$xlsx_track_a, sheet = 1, skip = 6, col_names = FALSE)
  ncol_read <- ncol(d)
  off <- if (ncol_read >= 22) 0L else -1L     # 0 if col A kept, -1 if dropped
  cY <- 2L + off; cM <- 3L + off; cD <- 4L + off; cT <- 5L + off; cTi <- 7L + off
  .info(paste0("Workbook read with ", ncol_read, " columns (offset ", off,
               "); Ti from column index ", cTi, "."))
  Year  <- d[[cY]]; Month <- d[[cM]]; Day <- d[[cD]]; Time <- d[[cT]]
  .ffill_vec <- function(x) { for (i in seq_along(x)) if (is.na(x[i]) && i > 1) x[i] <- x[i-1]; x }
  Year  <- .ffill_vec(Year); Month <- .ffill_vec(Month)
  hm <- format(Time, "%H:%M")
  bad <- is.na(hm) | hm == "NA"
  if (any(bad)) hm[bad] <- substr(as.character(Time[bad]), 12, 16)
  hour <- as.integer(substr(hm, 1, 2)); minute <- as.integer(substr(hm, 4, 5))
  ts <- as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00",
                           as.integer(Year), as.integer(Month), as.integer(Day),
                           hour, minute), tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  if (any(is.na(ts))) .die(paste0(sum(is.na(ts)),
                                  " timestamps failed to parse — inspect workbook cols B/C/D/E."))
  measured_T <- data.frame(timestamp = format(ts, "%Y-%m-%d %H:%M:%S"),
                           value = as.numeric(d[[cTi]]), stringsAsFactors = FALSE)
  measured_T <- measured_T[!is.na(measured_T$value), ]
  .ok(paste0("Measured T: ", nrow(measured_T), " points, ",
             format(min(ts)), " to ", format(max(ts)), "."))

  # --- build the grid (full by default; partial = pilot) ---
  gargs <- list()
  if (!is.null(crack_scale))  gargs$crack_scale  <- crack_scale
  if (!is.null(cap_mult))     gargs$cap_mult     <- cap_mult
  if (!is.null(heating_w_m2)) gargs$heating_w_m2 <- heating_w_m2
  grid <- if (length(gargs)) do.call(ce$TrackAGridV2, gargs) else ce$TrackAGridV2()
  ncomb <- tryCatch(
    length(py_to_r(grid$crack_scale)) *
    length(py_to_r(grid$cap_mult)) *
    length(py_to_r(grid$heating_w_m2)),
    error = function(e) NA_integer_)
  .info(paste0("Track A v2 sweep: ", ncomb,
               " grid points (crack_scale x cap_mult x heating). ",
               "Estimate total time = ncomb x (one-run time from phase3)."))

  resA <- ce$calibrate_track_a_v2(
    cfg = cfg, base_idf_path = PATHS$idf_base,
    measured_T = r_to_py(measured_T), grid = grid,
    p = p, eps_slope = eps_slope, radiant_fraction = radiant_fraction,
    normalize_wall3 = TRUE,
    results_csv = file.path(PATHS$runs_dir, "results_trackA_v2.csv"))

  # --- read results from the CSV the engine wrote (robust; no pandas coerce) ---
  resA_r <- read.csv(file.path(PATHS$runs_dir, "results_trackA_v2.csv"),
                     stringsAsFactors = FALSE)
  .ok(paste0("Track A v2 complete: ", nrow(resA_r), " grid points."))

  pass <- resA_r[resA_r$accept_triple_gate %in% c(TRUE, "TRUE", "True"), ]
  .ok(paste0(nrow(pass), " point(s) pass the triple gate ",
             "(|NMBE|<=10 & CVRMSE<=30 & |slope_err|<=", eps_slope, ")."))
  show_cols <- intersect(
    c("tag","crack_scale","cap_mult","heating_w_m2",
      "NMBE_hourly","CVRMSE_hourly","slope_abs_err_Cph","accept_triple_gate"),
    names(resA_r))
  print(utils::head(resA_r[, show_cols], 12))
  invisible(resA_r)
}


# =====================================================================
# End of Track A calibration orchestrator.
# Track B (CO2 mass-balance) is now a SEPARATE, self-contained pipeline
# (MASTER_trackB.R in the TRACK B project folder). It is no longer part
# of this file: the former phase5_track_b / phase6_loo_cv / phase7_export
# have been removed to keep Track A single-purpose.
# =====================================================================
cat("run_calibration.R (Track A, v3 + v2 engine) loaded. Functions defined, nothing executed.\n")
cat("Call: phase1_setup() -> phase2_check_inputs() -> phase3_test_run() -> phase4_track_a(...)\n")
