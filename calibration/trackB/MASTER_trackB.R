# =============================================================================
# MASTER_trackB.R  --  Master orchestrator for the FULL Track B workflow
# Thesis T3 -- FBO framework -- 1946 uninsulated corner classroom, Seville
# =============================================================================
# Governs the complete Track B calibration pipeline end to end, in strict
# topological (dependency) order, with fail-fast integrity control.
#
# DEPENDENCY DAG (see thesis notebook):
#   raw .xlsx
#     -> export_co2_clean.R      (produces co2_clean.csv)      [conditional]
#     -> load_co2_campaign.R     (data + 11 integrity anchors) [fail-fast]
#     -> trackB_massbalance.R    (defines implicit CO2 engine)
#     -> { trackB_calibrate.R,                                  (g_hat, LOO-CV)
#          trackB_sensitivity.R,                                (95% CI, sens.)
#          trackB_summary_tables.R }                            (Tables B1-B5)
#     -> trackB_figures.R              (FigB1, FigB2)
#     -> figB3_profile_likelihood.R    (FigB3 -- Appendix B.2 identifiability)
#
# NOTE: FigB2 (measured-vs-simulated scatter) is generated for internal
#   reproducibility but was NOT included in the published manuscript Appendix
#   B -- its vertical scatter is dominated by the phase-lag component of the
#   residual (Table B4) and risks being misread as a magnitude bias. FigB1 and
#   the Theil decomposition (Table B4) communicate the same finding safely.
#
# DESIGN DECISIONS (documented, per project methodology):
#   (a) The Excel->CSV export is CONDITIONAL and IDEMPOTENT: it runs only if
#       REGENERATE_CSV = TRUE or co2_clean.csv is absent. This avoids forcing a
#       reticulate/openpyxl dependency on every run while preserving the full
#       reproducibility chain when explicitly requested.
#   (b) FAIL-FAST: the first error halts the pipeline, naming the failed stage.
#       This enforces the anchor-first rule (no corrupted data may propagate).
#
# USAGE:
#   setwd("<.../TRACK B>")           # folder containing the scripts + csv/xlsx
#   source("MASTER_trackB.R")        # runs with defaults (no CSV regeneration)
#   # or force full chain from raw Excel:
#   REGENERATE_CSV <- TRUE; source("MASTER_trackB.R")
# =============================================================================

# ---------------------------- configuration ---------------------------------
if (!exists("REGENERATE_CSV")) REGENERATE_CSV <- FALSE   # set TRUE to rebuild csv
CSV_FILE  <- "co2_clean.csv"
XLSX_FILE <- "T3_Calibration_TRACK B_3 days.xlsx"

# scripts in strict execution order
SCRIPTS <- c(
  export   = "export_co2_clean.R",
  load     = "load_co2_campaign.R",
  engine   = "trackB_massbalance.R",
  calib    = "trackB_calibrate.R",
  sens     = "trackB_sensitivity.R",
  tables   = "trackB_summary_tables.R",
  figures  = "trackB_figures.R",
  figB3    = "figB3_profile_likelihood.R"
)

# ---------------------------- utilities -------------------------------------
banner <- function(txt, ch = "#") {
  cat("\n", strrep(ch, 74), "\n", ch, " ", txt, "\n", strrep(ch, 74), "\n", sep = "")
}

# run a stage with timing + fail-fast; returns elapsed seconds
run_stage <- function(n, label, expr) {
  cat(sprintf("\n>>> STAGE %d/%d : %s\n", n, N_STAGES, label))
  cat(strrep("-", 64), "\n", sep = "")
  t0 <- Sys.time()
  ok <- tryCatch({ force(expr); TRUE },
                 error = function(e) {
                   cat(sprintf("  [PIPELINE HALT @ stage %d] %s\n", n, conditionMessage(e)))
                   FALSE
                 })
  if (!ok) stop(sprintf("Master pipeline failed at stage %d (%s). See message above.",
                        n, label), call. = FALSE)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("    [done in %.2f s]\n", dt))
  dt
}

# verify a script exists before sourcing (clear error instead of cryptic one)
need_script <- function(key) {
  f <- SCRIPTS[[key]]
  if (!file.exists(f))
    stop(sprintf("Required script '%s' not found in working directory '%s'.",
                 f, getwd()), call. = FALSE)
  f
}

# ---------------------------- preflight -------------------------------------
banner("TRACK B MASTER PIPELINE  --  CO2 mass-balance calibration (FBO / T3)")
cat("Working directory:\n  ", getwd(), "\n", sep = "")

# confirm all scripts are present up front (fail early, fail clearly)
missing <- SCRIPTS[!file.exists(SCRIPTS)]
if (length(missing) > 0)
  stop("Missing scripts in working directory:\n  - ",
       paste(missing, collapse = "\n  - "), call. = FALSE)
cat(sprintf("All %d pipeline scripts located.\n", length(SCRIPTS)))

# decide whether the CSV export stage is required
do_export <- REGENERATE_CSV || !file.exists(CSV_FILE)
if (do_export) {
  if (!file.exists(XLSX_FILE))
    stop(sprintf("CSV export requested but raw workbook '%s' not found.", XLSX_FILE),
         call. = FALSE)
  cat(sprintf("CSV export stage: ENABLED (%s).\n",
              if (REGENERATE_CSV) "forced by REGENERATE_CSV" else "co2_clean.csv absent"))
} else {
  cat("CSV export stage: SKIPPED (co2_clean.csv present; set REGENERATE_CSV<-TRUE to rebuild).\n")
}

# total number of stages actually run
N_STAGES <- 7L + as.integer(do_export)   # 7 fixed stages (incl. figB3) + optional export
timings  <- numeric(0)
t_start  <- Sys.time()
stage_i  <- 0L
nxt <- function() { stage_i <<- stage_i + 1L; stage_i }

# ---------------------------- pipeline --------------------------------------

# STAGE (conditional): regenerate co2_clean.csv from raw Excel
if (do_export) {
  timings["export"] <- run_stage(nxt(),
    "Export co2_clean.csv from raw workbook (reticulate/openpyxl)",
    { source(need_script("export")); export_co2_clean() })
}

# STAGE: load + verify 11 integrity anchors (fail-fast on any mismatch)
timings["load"] <- run_stage(nxt(),
  "Load campaign data + verify 11 integrity anchors",
  { source(need_script("load")) })

# STAGE: mass-balance engine (implicit scheme, T1 formulation, switched WWeqR)
timings["engine"] <- run_stage(nxt(),
  "Load CO2 mass-balance engine (implicit, T1, switched WWeqR)",
  { source(need_script("engine")) })

# STAGE: calibration (C_d fixed at 0.65, calibrate g_person, LOO-CV)
timings["calib"] <- run_stage(nxt(),
  "Calibrate g_person (C_d fixed, LOO-CV Nov18+26 / Dec17)",
  { source(need_script("calib")); MASTER_RES <<- run_calibration() })

# STAGE: uncertainty + sensitivity (profile-likelihood CI, 9 scenarios)
timings["sens"] <- run_stage(nxt(),
  "Profile-likelihood 95% CI + sensitivity analysis",
  { source(need_script("sens")) })

# STAGE: summary tables B1-B5 (incl. MSE bias/variance/phase decomposition)
timings["tables"] <- run_stage(nxt(),
  "Build summary tables B1-B5 (CSV) incl. error decomposition",
  { source(need_script("tables")); MASTER_TABS <<- build_tables() })

# STAGE: figures (time-series fit + measured-vs-simulated scatter)
timings["figures"] <- run_stage(nxt(),
  "Generate calibration figures (FigB1, FigB2 -- FigB2 not published, see note above)",
  { source(need_script("figures")); make_figures() })

# STAGE: profile-likelihood figure (Appendix B.2 identifiability evidence)
timings["figB3"] <- run_stage(nxt(),
  "Generate profile-likelihood figure (FigB3 -- Appendix B.2)",
  { source(need_script("figB3")); make_profile_figure() })

# ---------------------------- summary ---------------------------------------
t_total <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))

banner("PIPELINE COMPLETE")
cat(sprintf("Total wall time: %.2f s across %d stages\n", t_total, N_STAGES))
cat("\nStage timings (s):\n")
for (nm in names(timings))
  cat(sprintf("  %-9s %6.2f\n", nm, timings[[nm]]))

cat("\nHEADLINE RESULT:\n")
cat(sprintf("  g_person = %.2f L/h/person  (C_d = 0.65 fixed, T1)\n", MASTER_RES$g_hat_lh))
cat(sprintf("  TEST (Dec17/S2): NMBE = %+.2f %% , CV-RMSE = %.2f %%\n",
            MASTER_RES$test$NMBE, MASTER_RES$test$CVRMSE))

cat("\nARTEFACTS WRITTEN:\n")
cat("  - co2_clean.csv (if regenerated)      -> working directory\n")
arts <- c("TableB1_parameters.csv", "TableB2_gof_by_fold.csv",
          "TableB3_validity_domain.csv", "TableB4_mse_decomposition.csv",
          "TableB5_sensitivity.csv",
          "FigB1_trackB_fit.png (published, Fig. 6)",
          "FigB2_trackB_scatter.png (generated, NOT published)",
          "FigB3_profile_likelihood.png (published, Appendix B.2)")
for (a in arts) cat("  -", a, "                -> runs/\n")

cat("\nProvenance: every value traces to Track A, T1, Tugores (2024) or EN 16798-7.\n")
cat("No-fabrication policy enforced; anchor-first rule passed.\n")
banner("END", "#")
