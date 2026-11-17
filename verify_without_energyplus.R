# =============================================================================
# verify_without_energyplus.R -- one-command reproducibility check
# -----------------------------------------------------------------------------
# Re-runs every stage of the pipeline that does NOT require EnergyPlus and
# compares the regenerated outputs with the released ones (anchor-first rule:
# no released number is trusted unless it is reproduced from the code).
#
# Usage (from the repository root):   Rscript verify_without_energyplus.R
# Runtime: ~2 min. The repository is copied to a temporary folder first, so
# no released file is modified.
#
# Required R packages: here, ggplot2, lubridate, readxl, patchwork, scales,
#                      DiceKriging, mco (DiceOptim only for re-running BO).
# =============================================================================

REPO <- normalizePath(getwd())
stopifnot("run from the repository root" = file.exists(file.path(REPO, "fbo", "FIS.Rproj")))
RSCRIPT <- file.path(R.home("bin"), "Rscript")
WORK <- file.path(tempdir(), "fbo_dnvp_verify")
unlink(WORK, recursive = TRUE); dir.create(WORK)
invisible(file.copy(list.files(REPO, full.names = TRUE, all.files = TRUE, no.. = TRUE),
                    WORK, recursive = TRUE))

results <- data.frame(check = character(), status = character(), detail = character())
add <- function(check, ok, detail = "") {
  results[nrow(results) + 1, ] <<- list(check, if (isTRUE(ok)) "PASS" else "FAIL", detail)
  cat(sprintf("[%s] %s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", check, detail))
}
run <- function(dir, expr) {
  out <- suppressWarnings(system2(RSCRIPT, c("-e", shQuote(expr)), stdout = TRUE, stderr = TRUE))
  st <- attr(out, "status"); list(ok = is.null(st) || st == 0, log = out)
}
in_dir <- function(dir, expr) { old <- setwd(dir); on.exit(setwd(old)); run(dir, expr) }
same_csv <- function(a, b, cols = NULL, tol = 1e-6) {
  x <- read.csv(a); y <- read.csv(b)
  if (!is.null(cols)) { x <- x[cols]; y <- y[cols] }
  if (!identical(dim(x), dim(y))) return(FALSE)
  num <- vapply(x, is.numeric, logical(1))
  all(abs(as.matrix(x[num]) - as.matrix(y[num])) <= tol, na.rm = TRUE) &&
    identical(as.matrix(x[!num]), as.matrix(y[!num]))
}

FBO <- file.path(WORK, "fbo"); TA <- file.path(WORK, "calibration", "trackA")
TB  <- file.path(WORK, "calibration", "trackB")

# ---- 1. Phase 0: FIS engine + indicators + control law (148 unit checks) ----
r <- in_dir(FBO, 'source(here::here("F0","MASTER_FIS.R"))')
n_pass <- sum(as.integer(sub(".*RESULT: ([0-9]+) passed, 0 failed.*", "\\1",
                             grep("RESULT: [0-9]+ passed, 0 failed", r$log, value = TRUE))))
add("Phase 0 unit tests (F0/MASTER_FIS.R)", r$ok && n_pass == 148, sprintf("%d/148 checks", n_pass))

# ---- 2. Phase 1: best_x rebuilt from the LHS seed + BO logs (Table 4) -------
ref_bx <- read.csv(file.path(REPO, "fbo", "F1", "bo", "best_x_9scen_v2.csv"))
r <- in_dir(FBO, 'source(here::here("F1","MASTER_FIS_F1.R"))')
bx <- read.csv(file.path(FBO, "F1", "bo", "best_x_9scen_v2.csv"))
m <- merge(ref_bx[c("scn", "U")], bx[c("scn", "U")], by = "scn")
add("Phase 1 optima U* (F1/MASTER_FIS_F1.R -> Table 4)", r$ok && nrow(m) == 9 &&
      all(abs(m$U.x - m$U.y) < 5e-5), sprintf("max|dU| = %.1e", max(abs(m$U.x - m$U.y))))

# ---- 3. Phase 2: surrogate gate, (b)/(c) metrics, direct-E+ re-scoring -------
r <- in_dir(FBO, 'source(here::here("F2","MASTER_FIS_F2.R"))')
add("Phase 2 orchestrator (F2/MASTER_FIS_F2.R)", r$ok)
r <- in_dir(FBO, 'source(here::here("F2","rescore_direct_ep.R"))')
add("Table 5 (rescore of direct-E+ NSGA-II logs)", r$ok &&
      same_csv(file.path(REPO, "fbo/F2/out/rescore_direct_summary.csv"),
               file.path(FBO, "F2/out/rescore_direct_summary.csv")))
r <- in_dir(FBO, 'source(here::here("F2","sens_g4_anchoring.R"))')
add("Table D1 (g4 anchoring sensitivity)", r$ok &&
      same_csv(file.path(REPO, "fbo/F2/out/sens_g4_anchoring.csv"),
               file.path(FBO, "F2/out/sens_g4_anchoring.csv")))

# ---- 4. Track B: CO2 mass-balance calibration (Tables B1-B5, A2) ------------
r <- in_dir(TB, 'source("MASTER_trackB.R")')
tabs <- list.files(file.path(REPO, "calibration/trackB/runs"), pattern = "^Table.*\\.csv$")
ok_tb <- vapply(tabs, function(f) same_csv(file.path(REPO, "calibration/trackB/runs", f),
                                           file.path(TB, "runs", f)), logical(1))
add("Track B calibration (MASTER_trackB.R)", r$ok && all(ok_tb),
    sprintf("%d/%d tables identical; G = 14.74 L/h", sum(ok_tb), length(ok_tb)))

# ---- 5. Track A: anchors from the archived calibrated EnergyPlus series ------
r <- in_dir(TA, paste(
  'source("unpack_runs.R");',
  'source("verify_trackA_geometry.R");',
  'verify_trackA_geometry(idf_v3="T1_aula_base_94_v3.idf", idf_v4="T1_aula_base_94_v4.idf",',
  '  htm=file.path("runs","v3_trackA","eplustbl.htm"), strict=TRUE);',
  'source("load_calibrated_series.R"); source("deviation_nocturnal.R"); source("figA3_tau.R")'))
add("Track A geometry audit (12 checks, UA = 99.3 W/K)", r$ok &&
      any(grepl("RESULT: 12 / 12 checks passed", r$log)))
add("Track A anchors v3 -5.41 % / 8.04 %, v6 +6.05 %", r$ok &&
      any(grepl("both series reproduce their validated anchors", r$log)))
add("Track A weekday/weekend split (-4.94 % / -6.03 %)", r$ok &&
      any(grepl("weekday nights  v3: NMBE -4.94%", r$log, fixed = TRUE)) &&
      any(grepl("weekend nights  v3: NMBE -6.03%", r$log, fixed = TRUE)))
add("Track A time constant (median tau = 38.09 h, 43 nights)", r$ok &&
      any(grepl("n=43 nights | median tau=38.09 h", r$log, fixed = TRUE)))

cat("\n================ SUMMARY ================\n")
print(results, row.names = FALSE, right = FALSE)
cat(sprintf("\n%d/%d checks passed. Work folder: %s\n",
            sum(results$status == "PASS"), nrow(results), WORK))
if (any(results$status != "PASS")) quit(status = 1)
