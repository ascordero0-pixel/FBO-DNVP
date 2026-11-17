# =============================================================================
# F2/run_classic800_baseline.R  --  Sensitivity baseline: classical DNVP at 800 ppm
#
#   PURPOSE. Second classical baseline for §3.4 (sensitivity case): the threshold
#   of the field-validated protocol of [11] (COVID-19 hygiene protocol, 800 ppm),
#   with every other set-point identical to the 1000-ppm baseline
#   (Δ_hys = 100, d_slot = 5, f_open = f_max(c), τ_AQ = 6). Same evaluator, same
#   physics, same SC7-corrected FIS, same EPW. 9 annual E+ simulations (~25 min).
#
#   The optimized-protocol values are NOT re-simulated: they are read from the
#   existing F2/out/classic_vs_opt.csv so that Table 10 stays anchored to one run.
#
# HOW TO RUN (RStudio, inside FIS.Rproj):
#   1) Pause OneDrive / antivirus on C:/EP_work. Make sure no other E+ run is active.
#   2) source(here::here("F2", "run_classic800_baseline.R"))
#   Outputs: F2/out/classic800_vs_opt.csv  +  printed comparison table & summary.
#   The 1000-ppm results (classic_vs_opt.csv) are never overwritten.
# =============================================================================

library(here)

# ---- baseline choices (only CL_TAU_CO2 differs from run_classic_baseline.R) ---
CL_TAU_CO2 <- 800      # threshold of the validated DNVP of [11]
CL_DELTA   <- 100
CL_DSLOT   <- 5
CL_FOPEN   <- "fmax"
CL_TAU_AQ  <- 6
WORK_ROOT  <- "C:/EP_work/FBO_classic800"

# ---- source evaluator + scenario helpers -----------------------------------
source(here::here("F2", "metrics_bc.R"))   # VARS, ALL9, scn_hi, SF
source(here::here("F1", "fbo_eval.R"))     # evaluate_U, fbo_setup

.NORTH <- c(NE = 0, N = 329, E = 59)
scn_north <- function(scn) .NORTH[[ sub("_.*", "", scn) ]]
fmax_of   <- function(scn) as.numeric(scn_hi(scn)[4])

x_classic <- function(scn) {
  fo <- if (identical(CL_FOPEN, "fmax")) fmax_of(scn) else as.numeric(CL_FOPEN)
  setNames(list(CL_TAU_CO2, CL_DELTA, CL_DSLOT, fo, CL_TAU_AQ), VARS)
}

# ---- pre-flight ------------------------------------------------------------
cat(sprintf("[pre-flight] cost(g3=30, g4=4) = %.3f  (SC7 -> expect >=0.6)\n", scoring_cost(SF, 30, 4)))
prev <- here::here("F2", "out", "classic_vs_opt.csv")
stopifnot("classic_vs_opt.csv (1000-ppm run) not found" = file.exists(prev))
opt1000 <- read.csv(prev)

# ---- resolve IDF / EPW (explicit production file, no list.files()[1]) -------
.idf <- c(here::here("F1","v4_FBO_annual_prod.idf"), "C:/EP_work/v4_FBO_annual_prod.idf")
base_idf <- .idf[file.exists(.idf)][1]
epw <- here::here("data", "ESP_Sevilla.083910_IWEC.epw")   # production EPW (anchor check 08-10-2026)
stopifnot("production IDF not found" = !is.na(base_idf) && file.exists(base_idf),
          "IWEC EPW not found"       = file.exists(epw))
fbo_setup(base_idf = base_idf, epw = epw, work_dir = file.path(WORK_ROOT, "w1"))

# ---- evaluate the 800-ppm classical protocol in the 9 scenarios ------------
cat("\n[run] classical 800-ppm baseline over the 9 scenarios (9 annual E+ sims) ...\n")
t0 <- Sys.time()
rows <- lapply(ALL9, function(scn) {
  rc <- evaluate_U(x_classic(scn), north = scn_north(scn), tag = paste0("cls800_", scn))
  o  <- opt1000[opt1000$scn == scn, ]
  cat(sprintf("  %-8s | U opt %.4f  cls800 %.4f  ΔU %+.4f | g1 opt %.0f cls800 %.0f | p_III opt %.2f cls800 %.2f | g3 opt %.1f cls800 %.1f | g4 opt %.1f cls800 %.1f\n",
              scn, o$U_opt, rc$U, o$U_opt - rc$U, o$g1_opt, rc$g1, o$pIII_opt, rc$g2[2], o$g3_opt, rc$g3, o$g4_opt, rc$g4))
  data.frame(scn = scn,
    U_opt = o$U_opt, U_cls800 = rc$U, dU = o$U_opt - rc$U,
    g1_opt = o$g1_opt, g1_cls800 = rc$g1,
    pII_opt = o$pII_opt, pII_cls800 = rc$g2[1],
    pIII_opt = o$pIII_opt, pIII_cls800 = rc$g2[2],
    g3_opt = o$g3_opt, g3_cls800 = rc$g3,
    g4_opt = o$g4_opt, g4_cls800 = rc$g4,
    g5_opt = o$g5_opt, g5_cls800 = rc$g5,
    U_cls1000 = o$U_cls, g1_cls1000 = o$g1_cls, pIII_cls1000 = o$pIII_cls,
    g3_cls1000 = o$g3_cls, g4_cls1000 = o$g4_cls)
})
out <- do.call(rbind, rows)
write.csv(out, here::here("F2", "out", "classic800_vs_opt.csv"), row.names = FALSE)
cat(sprintf("\n[done] %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ---- summary ---------------------------------------------------------------
pct <- function(o, c) ifelse(c == 0, NA, 100 * (c - o) / c)
cat("\n==========  SUMMARY: optimized vs classical 800-ppm DNVP (and vs 1000-ppm)  ==========\n")
cat(sprintf("  optimized U > cls800 U : %d / 9   (ties: %d)\n", sum(out$dU > 1e-4), sum(abs(out$dU) <= 1e-4)))
cat(sprintf("  ΔU (opt - cls800): mean %+.4f | range %+.4f .. %+.4f\n", mean(out$dU), min(out$dU), max(out$dU)))
cat(sprintf("  mean U : opt %.3f | cls800 %.3f | cls1000 %.3f\n", mean(out$U_opt), mean(out$U_cls800), mean(out$U_cls1000)))
cat(sprintf("  Cat III compliant (p_III<=5%%): opt %d/9 | cls800 %d/9 | cls1000 %d/9\n",
            sum(out$pIII_opt <= 5), sum(out$pIII_cls800 <= 5), sum(out$pIII_cls1000 <= 5)))
cat(sprintf("  mean per-scenario change opt vs cls800 | g1 %+.1f%%  g3 %+.1f%%  g4 %+.1f%%\n",
            -mean(pct(out$g1_opt, out$g1_cls800)), -mean(pct(out$g3_opt, out$g3_cls800)), -mean(pct(out$g4_opt, out$g4_cls800))))
cat(sprintf("  cls800 vs cls1000 (mean per-scenario) | g1 %+.1f%%  g3 %+.1f%%  g4 %+.1f%%\n",
            -mean(pct(out$g1_cls800, out$g1_cls1000)), -mean(pct(out$g3_cls800, out$g3_cls1000)), -mean(pct(out$g4_cls800, out$g4_cls1000))))
cat("  file: F2/out/classic800_vs_opt.csv\n")
cat("=======================================================================================\n")
invisible(out)
