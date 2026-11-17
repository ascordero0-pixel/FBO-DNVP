# =============================================================================
# F2/run_classic_baseline.R  --  Optimized DNVP vs classical fixed-set-point DNVP
#
#   PURPOSE. Close the §3 comparison promised in §2.9: the Bayesian-optimized
#   protocol vs the classical fixed-1000-ppm DNVP, scenario by scenario, on every
#   objective and on the utility U. Both are evaluated by the SAME annual-simulation
#   evaluator (F1/fbo_eval.R) under the SAME physics and scoring (SC7-corrected FIS),
#   so the only difference is the set-point vector. 9 scenarios × 2 protocols = 18
#   annual E+ simulations (cheap; sequential is fine).
#
#   CLASSICAL BASELINE (fixed, the protocol the framework improves upon):
#     τ_CO2 = 1000 ppm (regulatory / corpus median / [11]) ; Δ_hys = 100 ppm ;
#     d_slot = 5 min (reactive, no slot batching) ; f_open = f_max(c) (opens fully,
#     no modulation) ; τ_AQ = 6 (no outdoor-AQ gate — AQ-blind).
#   Adjust the three tunable choices below if you prefer a different baseline.
#
# HOW TO RUN (RStudio, inside FIS.Rproj):
#   1) Pause OneDrive / antivirus on C:/EP_work.
#   2) source(here::here("F2", "run_classic_baseline.R"))
#   Outputs: F2/out/classic_vs_opt.csv  +  a printed comparison table & summary.
# =============================================================================

library(here)

# ---- tunable classical-baseline choices ------------------------------------
CL_TAU_CO2 <- 1000     # fixed CO2 opening threshold (ppm)
CL_DELTA   <- 100      # fixed hysteresis (ppm)
CL_DSLOT   <- 5        # actuation slot (min); reactive protocol
CL_FOPEN   <- "fmax"   # "fmax" = open to the scenario ceiling (no modulation); or a number in (0,1]
CL_TAU_AQ  <- 6        # 6 = never gates (AQ-blind classical)
WORK_ROOT  <- "C:/EP_work/FBO_classic"

# ---- source evaluator + scenario helpers -----------------------------------
source(here::here("F2", "metrics_bc.R"))   # VARS, ALL9, scn_hi, f2_fbo_optimum, cost_transparent, SF, EPS_U
source(here::here("F1", "fbo_eval.R"))     # evaluate_U, fbo_setup

.NORTH <- c(NE = 0, N = 329, E = 59)
scn_north <- function(scn) .NORTH[[ sub("_.*", "", scn) ]]
fmax_of   <- function(scn) as.numeric(scn_hi(scn)[4])          # f_open ceiling = f_max(c)

x_classic <- function(scn) {
  fo <- if (identical(CL_FOPEN, "fmax")) fmax_of(scn) else as.numeric(CL_FOPEN)
  setNames(list(CL_TAU_CO2, CL_DELTA, CL_DSLOT, fo, CL_TAU_AQ), VARS)
}
x_opt <- function(scn) as.list(f2_fbo_optimum(scn)$x)          # optimized set-points (= Table 7)

# ---- pre-flight: FIS must be the SC7-corrected one --------------------------
cat(sprintf("[pre-flight] cost(g3=30, g4=4) = %.3f  (SC7 -> expect >=0.6)\n", scoring_cost(SF, 30, 4)))

# ---- resolve IDF/EPW (same idiom as the direct-E+ runner) ------------------
.idf <- c(here::here("F1","v4_FBO_annual_prod.idf"), "C:/EP_work/v4_FBO_annual_prod.idf")
base_idf <- .idf[file.exists(.idf)][1]
epw <- list.files(here::here("data"), pattern = "\\.epw$", full.names = TRUE)[1]
stopifnot("production IDF not found" = !is.na(base_idf) && file.exists(base_idf))
fbo_setup(base_idf = base_idf, epw = epw, work_dir = file.path(WORK_ROOT, "w1"))

# ---- evaluate optimized + classical in the 9 scenarios ---------------------
cat("\n[run] optimized vs classical over the 9 scenarios (18 annual E+ sims) ...\n")
rows <- lapply(ALL9, function(scn) {
  no <- scn_north(scn)
  ro <- evaluate_U(x_opt(scn),     north = no, tag = paste0("opt_", scn))
  rc <- evaluate_U(x_classic(scn), north = no, tag = paste0("cls_", scn))
  cat(sprintf("  %-8s | U opt %.4f  cls %.4f  ΔU %+.4f | p_III opt %.2f cls %.2f | g3 opt %.1f cls %.1f | g4 opt %.1f cls %.1f | g5 opt %.2f cls %.2f\n",
              scn, ro$U, rc$U, ro$U - rc$U, ro$g2[2], rc$g2[2], ro$g3, rc$g3, ro$g4, rc$g4, ro$g5, rc$g5))
  data.frame(scn = scn,
    U_opt = ro$U, U_cls = rc$U, dU = ro$U - rc$U,
    g1_opt = ro$g1, g1_cls = rc$g1,
    pII_opt = ro$g2[1], pII_cls = rc$g2[1],
    pIII_opt = ro$g2[2], pIII_cls = rc$g2[2],
    g3_opt = ro$g3, g3_cls = rc$g3,
    g4_opt = ro$g4, g4_cls = rc$g4,
    g5_opt = ro$g5, g5_cls = rc$g5)
})
out <- do.call(rbind, rows)
DE_OUT <- here::here("F2", "out"); if (!dir.exists(DE_OUT)) dir.create(DE_OUT, recursive = TRUE)
write.csv(out, file.path(DE_OUT, "classic_vs_opt.csv"), row.names = FALSE)

# ---- summary (per-objective improvement, compliance) -----------------------
pct <- function(o, c) ifelse(c == 0, NA, 100 * (c - o) / c)     # % reduction of a cost-type indicator (opt vs cls)
cat("\n================  SUMMARY: optimized vs classical DNVP  ================\n")
cat(sprintf("  scenarios where optimized U > classical U : %d / 9\n", sum(out$dU > 0)))
cat(sprintf("  ΔU (opt - cls): mean %+.4f | range %+.4f .. %+.4f\n", mean(out$dU), min(out$dU), max(out$dU)))
cat(sprintf("  Cat III compliant (p_III<=5%%): optimized %d/9 | classical %d/9\n",
            sum(out$pIII_opt <= 5), sum(out$pIII_cls <= 5)))
cat(sprintf("  mean reduction opt vs cls  |  g1(IAQ) %.1f%%  g3(energy) %.1f%%  g4(effort) %.1f%%  g5(import) %.1f%%\n",
            mean(pct(out$g1_opt, out$g1_cls), na.rm=TRUE), mean(pct(out$g3_opt, out$g3_cls), na.rm=TRUE),
            mean(pct(out$g4_opt, out$g4_cls), na.rm=TRUE), mean(pct(out$g5_opt, out$g5_cls), na.rm=TRUE)))
cat("  file: F2/out/classic_vs_opt.csv  (per-scenario, both protocols, all indicators)\n")
cat("=======================================================================\n")
invisible(out)
