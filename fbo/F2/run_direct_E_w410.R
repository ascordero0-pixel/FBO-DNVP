# =============================================================================
# F2/run_direct_E_w410.R  --  One-off: direct-E+ NSGA-II on the E_w410 scenario.
#
#   PURPOSE. Harden the weakest-surrogate flank. g3 in E_w410 has CV-RMSE ~22 %
#   (above the 20 % internal gate), so its surrogate (b)/(c) verdict is the one
#   most exposed to a "surrogate optimism" objection. This runs the classical
#   NSGA-II DIRECTLY on EnergyPlus for E_w410 (no surrogate), giving a physically
#   grounded (b)/(c) and an independent (a) sample-efficiency point, scored by the
#   SC7-corrected FIS. It closes E_w410 the same way N_w410 / NE_w216 are closed.
#
#   MACHINERY. Identical to the two headline runs: it just calls run_direct_ep()
#   for a different scenario. scn_box("E_w410") already applies f_max = 0.42 and
#   scn_north("E_w410") the E orientation (azimuth 59); nothing else changes.
#   Budget ≈ POP*(GEN+1) ≈ 1040 annual E+ simulations (a single extra scenario,
#   NOT the full 9 — that would contradict the sample-efficiency thesis).
#
# HOW TO RUN (RStudio, inside FIS.Rproj):
#   1) Pause OneDrive / antivirus sync on C:/EP_work (E+ writes many files).
#   2) Set NCORES below to the useful cores of your Ryzen (E+ is single-threaded:
#      one EnergyPlus worker per core).
#   3) source(here::here("F2", "run_direct_E_w410.R"))
#
#   Outputs: F2/out/de_log_E_w410.csv  and  F2/out/de_front_E_w410.csv,
#   plus a printed (a)/(b)/(c) summary ready to paste as the 3rd row of Table 9.
# =============================================================================

library(here)

# ---- config ----------------------------------------------------------------
NCORES    <- 8L                       # <-- adjust to your Ryzen (1 E+ worker per core)
POP       <- 40L                      # population (MUST be a multiple of 4 for mco::nsga2)
GEN       <- 25L                      # generations -> budget ≈ POP*(GEN+1) ≈ 1040 E+ runs
SEED      <- 1L                       # reproducibility (same as the headline runs)
WORK_ROOT <- "C:/EP_work/FBO_DE_E"    # own work root, isolated from prior DE runs
RUN_SMOKE <- TRUE                     # TRUE = cheap 16-eval coupling check first

# ---- source the direct-E+ pipeline (pulls metrics_bc.R + F1/fbo_eval.R) -----
source(here::here("F2", "nsga2_direct_ep.R"))

# ---- pre-flight: confirm the FIS in memory is the SC7-corrected one ----------
# Under SC7, (g3 High, g4 Low) -> Cost High: a high-energy / low-effort point must
# score in the HIGH cost region (≳0.6). If this prints ~0.5 you are still on the OLD
# FIS -> stop and replace core/fis_scoring.R with the corrected version before running.
.cost_HL <- scoring_cost(SF, g3 = 30, g4 = 4)
cat(sprintf("[pre-flight] cost(g3=30, g4=4) = %.3f  (SC7 corrected -> expect ≳0.6; old FIS -> ~0.5)\n", .cost_HL))
if (.cost_HL < 0.55)
  warning("cost(High,Low) looks LOW -> you may be running the pre-SC7 FIS. Check core/fis_scoring.R.")

# ---- optional smoke test (validates the E+ coupling cheaply) ----------------
if (isTRUE(RUN_SMOKE)) {
  cat("\n[smoke] E_w410 coupling check (16 E+ runs, sequential) ...\n")
  sm <- run_direct_ep("E_w410", popsize = 8L, generations = 1L, ncores = 1L,
                      work_root = WORK_ROOT, seed = SEED)
  cat("[smoke] OK -- coupling works. Proceeding to the real run.\n")
}

# ---- real run ---------------------------------------------------------------
cat(sprintf("\n[run] E_w410 direct-E+ NSGA-II  (pop=%d, gen=%d, ncores=%d) ...\n", POP, GEN, NCORES))
r_E <- run_direct_ep("E_w410", popsize = POP, generations = GEN, ncores = NCORES,
                     work_root = WORK_ROOT, seed = SEED)

# ---- one-line summary for Table 9 (3rd row) --------------------------------
with(r_E, {
  U_star <- U_target + EPS_U
  ratio  <- if (is.na(N_nsga)) sprintf(">%.1f", nrow(log) / N_fbo) else sprintf("%.1f", N_nsga / N_fbo)
  cat("\n================  E_w410  (paste into Table 9)  ================\n")
  cat(sprintf("  U(x*)          = %.4f\n", U_star))
  cat(sprintf("  (b) on front   = %s   (norm. dist %.3f)\n", bc$on_front, bc$dist_norm))
  cat(sprintf("  (c) ΔU         = %+.3f   (c_ok = %s)\n", bc$dU, bc$dU <= EPS_U))
  cat(sprintf("  N_FBO / N_NSGA = %d / %s\n", N_fbo, ifelse(is.na(N_nsga), paste0(">", nrow(log)), N_nsga)))
  cat(sprintf("  (a) ratio      = %s x\n", ratio))
  cat("  files: F2/out/de_log_E_w410.csv , F2/out/de_front_E_w410.csv\n")
  cat("===============================================================\n")
})

# ---- (optional) include E_w410 in the figures ------------------------------
# The paper's headline figures stay as the 2 bracketing scenarios (N_w410, NE_w216).
# For a supplementary/robustness figure that also shows E_w410:
#   source(here::here("F2", "figs.R"));  make_figs(c("N_w410", "NE_w216", "E_w410"))
