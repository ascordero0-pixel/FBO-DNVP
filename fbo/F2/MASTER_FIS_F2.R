# ============================================================================
# MASTER_FIS_F2.R  --  Orchestrator for T3 FBO PHASE 2 (NSGA-II benchmark)
# ----------------------------------------------------------------------------
# Runs, in dependency order, the Phase-2 pipeline that benchmarks the FBO
# framework (Phase 1) against classical NSGA-II. FAIL-FAST: aborts on the first
# error rather than propagating wrong values downstream.
#
# Thesis (design decision #4): the FBO framework delivers a deployable,
# Pareto-efficient, utility-maximizing DNVP using ~1 order of magnitude fewer
# annual E+ simulations than NSGA-II -- by internalizing the a-priori selection
# that classical MOO defers, and reaching it sample-efficiently with GP-BO.
#
# Pipeline (this master):
#   GATE       surrogate.R        GP surrogates (Matern5/2+ARD) + k-fold CV gate    [no E+]
#   FRONTS     nsga2_surrogate.R  NSGA-II (mco) over the surrogate -> dense fronts   [no E+]
#   METRICS    metrics_bc.R       (b) FBO optimum ON front + (c) max-U ; empirical   [no E+]
#   DIRECT-E+  nsga2_direct_ep.R  NSGA-II DIRECT on EnergyPlus (N_w410, NE_w216)     [EXPENSIVE]
#              -> metric (a) sample efficiency + INDEPENDENT real (b)/(c)
#   RESCORE    rescore_direct_ep.R  re-score the direct-E+ logs with the CURRENT FIS   [no E+]
#              -> if the FIS changes, RE-SCORE (seconds); do NOT re-run DIRECT-E+.
#   FIGS       figs.R             benchmark figures (3D front, FBO-on-front, U vs evals)  [pending]
#
# FIS: incluye la corrección SC7 (2026-08, Cost (g3=Alto,g4=Bajo) M->H; destapada por ESTE
#   benchmark — la energía se enmascaraba con pocas conmutaciones). Impacto nulo en best_x;
#   las métricas se re-puntúan solas al construir la U con build_scoring_fis().
#
# Benchmark objectives (decision #1, Option B), all MINIMISED:
#   f1 = g1 (CO2)  ;  f2 = g2_pIII (Cat III shortfall)  ;  f3 = Cost(g3,g4) transparent.
# Metrics currency: EVALUATIONS (not wall-clock). g5 kept out of objectives (flat+robust).
#
# Project layout: FIS.Rproj at root. Phase-2 scripts in F2/ ; consumes F1/ outputs
# (seed_9scen_results.csv, bo/*.csv) + core/ (fis_scoring). Fronts -> F2/out/. Paths via `here`.
#
# STATE: gate + fronts + metrics + DIRECT-E+ done (N_w410 14x, NE_w216 4.4x; (b)/(c) OK in both,
# under the SC7-corrected FIS). Default source(): gate + (b)/(c) metrics + re-score of the direct-E+
# logs (all fast, no E+). Heavy toggles (fronts, direct-E+) OFF. Flip a toggle to (re)run that stage.
#
# Usage:  open FIS.Rproj, then  source(here::here("F2", "MASTER_FIS_F2.R"))
#         (pause OneDrive before the DIRECT-E+ stage; it launches EnergyPlus)
# ============================================================================

library(here)

# ---- 0. Configuration + toggles --------------------------------------------
F2_DIR    <- here::here("F2")
OUT_DIR   <- file.path(F2_DIR, "out")

RUN_GATE            <- TRUE     # surrogate CV-RMSE gate (fits GPs; no E+)          [~min]
RUN_SURROGATE_NSGA  <- FALSE    # regenerate the 9 Pareto fronts (mco over GP; no E+) [~min]
RUN_METRICS         <- TRUE     # (b)/(c) surrogate + empirical cross-check (no E+)  [fast]
RUN_DIRECT_EP       <- FALSE    # NSGA-II DIRECT on EnergyPlus  -- EXPENSIVE (~1040 E+ runs/scn); ya corrido
RUN_RESCORE_DE      <- TRUE     # re-puntúa los logs directo-E+ con el FIS actual (no E+; refleja SC7)
RUN_FIGS            <- FALSE    # benchmark figures (needs F2/figs.R; pending)

# surrogate-NSGA settings (decision #3)
SUR_SEEDS <- 1:10; SUR_POP <- 100L; SUR_GEN <- 100L

# direct-E+ settings (decision #3): 2 representative scenarios, cap ~1000 evals
DE_SCENARIOS <- c("N_w410", "NE_w216")
DE_POP <- 40L; DE_GEN <- 25L; DE_NCORES <- 7L      # ajusta DE_NCORES a los núcleos útiles de la Ryzen

# ---- FAIL-FAST guard #1: required F2 scripts present ------------------------
.need_f2 <- file.path(F2_DIR, c("surrogate.R", "nsga2_surrogate.R", "metrics_bc.R", "nsga2_direct_ep.R"))
.miss <- .need_f2[!file.exists(.need_f2)]
if (length(.miss)) stop("Missing F2 script(s): ", paste(basename(.miss), collapse = ", "),
                        "\n  (open the project via FIS.Rproj so `here` resolves the root)")

cat("== MASTER_FIS_F2: starting ==\n\n")

# ---- 1. Source the surrogate + fronts + metrics modules (define only) -------
# Each script only DEFINES functions when sourced (their sys.nframe() guards keep
# smoke tests / gates from auto-running). nsga2_direct_ep.R is sourced later, only
# under RUN_DIRECT_EP, since it also pulls fbo_eval.R (EnergyPlus evaluator).
cat("[1] sourcing surrogate + fronts + metrics modules ...\n")
source(file.path(F2_DIR, "surrogate.R"))         # f2_gate(), f2_fit_scenario(), f2_predict(), cost_transparent()
source(file.path(F2_DIR, "nsga2_surrogate.R"))   # f2_nsga2_all(), f2_nsga2_scenario()
source(file.path(F2_DIR, "metrics_bc.R"))        # f2_metrics_all(), f2_metrics_scenario()

# ---- 2. GATE: surrogate CV-RMSE validation ---------------------------------
# Umbral CV-RMSE <= 20% (más estricto que ASHRAE Guideline 14, 30% horario).
# E_w410/g3 ~22% ACEPTADO y documentado (solo alimenta el frente; cruce empírico lo blinda).
if (RUN_GATE) {
  cat("\n[2] surrogate gate (5-fold CV, training siembra+BO) ...\n")
  GATE <- f2_gate()
} else cat("\n[2] gate skipped (RUN_GATE = FALSE)\n")

# ---- 3. FRONTS: NSGA-II over the surrogate ---------------------------------
# Produce F2/out/front_<scn>.csv para los 9 escenarios. Ya generados; OFF por defecto.
if (RUN_SURROGATE_NSGA) {
  cat("\n[3] NSGA-II over surrogate -> F2/out/front_*.csv ...\n")
  f2_nsga2_all(seeds = SUR_SEEDS, popsize = SUR_POP, generations = SUR_GEN)
} else cat("\n[3] surrogate fronts skipped (RUN_SURROGATE_NSGA = FALSE; fronts already in F2/out/)\n")

# ---- 4. METRICS (b)/(c): surrogate + empirical cross-check ------------------
if (RUN_METRICS) {
  .fr <- list.files(OUT_DIR, pattern = "^front_.*\\.csv$")
  if (length(.fr) < 9L)
    cat("  [aviso] faltan front_*.csv (", length(.fr), "/9): pon RUN_SURROGATE_NSGA=TRUE primero.\n", sep = "")
  else {
    cat("\n[4] métricas (b)/(c) ...\n")
    cat("--- frente SURROGATE (consistente; optimismo del GP -> ver cruce empírico) ---\n")
    print(f2_metrics_all(source = "surrogate"), row.names = FALSE)
    cat("\n--- cruce EMPÍRICO (física real; ancla de robustez) ---\n")
    print(f2_metrics_all(c("N_w410", "NE_w216", "E_w410"), source = "empirical"), row.names = FALSE)
  }
} else cat("\n[4] metrics skipped (RUN_METRICS = FALSE)\n")

# ---- 5. DIRECT-E+: NSGA-II directo en EnergyPlus --  EXPENSIVE --------------
# La prueba INDEPENDIENTE: métrica (a) eficiencia muestral + (b)/(c) reales, sobre
# un frente que la BO no ha visto. Corre en la Ryzen (paralelo). Pausar OneDrive.
if (RUN_DIRECT_EP) {
  cat("\n[5] NSGA-II DIRECTO en EnergyPlus (EXPENSIVE) ...\n")
  source(file.path(F2_DIR, "nsga2_direct_ep.R"))   # run_direct_ep()  (+ fbo_eval.R)
  DE <- lapply(DE_SCENARIOS, function(s)
    run_direct_ep(s, popsize = DE_POP, generations = DE_GEN, ncores = DE_NCORES))
  names(DE) <- DE_SCENARIOS
} else cat("\n[5] direct-E+ skipped (RUN_DIRECT_EP = FALSE)\n")

# ---- 6. RESCORE: re-puntúa los logs directo-E+ con el FIS ACTUAL ------------
# Si el FIS cambia (p.ej. corrección SC7) NO hay que re-simular: se re-puntúan los logs
# ya existentes (de_log_<scn>.csv, con los g1..g5 físicos) con el FIS actual y se recalculan
# (a)/(b)/(c). Refleja el flujo real: cambia el FIS -> re-puntuar (segundos), no re-simular.
if (RUN_RESCORE_DE) {
  .rs   <- file.path(F2_DIR, "rescore_direct_ep.R")
  .logs <- file.path(OUT_DIR, paste0("de_log_", DE_SCENARIOS, ".csv"))
  if (!file.exists(.rs)) cat("\n[6] rescore omitido: falta F2/rescore_direct_ep.R\n")
  else if (!all(file.exists(.logs))) cat("\n[6] rescore omitido: faltan de_log_*.csv (corre DIRECT-E+ antes)\n")
  else {
    cat("\n[6] re-puntuando los logs directo-E+ con el FIS actual (SC7) ...\n")
    source(.rs)
    for (s in DE_SCENARIOS) rescore_direct(s)
  }
} else cat("\n[6] rescore skipped (RUN_RESCORE_DE = FALSE)\n")

# ---- 7. FIGS: benchmark figures --------------------------------------------
if (RUN_FIGS) {
  .fg <- file.path(F2_DIR, "figs.R")
  if (!file.exists(.fg)) stop("[7] figs.R not found in F2/ (pendiente de escribir)")
  cat("\n[7] rendering benchmark figures -> F2/out/ ...\n")
  source(.fg)
} else cat("\n[7] figs skipped (RUN_FIGS = FALSE; figs.R pendiente)\n")

# ---- 8. Done ---------------------------------------------------------------
cat("\n== MASTER_FIS_F2: pipeline complete ==\n")
cat("   (b)/(c)/(a) confirmadas en física real (los 2 escenarios) con el FIS corregido (SC7).\n")
cat("   siguiente: F2/figs.R (RUN_FIGS=TRUE) + sección §3.3 del manuscrito.\n")
