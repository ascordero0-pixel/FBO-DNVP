# =============================================================================
# F2/nsga2_surrogate.R  --  T3 FBO Phase 2, step 2: NSGA-II over the SURROGATE
#
#   Runs classical NSGA-II (mco::nsga2, SBX crossover + polynomial mutation) on the
#   GP surrogates from F2/surrogate.R, per scenario and per random seed, to trace the
#   dense Pareto front in the 3 benchmark objectives (design decision #1, Option B):
#       f1 = g1 (CO2)  ·  f2 = g2_pIII (Cat III shortfall)  ·  f3 = Cost(g3,g4)
#   all MINIMISED. Each front point is also scored with the FIS utility U (from the
#   predicted indicators via core/fis_scoring.R) so that F2/metrics_bc.R can test:
#       (b) the FBO optimum lies ON this front, and (c) it is the MAX-U point of it.
#
#   This is the CHEAP track (no EnergyPlus). It does NOT feed the sample-efficiency
#   metric (a) — that runs NSGA-II directly on E+ (F2/nsga2_direct_ep.R), on the Ryzen.
#
#   Output: F2/out/front_<scn>.csv  with columns
#       scn, seed, tau_CO2..tau_AQ, g1..g5 (predicted), f1, f2, f3, U, pareto
#
# USAGE (RStudio):
#   source(here::here("F2","nsga2_surrogate.R"))
#   # SMOKE TEST primero (rápido, 1 escenario):
#   sm <- f2_nsga2_scenario("N_w410", seeds = 1, popsize = 40, generations = 30)
#   # corrida completa (multi-semilla, los 9 escenarios):
#   all <- f2_nsga2_all(seeds = 1:10, popsize = 100, generations = 100)
#
# Requires: mco  ->  install.packages("mco")   (además de DiceKriging, here).
#   mco::nsga2 usa exactamente los operadores canónicos NSGA-II (SBX + mutación polinómica).
#   Reproducible: set.seed(s) por semilla. La moneda del benchmark es nº de evaluaciones
#   (no wall-clock); aquí las evaluaciones son del surrogate (baratas), no de E+.
# =============================================================================

source(here::here("F2", "surrogate.R"))     # f2_fit_scenario, f2_predict, cost_transparent, SF, VARS, OBJS, ALL9, scn_box
suppressMessages(library(mco))

OUT_DIR <- here::here("F2", "out")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

F2_NSGA_VERSION <- "v3 (vectorized=TRUE, salida odim x popsize)"
message("[nsga2_surrogate] ", F2_NSGA_VERSION)   # confirma qué versión se ha cargado

# ---- objetivos del NSGA-II a partir de un fit de surrogate ------------------
#   Xraw: data.frame con columnas VARS (una o varias filas).
#   Devuelve lista: F = matriz (n x 3) de (f1,f2,f3) a MINIMIZAR, y pred = indicadores predichos.
f2_objectives <- function(fit, Xraw) {
  p  <- f2_predict(fit, Xraw, objs = OBJS)             # g1, g2_pII, g2_pIII, g3, g4, g5
  F  <- cbind(f1 = p$g1, f2 = p$g2_pIII, f3 = cost_transparent(p$g3, p$g4))
  list(F = F, pred = p)
}

# ---- NSGA-II sobre el surrogate de UN escenario, multi-semilla --------------
f2_nsga2_scenario <- function(scn, seeds = 1:10, popsize = 100, generations = 100,
                              cprob = 0.9, mprob = 1/length(VARS), write = TRUE) {
  fit <- f2_fit_scenario(scn); box <- fit$box
  # guardia: confirma que surrogate.R está ACTUALIZADO (f2_predict de 1 punto -> data.frame 1x6).
  # Si falla, el error de mco sería críptico; aquí damos el mensaje claro.
  .chk <- f2_predict(fit, as.data.frame(as.list(setNames(box$lo, VARS))))
  if (!all(c("g1", "g3", "g4") %in% names(.chk)) || nrow(.chk) != 1L)
    stop("f2_predict devuelve forma incorrecta para 1 punto -> reemplaza F2/surrogate.R por la versión corregida y re-sourcea.")
  # Función objetivo VECTORIZADA para mco (vectorized=TRUE). Contrato CONFIRMADO (código de
  # mco + mensaje de error de esta versión): la población entra como matriz popsize x idim
  # (individuos por FILA) y mco espera la salida como odim x popsize (objetivos por FILA,
  # individuos por COLUMNA) -> devolvemos t(F). Un solo predict GP por generación (rápido).
  fn <- function(X) {
    if (is.null(dim(X))) X <- matrix(X, nrow = 1L)
    Xr <- as.data.frame(X); names(Xr) <- VARS
    t(as.matrix(f2_objectives(fit, Xr)$F))          # odim x popsize  (3 x n)
  }
  out <- list()
  for (s in seeds) {
    set.seed(s)
    res <- nsga2(fn, idim = length(VARS), odim = 3L,
                 lower.bounds = box$lo, upper.bounds = box$hi,
                 popsize = popsize, generations = generations,
                 cprob = cprob, mprob = mprob, vectorized = TRUE)
    Xr <- as.data.frame(res$par); names(Xr) <- VARS
    ob <- f2_objectives(fit, Xr)
    U  <- vapply(seq_len(nrow(Xr)), function(i) scoring_utility(
            SF, g1 = ob$pred$g1[i], g2 = c(ob$pred$g2_pII[i], ob$pred$g2_pIII[i]),
            g3 = ob$pred$g3[i], g4 = ob$pred$g4[i], g5 = ob$pred$g5[i]), numeric(1))
    out[[length(out) + 1]] <- cbind(
      scn = scn, seed = s, Xr, ob$pred,
      f1 = ob$F[, "f1"], f2 = ob$F[, "f2"], f3 = ob$F[, "f3"],
      U = U, pareto = res$pareto.optimal)
  }
  D <- do.call(rbind, out)
  if (write) write.csv(D, file.path(OUT_DIR, paste0("front_", scn, ".csv")), row.names = FALSE)
  cat(sprintf("  %-8s | semillas %s | pop=%d gen=%d | pts Pareto (unión): %d | mejor U=%.4f\n",
              scn, paste(range(seeds), collapse="-"), popsize, generations,
              sum(D$pareto), max(D$U[D$pareto])))
  invisible(D)
}

# ---- todos los escenarios ---------------------------------------------------
f2_nsga2_all <- function(scns = ALL9, seeds = 1:10, popsize = 100, generations = 100, ...) {
  cat(sprintf("== NSGA-II sobre surrogate | %d escenarios | %d semillas | pop=%d gen=%d ==\n",
              length(scns), length(seeds), popsize, generations))
  res <- lapply(scns, function(s) f2_nsga2_scenario(s, seeds = seeds,
                popsize = popsize, generations = generations, ...))
  names(res) <- scns
  cat("  frentes escritos en:", OUT_DIR, "\n")
  invisible(res)
}

# ---- ejecución directa: SMOKE TEST (no corre los 9 por defecto) -------------
if (sys.nframe() == 0L) {
  cat("[nsga2_surrogate] smoke test: N_w410, 1 semilla, pop=40 gen=30\n")
  f2_nsga2_scenario("N_w410", seeds = 1, popsize = 40, generations = 30)
  cat("OK. Para la corrida completa: f2_nsga2_all(seeds=1:10, popsize=100, generations=100)\n")
}
