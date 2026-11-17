# =============================================================================
# F2/surrogate.R  --  T3 FBO Phase 2: Gaussian-process SURROGATES + validation gate
#
#   Builds one GP surrogate per (scenario, physical indicator) over the 5 operational
#   set-points x, using the SAME engine as Phase 1 (DiceKriging km, Matern 5/2 + ARD,
#   estimated nugget). Training data = LHS seed (50 pts/scenario) + the BO-augmented
#   points (F1/bo/bo_<scn>.csv). The augmentation is REQUIRED for the energy proxy g3:
#   with the seed alone g3 is under-resolved (CV-RMSE ~21%); adding the BO points that
#   densified the trade-off region brings it well under the gate (~12% on NE_w216).
#
#   The surrogates feed F2/nsga2_surrogate.R (dense Pareto fronts, metrics b/c). They
#   are NOT used for the sample-efficiency metric (a), which runs NSGA-II directly on
#   EnergyPlus (F2/nsga2_direct_ep.R).
#
#   Modelling decisions (baked in after cloud de-risking; the number of record is the
#   one THIS script prints on your machine):
#     * per-scenario GP, x normalised to [0,1]^5 within each scenario box (as bo_driver).
#     * log1p transform on the COST pair (g3 EDM, g4 switches): non-negative, skewed,
#       hinge-like -> log linearises them (g4 16.5->14.6%, g3 seed-only 21->20%; the
#       decisive lever for g3 is the BO augmentation, not the transform).
#     * linear on g1, g2_pII, g2_pIII, g5.
#     * g5 is ~constant within a scenario (empirically flat) -> handled as a constant
#       predictor; this also justifies its exclusion from the NSGA-II objectives (#1).
#
#   Objectives of the benchmark (design decision #1, Option B):
#     f1 = g1 (CO2) ; f2 = g2_pIII (Cat III shortfall) ; f3 = Cost(g3,g4) transparent.
#   The FIS utility U is reconstructed from the predicted indicators via core/fis_scoring.R
#   (needed for metrics b/c), so U is NOT surrogated directly.
#
# USAGE (RStudio):
#   source(here::here("F2","surrogate.R"))
#   gate <- f2_gate()                       # k-fold CV table + per-objective summary
#   fit  <- f2_fit_scenario("N_w410")       # fit all surrogates for one scenario
#   f2_predict(fit, data.frame(tau_CO2=900, delta_hys=100, d_slot=10, f_open=0.28, tau_AQ=5))
#
#   Rscript F2/surrogate.R                   # runs the gate on the 9 scenarios and prints it
#
# Requires: DiceKriging (already used in Phase 1), here.  Reproducible: set.seed in f2_gate.
# =============================================================================

suppressMessages({ library(here); library(DiceKriging) })

# ---- FIS de puntuación (base R, sin dependencias) --------------------------
source(here::here("core", "fis_scoring.R"))
SF <- build_scoring_fis(default_variables())     # terciles Cost-activo 5.40/12.02 | 6.44/8.44

# ---- rutas de datos --------------------------------------------------------
SEED_CSV <- here::here("F1", "seed_9scen_results.csv")   # 50 x 9 (re-scored in place by reprice.R)
BO_DIR   <- here::here("F1", "bo")                        # bo_<scn>.csv (puntos BO por escenario)

# ---- espacio de escenarios y de decisión (IDÉNTICO a bo_driver.R) ----------
ORIENT <- c(NE = 0, N = 329, E = 59)
FMAX   <- c(w154 = 0.16, w216 = 0.22, w410 = 0.42)
ALL9   <- as.vector(t(outer(names(ORIENT), names(FMAX), paste, sep = "_")))   # NE_w154, ...
VARS   <- c("tau_CO2", "delta_hys", "d_slot", "f_open", "tau_AQ")
LO     <- c(800, 50, 5, 0.10, 1)
scn_hi <- function(scn) { w <- sub("^[^_]+_", "", scn); c(1200, 200, 30, FMAX[[w]], 6) }
scn_box <- function(scn) list(lo = LO, hi = scn_hi(scn))

# ---- objetivos físicos y su transformación ---------------------------------
OBJS <- c("g1", "g2_pII", "g2_pIII", "g3", "g4", "g5")
TFWD <- list(g1 = identity, g2_pII = identity, g2_pIII = identity,
             g3 = log1p, g4 = log1p, g5 = identity)          # forward (fit space)
TINV <- list(g1 = identity, g2_pII = identity, g2_pIII = identity,
             g3 = expm1, g4 = expm1, g5 = identity)          # back-transform (physical)

# ---- normalización x -> [0,1]^5 --------------------------------------------
.nz <- function(M, box) sweep(sweep(M, 2, box$lo, "-"), 2, box$hi - box$lo, "/")

# ---- eje f3 = Cost transparente (ancla terciles FIS; monótono; SIN inferencia
#      difusa -> no contamina el baseline clásico). xmax por defecto = dominio FIS.
#      Ajustable a max observado si se quiere más spread en el frente. -----------
.T3T <- c(5.40, 12.02); .T4T <- c(6.44, 8.44)
.norm_tercile <- function(x, t1, t2, xmax) {
  pmin(1, pmax(0, ifelse(x <= t1, (1/3) * x / t1,
                  ifelse(x <= t2, 1/3 + (1/3) * (x - t1) / (t2 - t1),
                                  2/3 + (1/3) * (x - t2) / (xmax - t2)))))
}
cost_transparent <- function(g3, g4, g3max = 300, g4max = 60)
  0.5 * .norm_tercile(g3, .T3T[1], .T3T[2], g3max) +
  0.5 * .norm_tercile(g4, .T4T[1], .T4T[2], g4max)

# ---- carga de training: siembra + BO por escenario, U re-scored Cost-activo -
f2_load_training <- function(scn) {
  keep <- c("scn", VARS, OBJS)
  s <- read.csv(SEED_CSV, stringsAsFactors = FALSE)
  stopifnot("faltan columnas en el seed" = all(keep %in% names(s)))
  d <- s[s$scn == scn, keep]
  bo_path <- file.path(BO_DIR, paste0("bo_", scn, ".csv"))
  if (file.exists(bo_path)) {
    b <- read.csv(bo_path, stringsAsFactors = FALSE)
    if (all(keep %in% names(b))) d <- rbind(d, b[, keep])
  }
  d <- d[stats::complete.cases(d[, c(VARS, OBJS)]), , drop = FALSE]
  # U re-puntuada con Cost activo (el seed guardado es placeholder; los g's son físicos)
  d$U <- vapply(seq_len(nrow(d)), function(i) scoring_utility(
    SF, g1 = d$g1[i], g2 = c(d$g2_pII[i], d$g2_pIII[i]),
    g3 = d$g3[i], g4 = d$g4[i], g5 = d$g5[i]), numeric(1))
  d
}

# ---- ajuste de un GP (o predictor constante si el objetivo es plano) --------
.fit_one <- function(Xn, y, ob) {
  if (stats::sd(y) < 1e-9) return(list(type = "const", value = mean(y)))
  yt <- TFWD[[ob]](y)
  m <- tryCatch(km(~1, design = Xn, response = yt, covtype = "matern5_2",
                   nugget.estim = TRUE, control = list(trace = FALSE)),
                error = function(e) NULL)
  if (is.null(m)) list(type = "const", value = mean(y))
  else list(type = "km", km = m, inv = TINV[[ob]])
}

# ---- fit de TODOS los surrogates de un escenario ---------------------------
f2_fit_scenario <- function(scn, objs = OBJS) {
  d   <- f2_load_training(scn); box <- scn_box(scn)
  Xn  <- as.data.frame(.nz(as.matrix(d[, VARS]), box)); names(Xn) <- VARS
  mods <- setNames(lapply(objs, function(ob) .fit_one(Xn, d[[ob]], ob)), objs)
  list(scn = scn, box = box, models = mods, n = nrow(d), data = d)
}

# ---- predicción de indicadores en x crudos (back-transform aplicado) -------
#   Robusto para n=1 y n>1: usa lista con nombres -> data.frame con columnas = objs
#   (vapply colapsaba a 6x1 cuando nrow==1, rompiendo f2_objectives dentro de nsga2).
f2_predict <- function(fit, Xraw, objs = names(fit$models)) {
  Xn <- as.data.frame(.nz(as.matrix(Xraw[, VARS, drop = FALSE]), fit$box)); names(Xn) <- VARS
  n  <- nrow(Xn)
  out <- lapply(objs, function(ob) {
    mo <- fit$models[[ob]]
    if (mo$type == "const") rep(mo$value, n)
    else mo$inv(predict(mo$km, newdata = Xn, type = "UK", checkNames = FALSE)$mean)
  })
  names(out) <- objs
  as.data.frame(out, stringsAsFactors = FALSE)
}

# ---- GATE: validación cruzada k-fold por (escenario, objetivo) -------------
f2_gate <- function(k = 5, seed = 1, objs = OBJS, scns = ALL9, threshold = 20) {
  set.seed(seed)
  res <- list()
  for (scn in scns) {
    d <- f2_load_training(scn); box <- scn_box(scn)
    Xn <- as.data.frame(.nz(as.matrix(d[, VARS]), box)); names(Xn) <- VARS
    n <- nrow(d); fold <- sample(rep(seq_len(k), length.out = n))
    for (ob in objs) {
      y <- d[[ob]]
      if (stats::sd(y) < 1e-9) {
        res[[length(res) + 1]] <- data.frame(scn = scn, obj = ob, n = n,
          mean = mean(y), CVRMSE = 0, NMBE = 0); next
      }
      pred <- numeric(n)
      for (f in seq_len(k)) {
        tr <- fold != f; te <- !tr
        mo <- .fit_one(Xn[tr, , drop = FALSE], y[tr], ob)
        pred[te] <- if (mo$type == "const") mo$value
                    else mo$inv(predict(mo$km, newdata = Xn[te, , drop = FALSE],
                                        type = "UK", checkNames = FALSE)$mean)
      }
      rmse <- sqrt(mean((pred - y)^2)); mu <- mean(y)
      res[[length(res) + 1]] <- data.frame(scn = scn, obj = ob, n = n, mean = mu,
        CVRMSE = 100 * rmse / mu, NMBE = 100 * mean(pred - y) / mu)
    }
  }
  R <- do.call(rbind, res)
  cat(sprintf("\n== GATE surrogate  |  %d-fold CV  |  seed=%d  |  umbral CV-RMSE <= %g%% ==\n",
              k, seed, threshold))
  cat(sprintf("   (training = siembra + BO por escenario; log1p en g3,g4)\n\n"))
  for (ob in objs) {
    s <- R[R$obj == ob, ]
    cat(sprintf("  %-8s | CV-RMSE medio %5.1f%% | pasa %d/%d | rango [%4.1f, %5.1f]\n",
                ob, mean(s$CVRMSE), sum(s$CVRMSE <= threshold), nrow(s),
                min(s$CVRMSE), max(s$CVRMSE)))
  }
  bad <- R[R$CVRMSE > threshold, ]
  if (nrow(bad)) {
    cat("\n  escenarios que NO pasan el umbral:\n")
    for (i in seq_len(nrow(bad)))
      cat(sprintf("    %-8s %-8s  CV-RMSE=%.1f%%  (¿falta bo_%s.csv?)\n",
                  bad$scn[i], bad$obj[i], bad$CVRMSE[i], bad$scn[i]))
  }
  invisible(R)
}

# ---- ejecución directa: corre el gate ---------------------------------------
if (sys.nframe() == 0L) {
  R <- f2_gate()
  cat("\n== detalle por (escenario, objetivo) ==\n")
  print(R[order(R$obj, R$scn), ], row.names = FALSE, digits = 4)
}
