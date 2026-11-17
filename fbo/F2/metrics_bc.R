# =============================================================================
# F2/metrics_bc.R  --  T3 FBO Phase 2, step 3: métricas (b) y (c) del benchmark
#
#   Sobre los frentes de F2/nsga2_surrogate.R, evalúa las dos afirmaciones que
#   sostienen la tesis (decisión #2):
#     (b) el óptimo FBO CAE SOBRE el frente de Pareto  (no dominado + distancia).
#     (c) el óptimo FBO es el punto de MÁXIMA U del frente  (ΔU = maxU_front - U_FBO).
#   La métrica (a) (eficiencia muestral) va en F2/nsga2_direct_ep.R (E+, en la Ryzen).
#
#   El ÓPTIMO FBO se deriva de los datos: es el punto de MÁXIMA U entre todos los
#   evaluados (siembra + BO) de ese escenario -> idéntico a best_x por construcción,
#   pero sin depender de un fichero externo y siempre consistente con el surrogate.
#
#   Dos fuentes de frente:
#     * "surrogate": frente denso de mco sobre el GP. La U del óptimo FBO se recalcula
#       EN EL ESPACIO DEL SURROGATE (predecir indicadores en x_FBO y puntuar) para que
#       la comparación de (c) sea consistente (surrogate vs surrogate; evita el
#       "optimismo del surrogate" de comparar U-predicha del frente con U-real del FBO).
#     * "empirical": conjunto no-dominado de los puntos REALMENTE evaluados (siembra+BO).
#       Ancla de ROBUSTEZ, imprescindible en E_w410 (surrogate g3 ~22%): ahí (b)/(c) no
#       dependen del surrogate sino de física E+ real.
#
#   Objetivos (minimización): f1=g1(CO2), f2=g2_pIII(Cat III), f3=Cost(g3,g4) transparente.
#   Convención ΔU: ΔU <= EPS_U  => el óptimo FBO iguala o supera el mejor del frente => (c) OK.
#
# USAGE:
#   source(here::here("F2","metrics_bc.R"))
#   f2_metrics_all(source = "surrogate")
#   f2_metrics_all(c("N_w410","NE_w216","E_w410"), source = "empirical")
#   Rscript F2/metrics_bc.R
#
# Requires: F2/surrogate.R, y F2/out/front_<scn>.csv (de nsga2_surrogate). NO necesita best_x.
# =============================================================================

source(here::here("F2", "surrogate.R"))

FRONT_DIR <- here::here("F2", "out")
OBJC      <- c("f1", "f2", "f3")
EPS_U     <- 0.01                              # tolerancia de (c), coherente con métrica (a)

# ---- Pareto helpers (minimización) -----------------------------------------
.dominates <- function(a, b) all(a <= b) && any(a < b)          # a domina a b
.nd_mask <- function(M) {                                        # M: n x k -> lógico no-dominado
  n <- nrow(M); keep <- rep(TRUE, n)
  for (i in seq_len(n)) if (keep[i]) for (j in seq_len(n))
    if (i != j && keep[j] && .dominates(M[j, ], M[i, ])) { keep[i] <- FALSE; break }
  keep
}
.dominated_by_set <- function(p, M) any(apply(M, 1, function(r) .dominates(r, p)))

# objetivos (f1,f2,f3) desde un data.frame con g1, g2_pIII, g3, g4
.objs_from_g <- function(d) cbind(f1 = d$g1, f2 = d$g2_pIII, f3 = cost_transparent(d$g3, d$g4))

# óptimo FBO = punto de MÁXIMA U entre los evaluados (siembra + BO) del escenario
f2_fbo_optimum <- function(scn, d = NULL) {
  if (is.null(d)) d <- f2_load_training(scn)
  r <- d[which.max(d$U), , drop = FALSE]
  list(x = setNames(as.numeric(r[, VARS]), VARS),
       f = c(f1 = r$g1, f2 = r$g2_pIII, f3 = cost_transparent(r$g3, r$g4)),
       U = r$U)
}

# núcleo (b)+(c) contra un frente dado (matriz Ffront + U por punto)
.bc <- function(f_fbo, U_fbo_cmp, Ffront, Ufront) {
  on_front <- !.dominated_by_set(f_fbo, Ffront)
  rng <- apply(Ffront, 2, function(col) diff(range(col))); rng[rng <= 0] <- 1
  d <- sqrt(rowSums((sweep(Ffront, 2, f_fbo, "-") / matrix(rng, nrow(Ffront), 3, byrow = TRUE))^2))
  list(on_front = on_front, dist_norm = min(d),
       U_front_max = max(Ufront), dU = max(Ufront) - U_fbo_cmp)
}

# ---- métricas de un escenario ----------------------------------------------
f2_metrics_scenario <- function(scn, source = c("surrogate", "empirical")) {
  source <- match.arg(source)
  d   <- f2_load_training(scn)
  opt <- f2_fbo_optimum(scn, d)
  if (source == "surrogate") {
    fr <- read.csv(file.path(FRONT_DIR, paste0("front_", scn, ".csv")), stringsAsFactors = FALSE)
    po <- fr[fr$pareto == TRUE, , drop = FALSE]
    # frente conjunto (no-dominado) entre las semillas, para no inflar con duplicados
    ndp <- .nd_mask(as.matrix(po[, OBJC]))
    po  <- po[ndp, , drop = FALSE]
    Ff  <- as.matrix(po[, OBJC]); Uf <- po$U
    # CONSISTENCIA: el óptimo FBO se coloca en el MISMO espacio que el frente (predicho por
    # el surrogate), no en el real. Comparar f real vs frente predicho sesga (b) por el
    # optimismo del GP. Aquí: objetivos y U del óptimo FBO, ambos vía surrogate.
    fit <- f2_fit_scenario(scn)
    pf  <- f2_predict(fit, as.data.frame(as.list(opt$x)))
    f_cmp <- c(f1 = pf$g1, f2 = pf$g2_pIII, f3 = cost_transparent(pf$g3, pf$g4))
    U_cmp <- scoring_utility(SF, g1 = pf$g1, g2 = c(pf$g2_pII, pf$g2_pIII),
                             g3 = pf$g3, g4 = pf$g4, g5 = pf$g5)
  } else {
    Fall <- .objs_from_g(d); nd <- .nd_mask(Fall)
    Ff <- Fall[nd, , drop = FALSE]; Uf <- d$U[nd]
    f_cmp <- opt$f                                               # objetivos reales
    U_cmp <- opt$U                                               # U real del óptimo FBO
  }
  m <- .bc(f_cmp, U_cmp, Ff, Uf)
  data.frame(scn = scn, source = source, n_front = nrow(Ff),
             on_front = m$on_front, dist_norm = round(m$dist_norm, 4),
             U_fbo_true = round(opt$U, 4), U_fbo_cmp = round(U_cmp, 4),
             U_front_max = round(m$U_front_max, 4),
             dU = round(m$dU, 4), c_ok = (m$dU <= EPS_U), stringsAsFactors = FALSE)
}

f2_metrics_all <- function(scns = ALL9, source = "surrogate")
  do.call(rbind, lapply(scns, function(s) f2_metrics_scenario(s, source)))

# ---- ejecución directa ------------------------------------------------------
if (sys.nframe() == 0L) {
  cat("== (b)/(c) sobre el frente SURROGATE (9 escenarios) ==\n")
  print(f2_metrics_all(source = "surrogate"), row.names = FALSE)
  cat("\n== cruce EMPÍRICO (robustez; imprescindible en E_w410) ==\n")
  print(f2_metrics_all(c("N_w410", "NE_w216", "E_w410"), source = "empirical"), row.names = FALSE)
  cat("\n  on_front=TRUE => (b) OK.  c_ok=TRUE (ΔU<=EPS_U) => (c) OK.\n")
}
