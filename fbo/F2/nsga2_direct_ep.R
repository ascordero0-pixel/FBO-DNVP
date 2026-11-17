# =============================================================================
# F2/nsga2_direct_ep.R  --  T3 FBO Phase 2, step 4: NSGA-II DIRECTO en EnergyPlus
#
#   La prueba INDEPENDIENTE del benchmark. Corre NSGA-II clásico (mco) llamando a
#   evaluate_U (EnergyPlus real, F1/fbo_eval.R) por candidato, en los 2 escenarios
#   protagonistas (N_w410, NE_w216). Produce:
#     (a) EFICIENCIA MUESTRAL: nº de evaluaciones E+ hasta que el incumbente alcanza
#         U >= U*-eps (eps=0.01), comparado con el nº que necesitó el FBO (siembra+BO).
#         Moneda = EVALUACIONES E+ (no wall-clock). = el TITULAR del paper.
#     (b)/(c) sobre un frente REAL e independiente de la BO: ¿el óptimo FBO cae sobre
#         él y es el de máxima U? (aquí no hay optimismo de surrogate que valga).
#
#   Coste acotado: budget = popsize*generations (~1000). Se registra CADA evaluación,
#   y la métrica (a) se lee POST-HOC de la traza (no hace falta parar mco a mitad).
#   Paralelismo por población (cluster PSOCK, patrón de bo_driver.R): cada worker con
#   su propio work_dir de EnergyPlus. E+ monohilo -> escala con los núcleos de la Ryzen.
#
# USAGE (RStudio, dentro de FIS.Rproj):
#   source(here::here("F2","nsga2_direct_ep.R"))
#   # 1) SMOKE TEST barato (12 corridas E+, secuencial) para validar el acoplamiento:
#   sm <- run_direct_ep("N_w410", popsize = 6, generations = 2, ncores = 1)
#   # 2) corrida real (ajusta ncores al nº de núcleos útiles de la Ryzen):
#   r1 <- run_direct_ep("N_w410",  popsize = 40, generations = 25, ncores = 10)
#   r2 <- run_direct_ep("NE_w216", popsize = 40, generations = 25, ncores = 10)
#
# Requires: mco, F1/fbo_eval.R (evaluate_U, fbo_setup), F2/metrics_bc.R (helpers .bc,
#   .nd_mask, f2_fbo_optimum, f2_load_training). Reproducible: set.seed(seed).
# =============================================================================

source(here::here("F2", "metrics_bc.R"))       # trae surrogate.R + helpers (.nd_mask, .bc, f2_fbo_optimum, VARS, scn_box, ALL9, cost_transparent, SF, EPS_U)
source(here::here("F1", "fbo_eval.R"))         # evaluate_U, fbo_setup
suppressMessages({ library(mco); library(parallel) })

DE_OUT <- here::here("F2", "out")
if (!dir.exists(DE_OUT)) dir.create(DE_OUT, recursive = TRUE)
message("[nsga2_direct_ep] v1 (mco vectorized=TRUE, E+ real, paralelo por población)")

ORIENT_DE <- c(NE = 0, N = 329, E = 59)
scn_north <- function(scn) ORIENT_DE[[ sub("_.*", "", scn) ]]

# rutas E+ (idioma de bo_driver.R); resueltas en tiempo de ejecución
.resolve_idf <- function() {
  cand <- c(here::here("F1", "v4_FBO_annual_prod.idf"), "C:/EP_work/v4_FBO_annual_prod.idf")
  cand[file.exists(cand)][1]
}
.resolve_epw <- function() list.files(here::here("data"), pattern = "\\.epw$", full.names = TRUE)[1]

# columnas del vector devuelto por una evaluación
.DE_COLS <- c("f1","f2","f3","U","g1","g2_pII","g2_pIII","g3","g4","g5")
.PENALTY <- c(f1 = 1e6, f2 = 1e6, f3 = 1e6)     # objetivos de una corrida fallida (dominada por todo)

# ---- evaluación de UN candidato via EnergyPlus -----------------------------
.de_eval_one <- function(xvec, north, tag) {
  x <- setNames(as.list(as.numeric(xvec)), VARS)
  r <- tryCatch(evaluate_U(x, north = north, tag = tag), error = function(e) NULL)
  if (is.null(r) || is.na(r$U))
    return(c(f1 = NA, f2 = NA, f3 = NA, U = NA, g1 = NA, g2_pII = NA,
             g2_pIII = NA, g3 = NA, g4 = NA, g5 = NA))
  c(f1 = r$g1, f2 = r$g2[2], f3 = cost_transparent(r$g3, r$g4), U = r$U,
    g1 = r$g1, g2_pII = r$g2[1], g2_pIII = r$g2[2], g3 = r$g3, g4 = r$g4, g5 = r$g5)
}

# ---- driver principal -------------------------------------------------------
run_direct_ep <- function(scn, popsize = 40, generations = 25, ncores = 1L,
                          work_root = "C:/EP_work/FBO_DE", seed = 1,
                          base_idf = NULL, epw = NULL) {
  north <- scn_north(scn); box <- scn_box(scn)
  if (is.null(base_idf)) base_idf <- .resolve_idf()
  if (is.null(epw))      epw      <- .resolve_epw()
  stopifnot("IDF de producción no encontrado" = !is.na(base_idf) && file.exists(base_idf))

  # estado local: contador de evaluaciones + traza
  st <- new.env(); st$k <- 0L; st$rows <- list()

  # cluster (paralelo por población), cada worker con su propio work_dir E+
  cl <- NULL
  if (ncores > 1L) {
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, { suppressMessages({ library(here)
      source(here::here("F2", "surrogate.R")); source(here::here("F1", "fbo_eval.R")) }) })
    # cada worker: su propio work_dir E+ (args pasados EXPLÍCITOS, no por closure)
    parLapply(cl, seq_len(ncores),
      function(k, bi, ep, wr) fbo_setup(base_idf = bi, epw = ep,
                                        work_dir = file.path(wr, paste0("w", k))),
      base_idf, epw, work_root)
    clusterExport(cl, ".de_eval_one")            # desde .GlobalEnv (source lo dejó ahí)
  } else {
    fbo_setup(base_idf = base_idf, epw = epw, work_dir = file.path(work_root, "w1"))
  }

  # función objetivo VECTORIZADA para mco (salida odim x popsize = t(F))
  fn <- function(X) {
    if (is.null(dim(X))) X <- matrix(X, nrow = 1L)
    n <- nrow(X); idx <- st$k + seq_len(n); st$k <- st$k + n
    tags <- sprintf("de_%s_%05d", scn, idx)
    ev <- if (!is.null(cl))
      parLapply(cl, seq_len(n),
        function(i, Xm, no, tg) .de_eval_one(Xm[i, ], no, tg[i]), X, north, tags)
    else
      lapply(seq_len(n), function(i) .de_eval_one(X[i, ], north, tags[i]))
    M <- do.call(rbind, ev); colnames(M) <- .DE_COLS
    st$rows[[length(st$rows) + 1L]] <- cbind(eval = idx, setNames(as.data.frame(X), VARS),
                                             as.data.frame(M))
    Fm <- M[, c("f1", "f2", "f3"), drop = FALSE]
    bad <- !is.finite(Fm[, 1]); if (any(bad)) Fm[bad, ] <- matrix(.PENALTY, sum(bad), 3, byrow = TRUE)
    t(Fm)                                            # 3 x popsize
  }

  set.seed(seed)
  res <- nsga2(fn, idim = length(VARS), odim = 3L, lower.bounds = box$lo, upper.bounds = box$hi,
               popsize = popsize, generations = generations,
               cprob = 0.9, mprob = 1/length(VARS), vectorized = TRUE)
  if (!is.null(cl)) stopCluster(cl)

  log <- do.call(rbind, st$rows)
  write.csv(log, file.path(DE_OUT, paste0("de_log_", scn, ".csv")), row.names = FALSE)

  # ---- métricas ----
  opt  <- f2_fbo_optimum(scn)                        # óptimo FBO (real, de la BO)
  real <- log[is.finite(log$U), , drop = FALSE]
  Fr <- as.matrix(real[, c("f1", "f2", "f3")]); nd <- .nd_mask(Fr)
  Ff <- Fr[nd, , drop = FALSE]; Uf <- real$U[nd]
  write.csv(cbind(real[nd, ], pareto = TRUE),
            file.path(DE_OUT, paste0("de_front_", scn, ".csv")), row.names = FALSE)
  bc <- .bc(opt$f, opt$U, Ff, Uf)                    # (b)/(c) REALES

  # (a) eficiencia muestral: evals hasta alcanzar la CALIDAD del óptimo FBO (U_FBO - eps).
  #     Comparación justa: "¿cuántas corridas E+ necesita el MOO clásico para igualar el
  #     punto de operación que entrega el FBO?". Moneda = evaluaciones E+ (no wall-clock).
  U_target <- opt$U - EPS_U
  runU   <- cummax(ifelse(is.finite(log$U), log$U, -Inf))
  N_nsga <- which(runU >= U_target)[1]               # NA => NSGA no igualó al FBO en el budget
  d_bo   <- f2_load_training(scn)                    # siembra(50) + BO, en orden de evaluación
  N_fbo  <- which(cummax(d_bo$U) >= U_target)[1]
  maxU_nsga <- max(real$U)

  cat(sprintf("\n== DIRECTO-E+ %s ==\n", scn))
  cat(sprintf("  evaluaciones E+ totales: %d | frente real: %d pts\n", nrow(log), nrow(Ff)))
  cat(sprintf("  (b) óptimo FBO on_front: %s | dist_norm %.4f\n", bc$on_front, bc$dist_norm))
  cat(sprintf("  (c) U_FBO %.4f | maxU frente real %.4f | ΔU %.4f | c_ok %s%s\n",
              opt$U, bc$U_front_max, bc$dU, bc$dU <= EPS_U,
              if (maxU_nsga > opt$U + EPS_U) "  (¡NSGA supera al FBO: hueco real de la BO!)" else ""))
  cat(sprintf("  (a) objetivo U>=%.4f | N_FBO = %d evals | N_NSGA = %s evals | ratio ~ %s\n",
              U_target, N_fbo, ifelse(is.na(N_nsga), sprintf(">%d", nrow(log)), as.character(N_nsga)),
              ifelse(is.na(N_nsga), sprintf(">%.1fx", nrow(log)/N_fbo), sprintf("%.1fx", N_nsga/N_fbo))))
  invisible(list(scn = scn, res = res, log = log, front = cbind(Ff, U = Uf),
                 bc = bc, N_fbo = N_fbo, N_nsga = N_nsga, U_target = U_target, maxU_nsga = maxU_nsga))
}
