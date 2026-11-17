# =============================================================================
# F2/rescore_direct_ep.R  --  Re-puntúa los resultados directo-E+ ya simulados
#   con el FIS ACTUAL (corregido SC7), SIN re-simular. Recalcula (a)/(b)/(c) a
#   partir de F2/out/de_log_<scn>.csv (los g1..g5 físicos no cambian; solo la U).
#
# USO:
#   Rscript F2/rescore_direct_ep.R          # los tres escenarios + CSV de resumen
#   o, en sesión:
#     source(here::here("F2", "rescore_direct_ep.R"))
#     rescore_direct("N_w410")              # un escenario
#     rescore_direct_all()                  # los tres + CSV de resumen
#
#   N_w410 y NE_w216 se simularon ANTES de la corrección SC7 -> re-puntuado necesario.
#   E_w410 se simuló DESPUÉS (run_direct_E_w410.R ya usaba el FIS corregido) -> CONTROL:
#   debe reproducir N_FBO = 69, N_NSGA = 828, ratio = 12.0.
#
# SALIDA CANÓNICA (fuente de la Tabla 7 del manuscrito):
#   F2/out/rescore_direct_summary.csv  (precisión completa: 6 cifras significativas)
#
# Resultados de referencia (ejecución 2026-11-10):
#   N_w410  : on_front TRUE | dist 0.021521 | dU -0.005976 | N_FBO  51 | N_NSGA 716 | 14.0x
#   NE_w216 : on_front TRUE | dist 0.011480 | dU -0.005339 | N_FBO 118 | N_NSGA 525 |  4.4x
#   E_w410  : on_front TRUE | dist 0.020839 | dU -0.006486 | N_FBO  69 | N_NSGA 828 | 12.0x
# =============================================================================

source(here::here("F2", "metrics_bc.R"))   # trae surrogate.R + helpers (.nd_mask,.bc, f2_fbo_optimum,
                                           #  f2_load_training, cost_transparent, SF, VARS, EPS_U)
DE_OUT <- here::here("F2", "out")

rescore_direct <- function(scn) {
  f <- file.path(DE_OUT, paste0("de_log_", scn, ".csv"))
  stopifnot("no encuentro el de_log" = file.exists(f))
  log <- read.csv(f, stringsAsFactors = FALSE)

  # --- U re-puntuada con el FIS ACTUAL (SF se construyó al sourcear, ya corregido) ---
  log$U <- vapply(seq_len(nrow(log)), function(i) {
    if (!is.finite(log$g1[i])) return(NA_real_)
    scoring_utility(SF, g1 = log$g1[i], g2 = c(log$g2_pII[i], log$g2_pIII[i]),
                    g3 = log$g3[i], g4 = log$g4[i], g5 = log$g5[i])
  }, numeric(1))

  # --- óptimo FBO (de la BO), re-puntuado con el FIS actual ---
  opt <- f2_fbo_optimum(scn)                       # usa f2_load_training + SF corregido

  # --- frente real (no-dominado en f1,f2,f3 físicos; U re-puntuada) ---
  real <- log[is.finite(log$U), , drop = FALSE]
  Fr <- as.matrix(real[, c("f1", "f2", "f3")]); nd <- .nd_mask(Fr)
  Ff <- Fr[nd, , drop = FALSE]; Uf <- real$U[nd]

  # (b)/(c)
  bc <- .bc(opt$f, opt$U, Ff, Uf)

  # (a) evals hasta la calidad del óptimo FBO
  U_target <- opt$U - EPS_U
  runU   <- cummax(ifelse(is.finite(log$U), log$U, -Inf))
  N_nsga <- which(runU >= U_target)[1]
  d_bo   <- f2_load_training(scn)
  N_fbo  <- which(cummax(d_bo$U) >= U_target)[1]
  maxU_nsga <- max(real$U)

  cat(sprintf("\n== RE-PUNTUADO (FIS corregido) · %s ==\n", scn))
  cat(sprintf("  evals E+ = %d | frente real = %d pts | U_FBO = %.6f\n", nrow(log), nrow(Ff), opt$U))
  cat(sprintf("  (b) óptimo FBO on_front: %s | dist_norm %.6f\n", bc$on_front, bc$dist_norm))
  cat(sprintf("  (c) maxU frente real %.6f | ΔU %+.6f | c_ok %s%s\n",
              bc$U_front_max, bc$dU, bc$dU <= EPS_U,
              if (maxU_nsga > opt$U + EPS_U) "  (¡NSGA aún supera al FBO!)" else ""))
  cat(sprintf("  (a) objetivo U>=%.6f | N_FBO = %d | N_NSGA = %s | ratio ~ %s\n",
              U_target, N_fbo, ifelse(is.na(N_nsga), sprintf(">%d", nrow(log)), as.character(N_nsga)),
              ifelse(is.na(N_nsga), sprintf(">%.1fx", nrow(log)/N_fbo), sprintf("%.1fx", N_nsga/N_fbo))))
  invisible(list(scn = scn, bc = bc, U_fbo = opt$U, N_fbo = N_fbo, N_nsga = N_nsga))
}

# ---- los tres escenarios + resumen canónico (Tabla 7) ------------------------
rescore_direct_all <- function(scns = c("N_w410", "NE_w216", "E_w410"),
                               out = file.path(DE_OUT, "rescore_direct_summary.csv")) {
  res <- lapply(scns, rescore_direct)
  tab <- do.call(rbind, lapply(res, function(r) data.frame(
    scn         = r$scn,
    on_front    = r$bc$on_front,
    dist_norm   = signif(r$bc$dist_norm, 6),
    U_fbo       = signif(r$U_fbo, 6),
    U_front_max = signif(r$bc$U_front_max, 6),
    dU          = signif(r$bc$dU, 6),
    c_ok        = r$bc$dU <= EPS_U,
    N_FBO       = r$N_fbo,
    N_NSGA      = r$N_nsga,
    rho         = round(r$N_nsga / r$N_fbo, 2),
    stringsAsFactors = FALSE)))
  write.csv(tab, out, row.names = FALSE)
  cat("\nResumen escrito en:", out, "\n")
  print(tab, row.names = FALSE)
  invisible(tab)
}

# ---- ejecución directa (Rscript) ---------------------------------------------
if (sys.nframe() == 0L) rescore_direct_all()
