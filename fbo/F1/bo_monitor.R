# =============================================================================
# bo_monitor.R  --  FBO Phase 1: monitor del barrido de Optimizacion Bayesiana.
#   Lee los CSV incrementales F1/bo/bo_<scn>.csv (y la siembra) SIN tocar la
#   corrida en marcha: solo lectura. Resume por escenario cuantos pasos lleva,
#   el mejor U, la mejora sobre la siembra, si cumple Cat III y si ya se ha
#   estancado (plateau) -> util para decidir cuando esta hecho.
#
# USO (en otra sesion de R, con la corrida corriendo en la primera):
#   source(here::here("F1","bo_monitor.R"))
#   monitor_bo()                 # una foto del estado actual
#   monitor_bo(watch=TRUE)       # refresco cada REFRESH s hasta Ctrl-C
#   monitor_bo("NE_w216")        # detalle de un solo escenario
# =============================================================================

library(here)

BO_DIR   <- here::here("F1", "bo")
SEED_CSV <- here::here("F1", "seed_9scen_results.csv")
VARS     <- c("tau_CO2","delta_hys","d_slot","f_open","tau_AQ")
REFRESH  <- 60L        # segundos entre refrescos en modo watch
PLATEAU_TOL <- 1e-4    # mejora de U por debajo de esto = sin progreso
PLATEAU_N   <- 5L      # nº de ultimos pasos sin mejora para marcar "plateau"

# limites de x por escenario (para detectar variables pegadas a los topes) -----
.FMAX <- c(w154=0.16, w216=0.22, w410=0.42)
.hi_for <- function(scn) {
  w <- sub(".*_", "", scn)
  c(tau_CO2=1200, delta_hys=200, d_slot=30, f_open=unname(.FMAX[w]), tau_AQ=6)
}
.LO <- c(tau_CO2=800, delta_hys=50, d_slot=5, f_open=0.10, tau_AQ=1)

# mejor U de la siembra para ese escenario (referencia de la mejora) -----------
.seed_best <- function(scn) {
  if (!file.exists(SEED_CSV)) return(NA_real_)
  s <- read.csv(SEED_CSV, stringsAsFactors=FALSE)
  s <- s[s$scn == scn & !is.na(s$U), ]
  if (!nrow(s)) NA_real_ else max(s$U)
}

# resumen de UN escenario ------------------------------------------------------
.summ_one <- function(csv) {
  scn <- sub("^bo_", "", sub("\\.csv$", "", basename(csv)))
  d <- tryCatch(read.csv(csv, stringsAsFactors=FALSE), error=function(e) NULL)
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- d[!is.na(d$U), ]
  ib <- which.max(d$U); best <- d[ib, ]
  # racha final sin mejora (running max desde el final)
  rm_ <- cummax(d$U); stall <- 0L
  for (i in rev(seq_len(nrow(d)))) { if (d$U[i] >= rm_[nrow(d)] - PLATEAU_TOL) stall <- stall + 1L else break }
  last5 <- tail(d$U, 5); spread5 <- if (length(last5) > 1) diff(range(last5)) else 0
  seed <- .seed_best(scn)
  # variables pegadas a un tope (|dist normalizada al borde| < 1%)
  hi <- .hi_for(scn); un <- (unlist(best[VARS]) - .LO) / (hi - .LO)
  pinned <- VARS[un < 0.01 | un > 0.99]
  list(scn=scn, steps=nrow(d), bestU=best$U, seed=seed,
       dU=if (is.na(seed)) NA else best$U - seed,
       pIII=best$g2_pIII, cat3=isTRUE(best$g2_pIII < 5),
       g1=best$g1, spread5=spread5, stall=stall,
       plateau=isTRUE(stall >= PLATEAU_N),
       pinned=paste(sub("_.*","",pinned), collapse=","),
       best=best)
}

# banda cualitativa de U (segun el FIS: Poor / Aceptable / Good) ---------------
.band <- function(u) if (is.na(u)) "-" else if (u < 0.33) "Poor" else if (u < 0.67) "Acept." else "Good"

# ---- API: resumen de todos (o de uno) ---------------------------------------
monitor_bo <- function(scenario=NULL, watch=FALSE) {
  repeat {
    files <- list.files(BO_DIR, pattern="^bo_.*\\.csv$", full.names=TRUE)
    if (!is.null(scenario)) files <- files[grepl(paste0("bo_", scenario, "\\.csv$"), files)]
    if (!length(files)) { cat("(aun no hay CSV de BO en", BO_DIR, ")\n") }
    else {
      S <- Filter(Negate(is.null), lapply(files, .summ_one))
      S <- S[order(sapply(S, function(x) x$scn))]
      cat(sprintf("\n== MONITOR BO == %s | %d escenario(s) con datos ==\n",
                  format(Sys.time(), "%H:%M:%S"), length(S)))
      cat(sprintf("  %-9s %5s %8s %8s %8s %7s %6s %8s  %s\n",
                  "escenario","pasos","U_seed","U_best","dU","p_III","CatIII","banda","estado"))
      cat("  ", strrep("-", 78), "\n", sep="")
      for (x in S) {
        est <- if (x$plateau) sprintf("PLATEAU (%d)", x$stall) else "en curso"
        cat(sprintf("  %-9s %5d %8s %8.4f %+8.4f %6.2f%% %6s %8s  %s\n",
                    x$scn, x$steps,
                    ifelse(is.na(x$seed), "  -   ", sprintf("%.4f", x$seed)),
                    x$bestU, ifelse(is.na(x$dU), 0, x$dU),
                    x$pIII, ifelse(x$cat3, "SI", "no"),
                    .band(x$bestU), est))
      }
      # detalle si es un solo escenario
      if (!is.null(scenario) && length(S) == 1) {
        b <- S[[1]]$best
        cat("\n  Mejor punto (", S[[1]]$scn, "):\n", sep="")
        cat(sprintf("    tau_CO2=%.0f  delta_hys=%.0f  d_slot=%.1f  f_open=%.3f  tau_AQ=%.2f\n",
                    b$tau_CO2, b$delta_hys, b$d_slot, b$f_open, b$tau_AQ))
        cat(sprintf("    U=%.4f | CO2 medio g1=%.0f ppm | p_II=%.2f%% p_III=%.2f%% | g3=%.2f g4=%.2f g5=%.2f\n",
                    b$U, b$g1, b$g2_pII, b$g2_pIII, b$g3, b$g4, b$g5))
        if (nzchar(S[[1]]$pinned)) cat("    variables en el tope:", S[[1]]$pinned, "\n")
        cat(sprintf("    ultimos 5 pasos: spread(U)=%.5f | racha sin mejora=%d\n",
                    S[[1]]$spread5, S[[1]]$stall))
      }
      done <- sum(sapply(S, function(x) x$plateau))
      cat(sprintf("\n  %d/%d en plateau (candidatos a haber convergido). CatIII: %d/%d.\n",
                  done, length(S), sum(sapply(S, function(x) isTRUE(x$cat3))), length(S)))
    }
    if (!watch) break
    cat(sprintf("  [watch] refresco en %d s (Ctrl-C para salir) ...\n", REFRESH))
    Sys.sleep(REFRESH)
  }
  invisible(TRUE)
}

cat("[bo_monitor] listo. Uso:  monitor_bo()  |  monitor_bo(watch=TRUE)  |  monitor_bo(\"NE_w216\")\n")
