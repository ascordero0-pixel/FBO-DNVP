# =============================================================================
# seed_driver.R  --  FBO Phase 1: LHS SEEDING over the 9 scenarios.
#   Parallel (multi-core E+), RESUMABLE (skips done tags), incremental save.
#   Shared maximin LHS in [0,1]^5 mapped to each scenario (paired transferability
#   sweep). f_open scaled to [0.10, f_max(c)] per scenario.
#
# WORKFLOW:
#   1) source(here::here("F1","seed_driver.R"))
#   2) check_north_sign()        # 1 corrida rotada -> valida el signo del North Axis
#   3) run_seed()                # el barrido completo (resumible: re-llamar si se corta)
#
# Scenarios (Table 3): 3 orientations (North Axis) x 3 WWeqR_max (f_max).
#   orient : NE=0, N=329, E=59   |   f_max : 1.54%->0.16, 2.16%->0.22, 4.10%->0.42
# =============================================================================

library(here); library(parallel)

# ---- CONFIG (ajusta NCORES a tu maquina) ------------------------------------
NCORES   <- 7L   # workers E+ en paralelo; 1 = serie
N_LHS    <- 50L                                      # puntos por escenario
PROJECT  <- here::here()
P_IND    <- here::here("core", "indicators.R")
P_EVAL   <- here::here("F1",   "fbo_eval.R")
P_FIS    <- here::here("core", "fis_scoring.R")
P_MKAQI  <- here::here("F1",   "make_aqi_ext.R")
P_IDF    <- c(here::here("F1","v4_FBO_annual_prod.idf"), "C:/EP_work/v4_FBO_annual_prod.idf")
P_IDF    <- P_IDF[file.exists(P_IDF)][1]
P_EPW    <- list.files(here::here("data"), pattern="\\.epw$", full.names=TRUE)[1]
AQI_CSV  <- "C:/EP_work/AQI_ext_2024_5min.csv"
MASTER   <- here::here("F1", "seed_9scen_results.csv")
WORK_DIR <- "C:/EP_work/FBO"
V_ZONE   <- 126.7031
stopifnot("IDF prod no encontrado"=!is.na(P_IDF), "EPW no encontrado"=!is.na(P_EPW))

# ---- x-space (handoff sec.7) ------------------------------------------------
XB <- list(tau_CO2=c(800,1200), delta_hys=c(50,200), d_slot=c(5,30),
           f_open=c(0.10, NA), tau_AQ=c(1,6))

# ---- 9 escenarios -----------------------------------------------------------
ORIENT <- c(NE=0, N=329, E=59)
FMAX   <- c(w154=0.16, w216=0.22, w410=0.42)     # WWeqR 1.54/2.16/4.10 %
SCN <- do.call(rbind, lapply(names(ORIENT), function(o) do.call(rbind, lapply(names(FMAX), function(w)
  data.frame(name=sprintf("%s_%s", o, w), north=ORIENT[[o]], f_max=FMAX[[w]],
             orient=o, wweqr=w, stringsAsFactors=FALSE)))))
rownames(SCN) <- NULL

# ---- maximin LHS en [0,1]^5 (dependency-free; mejor min-distancia de N intentos) ----
maximin_lhs <- function(n, k, tries=300L, seed=20260728L) {
  set.seed(seed); best <- NULL; bd <- -Inf
  for (t in seq_len(tries)) {
    D <- sapply(seq_len(k), function(j) (sample.int(n) - runif(n)) / n)
    d <- min(dist(D)); if (d > bd) { bd <- d; best <- D }
  }
  best
}
D01 <- maximin_lhs(N_LHS, 5)
sc  <- function(u, lo, hi) lo + u*(hi-lo)

# ---- construir las 450 tareas (x, north, tag) -------------------------------
build_tasks <- function() {
  ts <- list()
  for (s in seq_len(nrow(SCN))) {
    fmx <- SCN$f_max[s]
    for (i in seq_len(N_LHS)) {
      x <- list(tau_CO2   = sc(D01[i,1], XB$tau_CO2[1],   XB$tau_CO2[2]),
                delta_hys = sc(D01[i,2], XB$delta_hys[1], XB$delta_hys[2]),
                d_slot    = sc(D01[i,3], XB$d_slot[1],    XB$d_slot[2]),
                f_open    = sc(D01[i,4], XB$f_open[1],    fmx),
                tau_AQ    = sc(D01[i,5], XB$tau_AQ[1],    XB$tau_AQ[2]))
      ts[[length(ts)+1]] <- list(tag=sprintf("seed_%s_%02d", SCN$name[s], i),
                                 north=SCN$north[s], scn=SCN$name[s], pt=i, x=x)
    }
  }
  ts
}

# ---- inicializar un worker (o el proceso principal en serie) ----------------
.init_env <- function() {
  setwd(PROJECT)                 # para que here::here() resuelva en el worker
  suppressMessages({ library(here) })
  source(P_IND); source(P_EVAL); source(P_FIS)
  fbo_setup(base_idf=P_IDF, epw=P_EPW, work_dir=WORK_DIR, V=V_ZONE)
  FBO$sf <- build_scoring_fis(default_variables())     # FIS re-anclado [5,10]
  invisible(TRUE)
}

# ---- evaluar una tarea -> fila de resultados --------------------------------
.eval_one <- function(task) {
  r <- tryCatch(evaluate_U(task$x, north=task$north, tag=task$tag),
                error=function(e) { message("[FAIL] ", task$tag, ": ", conditionMessage(e)); NULL })
  data.frame(scn=task$scn, pt=task$pt, tag=task$tag, north=task$north,
             tau_CO2=task$x$tau_CO2, delta_hys=task$x$delta_hys, d_slot=task$x$d_slot,
             f_open=task$x$f_open, tau_AQ=task$x$tau_AQ,
             U =if(is.null(r)) NA else r$U,  g1=if(is.null(r)) NA else r$g1,
             g2_pII =if(is.null(r)) NA else r$g2[1], g2_pIII=if(is.null(r)) NA else r$g2[2],
             g3=if(is.null(r)) NA else r$g3, g4=if(is.null(r)) NA else r$g4,
             g5=if(is.null(r)) NA else r$g5, stringsAsFactors=FALSE)
}

# ---- VERIFICACION DE SIGNO del North Axis (correr 1 vez antes del barrido) ---
check_north_sign <- function(north=329, target=0) {
  cat("[sign-check] corriendo 1 evaluacion a North Axis =", north, "...\n")
  .init_env()
  x <- list(tau_CO2=1000, delta_hys=100, d_slot=15, f_open=0.22, tau_AQ=4)
  evaluate_U(x, north=north, tag="signcheck")
  eio <- file.path(WORK_DIR, "signcheck", "eplusout.eio")
  if (!file.exists(eio)) { cat("  no encontre .eio; comprueba el azimut de Wall_3 a mano.\n"); return(invisible()) }
  L  <- readLines(eio, warn=FALSE)
  hd  <- grep("<HeatTransfer Surface>", L, value=TRUE)[1]              # cabecera
  wln <- grep("HeatTransfer Surface,.*WALL_3", L, value=TRUE, ignore.case=TRUE)[1]  # E+ pasa a MAYUS
  az <- NA
  if (!is.na(hd) && !is.na(wln)) {
    fields <- strsplit(sub(".*<HeatTransfer Surface>,?", "", hd), ",")[[1]]
    ai <- grep("Azimuth", fields, ignore.case=TRUE)[1]
    vals <- strsplit(wln, ",")[[1]]
    if (!is.na(ai) && length(vals) >= ai+1) az <- suppressWarnings(as.numeric(vals[ai+1]))
  }
  cat(sprintf("  Wall_3 (barlovento) azimut = %s  (esperado ~%d para North Axis %d)\n",
              ifelse(is.na(az),"?",round(az,1)), target, north))
  if (!is.na(az)) {
    if (abs(((az-target+180)%%360)-180) < 15) cat("  [OK] signo correcto -> usar ORIENT {0,329,59}\n")
    else cat("  [!!] signo INVERTIDO -> cambia ORIENT a c(NE=0, N=31, E=301) y re-source\n")
  }
  cat("  (linea .eio de Wall_3:)\n  ", if(length(wln)) wln else "no encontrada", "\n")
  invisible()
}

# ---- BARRIDO completo: paralelo/serie, resumible, guardado incremental -------
run_seed <- function() {
  # AQI_ext (autogenerar si falta)
  if (!file.exists(AQI_CSV)) { cat("[AQI_ext] generando ...\n"); source(P_MKAQI); make_aqi_ext(out_csv=AQI_CSV) }
  stopifnot("falta AQI_ext"=file.exists(AQI_CSV))

  tasks <- build_tasks()
  done <- character(0)
  if (file.exists(MASTER)) { prev <- read.csv(MASTER, stringsAsFactors=FALSE)
                             done <- prev$tag[!is.na(prev$U)] }
  todo <- Filter(function(t) !(t$tag %in% done), tasks)
  cat(sprintf("== seed 9 escenarios: %d tareas | %d hechas | %d por hacer | %d workers ==\n",
              length(tasks), length(done), length(todo), NCORES))
  if (!length(todo)) { cat("Todo hecho.\n"); return(invisible(read.csv(MASTER))) }

  append_rows <- function(df) {
    write.table(df, MASTER, sep=",", row.names=FALSE,
                col.names=!file.exists(MASTER), append=file.exists(MASTER))
  }

  t0 <- Sys.time(); n_done0 <- length(done)
  if (NCORES > 1L) {
    cl <- makeCluster(NCORES)
    on.exit(stopCluster(cl), add=TRUE)
    clusterExport(cl, c("PROJECT","P_IND","P_EVAL","P_FIS","P_IDF","P_EPW","WORK_DIR",
                        "V_ZONE",".init_env",".eval_one"), envir=.GlobalEnv)
    invisible(clusterCall(cl, .init_env))
    CHUNK <- 3L*NCORES
    idx <- split(seq_along(todo), ceiling(seq_along(todo)/CHUNK))
    for (k in seq_along(idx)) {
      res <- parLapplyLB(cl, todo[idx[[k]]], .eval_one)
      append_rows(do.call(rbind, res))
      nd <- sum(sapply(idx[seq_len(k)], length)); el <- as.numeric(difftime(Sys.time(),t0,units="mins"))
      cat(sprintf("  chunk %d/%d | %d/%d corridas | %.1f min | ETA %.1f min\n",
                  k, length(idx), nd, length(todo), el, el/nd*(length(todo)-nd)))
    }
  } else {                                  # serie
    .init_env()
    for (i in seq_along(todo)) {
      append_rows(.eval_one(todo[[i]]))
      el <- as.numeric(difftime(Sys.time(),t0,units="mins"))
      cat(sprintf("  [%d/%d] %s | %.1f min | ETA %.1f min\n",
                  i, length(todo), todo[[i]]$tag, el, el/i*(length(todo)-i)))
    }
  }

  DF <- read.csv(MASTER, stringsAsFactors=FALSE)
  ok <- !is.na(DF$U)
  cat(sprintf("\n== SEED COMPLETO ==\n  %d/%d OK | %.1f h total\n", sum(ok), nrow(DF),
              as.numeric(difftime(Sys.time(),t0,units="hours"))))
  cat("  U por escenario (min-max):\n")
  for (s in SCN$name) { d<-DF[DF$scn==s & ok,]; if(nrow(d)) cat(sprintf(
    "    %-8s U[%.3f, %.3f] p_III[%.1f, %.1f] g1[%.0f, %.0f]\n",
    s, min(d$U), max(d$U), min(d$g2_pIII), max(d$g2_pIII), min(d$g1), max(d$g1))) }
  cat("  master:", MASTER, "\n")
  invisible(DF)
}

cat("[seed_driver] listo. Ejecuta:  check_north_sign()  y luego  run_seed()\n")
cat(sprintf("  NCORES=%d | %d escenarios x %d LHS = %d corridas\n", NCORES, nrow(SCN), N_LHS, nrow(SCN)*N_LHS))
