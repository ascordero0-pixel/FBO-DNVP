# =============================================================================
# bo_driver.R  --  FBO Phase 1: Bayesian Optimization loop over the 9 scenarios.
#   Per scenario: GP surrogate (Matern 5/2 + ARD) on the LHS seed -> maximize EI
#   -> evaluate_U (EnergyPlus) -> augment -> refit. Stop when EI < 1% of the
#   observed U range, or after N_MAX extra evaluations.
#   Scenarios are INDEPENDENT -> run in parallel (one worker per scenario).
#   Each scenario is RESUMABLE: its BO points are appended to F1/bo/bo_<scn>.csv
#   (distinct file per scenario -> no write races), and re-running continues.
#
# Requiere: DiceKriging, DiceOptim   ->  install.packages(c("DiceKriging","DiceOptim"))
#
# WORKFLOW:
#   1) que exista F1/seed_9scen_results.csv (del seed_driver)
#   2) source(here::here("F1","bo_driver.R"))
#   3) run_bo(scenarios = "NE_w216")   # validar en 1 escenario primero (recomendado)
#   4) run_bo()                        # los 9 escenarios
# =============================================================================

library(here); library(parallel)

# ---- CONFIG -----------------------------------------------------------------
NCORES  <- 7L
N_MAX   <- 60L        # pasos BO maximos por escenario
EI_FRAC <- 0.01       # parada 1: EI < EI_FRAC * rango(U observado)
DUP_TOL <- 5e-3       # parada 2: EI apunta a un punto ya muestreado (dist. en [0,1]^5)
IMP_TOL <- 1e-4       # mejora minima de U que cuenta como "progreso"
STALL_MAX <- 5L       # parada 3: nº de pasos seguidos sin mejora > IMP_TOL
PROJECT <- here::here()
P_IND   <- here::here("core","indicators.R"); P_EVAL <- here::here("F1","fbo_eval.R")
P_FIS   <- here::here("core","fis_scoring.R")
P_IDF   <- c(here::here("F1","v4_FBO_annual_prod.idf"),"C:/EP_work/v4_FBO_annual_prod.idf")
P_IDF   <- P_IDF[file.exists(P_IDF)][1]
P_EPW   <- list.files(here::here("data"), pattern="\\.epw$", full.names=TRUE)[1]
WORK_DIR<- "C:/EP_work/FBO"; V_ZONE <- 126.7031
SEED_CSV<- here::here("F1","seed_9scen_results.csv")
BO_DIR  <- here::here("F1","bo"); if(!dir.exists(BO_DIR)) dir.create(BO_DIR, recursive=TRUE)
stopifnot("IDF prod no encontrado"=!is.na(P_IDF))   # [repo patch P5]
if (is.na(P_EPW)) message("[bo_driver] no EPW in data/: needed only by run_bo() (see data/README.md)")

# ---- escenarios + espacio x -------------------------------------------------
ORIENT <- c(NE=0, N=329, E=59); FMAX <- c(w154=0.16, w216=0.22, w410=0.42)
SCN <- do.call(rbind, lapply(names(ORIENT), function(o) do.call(rbind, lapply(names(FMAX), function(w)
  data.frame(name=sprintf("%s_%s",o,w), north=ORIENT[[o]], f_max=FMAX[[w]], stringsAsFactors=FALSE)))))
VARS <- c("tau_CO2","delta_hys","d_slot","f_open","tau_AQ")
LO   <- c(800, 50, 5, 0.10, 1)                       # limites inferiores (f_open sup = f_max(c))
COLS <- c("scn","tag","north","tau_CO2","delta_hys","d_slot","f_open","tau_AQ",
          "U","g1","g2_pII","g2_pIII","g3","g4","g5","phase")

# ---- init de un worker ------------------------------------------------------
.init_bo <- function() {
  setwd(PROJECT); suppressMessages({ library(here); library(DiceKriging); library(DiceOptim) })
  stopifnot("EPW no encontrado en data/" = !is.na(P_EPW))   # [repo patch P5]
  source(P_IND); source(P_EVAL); source(P_FIS)
  fbo_setup(base_idf=P_IDF, epw=P_EPW, work_dir=WORK_DIR, V=V_ZONE)
  FBO$sf <- build_scoring_fis(default_variables())
  invisible(TRUE)
}

# ---- fila de resultado con esquema comun ------------------------------------
.mk_row <- function(scn, tag, xr, r, phase) data.frame(
  scn=scn, tag=tag, north=NA, tau_CO2=xr[1], delta_hys=xr[2], d_slot=xr[3],
  f_open=xr[4], tau_AQ=xr[5], U=r$U, g1=r$g1, g2_pII=r$g2[1], g2_pIII=r$g2[2],
  g3=r$g3, g4=r$g4, g5=r$g5, phase=phase, stringsAsFactors=FALSE)

# ---- BO de UN escenario -----------------------------------------------------
bo_scenario <- function(scn, seed_df, n_max=N_MAX, ei_frac=EI_FRAC) {
  hi  <- c(1200, 200, 30, scn$f_max, 6)
  nz  <- function(df) { M<-as.matrix(df[,VARS,drop=FALSE]); sweep(sweep(M,2,LO,"-"),2,hi-LO,"/") }
  dz  <- function(u)  LO + as.numeric(u)*(hi-LO)
  csv <- file.path(BO_DIR, paste0("bo_", scn$name, ".csv"))

  base <- seed_df[!is.na(seed_df$U), , drop=FALSE]         # siembra (50 puntos)
  prev <- if (file.exists(csv)) read.csv(csv, stringsAsFactors=FALSE) else NULL
  hist <- base[, intersect(names(base), VARS)]             # placeholder; se re-arma abajo
  # tabla de trabajo: siembra + BO previa
  work <- rbind(base[, c(VARS,"U")], if(!is.null(prev)) prev[, c(VARS,"U")] else base[0,c(VARS,"U")])
  step <- if (is.null(prev)) 0L else nrow(prev)

  best   <- max(work$U); stall <- 0L; reason <- "n_max"    # guardianes de parada
  while (step < n_max) {
    Xn <- nz(work); y <- work$U; Ur <- diff(range(y)); if (Ur<=0) Ur <- 1
    m <- tryCatch(km(~1, design=as.data.frame(setNames(as.data.frame(Xn), VARS)),
                     response=y, covtype="matern5_2", nugget.estim=TRUE,
                     control=list(trace=FALSE)), error=function(e) NULL)
    if (is.null(m)) { reason <- "km_error"; break }
    ei <- tryCatch(max_EI(m, lower=setNames(rep(0,5),VARS), upper=setNames(rep(1,5),VARS),
                          minimization=FALSE, control=list(print.level=0)),
                   error=function(e) NULL)
    if (is.null(ei)) { reason <- "ei_error"; break }
    # --- parada 1: EI por debajo del 1% del rango de U observado -----------------
    if (as.numeric(ei$value) < ei_frac*Ur) { reason <- "ei_low"; break }
    # --- parada 2: EI apunta a un punto ya muestreado (vertice) -> convergido ----
    #     dmin = distancia euclidea minima en [0,1]^5 al diseno actual.
    un   <- as.numeric(ei$par)
    dmin <- min(sqrt(rowSums(sweep(nz(work), 2, un, "-")^2)))
    if (dmin < DUP_TOL) { reason <- "duplicado"; break }
    step <- step + 1L
    xr <- dz(ei$par); names(xr) <- VARS
    r  <- tryCatch(evaluate_U(as.list(xr), north=scn$north,
                              tag=sprintf("bo_%s_%03d", scn$name, step)), error=function(e) NULL)
    if (is.null(r)) next
    row <- .mk_row(scn$name, sprintf("bo_%s_%03d",scn$name,step), xr, r, "bo")
    write.table(row, csv, sep=",", row.names=FALSE,
                col.names=!file.exists(csv), append=file.exists(csv))   # guardado incremental
    work <- rbind(work, setNames(as.data.frame(as.list(c(xr, U=r$U))), c(VARS,"U")))
    # --- parada 3: sin mejora significativa en STALL_MAX pasos seguidos ----------
    if (r$U > best + IMP_TOL) { best <- r$U; stall <- 0L } else stall <- stall + 1L
    if (stall >= STALL_MAX) { reason <- "sin_mejora"; break }
  }
  list(scn=scn$name, n_bo=step, best_U=max(work$U), stop=reason,
       best_x=work[which.max(work$U), VARS])
}

# ---- funcion de tarea a NIVEL GLOBAL (evita serializar el cluster) ----------
.bo_task <- function(task) bo_scenario(task$scn, task$seed, task$n_max, task$ei_frac)

# ---- orquestador: paralelo entre escenarios ---------------------------------
run_bo <- function(scenarios=NULL, n_max=N_MAX, ei_frac=EI_FRAC) {
  stopifnot("falta el seed CSV"=file.exists(SEED_CSV))
  seed <- read.csv(SEED_CSV, stringsAsFactors=FALSE); seed <- seed[!is.na(seed$U), ]
  seed <- seed[!duplicated(seed$tag, fromLast=TRUE), ]        # por si hubo reintentos en la siembra
  sel  <- if (is.null(scenarios)) SCN$name else intersect(scenarios, SCN$name)
  # tareas autocontenidas: cada una lleva su escenario + su subconjunto de siembra
  tasks <- lapply(sel, function(nm) list(scn=SCN[SCN$name==nm, ], seed=seed[seed$scn==nm, ],
                                         n_max=n_max, ei_frac=ei_frac))
  K <- min(NCORES, length(tasks))
  cat(sprintf("== BO sobre %d escenario(s) | n_max=%d | EI_stop=%.0f%% | %d workers ==\n",
              length(tasks), n_max, 100*ei_frac, K))

  if (K > 1L) {
    cl <- makeCluster(K); on.exit(stopCluster(cl), add=TRUE)
    clusterExport(cl, c("PROJECT","P_IND","P_EVAL","P_FIS","P_IDF","P_EPW","WORK_DIR","V_ZONE",
                        "VARS","LO","BO_DIR","DUP_TOL","IMP_TOL","STALL_MAX",
                        ".init_bo",".mk_row","bo_scenario",".bo_task"),
                  envir=.GlobalEnv)
    invisible(clusterCall(cl, .init_bo))
    res <- parLapply(cl, tasks, .bo_task)
  } else {                          # 1 escenario o serie -> en el proceso principal
    .init_bo(); res <- lapply(tasks, .bo_task)
  }

  cat("\n== BO COMPLETO ==\n")
  for (r in res) cat(sprintf("  %-8s | %d pasos BO | mejor U = %.4f | parada: %s\n",
                             r$scn, r$n_bo, r$best_U, r$stop))
  cat("  puntos BO por escenario en:", BO_DIR, "\n")
  invisible(res)
}

cat("[bo_driver] listo. Requiere DiceKriging+DiceOptim. Uso:\n")
cat("  run_bo(scenarios=\"NE_w216\")   # validar 1 escenario primero\n")
cat("  run_bo()                        # los 9\n")
