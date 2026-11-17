# =============================================================================
# anchor_check_epw.R -- Prueba ancla: ¿qué EPW usaron las corridas de producción?
#
# Idea: re-evaluar HOY el óptimo registrado de N_w410 (Table 7: U = 0.5012) con
# un EPW elegido explícitamente. Si U reproduce el valor logueado, ese EPW es el
# que se usó entonces. Cada corrida anual tarda ~3 min.
#
# USO (RStudio, proyecto FIS.Rproj abierto):
#   1) Guarda este fichero como F1/anchor_check_epw.R
#   2) Ejecuta PASO 1 y lee en consola qué fichero de óptimos ha encontrado y
#      sus columnas. Si la fila N_w410 no se identifica sola, rellena x_manual.
#   3) Ejecuta PASO 2 (una corrida de E+). Compara la salida con lo esperado.
#   4) Solo si NO reproduce: ejecuta PASO 3 con los otros dos EPW.
# =============================================================================

library(here)
source(here("F1", "fbo_eval.R"))        # define fbo_setup() y evaluate_U()

# ------------------------------- PASO 1 --------------------------------------
# Localizar el CSV con los óptimos (best_x_9scen_v2.csv) y la fila N_w410
cand <- list.files(here(), pattern = "^best_x_9scen.*\\.csv$", recursive = TRUE, full.names = TRUE)
cat("Ficheros de óptimos encontrados:\n"); print(cand)
bx <- read.csv(cand[1], check.names = FALSE)
cat("\nColumnas:", paste(names(bx), collapse = " | "), "\n\n"); print(bx)

# Intento automático de identificar la fila N_w410 (orientación N, apertura 4.10 %)
row_N410 <- which(apply(bx, 1, function(r) any(grepl("N_w410", r)) ||
                        (any(grepl("^N$", trimws(r))) && any(grepl("410|4\\.10|0\\.42", r)))))
if (length(row_N410) == 1) {
  r <- bx[row_N410, ]
  pick <- function(nm) { j <- grep(nm, names(r), ignore.case = TRUE)[1]; as.numeric(r[[j]]) }
  x <- list(tau_CO2 = pick("tau_?co2"), delta_hys = pick("hys"), d_slot = pick("slot"),
            f_open = pick("f_?open"), tau_AQ = pick("tau_?aq"))
} else {
  # Si no se identifica sola, copia aquí los valores EXACTOS (todos los decimales)
  # de la fila N_w410 del CSV impreso arriba:
  x <- list(tau_CO2 = 872, delta_hys = 133, d_slot = 5.0, f_open = 0.286, tau_AQ = 5.29)  # <- x_manual
  warning("Fila N_w410 no identificada automáticamente: usando x_manual (valores redondeados de Table 7)")
}
cat("\nVector x a re-evaluar:\n"); str(x)

# ------------------------------- PASO 2 --------------------------------------
# Re-evaluación con el IWEC puro (el primero alfabético en data/, candidato principal)
NORTH_N <- 329            # North Axis del escenario N (NE = 0, E = 59; cierre Fase 1)
IDF     <- here("F1", "v4_FBO_annual_prod.idf")

run_check <- function(epw_name, tag) {
  epw <- here("data", epw_name); stopifnot(file.exists(epw))
  fbo_setup(base_idf = IDF, epw = epw, work_dir = "C:/EP_work/FBO_check", V = 126.7031)
  t0 <- Sys.time()
  r  <- evaluate_U(x, north = NORTH_N, tag = tag)
  cat(sprintf("\n[%s]  %.1f min\n", epw_name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  out <- c(U = r$U, g1 = r$g1, g3 = r$g3, g4 = r$g4, g5 = r$g5)
  print(round(out, 4))
  cat("Esperado (Table 7, N_w410): U 0.5012 | g1 565 | g3 22.3 | g4 19.3 | g5 3.00\n")
  invisible(r)
}

res_iwec <- run_check("ESP_Sevilla.083910_IWEC.epw", "N_w410_check_IWEC")

# ------------------------------- PASO 3 --------------------------------------
# SOLO si el PASO 2 no reproduce U = 0.5012: probar los otros dos ficheros
# res_hyb  <- run_check("Seville_measured_NovDec2021.epw", "N_w410_check_hybrid")
# res_swec <- run_check("ESP_Sevilla.083910_SWEC.epw",     "N_w410_check_SWEC")
