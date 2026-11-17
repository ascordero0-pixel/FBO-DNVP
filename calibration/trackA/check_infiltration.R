# =====================================================================
# check_infiltration.R  —  Verificacion fisica del crack_scale calibrado
# ---------------------------------------------------------------------
# Corre UNA simulacion en el punto optimo (crack=1.5, cap=15, heat=20)
# con los mismos parches de Track A, pide a EnergyPlus la infiltracion
# del AirflowNetwork y la expresa en renovaciones/hora (ren/h).
#
# Objetivo: comprobar que crack_scale=1.5 produce una infiltracion fisica
# coherente con un edificio de 1946 ('Pobre', nominal 0.700 ren/h).
#
# USO (tras phase1_setup()):
#   source("check_infiltration.R")
#   res <- check_infiltration_ach()
# =====================================================================

check_infiltration_ach <- function(crack = 1.5, cap = 15, heat = 20,
                                   zone = "ZXCASESTUDYCLASSROOM:CLASSROOM",
                                   zone_volume_m3 = 126.7031) {
  if (!exists(".STATE") || !exists("ce", envir = .STATE))
    stop("Ejecuta phase1_setup() primero.")
  ce  <- get("ce",  envir = .STATE)
  cfg <- get("cfg", envir = .STATE)

  cat(sprintf("Construyendo IDF del punto optimo (crack=%.1f, cap=%g, heat=%g) ",
              crack, cap, heat))
  cat("con salida de infiltracion AFN...\n")

  idf <- ce$load_idf(cfg, PATHS$idf_base)
  ce$set_timestep(idf, 12L)
  ce$scale_afn_crack_coeffs(idf, crack)
  ce$set_capacitance_multiplier(idf, cap)
  ce$add_heating_gain(idf, heat)            # estado calibrado (afecta al efecto chimenea)
  ce$close_windows_for_trackA(idf)          # ventanas NoVent + EMS off (solo cracks)
  ce$normalize_wall3_to_outdoors(idf)
  ce$add_output_variables(idf)

  # Pedir la infiltracion del AirflowNetwork (ACH directa + volumen de respaldo)
  idf$newidfobject("Output:Variable", Key_Value = "*",
                   Variable_Name = "AFN Zone Infiltration Air Change Rate",
                   Reporting_Frequency = "Timestep")
  idf$newidfobject("Output:Variable", Key_Value = "*",
                   Variable_Name = "AFN Zone Infiltration Volume",
                   Reporting_Frequency = "Timestep")

  run_idf <- ce$save_idf(idf, file.path(PATHS$runs_dir, "INFIL_check.idf"))
  out_dir <- ce$run_energyplus(cfg, run_idf, "INFIL_check")
  cat("Corrida EnergyPlus completada.\n")

  # --- lectura robusta (evita coaccionar el frame pandas) ---
  rd <- function(varname) {
    s <- ce$read_eplus_variable(out_dir, paste0(zone, ":", varname))
    if (inherits(s, "data.frame")) return(s)
    data.frame(timestamp = as.character(py_to_r(s$timestamp$astype("str")$tolist())),
               value     = as.numeric(py_to_r(s$value$tolist())),
               stringsAsFactors = FALSE)
  }

  # ACH directa del AFN; si no se generara, se deriva del volumen
  ach <- tryCatch(rd("AFN Zone Infiltration Air Change Rate"),
                  error = function(e) NULL)
  if (is.null(ach)) {
    cat("[aviso] Variable ACH no disponible; derivando de 'AFN Zone Infiltration Volume'.\n")
    vol <- rd("AFN Zone Infiltration Volume")          # m3 por timestep (5 min)
    # ren/h = (m3 por timestep) x (12 timesteps/h) / volumen de zona
    ach <- data.frame(timestamp = vol$timestamp,
                      value = vol$value * 12 / zone_volume_m3,
                      stringsAsFactors = FALSE)
  }

  ts   <- as.POSIXct(ach$timestamp, tz = "UTC")
  hh   <- as.integer(format(ts, "%H"))
  noct <- (hh >= 22 | hh < 7)

  mean_all  <- mean(ach$value, na.rm = TRUE)
  mean_noct <- mean(ach$value[noct], na.rm = TRUE)
  qs        <- stats::quantile(ach$value, c(.05, .5, .95), na.rm = TRUE)

  cat(sprintf("\n=== Infiltracion AirflowNetwork (crack_scale = %.1f) ===\n", crack))
  cat(sprintf("  Media global         : %.3f ren/h\n", mean_all))
  cat(sprintf("  Media nocturna (22-07): %.3f ren/h\n", mean_noct))
  cat(sprintf("  Percentiles 5/50/95  : %.3f / %.3f / %.3f ren/h\n", qs[1], qs[2], qs[3]))
  cat(sprintf("  --- referencias ---\n"))
  cat(sprintf("  Nominal plantilla    : 0.700 ren/h  (DesignBuilder, estanqueidad 'Pobre')\n"))
  cat(sprintf("  Estimacion lineal    : %.3f ren/h  (= 0.700 x %.1f)\n", 0.700 * crack, crack))
  cat(sprintf("  Rango esperable 1946 : ~0.5 - 1.5 ren/h\n"))

  verdict <- if (mean_all >= 0.5 && mean_all <= 1.5)
    "COHERENTE con un edificio de 1946 sin rehabilitar." else
    "FUERA del rango tipico — revisar."
  cat(sprintf("\n  Veredicto: %s\n", verdict))

  invisible(list(mean_all = mean_all, mean_noct = mean_noct,
                 quantiles = qs, series = ach, out_dir = out_dir))
}
