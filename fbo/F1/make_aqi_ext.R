# ============================================================================
# make_aqi_ext.R  --  build the outdoor worst-ICA-subindex series AQI_ext(t)
# for the FBO control gate (Schedule:File), full-year 5-min, from Torneo 2024.
#
# AQI_ext(t) = max_s  subindex( running_mean(C_out,s , ref_h_s) )   in 1..6
# using the SAME ICA logic as g5.R (running_mean, ica_2019_breaks, subindex),
# so the control gate (tau_AQ) and the g5 scoring speak the same units.
# Difference vs g5: applied to the OUTDOOR concentration directly (boundary
# condition, no indoor balance).
#
# Writes to C:/EP_work (ASCII path) so EnergyPlus can read it (non-ASCII paths
# under OneDrive break E+). One value per line, no header -> Schedule:File
# (Rows to Skip = 0). Feb 29 dropped -> 8760 h (105120 rows) for a non-leap sim.
#
# Run inside FIS (FIS.Rproj) so `here` finds core/ and data/.
# ============================================================================

library(here)
source(here::here("core", "g5.R"))            # running_mean, ica_2019_breaks, subindex, default_species_params
source(here::here("core", "load_outdoor.R"))  # load_outdoor()

make_aqi_ext <- function(out_csv  = "C:/EP_work/AQI_ext_2024_5min.csv",
                         data_dir = here::here("data"),
                         year = 2024, drop_leap = TRUE) {
  CSV <- list(PM25 = file.path(data_dir, "PM25_HH_2024_SEV_Torneo.csv"),
              NO2  = file.path(data_dir, "NO2_HH_2024_SEV_Torneo.csv"),
              O3   = file.path(data_dir, "O3_HH_2024_SEV_Torneo.csv"))
  od <- load_outdoor(CSV, year = year, five_min = TRUE, verbose = FALSE)$fivemin

  pr <- default_species_params(); br <- ica_2019_breaks(); dt_h <- 5/60
  w  <- function(ref_h) max(1L, as.integer(round(ref_h / dt_h)))

  aqi <- pmax(
    subindex(running_mean(od$Cout_PM25, w(pr$PM25$ref_h)), br$PM25),   # 24 h
    subindex(running_mean(od$Cout_NO2,  w(pr$NO2$ref_h)),  br$NO2),    # 1 h
    subindex(running_mean(od$Cout_O3,   w(pr$O3$ref_h)),   br$O3)      # 8 h
  )

  if (drop_leap) {                                   # non-leap sim (8760 h)
    keep <- format(od$timestamp, "%m-%d") != "02-29"
    aqi  <- aqi[keep]
  }

  if (!dir.exists(dirname(out_csv))) dir.create(dirname(out_csv), recursive = TRUE)
  writeLines(as.character(as.integer(aqi)), out_csv)

  cat(sprintf("[make_aqi_ext] %d rows (%.0f h) -> %s\n",
              length(aqi), length(aqi) * dt_h, out_csv))
  cat("  banda:conteo  ", paste(names(table(aqi)), as.integer(table(aqi)),
                                 sep = ":", collapse = "  "), "\n")
  invisible(aqi)
}

if (sys.nframe() == 0L) make_aqi_ext()
