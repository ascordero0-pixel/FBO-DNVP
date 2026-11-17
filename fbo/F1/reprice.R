# =============================================================================
# reprice.R  --  Re-puntúa (solo la columna U) de la siembra y de los CSV de BO
#   con los TERCILES CORREGIDOS de g3/g4, sin re-simular en EnergyPlus.
#   Los g1..g5 guardados son físicos (salidas E+) y se conservan; solo cambia U.
#   Motivo: la optimización corrió con los placeholders pre-LHS de g3/g4
#   (t=100/200, 20/40), con la dimensión Cost inactiva. Esto deja los CSV
#   consistentes con el FIS corregido para re-arrancar la BO en caliente.
#
# USO:
#   1) edita core/fis_scoring.R -> default_variables():  g3 t=5.40,12.02 | g4 t=6.44,8.44
#   2) source(here::here("F1","reprice.R"))     # crea backups *_placeholderU.csv y reescribe U
#   3) run_bo()                                  # re-optimiza con el objetivo corregido
# =============================================================================
library(here)
source(here::here("core", "fis_scoring.R"))

# terciles corregidos (siembra LHS, global sobre los 9 escenarios)
T3 <- c(5.40, 12.02)      # g3 (degC.min/day)
T4 <- c(6.44,  8.44)      # g4 (actions/session)

v      <- default_variables()
v$g3   <- tercile_partition("g3", 0, 300, T3[1], T3[2], "degC.min/day")
v$g4   <- tercile_partition("g4", 0,  60, T4[1], T4[2], "actions/session")
sf_fix <- build_scoring_fis(v)

reprice_file <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  need <- c("U","g1","g2_pII","g2_pIII","g3","g4","g5")
  if (!all(need %in% names(d))) { cat("  [skip]", basename(path), "(faltan columnas)\n"); return(invisible()) }
  bak <- sub("\\.csv$", "_placeholderU.csv", path)
  if (!file.exists(bak)) write.csv(d, bak, row.names = FALSE)   # backup una sola vez
  ok <- !is.na(d$U)
  Uold <- d$U
  d$U[ok] <- vapply(which(ok), function(i) scoring_utility(
    sf_fix, g1 = d$g1[i], g2 = c(d$g2_pII[i], d$g2_pIII[i]),
    g3 = d$g3[i], g4 = d$g4[i], g5 = d$g5[i]), numeric(1))
  write.csv(d, path, row.names = FALSE)
  cat(sprintf("  repriced %-28s | n=%4d | max|dU|=%.4f | best %0.4f -> %0.4f\n",
              basename(path), sum(ok), max(abs(d$U[ok]-Uold[ok])),
              max(Uold[ok]), max(d$U[ok])))
}

cat("== reprice: siembra ==\n")
reprice_file(here::here("F1", "seed_9scen_results.csv"))
cat("== reprice: BO por escenario ==\n")
for (f in list.files(here::here("F1","bo"), pattern="^bo_.*\\.csv$", full.names=TRUE)) reprice_file(f)
cat("\nHecho. Backups en *_placeholderU.csv. Ahora: run_bo()  (objetivo corregido, Cost activo).\n")
