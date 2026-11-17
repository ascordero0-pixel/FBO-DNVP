# Licence for data, models and evaluation logs

Except where noted below, all non-code files in this repository (EnergyPlus
input files `*.idf`, measured series `*.csv` / `*.xlsx`, calibrated simulation
outputs under `runs/`, and the evaluation logs under `fbo/F1/` and `fbo/F2/out/`)
are licensed under the **Creative Commons Attribution 4.0 International
licence (CC BY 4.0)**: <https://creativecommons.org/licenses/by/4.0/>.

Attribution: Sánchez Cordero, A., Gómez Melgar, S., & Andújar Márquez, J. M.
(2026). FBO-DNVP (Version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## Third-party data

- `fbo/data/*_HH_2024_SEV_Torneo.csv`: hourly PM2.5, NO2 and O3 of station
  41091009 (Torneo, Seville), year 2024. Source: Red de Vigilancia y Control de
  la Calidad del Aire de Andalucía, Junta de Andalucía open-data portal
  (<https://www.juntadeandalucia.es/datosabiertos/portal.html>), published under
  CC BY 4.0. Values are reproduced unchanged; only the station filter applies.
- Weather: the IWEC file `ESP_Sevilla.083910_IWEC.epw` (© 2001 ASHRAE,
  distributed by the U.S. DOE through energyplus.net) is **not** included.
  See README Section 3 for download and checksum instructions.
