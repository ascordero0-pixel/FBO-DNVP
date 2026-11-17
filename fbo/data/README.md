# fbo/data

| File | Content |
|---|---|
| `PM25_HH_2024_SEV_Torneo.csv`, `NO2_HH_2024_SEV_Torneo.csv`, `O3_HH_2024_SEV_Torneo.csv` | Hourly concentrations (µg m⁻³) of station 41091009 (Torneo, Seville), 2024. Red de Vigilancia y Control de la Calidad del Aire de Andalucía, Junta de Andalucía open data (CC BY 4.0). Read by `core/load_outdoor.R`, which keeps station 41091009, gap-fills and expands to 5 min (8784 h, leap year). |
| `ESP_Sevilla.083910_IWEC.epw` | **Not included.** Download from <https://energyplus.net/weather> (Spain → Sevilla, IWEC) and place it here. SHA-256 `c060cd44c9b8c5819979627992ccae82fc7623f43e7bc1a12967a8866f777b36`. It must be the only `.epw` file in this folder. |
