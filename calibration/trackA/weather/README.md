# calibration/trackA/weather

`build_measured_epw.py` builds `Seville_measured_NovDec2021.epw`, the weather
file of both calibration tracks. It splices the on-site hourly record of
1 November–21 December 2021 (from `../T3_Calibration_TRACK_A_50_nights.xlsx`)
into the IWEC Seville template. Steps: dropout cleaning, 5-min to hourly
aggregation, Magnus dew point, MSL-to-station pressure (z = 18 m), Erbs
decomposition of measured GHI into DNI + DHI (pvlib), vector-averaged wind
direction. `epw_build_report.json` records every transformation.

```bash
python build_measured_epw.py ../T3_Calibration_TRACK_A_50_nights.xlsx ESP_Sevilla_083910_IWEC.epw Seville_measured_NovDec2021.epw
```

Expected SHA-256 of the output:
`111fd24a618bf7fce3b2682f43c6e03e5753900b6a65304bfd47b8557ad9b8c7`
(input IWEC: `c060cd44c9b8c5819979627992ccae82fc7623f43e7bc1a12967a8866f777b36`).
