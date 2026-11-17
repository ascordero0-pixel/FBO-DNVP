# FBO-DNVP — Fuzzy–Bayesian Optimization of Dynamic Natural Ventilation Protocols

Code, input data and evaluation logs for:

> Sánchez Cordero, A., Gómez Melgar, S., & Andújar Márquez, J. M. (2026). *Balancing indoor air quality, thermal comfort and energy in naturally ventilated classrooms: A Fuzzy–Bayesian Optimization approach*. Submitted to *Building and Environment*.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

The framework tunes the five setpoints of a dynamic natural ventilation protocol (DNVP) for a case-study classroom in Seville (Spain). A calibrated EnergyPlus model evaluates each candidate over a school year. A hierarchical Mamdani fuzzy inference system (FIS) aggregates five annual indicators (g₁–g₅) into a utility *U* ∈ [0, 1]. A Gaussian-process Bayesian optimizer maximizes *U*. NSGA-II and two fixed-threshold DNVP serve as benchmarks.

---

## 1. Repository layout

```
FBO-DNVP/
├── verify_without_energyplus.R   one-command check of every stage that needs no EnergyPlus
├── fbo/                          optimization framework (RStudio project FIS.Rproj)
│   ├── core/                     FIS, indicators g1–g5, control law, outdoor-data loader
│   ├── data/                     Torneo 2024 outdoor series (PM2.5, NO2, O3); EPW goes here
│   ├── F0/                       Phase 0 orchestrator: unit tests + Fig. 3 + Fig. B2
│   ├── F1/                       Phase 1: production IDF, LHS seed, Bayesian optimization
│   │   └── bo/                   BO evaluation logs (9 scenarios) + optima (Table 4)
│   └── F2/                       Phase 2: NSGA-II benchmark, fixed-threshold baselines
│       └── out/                  fronts, direct-EnergyPlus logs, Tables 5, D1, E1
└── calibration/
    ├── trackA/                   passive thermal calibration (EnergyPlus 9.4 + eppy)
    │   ├── weather/              measured-EPW builder (Nov–Dec 2021)
    │   └── runs/                 calibrated v3/v6 series (gzip) + metric tables
    └── trackB/                   CO2 mass-balance calibration (pure R)
        └── runs/                 Tables B1–B5
```

Section 6 maps every script to the table or figure it produces.

## 2. Software environment

| Component | Version used in the study | Needed for |
|---|---|---|
| R | ≥ 4.3 (study runs on Windows; exact version in `environment/sessionInfo_author.txt`) | everything |
| R packages | `here`, `ggplot2`, `lubridate`, `readxl`, `patchwork`, `scales`, `parallel`, `DiceKriging`, `DiceOptim`, `mco`, `reticulate`, `rstudioapi` | see `environment/` |
| EnergyPlus | 9.4.0 (exported from DesignBuilder 7.3.1.003) | annual simulations, Track A |
| Python | 3.13 + `eppy` 0.5.69 (Track A); `pandas`, `numpy`, `pvlib` (EPW builder) | Track A, weather builder |

Exact versions of the run are recorded in `environment/sessionInfo_author.txt`. The verification environment (Linux, R 4.3.3) is recorded in `environment/sessionInfo_verification.txt`.

**Machine-specific settings.** The EnergyPlus-coupled scripts were run on Windows. They assume:

- EnergyPlus installed in `C:/EnergyPlusV9-4-0/` (`fbo/F1/fbo_eval.R`, `calibration/trackA/run_calibration.R`, `calibration_engine.py`).
- An ASCII-only work folder `C:/EP_work/`. EnergyPlus 9.4 fails on paths with non-ASCII characters. The production IDF reads its outdoor air-quality schedule from `C:/EP_work/AQI_ext_2024_5min.csv`, which `fbo/F1/make_aqi_ext.R` writes.
- The Python interpreter path in `calibration/trackA/.Rprofile` and `run_calibration.R` (`RETICULATE_PYTHON`).

Edit these lines for another machine. Cloud-sync clients (OneDrive, Dropbox) must be paused during EnergyPlus runs, because file locks corrupt the outputs.

## 3. Weather files (not redistributed)

The IWEC weather file for Seville is © ASHRAE and is distributed by the U.S. DOE through the EnergyPlus website. It is therefore not included here.

1. Download `ESP_Sevilla.083910_IWEC.epw` from <https://energyplus.net/weather> (Europe → Spain → Sevilla, IWEC).
2. Check its SHA-256: `c060cd44c9b8c5819979627992ccae82fc7623f43e7bc1a12967a8866f777b36`.
3. Copy it to `fbo/data/`. It must be the only `.epw` file there, because the drivers take the first `.epw` they find. It is the weather file of all optimization and benchmark runs (verified in `fbo/F1/anchor_check_epw.R`).
4. For Track A and B, build the calibration EPW from the IWEC template and the measured on-site record:

   ```bash
   cd calibration/trackA/weather
   cp ../../../fbo/data/ESP_Sevilla.083910_IWEC.epw ESP_Sevilla_083910_IWEC.epw
   python build_measured_epw.py ../T3_Calibration_TRACK_A_50_nights.xlsx ESP_Sevilla_083910_IWEC.epw Seville_measured_NovDec2021.epw
   ```

   The output must have SHA-256 `111fd24a618bf7fce3b2682f43c6e03e5753900b6a65304bfd47b8557ad9b8c7`. This byte-identical rebuild was verified with `pvlib` 0.15.1, `pandas` 3.0.2 and `numpy` 2.4.4.

## 4. Quick verification (no EnergyPlus)

```bash
Rscript verify_without_energyplus.R      # from the repository root, ~2 min
```

The script copies the repository to a temporary folder and re-runs all stages that need no EnergyPlus. It then compares the regenerated outputs with the released ones. Expected result: `10/10 checks passed`.

| Check | Reproduced value |
|---|---|
| Phase 0 unit tests (FIS, g1–g5, control law) | 148/148 |
| Optima U* of the nine scenarios (Table 4) | rebuilt from seed + BO logs, max abs(ΔU) = 0 |
| Surrogate gate, metrics (b)/(c) | `F2/MASTER_FIS_F2.R` completes (see note below) |
| Table 5: on-front distance, ΔU, N_FBO, N_NSGA, ρ | identical (14.0×, 12.0×, 4.4×) |
| Table D1: g₄ anchoring sensitivity | identical |
| Track B: G = 14.74 L h⁻¹, 95 % CI [12.48, 17.00]; Tables B1–B5 | identical |
| Track A: geometry audit, UA = 99.3 W K⁻¹ | 12/12 checks |
| Track A: v3 NMBE −5.41 % / CV-RMSE 8.04 %; v6 +6.05 % | reproduced |
| Track A: weekday/weekend nights −4.94 % / −6.03 % | reproduced |
| Track A: median τ = 38.09 h over 43 nights | reproduced |

**Note on the surrogate gate (§3.3.1).** The cross-validation errors depend on the maximum-likelihood fit of the GP hyperparameters, which differs slightly across platforms and `DiceKriging` builds. The study run (Windows) gave CV-RMSE = 22 % for g₃ in E-w410; the verification run (Linux, `DiceKriging` 1.6.1) gives 20.4 %. In both, E-w410/g₃ is the only case above the 20 % threshold and every other indicator passes; no reported result depends on the surrogate.

## 5. Full reproduction (EnergyPlus required)

Each orchestrator is fail-fast. Its expensive stages are switched off by default (`RUN_* <- FALSE`), so a plain `source()` only re-validates and rebuilds from the released logs. Set a toggle to `TRUE` to re-run that stage.

| Stage | Entry point (open the `.Rproj` first) | Cost |
|---|---|---|
| Track A calibration | `calibration/trackA/MASTER_trackA.R` (run `unpack_runs.R` first) | ≈ 3 runs of the Nov–Dec 2021 period |
| Track A parameter sweep | `calibration/trackA/run_calibration_order.R` (line by line) | 120 runs, ~50 min |
| Track B calibration | `calibration/trackB/MASTER_trackB.R` | ~5 s, no EnergyPlus |
| Phase 0 | `fbo/F0/MASTER_FIS.R` | ~10 s, no EnergyPlus |
| Phase 1 (seed + BO) | `fbo/F1/MASTER_FIS_F1.R`, `RUN_SEED`, `RUN_BO` | 450 + 270 = 720 annual runs |
| Phase 2 (NSGA-II direct) | `fbo/F2/MASTER_FIS_F2.R`, `RUN_DIRECT_EP`; `F2/run_direct_E_w410.R` | ≈ 1040 annual runs per scenario |
| Fixed-threshold baselines | `fbo/F2/run_classic_baseline.R`, `run_classic800_baseline.R` | 9 + 9 annual runs |

One annual run takes ≈ 2–2.5 min on a desktop workstation. The production IDF `fbo/F1/v4_FBO_annual_prod.idf` is regenerated from `v4_FBO_annual.idf` by `fbo/F1/make_prod_idf.R`; the regenerated file is object-for-object identical (336/336 objects). The published Track A IDF `T1_aula_base_94_v4.idf` is regenerated byte-identically from v3 by `make_v4.R`.

## 6. Script → manuscript map

| Manuscript element | Script | Output |
|---|---|---|
| Fig. 3 (membership functions) | `fbo/core/plot_mf_panel.R` (via `F0/MASTER_FIS.R`) | `fbo/runs/T3_Fig_MF_primaries_3panel.*` |
| Fig. 4 (Track A, representative night) | `calibration/trackA/fig5_representative_night.R` | `runs/deviation_figs/` |
| Fig. 5 (Track B, measured vs simulated CO₂) | `calibration/trackB/trackB_figures.R` | `runs/FigB1_trackB_fit.png` |
| Table 4 (optima) | `fbo/F1/MASTER_FIS_F1.R` → `build_best_x()` | `fbo/F1/bo/best_x_9scen_v2.csv` |
| Fig. 6 (U* by scenario) | `fbo/F1/fig_bo_utility.R` | `fbo/runs/Fig_BO_utility_by_scenario.*` |
| §3.3.1 (surrogate gate) | `fbo/F2/surrogate.R` (`f2_gate()`) | console |
| Figs. 7–8, Table 5 | `fbo/F2/figs.R`, `rescore_direct_ep.R` | `F2/fig/`, `F2/out/rescore_direct_summary.csv` |
| Fig. 9, Table E1 | `fbo/F2/run_classic_baseline.R`, `run_classic800_baseline.R`, `fig_classic3.R` | `F2/out/classic*_vs_opt.csv` |
| Table A1 (UA) | `calibration/trackA/UA_aula_v4_fixed.R` | console |
| Fig. A1 (τ per night) | `calibration/trackA/tau_per_night.R`, `figA3_tau.R` | `runs/tau_per_night.csv` |
| Fig. A2, Table A3 (v3 vs v6) | `close_trackA_v6.R`, `figA_calibration_v3v6.R`, `deviation_fullseries.R` | `runs/` |
| Weekend cross-validation (§3.1) | `calibration/trackA/deviation_nocturnal.R` | console |
| Fig. B1 (profile likelihood) | `calibration/trackB/figB1_profile_likelihood.R` | `runs/FigB1_profile_likelihood.*` |
| Tables B1–B4 | `calibration/trackB/trackB_summary_tables.R` | `runs/TableB*.csv` |
| Fig. B2 (g₅ sensitivity) | `fbo/core/sensitivity_g5.R` | `fbo/runs/T3_g5_sensitivity.*` |
| Table B5, Tables C1–C2 (FIS) | `fbo/core/fis_scoring.R` | rule bases in code |
| Table D1 | `fbo/F2/sens_g4_anchoring.R` | `F2/out/sens_g4_anchoring.csv` |

Items named in §2.9 of the paper: the wind-pressure-coefficient sets are the 16 `AirflowNetwork:MultiZone:WindPressureCoefficientValues` objects of the production IDF; the orientation scenarios rotate the building through `Building, North Axis` (NE = 0°, N = 329°, E = 59°). The school calendar is coded in `fbo/F1/fbo_eval.R` (`evaluate_U`). The LHS seed is generated by `fbo/F1/seed_driver.R` (maximin LHS, seed 20260728) and its evaluations are in `fbo/F1/seed_9scen_results.csv`.

## 7. Data and provenance

| File | Content | Source | Licence |
|---|---|---|---|
| `fbo/data/*_HH_2024_SEV_Torneo.csv` | hourly PM₂.₅, NO₂, O₃, station 41091009 (Torneo, Seville), 2024 | Red de Vigilancia y Control de la Calidad del Aire de Andalucía, Junta de Andalucía open-data portal | CC BY 4.0 (original publisher) |
| `calibration/trackA/meas_clean.csv`, `T3_Calibration_TRACK_A_50_nights.xlsx` | 5-min indoor/outdoor temperature and on-site weather, 1 Nov–20 Dec 2021 | authors' monitoring campaign | CC BY 4.0 |
| `calibration/trackB/co2_clean.csv`, `T3_Calibration_TRACK B_3 days.xlsx` | 5-min indoor CO₂, occupancy and opening log, 18 Nov, 26 Nov, 17 Dec 2021 | authors' monitoring campaign | CC BY 4.0 |
| `*.idf`, `runs/`, `F1/bo/`, `F2/out/` | models, calibrated series and evaluation logs | this study | CC BY 4.0 |

The monitoring data contain no personal information. Occupancy is an aggregate head count. The school is not identified; identifying names were removed from non-functional IDF fields (location name, comments, construction and schedule labels) without any effect on the simulation.

## 8. Licence

Code (`*.R`, `*.py`): MIT (`LICENSE`). Data, models and logs: CC BY 4.0 (`LICENSE-DATA.md`). Third-party data keep their original terms (Section 7).

## 9. How to cite

Cite the article and the archived software version (see `CITATION.cff`, or the "Cite this repository" button on GitHub):

> Sánchez Cordero, A., Gómez Melgar, S., & Andújar Márquez, J. M. (2026). *FBO-DNVP: Fuzzy–Bayesian Optimization of dynamic natural ventilation protocols* (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX

## 10. Contact

Antonio Sánchez Cordero — antonio.sanchez443@alu.uhu.es
Centro de Investigación en Tecnología, Energía y Sostenibilidad (CITES), Universidad de Huelva, Spain.

Changes made to the code for this public release are listed in `CHANGELOG.md`. Generative-AI assistance in drafting parts of the R scripts is declared in the article; every module was verified by execution against the anchor checks and unit tests.
