# Changelog

## [1.0.0] – 2026-11-17 — release accompanying the Building and Environment submission

First public release. The files are the authors' production scripts and logs,
with the changes below. None of them alters a reported number:
`verify_without_energyplus.R` reproduces all released outputs (10/10 checks).

### Portability patches (marked `[repo patch Pn]` in the code)

| Id | File | Change | Reason |
|---|---|---|---|
| P1 | `fbo/F1/MASTER_FIS_F1.R` | `source(core/control_law.R)` added before the validation gate | `run_tests_control()` was undefined in a fresh R session |
| P2 | `calibration/trackA/weather/` | folder renamed from `Weather/` to lower case | `run_calibration.R` and `calibration_engine.py` read `weather/…`; case-sensitive file systems failed |
| P3 | `calibration/trackB/figB1_profile_likelihood.R` | fallback working directory = `getwd()` | the fallback was an absolute path of the authors' machine |
| P4 | `fbo/F1/MASTER_FIS_F1.R` | `build_best_x()` writes to `F1/bo/best_x_9scen_v2.csv` | canonical location read by `F2/` and `anchor_check_epw.R` |
| P5 | `fbo/F1/bo_driver.R` | EPW presence checked inside `run_bo()` instead of at `source()` | the IWEC EPW is not redistributed; sourcing must work without it |
| P6 | `fbo/F2/fig_classic3.R` | invisible `geom_blank()` layer added first | ggplot2 < 3.5 failed ("Discrete value supplied to continuous scale"); figure unchanged |

### Added

- `fbo/F1/fig_bo_utility.R`: R script for Fig. 6, replacing the original matplotlib figure for which no script had been kept. Same data and layout; utilities labelled from full-precision U* (N-w216 0.410, E-w216 0.443).
- `verify_without_energyplus.R`, `calibration/trackA/unpack_runs.R`, `environment/capture_session_info.R`.

### Anonymization

Identifying names of the school and local file paths were replaced in
non-functional fields of the IDF and `eplustbl.htm` files (location name,
design-day names, construction and schedule labels, header comments). Object
references were renamed consistently. `make_v4.R` regenerates the released v4
IDF byte-identically, and `make_prod_idf.R` regenerates the production IDF
object-for-object (336/336 objects).

### Packaging

- Track A calibrated series stored as `runs/v3_trackA/eplusout.csv.gz` and
  `runs/v6_trackA/eplusout.csv.gz`; `unpack_runs.R` restores them.
- `.here` marker in `fbo/` so `here::here()` resolves to the FIS project root.
- RStudio projects renamed `trackA.Rproj` / `trackB.Rproj` (no spaces).

### Not included (retained in the thesis archive)

Laboratory notebooks and superseded scripts: `_archive/` (incl. `PRUEBAS
250724.R`, `g2_transfer.R`), `F1/_archive/` (G0 smoke tests), `lhs_pilot.R` and
the 12-point pilot design, `fig_classic.R` (superseded by `fig_classic3.R`),
`closure_trackA_v3_v6.R` and `UA_aula_v3.R` (retired; see `MASTER_trackA.R`),
`run_trackA_pipeline.R` (superseded by `MASTER_trackA.R`),
`weekend_freerun_verification.R` (superseded by `deviation_nocturnal.R`),
`plot_trackA_verification.R`, `g0_anchor.R` (Track B to Phase-1 hand-over test),
the Python port `fis_scoring.py` (pre-dates the final FIS), intermediate IDFs,
the raw 120-run sweep folders and internal hand-over notes.
