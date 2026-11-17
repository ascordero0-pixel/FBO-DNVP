# -*- coding: utf-8 -*-
"""
calibration_engine.py
=====================================================================
T3 digital-twin calibration engine — Seville classroom (T1 model).
Tool stack (FIXED): EnergyPlus 9.4.0 / DesignBuilder 7.3.1.003 /
                    eppy 0.5.69 / Python 3.13 / R 4.5.1 + reticulate.

Design contract (NO-FABRICATED-DATA POLICY):
  - This module ONLY computes NMBE / CV-RMSE from measured vs. simulated
    series. It never returns hard-coded metric values.
  - Every numeric metric that reaches the manuscript must originate from
    `compute_metrics()` applied to a real EnergyPlus run aligned with a
    real measured file. Placeholders stay labelled until then.

Two calibration tracks (see protocol §1, §5):
  Track A — passive thermal : calibrate infiltration crack multiplier
            and ZoneCapacitanceMultiplier C against measured indoor T.
  Track B — CO2 mass balance: calibrate g_person and C_d against
            measured indoor CO2, with Track-A infiltration FIXED.

Author: A. Sánchez Cordero (ASC), TEP192, Universidad de Huelva, 2026.
=====================================================================
"""

from __future__ import annotations

import os
import csv
import shutil
import subprocess
import itertools
from dataclasses import dataclass, field, asdict
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

# eppy is imported lazily inside functions so that the metrics layer can be
# unit-tested without an EnergyPlus install present.


# ---------------------------------------------------------------------------
# 0. Global configuration (edit paths to the local machine ONCE)
# ---------------------------------------------------------------------------
@dataclass
class EPConfig:
    """Absolute paths and fixed engine identifiers (protocol §2)."""
    idd_path: str = r"C:/EnergyPlusV9-4-0/Energy+.idd"
    eplus_exe: str = r"C:/EnergyPlusV9-4-0/energyplus.exe"
    epw_path: str = r"./weather/Seville_measured_NovDec2021.epw"
    work_dir: str = r"./runs"
    engine_version: str = "9.4.0.002"          # confirmed from IDF (§8)
    designbuilder_version: str = "7.3.1.003"   # documented for manuscript

    def ensure_dirs(self) -> None:
        os.makedirs(self.work_dir, exist_ok=True)


# ---------------------------------------------------------------------------
# 1. ASHRAE Guideline 14 metrics
# ---------------------------------------------------------------------------
# Definitions (hourly basis is the G14 comparison basis; 5-min reported as
# the conservative companion — protocol §6).
#
#   NMBE (%) = 100 * ( Σ (s_i - m_i) ) / ( (n - p) * mean(m) )
#   CV-RMSE (%) = 100 * sqrt( Σ (m_i - s_i)^2 / (n - p) ) / mean(m)
#
# where m_i = measured, s_i = simulated, n = number of paired points,
# p = number of adjustable model parameters in the calibration
# (G14 "degrees of freedom" deflation). For a single held-out validation
# period p is conventionally set to 0 or 1; this is a documented modelling
# choice — see the technical note. We expose `p` explicitly so the manuscript
# can state the value used and run sensitivity on it.
#
# Sign convention (protocol §6):  Σ (m_i - s_i) < 0  ⇔  model OVER-predicts.
# NMBE here is built on (s_i - m_i); therefore NMBE > 0 ⇔ over-prediction.
# This is declared explicitly to avoid the recurrent sign ambiguity.

def _validate_pair(measured: np.ndarray, simulated: np.ndarray) -> None:
    if measured.shape != simulated.shape:
        raise ValueError(
            f"Shape mismatch: measured {measured.shape} vs simulated "
            f"{simulated.shape}. Align series before computing metrics."
        )
    if measured.size == 0:
        raise ValueError("Empty series passed to metrics.")
    if np.any(~np.isfinite(measured)) or np.any(~np.isfinite(simulated)):
        raise ValueError("NaN/Inf present. Clean/align series first.")


def nmbe(measured: np.ndarray, simulated: np.ndarray, p: int = 0) -> float:
    """Normalised Mean Bias Error (%), sign: >0 ⇔ model over-predicts."""
    measured = np.asarray(measured, dtype=float)
    simulated = np.asarray(simulated, dtype=float)
    _validate_pair(measured, simulated)
    n = measured.size
    denom_dof = (n - p)
    if denom_dof <= 0:
        raise ValueError(f"(n - p) = {denom_dof} ≤ 0. Reduce p or add data.")
    mbar = measured.mean()
    if np.isclose(mbar, 0.0):
        raise ValueError("Mean of measured ≈ 0; NMBE undefined.")
    return 100.0 * np.sum(simulated - measured) / (denom_dof * mbar)


def cv_rmse(measured: np.ndarray, simulated: np.ndarray, p: int = 0) -> float:
    """Coefficient of Variation of the RMSE (%)."""
    measured = np.asarray(measured, dtype=float)
    simulated = np.asarray(simulated, dtype=float)
    _validate_pair(measured, simulated)
    n = measured.size
    denom_dof = (n - p)
    if denom_dof <= 0:
        raise ValueError(f"(n - p) = {denom_dof} ≤ 0. Reduce p or add data.")
    mbar = measured.mean()
    if np.isclose(mbar, 0.0):
        raise ValueError("Mean of measured ≈ 0; CV-RMSE undefined.")
    rmse = np.sqrt(np.sum((measured - simulated) ** 2) / denom_dof)
    return 100.0 * rmse / mbar


# ASHRAE Guideline 14 hourly acceptance thresholds (protocol §1).
G14 = {"NMBE_abs_max": 10.0, "CVRMSE_max": 30.0}


def passes_g14(nmbe_val: float, cvrmse_val: float) -> bool:
    return (abs(nmbe_val) <= G14["NMBE_abs_max"]) and (cvrmse_val <= G14["CVRMSE_max"])


def compute_metrics(measured: np.ndarray, simulated: np.ndarray,
                    p: int = 0, label: str = "") -> dict:
    """Return both metrics + G14 verdict for a single aligned pair."""
    nv = nmbe(measured, simulated, p=p)
    cv = cv_rmse(measured, simulated, p=p)
    return {
        "label": label,
        "n": int(np.asarray(measured).size),
        "p": p,
        "NMBE_pct": round(nv, 4),
        "CVRMSE_pct": round(cv, 4),
        "mean_measured": round(float(np.asarray(measured, float).mean()), 4),
        "G14_pass": passes_g14(nv, cv),
    }


# ---------------------------------------------------------------------------
# 2. Series alignment (5-min native → hourly aggregation)
# ---------------------------------------------------------------------------
def align_series(meas: pd.DataFrame, sim: pd.DataFrame,
                 ts_col: str = "timestamp", val_col: str = "value",
                 how: str = "inner") -> pd.DataFrame:
    """
    Inner-join measured and simulated on a shared timestamp index.
    Both frames must carry [ts_col, val_col]. Returns columns
    [timestamp, measured, simulated]. Strict: no interpolation across gaps —
    a missing measured point drops the pair (protects metric integrity).
    """
    m = meas[[ts_col, val_col]].rename(columns={val_col: "measured"}).copy()
    s = sim[[ts_col, val_col]].rename(columns={val_col: "simulated"}).copy()
    m[ts_col] = pd.to_datetime(m[ts_col])
    s[ts_col] = pd.to_datetime(s[ts_col])
    out = pd.merge(m, s, on=ts_col, how=how).sort_values(ts_col)
    out = out.dropna(subset=["measured", "simulated"]).reset_index(drop=True)
    return out


def to_hourly(aligned: pd.DataFrame, ts_col: str = "timestamp") -> pd.DataFrame:
    """Mean-aggregate an aligned 5-min frame to the hourly G14 basis."""
    a = aligned.copy()
    a[ts_col] = pd.to_datetime(a[ts_col])
    a = a.set_index(ts_col)
    hourly = a.resample("1h").mean().dropna().reset_index()
    return hourly


# ---------------------------------------------------------------------------
# 3. eppy IDF editing helpers (protocol §3 Phase 4, §5)
# ---------------------------------------------------------------------------
def load_idf(cfg: EPConfig, idf_path: str):
    """Load an IDF with the fixed IDD (eppy)."""
    from eppy.modeleditor import IDF
    IDF.setiddname(cfg.idd_path)
    return IDF(idf_path)


def set_timestep(idf, intervals_per_hour: int = 12) -> None:
    """Timestep 1 (hourly) → 12 (5-min). Protocol §4/§8."""
    ts = idf.idfobjects["TIMESTEP"]
    if ts:
        ts[0].Number_of_Timesteps_per_Hour = intervals_per_hour
    else:
        idf.newidfobject("TIMESTEP", Number_of_Timesteps_per_Hour=intervals_per_hour)


def set_capacitance_multiplier(idf, c_mult: float) -> None:
    """ZoneCapacitanceMultiplier:ResearchSpecial — Track A parameter C."""
    objs = idf.idfobjects.get("ZONECAPACITANCEMULTIPLIER:RESEARCHSPECIAL", [])
    if not objs:
        raise RuntimeError("ZoneCapacitanceMultiplier:ResearchSpecial absent; "
                           "expected present at 1.0 per §8.")
    # Temperature capacitance multiplier field (sensible heat).
    objs[0].Temperature_Capacity_Multiplier = c_mult


def scale_afn_crack_coeffs(idf, scale: float) -> int:
    """
    Scale AirflowNetwork:MultiZone:Surface:Crack flow coefficients
    (infiltration). 11 cracks present per §8. Returns count modified.
    Calibration acts on the multiplier, not on a single ACH value.
    """
    cracks = idf.idfobjects.get("AIRFLOWNETWORK:MULTIZONE:SURFACE:CRACK", [])
    if not cracks:
        raise RuntimeError("No AFN Surface:Crack objects found; check IDF.")
    for cr in cracks:
        base = float(cr.Air_Mass_Flow_Coefficient_at_Reference_Conditions)
        cr.Air_Mass_Flow_Coefficient_at_Reference_Conditions = base * scale
    return len(cracks)


def close_operable_openings(idf) -> int:
    """
    Force all OPERABLE AFN openings to NoVent (closed) for the Track A
    thermal calibration. Rationale: Track A targets the free-running NIGHT
    behaviour with the classroom CLOSED (no purpose-provided ventilation),
    so window/door openings must not vent — only infiltration through the
    fabric cracks remains active (those are what Track A calibrates).

    Distinguishes openings from cracks by Ventilation_Control_Mode: surfaces
    with a real control mode (e.g. 'Temperature') are operable openings and
    are set to 'NoVent'; crack/leak surfaces (blank control) are left intact.

    Track B does NOT call this — there the openings operate per the DNVP
    schedules during class hours, since ventilation is the studied phenomenon.

    Returns the number of openings closed.
    """
    surfs = idf.idfobjects.get("AIRFLOWNETWORK:MULTIZONE:SURFACE", [])
    closed = 0
    for s in surfs:
        ctrl = str(getattr(s, "Ventilation_Control_Mode", "") or "").strip()
        # Operable openings carry an active control mode (Temperature, Constant,
        # Enthalpy, ASHRAE55Adaptive, ...). NoVent / blank are already closed.
        if ctrl and ctrl.lower() not in ("novent", ""):
            s.Ventilation_Control_Mode = "NoVent"
            closed += 1
    return closed


def add_output_variables(idf, names: Sequence[str] = (
        "Zone Mean Air Temperature",
        "Site Outdoor Air Drybulb Temperature",
        "Zone Air CO2 Concentration",
), freq: str = "Timestep") -> None:
    """Add Output:Variable rows (§4: Ti interior + Ot exterior at 5-min)."""
    existing = {(o.Variable_Name, o.Key_Value)
                for o in idf.idfobjects.get("OUTPUT:VARIABLE", [])}
    for nm in names:
        if (nm, "*") not in existing:
            idf.newidfobject("OUTPUT:VARIABLE", Key_Value="*",
                             Variable_Name=nm, Reporting_Frequency=freq)


def save_idf(idf, out_path: str) -> str:
    idf.saveas(out_path, encoding='utf-8')
    return out_path


# ---------------------------------------------------------------------------
# 4. EnergyPlus run wrapper (protocol §7 step 2)
# ---------------------------------------------------------------------------
def run_energyplus(cfg: EPConfig, idf_path: str, run_tag: str) -> str:
    """
    Invoke EnergyPlus via subprocess. Returns the run output directory.
    Raises on non-zero exit so a failed run never silently yields metrics.
    """
    cfg.ensure_dirs()
    out_dir = os.path.join(cfg.work_dir, run_tag)
    os.makedirs(out_dir, exist_ok=True)
    cmd = [cfg.eplus_exe, "-w", cfg.epw_path, "-d", out_dir,
           "-r", "-x", idf_path]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"EnergyPlus run '{run_tag}' failed (exit {proc.returncode}).\n"
            f"STDERR tail:\n{proc.stderr[-2000:]}"
        )
    return out_dir


def read_eplus_variable(out_dir: str, column_contains: str,
                        csv_name: str = "eplusout.csv") -> pd.DataFrame:
    """
    Read a single time-series column from eplusout.csv by substring match.
    Returns [timestamp, value]. EnergyPlus 'Date/Time' is parsed to a
    naive datetime on the run year.
    """
    path = os.path.join(out_dir, csv_name)
    df = pd.read_csv(path)
    time_col = df.columns[0]  # 'Date/Time'
    matches = [c for c in df.columns if column_contains.lower() in c.lower()]
    if not matches:
        raise KeyError(f"No column matching '{column_contains}' in {path}. "
                       f"Available: {list(df.columns)[:6]} ...")
    col = matches[0]
    ts = _parse_eplus_datetime(df[time_col])
    return pd.DataFrame({"timestamp": ts, "value": df[col].astype(float)})


def _parse_eplus_datetime(series: pd.Series, year: int = 2021) -> pd.Series:
    """
    Parse EnergyPlus ' MM/DD  HH:MM:SS' strings. EnergyPlus writes hour 24:00
    for end-of-day; remap to 00:00 of the next day to keep pandas happy.
    """
    out = []
    for raw in series.astype(str):
        s = raw.strip()
        date_part, time_part = s.split("  ") if "  " in s else s.split(maxsplit=1)
        mm, dd = date_part.split("/")
        hh, mi, ss = time_part.split(":")
        roll = False
        if hh == "24":
            hh = "00"
            roll = True
        dt = pd.Timestamp(year=year, month=int(mm), day=int(dd),
                          hour=int(hh), minute=int(mi), second=int(ss))
        if roll:
            dt += pd.Timedelta(days=1)
        out.append(dt)
    return pd.Series(out)


# ---------------------------------------------------------------------------
# 5. Parameter grids and the calibration loop (protocol §7)
# ---------------------------------------------------------------------------
@dataclass
class TrackAGrid:
    """Track A search space: infiltration crack scale × capacitance C."""
    crack_scale: Sequence[float] = field(
        default_factory=lambda: [0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    cap_mult: Sequence[float] = field(
        default_factory=lambda: [1.0, 3.0, 6.0, 10.0, 15.0])

    def combinations(self) -> Iterable[tuple[float, float]]:
        return itertools.product(self.crack_scale, self.cap_mult)


@dataclass
class TrackBGrid:
    """Track B search space: g_person × C_d (infiltration fixed from A)."""
    g_person: Sequence[float] = field(
        default_factory=lambda: [0.0030, 0.0033, 0.0036])  # L/s/person, Persily
    c_d: Sequence[float] = field(
        default_factory=lambda: [0.55, 0.60, 0.65, 0.70])   # T1 prior 0.65

    def combinations(self) -> Iterable[tuple[float, float]]:
        return itertools.product(self.g_person, self.c_d)


def calibrate_track_a(cfg: EPConfig, base_idf_path: str,
                      measured_T: pd.DataFrame, grid: TrackAGrid,
                      p: int = 1, results_csv: str = "results_trackA.csv") -> pd.DataFrame:
    """
    Sweep (crack_scale, cap_mult); for each, edit IDF → run E+ → read Ti →
    align with measured → compute 5-min and hourly metrics. Returns a frame
    sorted by hourly CV-RMSE ascending. Writes results_csv.

    `measured_T` must be [timestamp, value] in °C at 5-min resolution.
    """
    rows = []
    for scale, cmult in grid.combinations():
        tag = f"A_s{scale}_c{cmult}"
        idf = load_idf(cfg, base_idf_path)
        set_timestep(idf, 12)
        n_cracks = scale_afn_crack_coeffs(idf, scale)
        n_closed = close_operable_openings(idf)   # Track A: classroom CLOSED at night
        set_capacitance_multiplier(idf, cmult)
        add_output_variables(idf)
        run_idf = save_idf(idf, os.path.join(cfg.work_dir, f"{tag}.idf"))
        out_dir = run_energyplus(cfg, run_idf, tag)
        sim = read_eplus_variable(out_dir, "Zone Mean Air Temperature")

        aligned = align_series(measured_T, sim)
        m5 = compute_metrics(aligned["measured"].values,
                             aligned["simulated"].values, p=p, label=f"{tag}_5min")
        hourly = to_hourly(aligned)
        mh = compute_metrics(hourly["measured"].values,
                             hourly["simulated"].values, p=p, label=f"{tag}_hourly")
        rows.append({
            "tag": tag, "crack_scale": scale, "cap_mult": cmult,
            "n_cracks": n_cracks,
            "NMBE_5min": m5["NMBE_pct"], "CVRMSE_5min": m5["CVRMSE_pct"],
            "NMBE_hourly": mh["NMBE_pct"], "CVRMSE_hourly": mh["CVRMSE_pct"],
            "G14_hourly_pass": mh["G14_pass"], "n_hourly": mh["n"],
        })
    res = pd.DataFrame(rows).sort_values("CVRMSE_hourly").reset_index(drop=True)
    res.to_csv(results_csv, index=False)
    return res


def calibrate_track_b(cfg: EPConfig, base_idf_path: str,
                      measured_CO2: pd.DataFrame, grid: TrackBGrid,
                      fixed_crack_scale: float, p: int = 1,
                      results_csv: str = "results_trackB.csv") -> pd.DataFrame:
    """
    Track B: infiltration crack scale FIXED to Track-A optimum; sweep
    (g_person, C_d) against measured indoor CO2. g_person is applied via the
    occupancy CO2 generation; C_d via AFN:MultiZone:Opening discharge coeff.
    Returns frame sorted by hourly CV-RMSE; writes results_csv.

    NOTE on identifiability (open issue in project notes): C_d and ΔC_p can
    trade off. Here infiltration is locked from Track A to reduce the
    parameter space; the residual C_d/ΔC_p coupling is examined in the
    technical note, not resolved by brute force.
    """
    rows = []
    for g_p, cd in grid.combinations():
        tag = f"B_g{g_p}_cd{cd}"
        idf = load_idf(cfg, base_idf_path)
        set_timestep(idf, 12)
        scale_afn_crack_coeffs(idf, fixed_crack_scale)   # locked from Track A
        _set_co2_generation(idf, g_p)
        _set_opening_discharge(idf, cd)
        add_output_variables(idf)
        run_idf = save_idf(idf, os.path.join(cfg.work_dir, f"{tag}.idf"))
        out_dir = run_energyplus(cfg, run_idf, tag)
        sim = read_eplus_variable(out_dir, "Zone Air CO2 Concentration")

        aligned = align_series(measured_CO2, sim)
        m5 = compute_metrics(aligned["measured"].values,
                             aligned["simulated"].values, p=p, label=f"{tag}_5min")
        hourly = to_hourly(aligned)
        mh = compute_metrics(hourly["measured"].values,
                             hourly["simulated"].values, p=p, label=f"{tag}_hourly")
        rows.append({
            "tag": tag, "g_person": g_p, "c_d": cd,
            "fixed_crack_scale": fixed_crack_scale,
            "NMBE_5min": m5["NMBE_pct"], "CVRMSE_5min": m5["CVRMSE_pct"],
            "NMBE_hourly": mh["NMBE_pct"], "CVRMSE_hourly": mh["CVRMSE_pct"],
            "G14_hourly_pass": mh["G14_pass"], "n_hourly": mh["n"],
        })
    res = pd.DataFrame(rows).sort_values("CVRMSE_hourly").reset_index(drop=True)
    res.to_csv(results_csv, index=False)
    return res


def _set_co2_generation(idf, g_person_l_s: float) -> None:
    """
    Apply per-person CO2 generation rate. EnergyPlus drives People CO2 via
    Carbon Dioxide Generation Rate [m3/s-person] on People objects (with
    ZoneAirContaminantBalance enabled). Convert L/s → m3/s (×1e-3).
    """
    rate_m3s = g_person_l_s * 1e-3
    people = idf.idfobjects.get("PEOPLE", [])
    if not people:
        raise RuntimeError("No People objects; cannot set CO2 generation.")
    for pe in people:
        pe.Carbon_Dioxide_Generation_Rate = rate_m3s


def _set_opening_discharge(idf, c_d: float) -> None:
    """Set discharge coefficient on AFN:MultiZone:Component:DetailedOpening
    / SimpleOpening objects (Track B parameter C_d)."""
    set_any = False
    for key in ("AIRFLOWNETWORK:MULTIZONE:COMPONENT:DETAILEDOPENING",
                "AIRFLOWNETWORK:MULTIZONE:COMPONENT:SIMPLEOPENING"):
        for op in idf.idfobjects.get(key, []):
            if hasattr(op, "Discharge_Coefficient"):
                op.Discharge_Coefficient = c_d
                set_any = True
            # DetailedOpening uses per-factor discharge coeffs:
            for fld in ("Discharge_Coefficient_for_Opening_Factor_1",
                        "Discharge_Coefficient_for_Opening_Factor_2",
                        "Discharge_Coefficient_for_Opening_Factor_3",
                        "Discharge_Coefficient_for_Opening_Factor_4"):
                if hasattr(op, fld):
                    setattr(op, fld, c_d)
                    set_any = True
    if not set_any:
        raise RuntimeError("No AFN opening discharge coefficient field found.")


# ---------------------------------------------------------------------------
# 6. Validation folds (protocol §6)
# ---------------------------------------------------------------------------
def loo_cv_co2(cfg: EPConfig, base_idf_path: str,
               days: dict[str, pd.DataFrame], best_params: dict,
               fixed_crack_scale: float, p: int = 1) -> pd.DataFrame:
    """
    3-fold leave-one-out over the 3 monitored CO2 days. Headline fold:
    train {18 Nov, 26 Nov} → validate {17 Dec}. `days` maps a day label to
    its measured [timestamp, value] CO2 frame. Parameters are FIXED to
    `best_params` (selected on training folds); this routine reports the
    held-out validation metrics per fold without re-fitting (clean LOO).
    """
    labels = list(days.keys())
    rows = []
    for holdout in labels:
        idf = load_idf(cfg, base_idf_path)
        set_timestep(idf, 12)
        scale_afn_crack_coeffs(idf, fixed_crack_scale)
        _set_co2_generation(idf, best_params["g_person"])
        _set_opening_discharge(idf, best_params["c_d"])
        add_output_variables(idf)
        tag = f"LOO_val_{holdout}"
        run_idf = save_idf(idf, os.path.join(cfg.work_dir, f"{tag}.idf"))
        out_dir = run_energyplus(cfg, run_idf, tag)
        sim = read_eplus_variable(out_dir, "Zone Air CO2 Concentration")
        aligned = align_series(days[holdout], sim)
        hourly = to_hourly(aligned)
        mh = compute_metrics(hourly["measured"].values,
                             hourly["simulated"].values, p=p,
                             label=f"holdout_{holdout}")
        rows.append({"validation_day": holdout, **{k: mh[k] for k in
                     ("n", "NMBE_pct", "CVRMSE_pct", "G14_pass")}})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 7. Self-test of the metrics layer (runs without EnergyPlus)
# ---------------------------------------------------------------------------
def _selftest() -> None:
    rng = np.random.default_rng(42)
    n = 200
    measured = 18.0 + 2.0 * np.sin(np.linspace(0, 6, n))
    # simulated = measured with small bias + noise → known-sign metrics
    simulated = measured + 0.2 + rng.normal(0, 0.15, n)  # +0.2 ⇒ over-predict
    m = compute_metrics(measured, simulated, p=1, label="selftest")
    assert m["NMBE_pct"] > 0, "Sign convention broken: +bias must give NMBE>0"
    assert m["CVRMSE_pct"] > 0
    # zero-error case
    m0 = compute_metrics(measured, measured.copy(), p=0)
    assert abs(m0["NMBE_pct"]) < 1e-9 and m0["CVRMSE_pct"] < 1e-9
    print("[selftest] metrics layer OK:", m)


if __name__ == "__main__":
    _selftest()
# =====================================================================
# >>> APPEND THIS BLOCK TO THE END OF calibration_engine.py <<<
#     (paste AFTER calibrate_track_b / loo_cv_co2, BEFORE the
#      `def _selftest()` / `if __name__ == "__main__"` block)
# ---------------------------------------------------------------------
# Track A — v2 extension for the rebuilt IDF (T1_aula_base_94_v2.idf).
#
# What changed vs the original Track A loop, and WHY:
#   * v2 already contains People / Lights / OtherEquipment -> we DO NOT
#     add occupancy; we ONLY add heating (the 3rd calibration parameter).
#   * v2's window control is an EMS (blank AFN control mode) -> the legacy
#     close_operable_openings() is a NO-OP on v2; we remove the EMS
#     ProgramCallingManager AND force the window surfaces to NoVent.
#   * v2's Wall_3 exported a spurious GapConvectionRadiation cavity ->
#     normalised to Outdoors (sibling-consistent), shading preserved.
#   * v2 has 2 zones -> read the CLASSROOM Zone Mean Air Temperature by
#     an explicit zone-qualified column, never matches[0].
#   * 3rd acceptance gate: nocturnal cooling-rate error |slope_sim-slope_meas|
#     so heating power cannot silently absorb the mass/coupling deficit.
#
# Reuses the existing engine helpers (load_idf, set_timestep,
# scale_afn_crack_coeffs, set_capacitance_multiplier, add_output_variables,
# save_idf, run_energyplus, read_eplus_variable, align_series, to_hourly,
# compute_metrics) — all defined above in this same module.
# =====================================================================

CLASSROOM_ZONE = "ZXCASESTUDYCLASSROOM:CLASSROOM"
CLASSROOM_FLOOR_AREA_M2 = 45.2511   # E+ zone floor area, v3 (doc only; not used in calc)


# ---------------------------------------------------------------------
# A. IDF patches (Track A, v2)
# ---------------------------------------------------------------------
def add_heating_gain(idf, w_per_m2, *,
                     sched_name="T3_TrackA_Heating_Sched",
                     gain_name="T3_TrackA_Heating_Gain",
                     radiant_fraction=0.30, start_h=7, end_h=18):
    """Add/replace an ElectricEquipment heating gain (W/m2) on WEEKDAYS
    start_h..end_h. Idempotent. radiant_fraction (0.2-0.4 physical) is also
    a secondary lever on the nocturnal slope; keep it fixed during the P
    sweep. Holidays 6/8 Dec are treated as weekdays (RunPeriod holidays=No);
    add a 'For: Holidays' line returning 0 if you later flag school breaks.
    """
    stl = "Fraction"
    if not any(o.Name == stl for o in idf.idfobjects["ScheduleTypeLimits"]):
        idf.newidfobject("ScheduleTypeLimits", Name=stl,
                         Lower_Limit_Value=0, Upper_Limit_Value=1,
                         Numeric_Type="Continuous")
    for o in list(idf.idfobjects["Schedule:Compact"]):
        if o.Name == sched_name:
            idf.removeidfobject(o)
    sch = idf.newidfobject("Schedule:Compact", Name=sched_name)
    sch.Schedule_Type_Limits_Name = stl
    fields = ["Through: 12/31", "For: Weekdays",
              f"Until: {start_h:02d}:00", "0",
              f"Until: {end_h:02d}:00", "1",
              "Until: 24:00", "0",
              "For: AllOtherDays", "Until: 24:00", "0"]
    for i, val in enumerate(fields, start=1):
        setattr(sch, f"Field_{i}", val)

    for o in list(idf.idfobjects["ElectricEquipment"]):
        if o.Name == gain_name:
            idf.removeidfobject(o)
    eq = idf.newidfobject("ElectricEquipment", Name=gain_name)
    eq.Zone_or_ZoneList_Name = CLASSROOM_ZONE
    eq.Schedule_Name = sched_name
    eq.Design_Level_Calculation_Method = "Watts/Area"
    eq.Watts_per_Zone_Floor_Area = float(w_per_m2)
    eq.Fraction_Latent = 0.0
    eq.Fraction_Radiant = float(radiant_fraction)
    eq.Fraction_Lost = 0.0


def close_windows_for_trackA(idf):
    """Hard-close operable windows for free-running nights (supersedes
    close_operable_openings on v2, where AFN control modes are blank and an
    EMS drives the windows). (a) remove EMS ProgramCallingManager(s) so the
    Venting_Opening_Factor actuators never fire; (b) set the window AFN
    surfaces to NoVent. Cracks (Surface:Crack) untouched. Returns a report.
    """
    removed = 0
    for o in list(idf.idfobjects["EnergyManagementSystem:ProgramCallingManager"]):
        idf.removeidfobject(o)
        removed += 1
    window_names = {w.Name for w in idf.idfobjects["FenestrationSurface:Detailed"]
                    if str(w.Surface_Type).lower() == "window"}
    novent = 0
    for s in idf.idfobjects["AirflowNetwork:MultiZone:Surface"]:
        if s.Surface_Name in window_names:
            s.Ventilation_Control_Mode = "NoVent"
            novent += 1
    return {"ems_callmanagers_removed": removed, "afn_windows_set_novent": novent}


def normalize_wall3_to_outdoors(
        idf, surface=CLASSROOM_ZONE + "_Wall_3_0_0",
        oscm_name=CLASSROOM_ZONE + "_Wall_3_0_0 BoundaryObject"):
    """Make the 31 deg facade a plain exterior wall (sibling-consistent),
    dropping the spurious GapConvectionRadiation cavity. Shading untouched.
    """
    surf = next((s for s in idf.idfobjects["BuildingSurface:Detailed"]
                 if s.Name == surface), None)
    if surf is None:
        raise ValueError(f"Surface not found: {surface}")
    before = surf.Outside_Boundary_Condition
    surf.Outside_Boundary_Condition = "Outdoors"
    surf.Outside_Boundary_Condition_Object = ""
    surf.Sun_Exposure = "SunExposed"
    surf.Wind_Exposure = "WindExposed"
    surf.View_Factor_to_Ground = "autocalculate"
    dropped = 0
    for o in list(idf.idfobjects["SurfaceProperty:OtherSideConditionsModel"]):
        if o.Name == oscm_name:
            idf.removeidfobject(o)
            dropped += 1
    return {"wall3_obc_before": before, "wall3_obc_after": "Outdoors",
            "oscm_dropped": dropped}


# ---------------------------------------------------------------------
# B. Nocturnal cooling-rate metric (3rd acceptance gate)
# ---------------------------------------------------------------------
def _nightly_slope(aligned: pd.DataFrame, value_col: str,
                   night_start_h: int = 22, night_end_h: int = 7,
                   min_pts: int = 6) -> float:
    """Mean nocturnal cooling slope (deg C / h, negative) over every
    22:00->07:00 window of an aligned frame, each night fit by OLS.
    `aligned` carries [timestamp, measured, simulated]."""
    s = aligned[["timestamp", value_col]].copy()
    s["ts"] = pd.to_datetime(s["timestamp"])
    s = s.sort_values("ts")
    h = s["ts"].dt.hour + s["ts"].dt.minute / 60.0
    s = s[(h >= night_start_h) | (h < night_end_h)].copy()
    anchor = s["ts"].dt.normalize()
    anchor = anchor.where(s["ts"].dt.hour >= night_start_h,
                          anchor - pd.Timedelta(days=1))
    s["night"] = anchor
    slopes = []
    for _, g in s.groupby("night"):
        if len(g) < min_pts:
            continue
        x = (g["ts"] - g["ts"].iloc[0]).dt.total_seconds().to_numpy() / 3600.0
        y = g[value_col].to_numpy(dtype=float)
        slopes.append(float(np.polyfit(x, y, 1)[0]))
    return float(np.mean(slopes)) if slopes else float("nan")


# ---------------------------------------------------------------------
# C. 3-D parameter grid: crack_scale x cap_mult x heating
# ---------------------------------------------------------------------
@dataclass
class TrackAGridV2:
    """Track A v2 search space. heating_w_m2=0.0 is the gains-but-no-heating
    reference (v2 still has occupancy+lights+equipment). Physical radiator
    band is [20, 40] W/m2 for the rebuilt envelope."""
    crack_scale: Sequence[float] = field(
        default_factory=lambda: [0.5, 0.75, 1.0, 1.5, 2.0, 3.0])
    cap_mult: Sequence[float] = field(
        default_factory=lambda: [1.0, 3.0, 6.0, 10.0, 15.0])
    heating_w_m2: Sequence[float] = field(
        default_factory=lambda: [0.0, 20.0, 30.0, 40.0])

    def combinations(self):
        return itertools.product(self.crack_scale, self.cap_mult,
                                 self.heating_w_m2)


# ---------------------------------------------------------------------
# D. Track A v2 calibration loop (3 params, triple acceptance gate)
# ---------------------------------------------------------------------
def calibrate_track_a_v2(cfg: EPConfig, base_idf_path: str,
                         measured_T: pd.DataFrame, grid: TrackAGridV2,
                         p: int = 1, eps_slope: float = 0.06,
                         radiant_fraction: float = 0.30,
                         normalize_wall3: bool = True,
                         results_csv: str = "results_trackA_v2.csv") -> pd.DataFrame:
    """Sweep (crack_scale, cap_mult, heating_w_m2) on the v2 IDF.
    For each point: scale cracks, set capacitance, add heating, close
    windows/EMS, (optionally) normalise Wall_3 -> run E+ -> read CLASSROOM
    Ti (zone-qualified) -> align -> NMBE/CV-RMSE + nocturnal slope.

    Acceptance (triple gate): |NMBE_hourly|<=10 AND CVRMSE_hourly<=30 AND
    slope_abs_err<=eps_slope. Sorted by CVRMSE_hourly. Writes results_csv.

    `measured_T` = [timestamp, value] in degC at 5-min, nocturnal windows.
    """
    temp_col = f"{CLASSROOM_ZONE}:Zone Mean Air Temperature"  # 2-zone safe
    rows = []
    for scale, cmult, w_heat in grid.combinations():
        tag = f"Av2_s{scale}_c{cmult}_h{w_heat}"
        idf = load_idf(cfg, base_idf_path)
        set_timestep(idf, 12)
        n_cracks = scale_afn_crack_coeffs(idf, scale)
        set_capacitance_multiplier(idf, cmult)
        add_heating_gain(idf, w_heat, radiant_fraction=radiant_fraction)
        win = close_windows_for_trackA(idf)
        if normalize_wall3:
            normalize_wall3_to_outdoors(idf)
        add_output_variables(idf)
        run_idf = save_idf(idf, os.path.join(cfg.work_dir, f"{tag}.idf"))
        out_dir = run_energyplus(cfg, run_idf, tag)

        sim = read_eplus_variable(out_dir, temp_col)   # CLASSROOM zone only
        aligned = align_series(measured_T, sim)
        m5 = compute_metrics(aligned["measured"].values,
                             aligned["simulated"].values, p=p, label=f"{tag}_5min")
        hourly = to_hourly(aligned)
        mh = compute_metrics(hourly["measured"].values,
                             hourly["simulated"].values, p=p, label=f"{tag}_hourly")

        slope_meas = _nightly_slope(aligned, "measured")
        slope_sim = _nightly_slope(aligned, "simulated")
        slope_err = abs(slope_sim - slope_meas)

        accept = (abs(mh["NMBE_pct"]) <= G14["NMBE_abs_max"]
                  and mh["CVRMSE_pct"] <= G14["CVRMSE_max"]
                  and slope_err <= eps_slope)

        rows.append({
            "tag": tag, "crack_scale": scale, "cap_mult": cmult,
            "heating_w_m2": w_heat, "n_cracks": n_cracks,
            "ems_removed": win["ems_callmanagers_removed"],
            "windows_novent": win["afn_windows_set_novent"],
            "NMBE_5min": m5["NMBE_pct"], "CVRMSE_5min": m5["CVRMSE_pct"],
            "NMBE_hourly": mh["NMBE_pct"], "CVRMSE_hourly": mh["CVRMSE_pct"],
            "slope_sim_Cph": round(slope_sim, 4),
            "slope_meas_Cph": round(slope_meas, 4),
            "slope_abs_err_Cph": round(slope_err, 4),
            "G14_hourly_pass": mh["G14_pass"],
            "accept_triple_gate": bool(accept),
            "n_hourly": mh["n"],
        })
    res = pd.DataFrame(rows).sort_values("CVRMSE_hourly").reset_index(drop=True)
    res.to_csv(results_csv, index=False)
    return res
# =====================================================================
# patch_normalize_wall3.py
# Replacement for normalize_wall3_to_outdoors() in calibration_engine.py
# Fixes the fatal 'is not an OtherSideConditionedModel surface' error by
# removing the FULL cavity chain of Wall_3, not just the OSCM:
#
#   ExteriorNaturalVentedCavity --(Boundary Conditions Model Name)--> OSCM
#         |                                                            |
#         +--(Surface 1 Name)--> Wall_3 <--(Outside BC Object)---------+
#
# We (1) repoint Wall_3 to Outdoors, (2) drop the OSCM, (3) drop the
# ExteriorNaturalVentedCavity that referenced this surface. Real solar
# SHADING surfaces are NOT touched.
#
# This file is applied to calibration_engine.py via an R replacement
# script (apply_wall3_patch.R) — do NOT paste it by hand.
# =====================================================================

def normalize_wall3_to_outdoors(
        idf, surface=CLASSROOM_ZONE + "_Wall_3_0_0",
        oscm_name=CLASSROOM_ZONE + "_Wall_3_0_0 BoundaryObject"):
    """Make the 31 deg facade a plain exterior wall (sibling-consistent),
    removing the spurious ExteriorNaturalVentedCavity chain. E+ aborts if
    the cavity (or its OSCM) is left pointing at a surface that is no longer
    an OtherSideConditionedModel surface, so all three links are cleaned.
    Returns a report dict. Shading surfaces are preserved.
    """
    surf = next((s for s in idf.idfobjects["BuildingSurface:Detailed"]
                 if s.Name == surface), None)
    if surf is None:
        raise ValueError(f"Surface not found: {surface}")

    before = surf.Outside_Boundary_Condition
    surf.Outside_Boundary_Condition = "Outdoors"
    surf.Outside_Boundary_Condition_Object = ""
    surf.Sun_Exposure = "SunExposed"
    surf.Wind_Exposure = "WindExposed"
    surf.View_Factor_to_Ground = "autocalculate"

    dropped_cav = 0
    cav_key = "SurfaceProperty:ExteriorNaturalVentedCavity"
    for o in list(idf.idfobjects[cav_key]):
        # match by the surface it sits on OR the OSCM it couples through
        if (str(o.Surface_1_Name) == surface
                or str(o.Boundary_Conditions_Model_Name) == oscm_name):
            idf.removeidfobject(o)
            dropped_cav += 1

    dropped_oscm = 0
    for o in list(idf.idfobjects["SurfaceProperty:OtherSideConditionsModel"]):
        if o.Name == oscm_name:
            idf.removeidfobject(o)
            dropped_oscm += 1

    return {"wall3_obc_before": before, "wall3_obc_after": "Outdoors",
            "natventcavity_dropped": dropped_cav, "oscm_dropped": dropped_oscm}
