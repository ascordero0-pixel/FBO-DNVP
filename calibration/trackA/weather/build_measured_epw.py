# -*- coding: utf-8 -*-
"""
build_measured_epw.py
=====================================================================
Reproducible builder of a MEASURED-data EPW for the T3 calibration of
the Seville classroom digital twin.

What it does (fully documented, replicable by reviewers):
  1. Reads measured 5-min meteorology from the Track A workbook
     (columns J,K,L,M = outdoor T, RH, wind dir, wind speed;
      column O = solar radiation already ÷10-corrected;
      column Q = atmospheric pressure, MSL-reduced, in hPa).
  2. Cleans sensor dropouts (zeros) and aggregates 5-min → hourly.
  3. Derives dew point (Magnus) from T and RH.
  4. Corrects pressure: hPa → Pa and MSL → station (z = 18 m).
  5. Decomposes measured Global Horizontal Irradiance into Direct
     Normal (DNI) and Diffuse Horizontal (DHI) via the Erbs model
     (pvlib), using the true solar geometry of Seville for 2021.
     >>> This is required: EnergyPlus reconstructs solar gains from
         DNI + DHI, NOT from the GHI field. Feeding only GHI yields a
         sunless model. <<<
  6. Splices the measured block (1 Nov – 21 Dec 2021) into the IWEC
     Seville template, leaving the rest of the year untouched (it is
     not simulated; the IDF RunPeriod is Nov 1 – Dec 20).
  7. Writes the new EPW and a transformation report.

Provenance policy (no-fabricated-data): every field is either measured,
or derived from measured values by a named, citable transformation.
Nothing is invented. Original "wip" columns are preserved in the source
workbook for audit.

Author: A. Sánchez Cordero (ASC), TEP192, Universidad de Huelva, 2026.
Requires: pandas, numpy, pvlib.
=====================================================================
"""
from __future__ import annotations

import os
import json
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd

# --------------------------------------------------------------------
# Site metadata (Seville classroom). Latitude/longitude from the IWEC
# header and T1; elevations supplied by ASC.
# --------------------------------------------------------------------
@dataclass
class Site:
    latitude: float = 37.42       # IWEC header / T1
    longitude: float = -5.90      # IWEC header / T1
    tz_utc_offset: float = 1.0    # UTC+1 (winter, standard time)
    elev_building_m: float = 12.0 # building ground level → EPW Elevation
    elev_station_m: float = 18.0  # met-station barometer height
    year: int = 2021              # measurement year


# EPW field layout (0-based index after comma split). 35 fields total.
EPW_FIELDS = [
    "Year", "Month", "Day", "Hour", "Minute", "DataSourceFlags",
    "DryBulb", "DewPoint", "RelHum", "AtmPressure",
    "ExtHorizRad", "ExtDirNormRad", "HorizIRSky",
    "GlobHorizRad", "DirNormRad", "DiffHorizRad",
    "GlobHorizIllum", "DirNormIllum", "DiffHorizIllum", "ZenithLum",
    "WindDir", "WindSpeed", "TotalSkyCover", "OpaqueSkyCover",
    "Visibility", "CeilingHeight", "PresWeathObs", "PresWeathCodes",
    "PrecipitableWater", "AerosolOpticalDepth", "SnowDepth",
    "DaysSinceLastSnow", "Albedo", "LiquidPrecipDepth", "LiquidPrecipRate",
]
IDX = {name: i for i, name in enumerate(EPW_FIELDS)}


# --------------------------------------------------------------------
# 1. Read and clean the measured 5-min data
# --------------------------------------------------------------------
def read_measured(xlsx_path: str, site: Site) -> pd.DataFrame:
    """Return a clean 5-min DataFrame indexed by tz-aware datetime."""
    raw = pd.read_excel(xlsx_path, sheet_name=0, header=None, skiprows=6)
    # Column letters → 0-based positions (B=1 ... V=21). See header Row6.
    colmap = {
        "Year": 1, "Month": 2, "Day": 3, "Time": 4,
        "To": 9,        # J  outdoor dry-bulb (°C)
        "Ho": 10,       # K  outdoor RH (%)
        "WindD": 11,    # L  wind direction (°)
        "WindS": 12,    # M  wind speed (m/s)
        "Rsolar": 14,   # O  GHI corrected (W/m²)
        "Pressure": 16, # Q  pressure MSL-reduced (hPa)
    }
    df = pd.DataFrame({k: raw[v] for k, v in colmap.items()})
    df["Year"] = df["Year"].ffill().astype(int)
    df["Month"] = df["Month"].ffill().astype(int)
    df["Day"] = pd.to_numeric(df["Day"], errors="coerce").ffill().astype(int)

    # Build timestamp from Year/Month/Day + Time (datetime.time or string)
    times = pd.to_datetime(df["Time"].astype(str), errors="coerce")
    df["hour"] = times.dt.hour.fillna(0).astype(int)
    df["minute"] = times.dt.minute.fillna(0).astype(int)
    ts = pd.to_datetime(dict(year=df["Year"], month=df["Month"], day=df["Day"],
                             hour=df["hour"], minute=df["minute"]))
    df = df.set_index(ts)

    for c in ["To", "Ho", "WindD", "WindS", "Rsolar", "Pressure"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Sensor dropouts encoded as 0 in pressure (physically impossible) and
    # occasional 0 in wind/RH. Pressure: replace 0 → NaN then interpolate.
    df.loc[df["Pressure"] <= 900, "Pressure"] = np.nan  # <900 hPa impossible here
    df["Pressure"] = df["Pressure"].interpolate(limit_direction="both")
    # RH must be within (0,100]; 0 is a dropout for outdoor winter air.
    df.loc[df["Ho"] <= 0, "Ho"] = np.nan
    df["Ho"] = df["Ho"].interpolate(limit_direction="both").clip(1, 100)
    # Wind speed 0 is physically valid (calm) → keep. Direction 0 valid (N).
    # Solar 0 valid at night → keep.
    return df[["To", "Ho", "WindD", "WindS", "Rsolar", "Pressure"]]


# --------------------------------------------------------------------
# 2. Aggregate 5-min → hourly (EPW convention: hour h = interval ending h:00)
# --------------------------------------------------------------------
def aggregate_hourly(df5: pd.DataFrame) -> pd.DataFrame:
    """
    EPW hour H represents the interval (H-1:00 , H:00]. We therefore
    label each hourly mean by its RIGHT edge. Wind direction is averaged
    vectorially to avoid the 359°/1° wrap artefact.
    """
    # Vector mean for wind direction
    wd = np.deg2rad(df5["WindD"])
    u = np.cos(wd); v = np.sin(wd)
    base = df5.drop(columns=["WindD"]).copy()
    base["u"] = u; base["v"] = v

    hourly = base.resample("1h", label="right", closed="right").mean()
    # Recover direction from averaged components
    ang = np.rad2deg(np.arctan2(hourly["v"], hourly["u"])) % 360.0
    hourly["WindD"] = ang
    hourly = hourly.drop(columns=["u", "v"])
    hourly = hourly.dropna(subset=["To", "Ho", "WindS", "Pressure"])
    return hourly


# --------------------------------------------------------------------
# 3. Dew point from T and RH (Magnus formula)
# --------------------------------------------------------------------
def dew_point_magnus(T: np.ndarray, RH: np.ndarray) -> np.ndarray:
    """Magnus-Tetens, valid for the temperature range of interest."""
    a, b = 17.62, 243.12  # °C constants (WMO)
    gamma = np.log(RH / 100.0) + (a * T) / (b + T)
    return (b * gamma) / (a - gamma)


# --------------------------------------------------------------------
# 4. Pressure correction: hPa → Pa, MSL → station
# --------------------------------------------------------------------
def correct_pressure(p_msl_hpa: np.ndarray, station_h_m: float) -> np.ndarray:
    """
    Convert MSL-reduced pressure (hPa) to station pressure (Pa) at the
    barometer height using the standard barometric formula.
    """
    factor = (1 - 0.0065 * station_h_m / 288.15) ** 5.255
    return p_msl_hpa * factor * 100.0  # hPa → Pa


# --------------------------------------------------------------------
# 5. Erbs decomposition: GHI → DNI + DHI
# --------------------------------------------------------------------
def decompose_ghi(hourly: pd.DataFrame, site: Site) -> pd.DataFrame:
    """
    Use pvlib solar position + Erbs model. Timestamps are localized to the
    site's fixed UTC offset (standard winter time, no DST in the data).
    Returns hourly frame with added DNI, DHI columns (W/m²).
    """
    import pvlib
    # Localize index to a fixed offset matching the data's clock.
    tz = f"Etc/GMT{'-' if site.tz_utc_offset >= 0 else '+'}{abs(int(site.tz_utc_offset))}"
    idx = hourly.index.tz_localize(tz)
    solpos = pvlib.solarposition.get_solarposition(
        idx, site.latitude, site.longitude,
        altitude=site.elev_building_m)
    ghi = hourly["Rsolar"].clip(lower=0).values
    zenith = solpos["apparent_zenith"].values
    doy = idx.dayofyear
    erbs = pvlib.irradiance.erbs(ghi, zenith, doy)
    out = hourly.copy()
    out["DNI"] = np.nan_to_num(erbs["dni"], nan=0.0)
    out["DHI"] = np.nan_to_num(erbs["dhi"], nan=0.0)
    out["GHI"] = ghi
    # Guard: zero out any spurious irradiance when sun is below horizon.
    night = zenith >= 90.0
    out.loc[night, ["DNI", "DHI", "GHI"]] = 0.0
    return out


# --------------------------------------------------------------------
# 6. Splice measured block into the IWEC template
# --------------------------------------------------------------------
def build_epw(measured_hourly: pd.DataFrame, template_epw: str,
              out_epw: str, site: Site,
              start=(11, 1), end=(12, 21)) -> dict:
    """
    Replace template rows in [start, end] (inclusive) with measured values.
    Untouched rows keep IWEC values. Returns a provenance dict.
    """
    with open(template_epw, encoding="latin-1") as f:
        lines = f.readlines()
    header, data = lines[:8], lines[8:]

    # Index measured frame by (month, day, hour) with EPW hour 1..24.
    mh = measured_hourly.copy()
    mh["epw_hour"] = mh.index.hour
    # EPW uses hour 1..24; our right-labelled means put 00:00 as next day's
    # hour 24 of the previous day. Remap: hour 0 → 24 of previous day.
    epw_key = {}
    for ts, row in mh.iterrows():
        h = ts.hour
        m, d = ts.month, ts.day
        if h == 0:
            # belongs to hour 24 of the previous calendar day
            prev = ts - pd.Timedelta(hours=1)
            m, d, h = prev.month, prev.day, 24
        epw_key[(m, d, h)] = row

    def in_window(mo, dy):
        return (start <= (mo, dy) <= end)

    replaced, missing = 0, 0
    for i, ln in enumerate(data):
        p = ln.rstrip("\n").split(",")
        mo, dy, hr = int(p[1]), int(p[2]), int(p[3])
        if not in_window(mo, dy):
            continue
        key = (mo, dy, hr)
        if key not in epw_key:
            missing += 1
            continue
        r = epw_key[key]
        p[IDX["Year"]] = str(site.year)
        p[IDX["DryBulb"]] = f"{r['To']:.1f}"
        p[IDX["DewPoint"]] = f"{dew_point_magnus(r['To'], r['Ho']):.1f}"
        p[IDX["RelHum"]] = f"{int(round(r['Ho']))}"
        p[IDX["AtmPressure"]] = f"{int(round(r['P_station_Pa']))}"
        p[IDX["GlobHorizRad"]] = f"{int(round(r['GHI']))}"
        p[IDX["DirNormRad"]] = f"{int(round(r['DNI']))}"
        p[IDX["DiffHorizRad"]] = f"{int(round(r['DHI']))}"
        p[IDX["WindDir"]] = f"{int(round(r['WindD']))}"
        p[IDX["WindSpeed"]] = f"{r['WindS']:.1f}"
        p[IDX["DataSourceFlags"]] = "MEASURED_ASC_2021"
        data[i] = ",".join(p) + "\n"
        replaced += 1

    with open(out_epw, "w", encoding="latin-1") as f:
        f.writelines(header + data)

    return {"rows_replaced": replaced, "rows_missing_in_measured": missing,
            "window": {"start": start, "end": end}, "year": site.year}


# --------------------------------------------------------------------
# 7. Orchestration
# --------------------------------------------------------------------
def generate(xlsx_path: str, template_epw: str, out_epw: str,
             report_path: str, site: Optional[Site] = None) -> dict:
    site = site or Site()
    df5 = read_measured(xlsx_path, site)
    hourly = aggregate_hourly(df5)
    hourly["P_station_Pa"] = correct_pressure(
        hourly["Pressure"].values, site.elev_station_m)
    hourly = decompose_ghi(hourly, site)
    prov = build_epw(hourly, template_epw, out_epw, site)

    report = {
        "site": site.__dict__,
        "splice": prov,
        "transformations": {
            "dry_bulb": "col J, 5-min mean → hourly, direct",
            "rel_hum": "col K, dropouts→interp, clip[1,100], hourly mean",
            "wind_dir": "col L, vector (u,v) hourly mean to avoid wrap",
            "wind_speed": "col M, hourly mean",
            "pressure": ("col Q (MSL hPa) → zeros interpolated → "
                         f"×factor(z={site.elev_station_m} m)→station → ×100 → Pa"),
            "dew_point": "derived from T,RH via Magnus (WMO a=17.62,b=243.12)",
            "solar": ("col O (GHI, ÷10-corrected) → Erbs decomposition "
                      "(pvlib) → DNI+DHI; REQUIRED for EnergyPlus solar gains"),
            "hour_convention": "EPW hour h = interval ending h:00 (right edge)",
        },
        "policy": "All EPW fields measured or derived by named transforms; none invented.",
        "tools": {"pvlib": __import__("pvlib").__version__,
                  "pandas": pd.__version__, "numpy": np.__version__},
    }
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    return report


if __name__ == "__main__":
    import sys
    xlsx = sys.argv[1] if len(sys.argv) > 1 else "T3_Calibration_TRACK_A_50_nights.xlsx"
    tmpl = sys.argv[2] if len(sys.argv) > 2 else "ESP_Sevilla_083910_IWEC.epw"
    out = sys.argv[3] if len(sys.argv) > 3 else "Seville_measured_NovDec2021.epw"
    rep = "epw_build_report.json"
    r = generate(xlsx, tmpl, out, rep)
    print(json.dumps(r["splice"], indent=2))
    print("Report:", rep)
