#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
from pathlib import Path
import pandas as pd

ERA = "2025G"

YIELDS_CSV = Path("../yield/results/2025G/yields.csv")
ACCEPTANCE_CSV = Path("../acceptance/results/acceptance.csv")
EFFICIENCY_CSV = Path("../efficiency/mc_efficiency/results/efficiency.csv")
LUMINOSITY_CSV = Path(f"../luminosity/results/{ERA}.csv")

RESULTS_DIR = Path("results")
OUT_CSV = RESULTS_DIR / "cross_section.csv"


def safe_div(a: float, b: float) -> float:
    return a / b if b != 0.0 else float("nan")


def main() -> int:
    yields = pd.read_csv(YIELDS_CSV)
    acc = pd.read_csv(ACCEPTANCE_CSV)
    eff = pd.read_csv(EFFICIENCY_CSV)
    lumi_df = pd.read_csv(LUMINOSITY_CSV)

    keys = ["pt_min", "pt_max", "y_abs_min", "y_abs_max"]
    for df_name, df in [("yields", yields), ("acceptance", acc), ("efficiency", eff)]:
        for c in keys:
            if c not in df.columns:
                raise KeyError(f"Missing column '{c}' in {df_name}")

    needed_y = [
        "N_1S", "N_1S_err",
        "N_2S", "N_2S_err",
        "N_3S", "N_3S_err",
    ]
    for c in needed_y:
        if c not in yields.columns:
            raise KeyError(f"Missing column '{c}' in yields.csv")

    if "acceptance" not in acc.columns or "acceptance_error" not in acc.columns:
        raise KeyError("Missing 'acceptance' or 'acceptance_error' in acceptance.csv")

    if "efficiency" not in eff.columns or "efficiency_error" not in eff.columns:
        raise KeyError("Missing 'efficiency' or 'efficiency_error' in efficiency.csv")

    if "lumi" not in lumi_df.columns:
        raise KeyError("Missing 'lumi' in luminosity csv")

    lumi_fb = float(lumi_df["lumi"].sum())
    lumi_pb = lumi_fb * 1000.0  # 1 fb = 1000 pb

    df = yields[keys + needed_y].merge(
        acc[keys + ["acceptance", "acceptance_error"]], on=keys, how="inner"
    ).merge(
        eff[keys + ["efficiency", "efficiency_error"]], on=keys, how="inner"
    )

    if df.empty:
        raise RuntimeError("No matched bins after merging yields/acceptance/efficiency. Check bin definitions.")

    def calc_state(N: float, dN: float, A: float, dA: float, E: float, dE: float, denom: float) -> tuple[float, float, float, float]:
        cs = safe_div(N, denom)  # pb / GeV
        dcs_stat = abs(cs) * abs(safe_div(dN, N)) if (N != 0.0 and math.isfinite(cs)) else float("nan")
        dcs_syst = abs(cs) * math.sqrt((safe_div(dA, A) ** 2) + (safe_div(dE, E) ** 2)) if (A != 0.0 and E != 0.0 and math.isfinite(cs)) else float("nan")
        dcs_tot = math.sqrt(dcs_stat * dcs_stat + dcs_syst * dcs_syst) if (math.isfinite(dcs_stat) and math.isfinite(dcs_syst)) else float("nan")
        return cs, dcs_tot, dcs_stat, dcs_syst

    def bin_integral(cs_pb_per_gev: float, pt_min: float, pt_max: float, y_abs_min: float, y_abs_max: float) -> float:
        pt_w = pt_max - pt_min
        y_w = 2.0 * (y_abs_max - y_abs_min)
        return cs_pb_per_gev * pt_w * y_w

    out_rows = []
    for _, r in df.iterrows():
        pt_min = float(r["pt_min"])
        pt_max = float(r["pt_max"])
        y0 = float(r["y_abs_min"])
        y1 = float(r["y_abs_max"])

        pt_w = pt_max - pt_min
        y_w = 2.0 * (y1 - y0)  # abs(y) bin width -> full rapidity width

        A = float(r["acceptance"])
        dA = float(r["acceptance_error"])
        E = float(r["efficiency"])
        dE = float(r["efficiency_error"])

        denom = E * A * lumi_pb * pt_w * y_w

        cs1, e1, s1, y1s = calc_state(float(r["N_1S"]), float(r["N_1S_err"]), A, dA, E, dE, denom)
        cs2, e2, s2, y2s = calc_state(float(r["N_2S"]), float(r["N_2S_err"]), A, dA, E, dE, denom)
        cs3, e3, s3, y3s = calc_state(float(r["N_3S"]), float(r["N_3S_err"]), A, dA, E, dE, denom)

        out_rows.append(
            {
                "pt_min": pt_min,
                "pt_max": pt_max,
                "y_abs_min": y0,
                "y_abs_max": y1,

                "1S_cross_section": cs1,
                "1S_cross_section_error": abs(e1) if math.isfinite(e1) else e1,
                "1S_cross_section_stat_error": abs(s1) if math.isfinite(s1) else s1,
                "1S_cross_section_syst_error": abs(y1s) if math.isfinite(y1s) else y1s,

                "2S_cross_section": cs2,
                "2S_cross_section_error": abs(e2) if math.isfinite(e2) else e2,
                "2S_cross_section_stat_error": abs(s2) if math.isfinite(s2) else s2,
                "2S_cross_section_syst_error": abs(y2s) if math.isfinite(y2s) else y2s,

                "3S_cross_section": cs3,
                "3S_cross_section_error": abs(e3) if math.isfinite(e3) else e3,
                "3S_cross_section_stat_error": abs(s3) if math.isfinite(s3) else s3,
                "3S_cross_section_syst_error": abs(y3s) if math.isfinite(y3s) else y3s,
            }
        )

    out_cols = [
        "pt_min", "pt_max", "y_abs_min", "y_abs_max",
        "1S_cross_section", "1S_cross_section_error", "1S_cross_section_stat_error", "1S_cross_section_syst_error",
        "2S_cross_section", "2S_cross_section_error", "2S_cross_section_stat_error", "2S_cross_section_syst_error",
        "3S_cross_section", "3S_cross_section_error", "3S_cross_section_stat_error", "3S_cross_section_syst_error",
    ]

    out = pd.DataFrame(out_rows)[out_cols].sort_values(["y_abs_min", "y_abs_max", "pt_min", "pt_max"]).reset_index(drop=True)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False, float_format="%.10g")

    out_int = out[(out["pt_max"] <= 130.0) & (out["y_abs_max"] <= 2.4)].copy()
    if out_int.empty:
        raise RuntimeError("No bins found within |y|<2.4 and pT<130 for integration.")

    def integrate_total(state: int) -> tuple[float, float, float, float]:
        cs = out_int[f"{state}S_cross_section"].astype(float).values
        et = out_int[f"{state}S_cross_section_error"].astype(float).values
        es = out_int[f"{state}S_cross_section_stat_error"].astype(float).values
        ey = out_int[f"{state}S_cross_section_syst_error"].astype(float).values

        pt_min = out_int["pt_min"].astype(float).values
        pt_max = out_int["pt_max"].astype(float).values
        y0 = out_int["y_abs_min"].astype(float).values
        y1 = out_int["y_abs_max"].astype(float).values

        pt_w = pt_max - pt_min
        y_w = 2.0 * (y1 - y0)
        w = pt_w * y_w

        sigma = float(np.nansum(cs * w))
        sigma_stat = float(math.sqrt(np.nansum((es * w) ** 2)))
        sigma_syst = float(math.sqrt(np.nansum((ey * w) ** 2)))
        sigma_tot = float(math.sqrt(np.nansum((et * w) ** 2)))
        return sigma, sigma_tot, sigma_stat, sigma_syst

    import numpy as np  # keep scope local impact minimal

    s1, s1_tot, s1_stat, s1_syst = integrate_total(1)
    s2, s2_tot, s2_stat, s2_syst = integrate_total(2)
    s3, s3_tot, s3_stat, s3_syst = integrate_total(3)

    print(f"[OK] ERA={ERA}, lumi={lumi_fb:.6f} fb = {lumi_pb:.6f} pb")
    print(f"[OK] Wrote: {OUT_CSV}")
    print("[INT] Integrated over |y|<2.4 and pT<130 (sigma*BR) in pb:")
    print(f"      1S: {s1:.6g} ± {s1_tot:.6g} (stat {s1_stat:.6g}, syst {s1_syst:.6g})")
    print(f"      2S: {s2:.6g} ± {s2_tot:.6g} (stat {s2_stat:.6g}, syst {s2_syst:.6g})")
    print(f"      3S: {s3:.6g} ± {s3_tot:.6g} (stat {s3_stat:.6g}, syst {s3_syst:.6g})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())