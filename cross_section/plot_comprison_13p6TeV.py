#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.legend_handler import HandlerTuple
from matplotlib.ticker import ScalarFormatter, NullLocator

CROSS_2025_CSV = Path("results/cross_section.csv")
LUMI_2025_CSV = Path("../luminosity/results/2025G.csv")

DATA2022_DIR = Path("/afs/ihep.ac.cn/users/y/yiyangzhao/Research/CMS_Run3_Upsilon/CMS_Plot")
DATA2022_FILES = {
    (1, 1): DATA2022_DIR / "data11.txt",
    (1, 2): DATA2022_DIR / "data12.txt",
    (2, 1): DATA2022_DIR / "data21.txt",
    (2, 2): DATA2022_DIR / "data22.txt",
    (3, 1): DATA2022_DIR / "data31.txt",
    (3, 2): DATA2022_DIR / "data32.txt",
}

OUT_PDF = Path("results/2025vs2022.pdf")

PT_MIN_PLOT = 20.0
PT_MAX_PLOT = 130.0

COLORS = {1: "#e42536", 2: "#3a981dff", 3: "#5790fc"}
STATE_LABELS = {1: r"$\Upsilon(1S)$", 2: r"$\Upsilon(2S)$", 3: r"$\Upsilon(3S)$"}

YBIN_TO_RAPIDX = {(0.0, 0.6): 1, (0.6, 1.2): 2}
RAP_LABELS = {1: r"$|\mathit{y}|<0.6$", 2: r"$0.6<|\mathit{y}|<1.2$"}

LOG10_OFFSET = 0.006  # half of previous (was 0.012)

CENTER_MATCH_ATOL = 1e-6
CENTER_MATCH_RTOL = 0.0

POINT_MS = 9.0
PROXY_MS = 9.0

TICK_LABEL_SIZE = 26
AXLABEL_SIZE = 28
LEGEND_SIZE = 28
YLABEL_SIZE = 28

XTICKS = [20.0, 40.0, 60.0, 80.0, 100.0]
XTICKLABELS = ["20", "40", "60", "80", "100"]
MAJOR_TICK_LEN = 8.0


def errorbar_proxy(ax, marker: str, color: str, filled: bool):
    mfc = color if filled else "none"
    mec = color
    return ax.errorbar(
        [np.nan], [np.nan],
        yerr=[1.0],
        fmt=marker,
        ms=PROXY_MS,
        mfc=mfc,
        mec=mec,
        mew=2.0,
        ecolor=color,
        elinewidth=2.0,
        capsize=0,
        linestyle="none",
    )


def read_2022_file(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+", header=None, names=["pt_center", "cs", "err"], engine="python")
    df["pt_center"] = df["pt_center"].astype(float)
    df["cs"] = df["cs"].astype(float)
    df["err"] = df["err"].astype(float).abs()
    return df.sort_values("pt_center").reset_index(drop=True)


def apply_axis_style(ax):
    ax.set_xscale("log")
    ax.set_xlim(PT_MIN_PLOT, PT_MAX_PLOT)
    ax.set_ylim(0.3, 1.15)

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=2.0, zorder=1)

    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    xfmt = ScalarFormatter(useOffset=False)
    xfmt.set_scientific(False)
    ax.xaxis.set_major_formatter(xfmt)

    yfmt = ScalarFormatter(useOffset=False)
    yfmt.set_scientific(False)
    ax.yaxis.set_major_formatter(yfmt)

    ax.tick_params(
        axis="both", which="major",
        labelsize=TICK_LABEL_SIZE,
        top=False, right=False,
        length=MAJOR_TICK_LEN
    )

def read_lumi_fb(path: Path) -> float:
    df = pd.read_csv(path)
    if "lumi" not in df.columns:
        raise KeyError(f"Missing 'lumi' in {path}")
    return float(df["lumi"].sum())

def plot_panel(ax, df25: pd.DataFrame, y0: float, y1: float, rapidx: int):
    sel = df25[np.isclose(df25["y_abs_min"], y0) & np.isclose(df25["y_abs_max"], y1)].copy()
    if sel.empty:
        return

    for state in (1, 2, 3):
        f22 = DATA2022_FILES[(state, rapidx)]
        df22 = read_2022_file(f22)

        cs25_col = f"{state}S_cross_section"
        e25_col = f"{state}S_cross_section_error"

        x_all, y_all, yerr_all = [], [], []
        for _, r22 in df22.iterrows():
            x = float(r22["pt_center"])
            if x < PT_MIN_PLOT or x > PT_MAX_PLOT:
                continue

            m = sel[np.isclose(sel["pt_center"].values, x, atol=CENTER_MATCH_ATOL, rtol=CENTER_MATCH_RTOL)]
            if len(m) != 1:
                continue

            r25 = m.iloc[0]
            cs25 = float(r25[cs25_col])
            e25 = float(r25[e25_col])
            cs22 = float(r22["cs"])
            e22 = float(r22["err"])

            if not (math.isfinite(cs25) and math.isfinite(e25) and math.isfinite(cs22) and math.isfinite(e22)):
                continue
            if cs22 == 0.0 or cs25 == 0.0:
                continue

            ratio = cs25 / cs22
            ratio_err = abs(ratio) * math.sqrt((e25 / cs25) ** 2 + (e22 / cs22) ** 2)

            x_all.append(x)
            y_all.append(ratio)
            yerr_all.append(abs(ratio_err))

        if not x_all:
            continue

        x_all = np.array(x_all, dtype=float)
        y_all = np.array(y_all, dtype=float)
        yerr_all = np.array(yerr_all, dtype=float)

        if state == 1:
            x_plot = x_all
        elif state == 2:
            x_plot = x_all * (10.0 ** (-LOG10_OFFSET))
        else:
            x_plot = x_all * (10.0 ** (LOG10_OFFSET))

        color = COLORS[state]
        if rapidx == 1:
            marker = "o"
            mfc = "none"
            mec = color
        else:
            marker = "s"
            mfc = "none"
            mec = color

        ax.errorbar(
            x_plot, y_all,
            yerr=yerr_all,
            fmt=marker,
            linestyle="none",
            ms=POINT_MS,
            mfc=mfc,
            mec=mec,
            mew=2.0,
            ecolor=color,
            elinewidth=2.0,
            capsize=0,
            zorder=3,
        )


def main() -> int:
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

    df25 = pd.read_csv(CROSS_2025_CSV)

    for c in ["pt_min", "pt_max", "y_abs_min", "y_abs_max"]:
        if c not in df25.columns:
            raise KeyError(f"Missing '{c}' in {CROSS_2025_CSV}")

    for s in (1, 2, 3):
        for suffix in ["cross_section", "cross_section_error"]:
            col = f"{s}S_{suffix}"
            if col not in df25.columns:
                raise KeyError(f"Missing '{col}' in {CROSS_2025_CSV}")

    df25["pt_center"] = 0.5 * (df25["pt_min"].astype(float) + df25["pt_max"].astype(float))
    df25["y_abs_min"] = df25["y_abs_min"].astype(float)
    df25["y_abs_max"] = df25["y_abs_max"].astype(float)
    df25 = df25[(df25["pt_center"] >= PT_MIN_PLOT) & (df25["pt_center"] <= PT_MAX_PLOT)].copy()

    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["mathtext.rm"] = "serif"
    plt.style.use(hep.style.CMS)

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(12, 12), sharex=True,
        gridspec_kw={"hspace": 0.0, "height_ratios": [1, 1]}
    )

    lumi_2025_fb = read_lumi_fb(LUMI_2025_CSV)
    hep.cms.label(ax=ax_top, data=True, com=13.6, fontsize=28, label="Preliminary")

    apply_axis_style(ax_top)
    apply_axis_style(ax_bot)

    ax_top.set_ylabel("")
    ax_bot.set_ylabel("")

    ax_top.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
    ax_top.set_xlabel("")

    ax_bot.tick_params(axis="x", which="both", top=False)
    ax_bot.set_xlabel(r"$\mathit{p}_{T}\;\mathrm{(GeV)}$", fontsize=AXLABEL_SIZE)

    ax_bot.set_xticks(XTICKS)
    ax_bot.set_xticklabels(XTICKLABELS, fontsize=TICK_LABEL_SIZE)
    ax_top.set_xticks(XTICKS)

    plot_panel(ax_top, df25, 0.0, 0.6, rapidx=1)
    plot_panel(ax_bot, df25, 0.6, 1.2, rapidx=2)

    fig.text(0.005, 0.5, r"$2025\;/\;2022$", va="center", rotation="vertical", fontsize=YLABEL_SIZE)

    colors = [COLORS[1], COLORS[2], COLORS[3]]

    state_labels = [STATE_LABELS[1], STATE_LABELS[2], STATE_LABELS[3]]
    state_handles = []
    for c in colors:
        tup = (errorbar_proxy(ax_top, "o", c, filled=False),
               errorbar_proxy(ax_top, "s", c, filled=False))
        state_handles.append(tup)

    rap_labels = [RAP_LABELS[1], RAP_LABELS[2]]
    rap_handles = [
        tuple(errorbar_proxy(ax_bot, "o", c, filled=False) for c in colors),
        tuple(errorbar_proxy(ax_bot, "s", c, filled=False) for c in colors),
    ]

    leg_states = ax_top.legend(
        state_handles, state_labels,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        frameon=False, fontsize=LEGEND_SIZE, borderpad=0.3, handlelength=2.8,
        handletextpad=0.8, markerscale=1.4
    )

    leg_rap = ax_bot.legend(
        rap_handles, rap_labels,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        loc="lower left", bbox_to_anchor=(0.02, 0.02),
        frameon=False, fontsize=LEGEND_SIZE, borderpad=0.3, handlelength=2.8,
        handletextpad=0.8, markerscale=1.4
    )

    fig.tight_layout()
    fig.savefig(OUT_PDF)
    print(f"[OK] Wrote: {OUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())