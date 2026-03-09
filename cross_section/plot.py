#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as hep
import matplotlib as mpl
from matplotlib.ticker import ScalarFormatter, NullLocator
from matplotlib.legend_handler import HandlerTuple


CROSS_2025_CSV = Path("results/cross_section.csv")
LUMI_2025_CSV = Path("../luminosity/results/2025G.csv")

OUT_PDF = Path("results/cross_section.pdf")

PT_MIN_PLOT = 0.0
PT_MAX_PLOT = 130.0

YBINS = [(0.0, 0.6), (0.6, 1.2), (1.2, 1.8), (1.8, 2.4)]
YBIN_LABELS = {
    (0.0, 0.6): r"$|\mathit{y}|<0.6$",
    (0.6, 1.2): r"$0.6<|\mathit{y}|<1.2$",
    (1.2, 1.8): r"$1.2<|\mathit{y}|<1.8$",
    (1.8, 2.4): r"$1.8<|\mathit{y}|<2.4$",
}

COLORS = {1: "#e42536", 2: "#3a981dff", 3: "#5790fc"}  # 1S, 2S, 3S
STATE_LABELS = {
    1: r"$\Upsilon(1S)\times100$",
    2: r"$\Upsilon(2S)\times10$",
    3: r"$\Upsilon(3S)$",
}
STATE_SCALE = {1: 100.0, 2: 10.0, 3: 1.0}

# marker per |y| bin: y0 filled o, y1 filled s, y2 open o, y3 open s
YBIN_MARKER = {
    (0.0, 0.6): ("o", False),
    (0.6, 1.2): ("s", False),
    (1.2, 1.8): ("^", False),
    (1.8, 2.4): ("v", False),
}

POINT_MS = 8
PROXY_MS = 10
MARKER_EDGE_W = 1.5
ELINE_W = 1.5
CAPSIZE = 0

AXLABEL_SIZE = 28
LEGEND_SIZE = 28
TICK_LABEL_SIZE = 28
MAJOR_TICK_LEN = 8.0

XTICKS = [0, 20, 40, 60, 80, 100]
YTICKS = [1e-4, 1e-2, 1e0, 1e2, 1e4]
YMIN, YMAX = 6e-5, 5e4

CMS_FONT = 30
COM_ENERGY = 13.6

YABS_MATCH_ATOL = 1e-8
YABS_MATCH_RTOL = 0.0


def read_lumi_fb(path: Path) -> float:
    df = pd.read_csv(path)
    if "lumi" not in df.columns:
        raise KeyError(f"Missing 'lumi' in {path}")
    return float(df["lumi"].sum())


def errorbar_proxy(ax, marker: str, color: str, filled: bool, ms: float = PROXY_MS):
    mfc = color if filled else "none"
    mec = color
    return ax.errorbar(
        [np.nan], [np.nan],
        yerr=[1], #xerr=[1],
        fmt=marker,
        markersize=ms,
        markerfacecolor=mfc,
        markeredgewidth=MARKER_EDGE_W,
        color=color,
        elinewidth=ELINE_W,
        capsize=CAPSIZE,
        linestyle="none",
    )


def main() -> int:
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

    lumi_2025_fb = read_lumi_fb(LUMI_2025_CSV)

    df = pd.read_csv(CROSS_2025_CSV)

    needed_base = ["pt_min", "pt_max", "y_abs_min", "y_abs_max"]
    for c in needed_base:
        if c not in df.columns:
            raise KeyError(f"Missing '{c}' in {CROSS_2025_CSV}")

    for s in (1, 2, 3):
        for suff in ["cross_section", "cross_section_error"]:
            c = f"{s}S_{suff}"
            if c not in df.columns:
                raise KeyError(f"Missing '{c}' in {CROSS_2025_CSV}")

    df["pt_min"] = df["pt_min"].astype(float)
    df["pt_max"] = df["pt_max"].astype(float)
    df["y_abs_min"] = df["y_abs_min"].astype(float)
    df["y_abs_max"] = df["y_abs_max"].astype(float)
    df["pt_center"] = 0.5 * (df["pt_min"] + df["pt_max"])

    df = df[(df["pt_center"] >= PT_MIN_PLOT) & (df["pt_center"] <= PT_MAX_PLOT)].copy()

    plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams["mathtext.rm"] = "serif"
    plt.style.use(hep.style.CMS)

    mpl.rcParams.update({
        "ytick.minor.visible": False,
        "xtick.minor.visible": False,
    })

    fig, ax = plt.subplots(1, 1, figsize=(14, 12))

    ax.set_yscale("log")
    ax.set_xlim(PT_MIN_PLOT, PT_MAX_PLOT)
    ax.set_ylim(YMIN, YMAX)

    xfmt = ScalarFormatter()
    xfmt.set_scientific(False)
    xfmt.set_useOffset(False)
    ax.xaxis.set_major_formatter(xfmt)

    ax.xaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_minor_locator(NullLocator())

    ax.tick_params(axis="both", which="both", top=False, right=False)
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_SIZE, pad=10, length=MAJOR_TICK_LEN)

    ax.set_xticks(XTICKS)
    ax.set_yticks(YTICKS)
    ax.set_xticklabels([str(int(x)) for x in XTICKS], fontsize=TICK_LABEL_SIZE)
    ax.set_yticklabels([rf"$10^{{{int(math.log10(y))}}}$" for y in YTICKS], fontsize=TICK_LABEL_SIZE)

    for (y0, y1) in YBINS:
        sel = df[
            np.isclose(df["y_abs_min"].values, y0, atol=YABS_MATCH_ATOL, rtol=YABS_MATCH_RTOL)
            & np.isclose(df["y_abs_max"].values, y1, atol=YABS_MATCH_ATOL, rtol=YABS_MATCH_RTOL)
        ].copy()
        if sel.empty:
            continue
        sel = sel.sort_values("pt_center")

        marker, filled = YBIN_MARKER[(y0, y1)]

        for state in (1, 2, 3):
            col = COLORS[state]
            scale = STATE_SCALE[state]

            y = sel[f"{state}S_cross_section"].astype(float).values * scale
            yerr = sel[f"{state}S_cross_section_error"].astype(float).abs().values * scale
            x = sel["pt_center"].astype(float).values

            ax.errorbar(
                x, y, yerr=yerr,
                fmt=marker,
                markersize=POINT_MS,
                markerfacecolor=(col if filled else "none"),
                markeredgewidth=MARKER_EDGE_W,
                color=col,
                linewidth=ELINE_W,
                elinewidth=ELINE_W,
                capsize=CAPSIZE,
                linestyle="none",
            )

    ax.set_ylabel(
        r"$\mathit{B}\,\frac{d^{2}\sigma}{d\mathit{p}_{T}d\mathit{y}}\;\mathrm{(pb \, / \, GeV)}$",
        size=34
    )
    #ax.yaxis.set_label_coords(-0.05, 1.0)
    #ax.xaxis.set_label_coords(0.9, -0.07)
    ax.set_xlabel(r"$\mathit{p}_{T}\;\mathrm{(GeV)}$", size=34)

    colors = [COLORS[1], COLORS[2], COLORS[3]]

    state_labels = [STATE_LABELS[1], STATE_LABELS[2], STATE_LABELS[3]]
    state_handles = []
    for c in colors:
        tup = tuple(
            errorbar_proxy(ax, YBIN_MARKER[ybin][0], c, filled=YBIN_MARKER[ybin][1])
            for ybin in YBINS
        )
        state_handles.append(tup)

    rap_labels = [YBIN_LABELS[ybin] for ybin in YBINS]
    rap_handles = [
        tuple(errorbar_proxy(ax, YBIN_MARKER[ybin][0], c, filled=YBIN_MARKER[ybin][1]) for c in colors)
        for ybin in YBINS
    ]

    leg_states = ax.legend(
        state_handles, state_labels,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        loc="upper right", bbox_to_anchor=(0.98, 1),
        frameon=False, fontsize=LEGEND_SIZE, borderpad=0.3, handlelength=2.8,
        handletextpad=0.8, markerscale=1.5
    )

    leg_rap = ax.legend(
        rap_handles, rap_labels,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        loc="upper right", bbox_to_anchor=(0.63, 1),
        frameon=False, fontsize=LEGEND_SIZE, borderpad=0.3, handlelength=2.8,
        handletextpad=0.8, markerscale=1.5
    )

    ax.add_artist(leg_states)
    ax.add_artist(leg_rap)

    hep.cms.label(
        ax=ax, data=True,
        com=COM_ENERGY,
        fontsize=CMS_FONT,
        lumi=f"{lumi_2025_fb:.1f}"
    )

    plt.tight_layout()
    plt.savefig(OUT_PDF)
    plt.close(fig)
    print(f"[OK] Wrote: {OUT_PDF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())