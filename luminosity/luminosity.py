#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

JSON_FILENAME = "Cert_Collisions2025_391658_398903_Muon.json"
ERAS: List[str] = ["2025G"]
RESULTS_DIR = Path("results")

ERA_RUN_RANGES: Dict[str, Tuple[int, int]] = {
    "2025Av1": (390735, 390936),
    "2025Av2": (390937, 391530),
    "2025B":   (391531, 392158),
    "2025Cv1": (392159, 393108),
    "2025Cv2": (393111, 393609),
    "2025D":   (394286, 395967),
    "2025E":   (395968, 396597),
    "2025F":   (396598, 397853),
    "2025G":   (397854, 398903),
}

BRILCALC_EXE = "brilcalc"
BRILCALC_CONNECT = "web"
NORMTAG: Optional[str] = None
LUMI_UNIT = "/fb"

@dataclass(frozen=True)
class RunRange:
    start: int
    end: int


def ensure_brilcalc_exists(exe: str) -> None:
    if shutil.which(exe) is None:
        raise RuntimeError(
            f'Cannot find "{exe}" in PATH. '
            f"Please set up brilcalc (CMS environment) or set BRILCALC_EXE to the full path."
        )


def load_cert_json(path: Path) -> Dict[str, List[List[int]]]:
    if not path.exists():
        raise FileNotFoundError(f"Certification JSON not found: {path}")
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError("Certification JSON root must be an object/dict.")
    return data


def validate_eras(eras: List[str], mapping: Dict[str, Tuple[int, int]]) -> List[str]:
    missing = [e for e in eras if e not in mapping]
    if missing:
        raise ValueError(f"Unknown era(s): {missing}. Known eras: {sorted(mapping.keys())}")
    return eras


def combined_run_range(eras: List[str], mapping: Dict[str, Tuple[int, int]]) -> RunRange:
    starts = []
    ends = []
    for e in eras:
        s, t = mapping[e]
        starts.append(s)
        ends.append(t)
    return RunRange(start=min(starts), end=max(ends))


def run_brilcalc_lumi(
    exe: str,
    connect: str,
    cert_json: Path,
    run_range: RunRange,
    unit: str,
    normtag: Optional[str] = None,
) -> str:
    cmd = [
        exe, "lumi",
        "-c", connect,
        "--begin", str(run_range.start),
        "--end", str(run_range.end),
        "-i", str(cert_json),
        "-u", unit,
        "--output-style", "csv",
    ]
    if normtag:
        cmd.extend(["--normtag", str(normtag)])

    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "brilcalc failed.\n"
            f"Command: {' '.join(cmd)}\n"
            f"Return code: {proc.returncode}\n"
            f"STDERR:\n{proc.stderr}\n"
            f"STDOUT:\n{proc.stdout}\n"
        )
    return proc.stdout


def parse_brilcalc_csv_total_recorded(stdout_csv: str) -> float:
    lines = [ln.strip() for ln in stdout_csv.splitlines() if ln.strip()]
    if not lines:
        return 0.0

    header_line = None
    for ln in lines:
        if ln.startswith("#") and "recorded" in ln:
            header_line = ln.lstrip("#").strip()
            break
    if header_line is None:
        raise ValueError("Could not find a CSV header line containing 'recorded' in brilcalc output.")

    header = [h.strip() for h in header_line.split(",")]
    rec_idx = None
    for i, col in enumerate(header):
        if col.startswith("recorded(") or col == "recorded":
            rec_idx = i
            break
    if rec_idx is None:
        for i, col in enumerate(header):
            if "recorded" in col:
                rec_idx = i
                break
    if rec_idx is None:
        raise ValueError(f"Could not locate 'recorded' column in header: {header}")

    total = 0.0
    for ln in lines:
        if ln.startswith("#"):
            continue
        row = [x.strip() for x in ln.split(",")]
        if len(row) <= rec_idx:
            continue
        val = row[rec_idx]
        if val == "" or val.lower() == "nan":
            continue
        try:
            total += float(val)
        except ValueError:
            continue

    return total


def write_result_csv(out_path: Path, run_range: RunRange, lumi_value: float) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["start_run", "end_run", "lumi"])
        w.writerow([run_range.start, run_range.end, f"{lumi_value:.6f}"])


def main() -> int:
    cert_path = Path(JSON_FILENAME).expanduser().resolve()

    ensure_brilcalc_exists(BRILCALC_EXE)
    _ = load_cert_json(cert_path)

    eras = validate_eras(ERAS, ERA_RUN_RANGES)
    rr = combined_run_range(eras, ERA_RUN_RANGES)

    stdout = run_brilcalc_lumi(
        exe=BRILCALC_EXE,
        connect=BRILCALC_CONNECT,
        cert_json=cert_path,
        run_range=rr,
        unit=LUMI_UNIT,
        normtag=NORMTAG,
    )
    total_recorded = parse_brilcalc_csv_total_recorded(stdout)

    out_name = "_".join(eras) + ".csv"
    out_path = RESULTS_DIR / out_name
    write_result_csv(out_path, rr, total_recorded)

    print(f"[OK] ERAS={eras}")
    print(f"[OK] Run range: {rr.start} - {rr.end}")
    print(f"[OK] Integrated recorded lumi: {total_recorded:.6f} ({LUMI_UNIT})")
    print(f"[OK] Wrote: {out_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        raise