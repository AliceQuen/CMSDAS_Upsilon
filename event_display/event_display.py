#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple


@dataclass(frozen=True)
class EventID:
    run: int
    lumi: int
    event: int
    fill: str = "N/A"

    def as_pick_line(self) -> str:
        return f"{self.run}:{self.lumi}:{self.event}"


def parse_event_line(line: str) -> EventID | None:
    text = line.strip()
    if not text or text.startswith("#"):
        return None
    fields = re.split(r"[\s,]+", text)
    core = fields[0]
    parts = core.split(":")
    if len(parts) != 3:
        raise ValueError(f"Invalid event entry: {line.rstrip()} (expected run:lumi:event)")
    run, lumi, event = (int(parts[0]), int(parts[1]), int(parts[2]))
    fill = fields[1] if len(fields) > 1 else "N/A"
    return EventID(run=run, lumi=lumi, event=event, fill=fill)


def read_event_list(path: Path) -> List[EventID]:
    events: List[EventID] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parsed = parse_event_line(line)
            if parsed is not None:
                events.append(parsed)
    return events


def unique_preserve_order(events: Sequence[EventID]) -> List[EventID]:
    seen = set()
    out: List[EventID] = []
    for evt in events:
        key = (evt.run, evt.lumi, evt.event)
        if key in seen:
            continue
        seen.add(key)
        out.append(evt)
    return out


def select_events_from_ntuple(ntuple_pattern: str, tree: str, max_events: int) -> List[EventID]:
    try:
        import numpy as np
        import uproot
    except ImportError as exc:
        raise RuntimeError(
            "Selecting events from ntuples requires uproot and numpy. "
            "Install them first in your Python environment."
        ) from exc

    files = sorted(glob.glob(ntuple_pattern))
    if not files and Path(ntuple_pattern).exists():
        files = [ntuple_pattern]
    if not files:
        raise FileNotFoundError(f"No files match ntuple pattern: {ntuple_pattern}")

    branches = ["run", "lumiblock", "event", "trigger", "charge", "vProb"]
    arrays = uproot.concatenate([f"{f}:{tree}" for f in files], branches, library="np")

    mask = np.ones_like(arrays["run"], dtype=bool)
    if "trigger" in arrays:
        mask &= arrays["trigger"] > 0
    if "charge" in arrays:
        mask &= arrays["charge"] == 0
    if "vProb" in arrays:
        mask &= arrays["vProb"] > 0.01

    runs = arrays["run"][mask]
    lumis = arrays["lumiblock"][mask]
    events = arrays["event"][mask]

    selected: List[EventID] = []
    seen = set()
    for run, lumi, evt in zip(runs, lumis, events):
        key = (int(run), int(lumi), int(evt))
        if key in seen:
            continue
        selected.append(EventID(run=key[0], lumi=key[1], event=key[2], fill="N/A"))
        seen.add(key)
        if len(selected) >= max_events:
            break

    if not selected:
        raise RuntimeError("No events selected from ntuple with the current baseline selection.")
    return selected


def write_event_files(events: Sequence[EventID], txt_path: Path, csv_path: Path) -> None:
    txt_path.parent.mkdir(parents=True, exist_ok=True)

    with txt_path.open("w", encoding="utf-8") as handle:
        for evt in events:
            handle.write(evt.as_pick_line() + "\n")

    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["run", "lumi", "event", "fill"])
        for evt in events:
            writer.writerow([evt.run, evt.lumi, evt.event, evt.fill])


def extract_copy_command(text: str) -> str:
    lines = text.splitlines()
    start = None
    for i, line in enumerate(lines):
        if line.strip().startswith("edmCopyPickMerge"):
            start = i
            break
    if start is None:
        raise RuntimeError(
            "Could not find edmCopyPickMerge command in edmPickEvents output. "
            "Please inspect stdout/stderr manually."
        )

    command_parts: List[str] = []
    i = start
    while i < len(lines):
        part = lines[i].strip()
        if not part:
            break
        command_parts.append(part.rstrip("\\").strip())
        if not lines[i].rstrip().endswith("\\"):
            break
        i += 1

    command = " ".join(command_parts)
    if not command.startswith("edmCopyPickMerge"):
        raise RuntimeError("Parsed command does not start with edmCopyPickMerge.")
    return command


def run_pick_commands(dataset: str, event_txt: Path, pick_root: Path) -> None:
    pick_cmd = ["edmPickEvents.py", dataset, str(event_txt)]
    pick_proc = subprocess.run(
        pick_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if pick_proc.returncode != 0:
        raise RuntimeError(
            "edmPickEvents.py failed\n"
            f"Command: {' '.join(pick_cmd)}\n"
            f"STDOUT:\n{pick_proc.stdout}\n"
            f"STDERR:\n{pick_proc.stderr}"
        )

    combined = pick_proc.stdout + "\n" + pick_proc.stderr
    copy_command = extract_copy_command(combined)

    if "outputFile=" in copy_command:
        copy_command = re.sub(r"outputFile=[^\s]+", f"outputFile={pick_root}", copy_command)
    else:
        copy_command += f" outputFile={pick_root}"

    copy_proc = subprocess.run(
        copy_command,
        shell=True,
        executable="/bin/bash",
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if copy_proc.returncode != 0:
        raise RuntimeError(
            "edmCopyPickMerge failed\n"
            f"Command: {copy_command}\n"
            f"STDOUT:\n{copy_proc.stdout}\n"
            f"STDERR:\n{copy_proc.stderr}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Create pickevents.root from selected events.")
    parser.add_argument("--dataset", required=True, help="CMS dataset path for edmPickEvents.py")
    parser.add_argument("--event-list", help="Existing text file with run:lumi:event entries")
    parser.add_argument("--ntuple", help="Ntuple file or glob pattern used to select events")
    parser.add_argument("--tree", default="rootuple/mm_tree", help="Tree path inside ntuple ROOT files")
    parser.add_argument("--max-events", type=int, default=5, help="Maximum number of events to keep")
    parser.add_argument("--output-dir", default="event_display/output", help="Output directory")
    parser.add_argument("--event-txt-name", default="event_display.txt", help="Name of run:lumi:event text file")
    parser.add_argument("--pick-root-name", default="pickevents.root", help="Output picked ROOT filename")
    return parser


def main() -> int:
    args = build_parser().parse_args()

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    event_txt = out_dir / args.event_txt_name
    event_csv = out_dir / "selected_events.csv"
    pick_root = out_dir / args.pick_root_name

    if args.event_list:
        events = read_event_list(Path(args.event_list))
    else:
        if not args.ntuple:
            raise RuntimeError("Provide either --event-list or --ntuple.")
        events = select_events_from_ntuple(args.ntuple, args.tree, args.max_events)

    events = unique_preserve_order(events)
    if len(events) > args.max_events:
        events = events[: args.max_events]
    if not events:
        raise RuntimeError("No events available after parsing/selection.")

    write_event_files(events, event_txt, event_csv)

    print(f"[OK] Wrote event list: {event_txt}")
    print(f"[OK] Wrote event table: {event_csv}")
    print("[INFO] Events to pick:")
    for evt in events:
        print(f"       {evt.run}:{evt.lumi}:{evt.event} (fill={evt.fill})")

    run_pick_commands(args.dataset, event_txt, pick_root)

    if not pick_root.exists():
        raise RuntimeError(f"pickevents output not found: {pick_root}")

    print(f"[OK] Built picked file: {pick_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
