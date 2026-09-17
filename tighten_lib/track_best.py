#!/usr/bin/env python3
"""Record where each arm's best ropelength and best residual actually occurred.

A run's best configuration is only recoverable from a SNAPSHOT, and snapshots
are written every --snapshot-interval steps. If the best step falls between two
of them, the geometry is gone: the 4_LC mirror run touched 113.58477 at step
19405 with snapshots at 18000 and 20000, and the best file on disk was 113.5896.

This does not change that, but it records what was lost and where, so a missed
optimum is a known quantity rather than a silent one. Run it any time, during or
after a run. It reads only the .dat traces and never touches the run.

For each arm it reports the best step, the nearest snapshot on each side, and
the penalty -- how much worse the recoverable snapshot is than the trace best.
A penalty small against the differences being compared is fine; one comparable
to them means the snapshot interval was too coarse for the claim being made.
"""
from __future__ import annotations

import sys
from pathlib import Path

W = Path(__file__).resolve().parent
ARMS = ["3_OC", "4BL_wider", "4_LC", "2_BC"]


def trace(arm: str, name: str, col: int = 1):
    g = list(W.glob(f"rr_{arm}/*.rr/logfiles/{name}.dat"))
    if not g:
        return []
    out = []
    for line in g[0].read_text().split("\n")[:-1]:      # drop a half-written last line
        p = line.split()
        if len(p) > col:
            try:
                out.append((int(p[0]), float(p[col])))
            except ValueError:
                pass
    return out


def snapshots(arm: str) -> list[int]:
    steps = []
    for p in W.glob(f"rr_{arm}/*.rr/snapshots/*.vect"):
        stem = p.name.rsplit(".", 2)
        if len(stem) == 3 and stem[1].isdigit():
            steps.append(int(stem[1]))
    return sorted(set(steps))


def report(arm: str) -> None:
    rope, res = trace(arm, "ropelength"), trace(arm, "residual")
    if not rope:
        print(f"  {arm:<11} no trace yet")
        return
    snaps = snapshots(arm)
    print(f"  {arm:<11} {len(rope)} logged steps, {len(snaps)} snapshots")
    for label, series, better in (("ropelength", rope, min), ("residual", res, min)):
        if not series:
            continue
        step, val = better(series, key=lambda r: r[1])
        near = [s for s in snaps if s <= step] or [None]
        above = [s for s in snaps if s >= step] or [None]
        lo, hi = near[-1], above[0]
        at = lambda s: (min(series, key=lambda r: abs(r[0] - s))[1] if s is not None else None)
        best_snap = min([v for v in (at(lo), at(hi)) if v is not None], default=None)
        pen = (best_snap - val) if best_snap is not None else None
        print(f"      best {label:<10} {val:.5f} at step {step}"
              f"   snapshots {lo}/{hi}"
              + (f"   recoverable {best_snap:.5f}, penalty {pen:+.5f}" if pen is not None else ""))


def main() -> int:
    print(f"best-so-far tracking, {W}\n")
    for arm in ARMS:
        report(arm)
    print("\n  penalty = how much worse the best RECOVERABLE snapshot is than the")
    print("  best the trace ever saw. Judge it against the differences being compared.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
