#!/usr/bin/env python3
"""Per-2000-step ropelength record for every arm of a run set.

Reads each arm's logfiles/ropelength.dat directly rather than its snapshots, so
it works on runs still in progress. Handles the two parsing traps: a live run
leaves its last line half-written, and the logs are decimated as a run grows so
an exact step number often is not present -- the nearest logged step is used.

Usage:  python3 blocks2k.py <run-set-dir> [block-size]
        python3 blocks2k.py lc_fromraw_sq 2000
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

TMP = Path("/Users/diliu/Documents/Claude_tmp")


def read(run_dir: Path, name: str, col: int = 1) -> dict[int, float]:
    path = run_dir / "logfiles" / f"{name}.dat"
    if not path.is_file():
        return {}
    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    if lines and not text.endswith("\n"):
        lines.pop()
    out, last = {}, -1
    for line in lines:
        f = line.split()
        if len(f) < col + 1:
            continue
        try:
            step, value = int(float(f[0])), float(f[col])
        except ValueError:
            continue
        if step <= last:
            continue
        last = step
        out[step] = value
    return out


def main() -> int:
    root = TMP / sys.argv[1] if len(sys.argv) > 1 else None
    width = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    if root is None or not root.is_dir():
        print(__doc__)
        return 1

    arms = sorted(root.glob("rr_*"))
    for arm in arms:
        run = next(arm.glob("*.rr"), None)
        if run is None:
            continue
        rope = read(run, "ropelength")
        if not rope:
            continue
        res, st = read(run, "residual"), read(run, "strutcount", 1)
        near = lambda d, k: d[min(d, key=lambda x: abs(x - k))] if d else float("nan")
        last = max(rope)
        print(f"\n=== {arm.name[3:]}   {last} steps logged")
        print(f"  {'block':>14}{'rope end':>10}{'change':>9}{'%':>8}{'resid':>8}{'struts':>8}")
        prev = near(rope, 0)
        for a in range(0, last, width):
            b = min(a + width, last)
            rb = near(rope, b)
            tail = sorted(v for k, v in res.items() if a < k <= b) or [float("nan")]
            stm = np.median([v for k, v in st.items() if a < k <= b]) if st else float("nan")
            print(f"  {a:>5}-{b:<6}{rb:>10.2f}{rb - prev:>9.2f}"
                  f"{100 * (rb / prev - 1):>7.2f}%{tail[len(tail) // 10]:>8.3f}{stm:>8.0f}")
            prev = rb
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
