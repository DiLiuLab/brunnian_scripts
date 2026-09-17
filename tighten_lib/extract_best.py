#!/usr/bin/env python3
"""Pull out BOTH the lowest-ropelength and the most nearly critical snapshot.

A run does not optimise one thing. Ropelength and residual can diverge, and on
the 2_BC continuation they did: residual bottomed near step 4000 at 0.087 and
then climbed monotonically to 0.151 over the next 16000 steps while ropelength
kept falling. The run was buying length with criticality.

ridgerunner's wrapper selects by DOMINANCE -- it keeps a configuration only if
it is no worse on ropelength, residual and struts, and better on one. When the
two measures diverge neither configuration dominates, so the rule cannot choose
and the run's own pick is whichever it happened to land on. That is the right
conservative behaviour and the wrong thing to report, because the lowest bound
and the best-conditioned configuration are different scientific claims:

  lowest ropelength  -> the tightest embedded curve found: Rop(L) <= this
  lowest residual    -> the closest to ropelength-critical, i.e. the one whose
                        geometry is most nearly an actual critical point

This finds the snapshot nearest each optimum, writes it as .xyz, and measures
it. Recovery is limited by --SnapshotInterval: the trace optimum usually falls
between two snapshots, and the penalty for that is reported so a missed optimum
is a known quantity rather than a silent one.

Usage:  python3 extract_best.py <run-dir> [outdir]
        run-dir is the --work-dir of the run, holding a *.rr directory.
"""
from __future__ import annotations

import re
import subprocess
import sys
import tempfile
from pathlib import Path

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked
# out -- these scripts previously carried an absolute path to one machine.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_vect, write_xyz, write_vect  # noqa: E402


def trace(logs: Path, name: str) -> list[tuple[int, float]]:
    out = []
    for line in (logs / f"{name}.dat").read_text().split("\n")[:-1]:
        p = line.split()
        if len(p) > 1:
            try:
                out.append((int(p[0]), float(p[1])))
            except ValueError:
                pass
    return out


def measure(path: Path) -> str:
    comps = read_vect(path) if path.suffix == ".vect" else None
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / "x.vect"
        write_vect(vect, comps)
        out = subprocess.run(["residual", str(vect)], capture_output=True,
                             text=True, timeout=1800).stdout
    g = lambda pat: (re.search(pat, out).group(1) if re.search(pat, out) else "?")
    return (f"Rop {g(r'Rop: ([0-9.]+)')}  Thi {g(r'Thi: ([0-9.]+)')}  "
            f"struts {g(r'Struts: ([0-9]+)')}+{g(r'MrStruts: ([0-9]+)')}  "
            f"resid {g(r'Residual: ([0-9.eE+-]+)')}")


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__.strip().splitlines()[-2], file=sys.stderr)
        return 2
    run = Path(sys.argv[1])
    outdir = Path(sys.argv[2]) if len(sys.argv) > 2 else run
    rr = next(iter(run.glob("*.rr")), None)
    if rr is None:
        print(f"no *.rr directory under {run}", file=sys.stderr)
        return 1
    logs, snaps = rr / "logfiles", rr / "snapshots"
    rope, res = trace(logs, "ropelength"), trace(logs, "residual")

    available = {}
    for p in snaps.glob("*.vect"):
        parts = p.name.rsplit(".", 2)
        if len(parts) == 3 and parts[1].isdigit() and "struts" not in p.name:
            available[int(parts[1])] = p
    if not available:
        print("no numbered snapshots found", file=sys.stderr)
        return 1

    # Ropelength is smooth enough to take its single-step minimum. RESIDUAL IS
    # NOT: it spikes to 1.0 whenever the strut set is rebuilt and momentarily
    # empties, so its single-step minimum is a spike, not a state. Choosing a
    # snapshot that way is worthless -- on the 2_BC run the best single step was
    # 0.075 while the trace at the nearest snapshot read 0.307. Each snapshot is
    # therefore scored by the MEDIAN residual in a window centred on it, which
    # is the same rule the report uses for every residual it quotes.
    import statistics
    HALF = 250

    def windowed(step: int) -> float:
        w = [v for s, v in res if abs(s - step) <= HALF]
        return statistics.median(w) if w else float("inf")

    for label, series in (("ropelength", rope), ("residual", res)):
        if not series:
            continue
        if label == "residual":
            near = min(available, key=windowed)
            val = windowed(near)
            step = near
            got = val
        else:
            step, val = min(series, key=lambda r: r[1])
            near = min(available, key=lambda s: abs(s - step))
            got = min(series, key=lambda r: abs(r[0] - near))[1]
        dest = outdir / f"best_{label}_step{near}.xyz"
        write_xyz(dest, read_vect(available[near]), 9)
        tag = "windowed median" if label == "residual" else "single step"
        print(f"  best {label:<11} {val:.5f} at step {step}   ({tag})")
        print(f"      nearest snapshot step {near} -> {dest.name}")
        print(f"      trace there {got:.5f}, penalty {got - val:+.5f}")
        print(f"      measured: {measure(available[near])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
