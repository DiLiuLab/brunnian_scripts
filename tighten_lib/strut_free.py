#!/usr/bin/env python3
"""Measure how much of a curve is strut-free -- i.e. whether --Timewarp helps.

--Timewarp rescales the length gradient on sections of the curve that carry no
strut and no kink. Such a section straightens asymptotically slowly, because the
length gradient at a vertex is 2*sin(theta/2) and so vanishes as the section
flattens; Timewarp compensates. The ridgerunner(1) man page adds the other half:
"There is a performance cost in turning this on if there are no such segments."

So --Timewarp is not a default to leave on. It is worth it exactly when the
curve has free stretches, and costs throughput when it does not. This measures
that directly from a file, with no run required.

A vertex counts as in contact if some point of the link, excluding its own
along-curve neighbourhood, lies within a chosen multiple of the tube diameter.
The exclusion is by ARCLENGTH, not by vertex index, so it does not depend on the
resolution: neighbours within 1.2 diameters along the same component are skipped.

The number that decides the flag is the longest unbroken free stretch, not the
total free fraction. Scattered free vertices between contacts have no slow mode
to accelerate; a long free arc does.

Usage:  python3 strut_free.py file.xyz [file2.xyz ...]
"""
from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked
# out -- these scripts previously carried an absolute path to one machine.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_vect  # noqa: E402


def thickness(comps) -> float:
    """Tube RADIUS from octrope, so the diameter is 2r. Falls back to 0.5."""
    try:
        with tempfile.TemporaryDirectory() as tmp:
            vect = Path(tmp) / "x.vect"
            write_vect(vect, comps)
            out = subprocess.run(["residual", str(vect)], capture_output=True,
                                 text=True, timeout=1800).stdout
        import re
        m = re.search(r"Thi: ([\d.]+)", out)
        return float(m.group(1)) if m else 0.5
    except Exception:
        return 0.5


def analyse(path: Path, loose: float = 1.15, search: float = 1.30) -> dict:
    comps = [np.asarray(c, float) for c in read_xyz(path)]
    allp = np.vstack(comps)
    n = len(allp)
    lens = [len(c) for c in comps]
    lab = np.concatenate([np.full(m, i) for i, m in enumerate(lens)])
    idx = np.concatenate([np.arange(m) for m in lens])
    D = 2 * thickness(comps)
    edge = np.concatenate([np.linalg.norm(np.roll(c, -1, 0) - c, axis=1) for c in comps])
    mean_edge = float(edge.mean())
    skip = int(np.ceil(1.2 * D / mean_edge))          # arclength exclusion
    tree = cKDTree(allp)
    mind = np.full(n, np.inf)
    for i in range(n):
        for j in tree.query_ball_point(allp[i], search * D):
            if lab[j] == lab[i]:
                m = lens[lab[i]]
                d = abs(idx[j] - idx[i])
                if min(d, m - d) <= skip:
                    continue
            v = float(np.linalg.norm(allp[i] - allp[j]))
            if v < mind[i]:
                mind[i] = v
    free = (mind > loose * D) | np.isinf(mind)
    longest, off = 0, 0
    for m in lens:                                     # longest run, wrapping
        f = free[off:off + m]
        off += m
        run = best = 0
        for k in range(2 * m):
            if f[k % m]:
                run += 1
                best = max(best, run)
            else:
                run = 0
        longest = max(longest, min(best, m))
    return dict(n=n, D=D, touching=int((mind < 1.02 * D).sum()),
                near=int((mind < 1.05 * D).sum()), free=int(free.sum()),
                longest=longest, arc=longest * mean_edge / D)


def main() -> int:
    paths = [Path(a) for a in sys.argv[1:]]
    if not paths:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    print(f"{'file':<42}{'verts':>6}{'touch%':>8}{'free':>6}{'longest':>9}{'arc/D':>8}  verdict")
    for p in paths:
        try:
            r = analyse(p)
        except Exception as exc:
            print(f"{p.name:<42} failed: {str(exc)[:40]}")
            continue
        verdict = ("Timewarp WORTH IT" if r["arc"] >= 2.0 else
                   "marginal" if r["longest"] else "Timewarp is pure overhead")
        print(f"{p.name:<42}{r['n']:>6}{r['touching']/r['n']*100:>7.1f}%"
              f"{r['free']:>6}{r['longest']:>9}{r['arc']:>8.1f}  {verdict}")
    print("\n  free = no other strand within 1.15 tube diameters (along-curve neighbours excluded)")
    print("  longest = longest unbroken free run of vertices; arc/D is its length in tube diameters")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
