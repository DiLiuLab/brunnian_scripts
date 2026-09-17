#!/usr/bin/env python3
"""Make a link exactly mirror-symmetric across z=0, so ridgerunner can constrain it.

ridgerunner's --Symmetry=D2 calls plCurve's plc_build_symmetry, which matches
each vertex to the NEAREST vertex under the reflection and fails outright if two
vertices claim the same target. Its docstring: "This can fail if the curve is
not symmetric to begin with to within < (1/2) an edgelength. In this case, you
should build the symmetry map yourself and then symmetrize the results."

Measured on the 4_LC fork at 130.65: the best mirror had mean deviation 0.028 of
the mean radius (~0.08 D) but a maximum of 0.051 (~0.15 D) against a half-edge
of roughly 0.1 D, so the build failed.

This constructs the correspondence explicitly rather than by nearest-vertex:
for each component it finds which component the reflection sends it to, then the
cyclic offset and orientation that best align them, then averages each vertex
with the reflection of its partner. The result is symmetric to machine
precision, so the nearest-vertex matching downstream is unambiguous.

Usage:  python3 symmetrize_mirror.py in.xyz out.xyz
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked
# out -- these scripts previously carried an absolute path to one machine.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_xyz  # noqa: E402

MIRROR = np.diag([1.0, 1.0, -1.0])      # reflection in the plane z = 0


def best_alignment(target, source):
    """Cyclic offset and direction aligning `source` to `target`, both closed."""
    n = len(target)
    best = None
    for flip in (False, True):
        s = source[::-1] if flip else source
        for off in range(len(s)):
            rolled = np.roll(s, -off, axis=0)
            if len(rolled) != n:
                idx = np.linspace(0, len(rolled), n, endpoint=False).astype(int)
                rolled = rolled[idx]
            cost = float(np.linalg.norm(rolled - target, axis=1).mean())
            if best is None or cost < best[0]:
                best = (cost, flip, off)
    return best


def main() -> int:
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    comps = [np.array(c) for c in read_xyz(Path(sys.argv[1]))]
    centre = np.vstack(comps).mean(0)
    comps = [c - centre for c in comps]
    mirrored = [c @ MIRROR.T for c in comps]

    # which component does each map to?
    pairing = {}
    for i, m in enumerate(mirrored):
        costs = [(best_alignment(comps[j], m)[0], j) for j in range(len(comps))]
        costs.sort()
        pairing[i] = costs[0][1]
        print(f"  component {i} -> {costs[0][1]}   mean mismatch {costs[0][0]:.4f}")

    if sorted(pairing.values()) != list(range(len(comps))):
        print(f"error: reflection is not a permutation of components ({pairing})",
              file=sys.stderr)
        return 1

    out = [None] * len(comps)
    done = set()
    for i, j in pairing.items():
        if i in done:
            continue
        _, flip, off = best_alignment(comps[j], mirrored[i])
        src = mirrored[i][::-1] if flip else mirrored[i]
        src = np.roll(src, -off, axis=0)
        if len(src) != len(comps[j]):
            idx = np.linspace(0, len(src), len(comps[j]), endpoint=False).astype(int)
            src = src[idx]
        avg = 0.5 * (comps[j] + src)                 # component j, symmetrized
        out[j] = avg
        if i != j:                                   # and its partner, by reflection
            back = avg @ MIRROR.T
            back = np.roll(back, off, axis=0)
            out[i] = back[::-1] if flip else back
            done.add(i)
        done.add(j)

    out = [c + centre for c in out]
    write_xyz(Path(sys.argv[2]), [c.tolist() for c in out], 9)

    # report how exact the result is
    P = np.vstack([c - centre for c in out])
    from scipy.spatial import cKDTree
    d, _ = cKDTree(P).query(P @ MIRROR.T)
    scale = float(np.linalg.norm(P, axis=1).mean())
    print(f"\nmirror residual after symmetrizing: mean {d.mean()/scale:.2e}  "
          f"max {d.max()/scale:.2e}   (was 0.028 / 0.051)")
    print(f"wrote {sys.argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
