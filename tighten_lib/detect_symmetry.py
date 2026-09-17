#!/usr/bin/env python3
"""Find the symmetries of a link that ridgerunner could be asked to enforce.

ridgerunner offers exactly three (ridgerunner_main.c:843-893):
  --Symmetry=Z/pZ   plc_rotation_group(link, zaxis, p)  -- rotation by 2pi/p about z
  --Symmetry=D2     plc_reflection_group(link, zaxis)   -- ONE mirror, the plane
                                                           normal to z (symmetries.c:377)
  --Symmetry=cplanes  reflections over all three coordinate planes

So the useful question per structure is: is there an axis about which it is
nearly C_p symmetric, or a plane it nearly mirrors across, and how near?

The threshold that matters is set by plc_build_symmetry, which matches each
vertex to its NEAREST image and fails if two vertices claim one target. Its
docstring: this needs symmetry "to within < (1/2) an edgelength". Deviations
are therefore reported in units of the mean edge length, not of the radius --
under about 0.5 the constraint builds directly; above it the structure must be
symmetrized first (symmetrize_mirror.py), which is cheap if the deviation is
small and distorting if it is not.

Usage:  python3 detect_symmetry.py file.xyz [file2.xyz ...]
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


def octrope(path):
    with tempfile.TemporaryDirectory() as t:
        v = Path(t) / "x.vect"
        write_vect(v, read_xyz(path))
        return float(subprocess.run(["ropelength", "-q", str(v)],
                                    capture_output=True, text=True, timeout=300).stdout)


def directions(n=2500):
    """Fibonacci sphere -- even coverage without a preferred axis."""
    for i in range(n):
        z = 1 - 2 * (i + 0.5) / n
        th = np.pi * (1 + 5 ** 0.5) * i
        r = np.sqrt(max(0.0, 1 - z * z))
        yield np.array([r * np.cos(th), r * np.sin(th), z])


def scan(path):
    comps = [np.array(c) for c in read_xyz(path)]
    P = np.vstack(comps)
    edge = float(np.mean([np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).mean()
                          for c in comps]))
    X = P - P.mean(0)
    tree = cKDTree(X)

    def dev(M):
        d, _ = tree.query(X @ M.T)
        return float(d.mean() / edge), float(d.max() / edge)

    best = {}
    for n in directions():
        M = np.eye(3) - 2 * np.outer(n, n)                    # mirror
        m = dev(M)
        if "mirror" not in best or m[0] < best["mirror"][0]:
            best["mirror"] = (*m, n)
        for p in (2, 3, 4):
            th = 2 * np.pi / p
            K = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
            R = np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)
            m = dev(R)
            k = f"C{p}"
            if k not in best or m[0] < best[k][0]:
                best[k] = (*m, n)
    return best, edge


def main() -> int:
    print(f"{'structure':<30}{'rope':>8}{'edge':>7}   "
          f"{'best symmetry':<10}{'mean dev':>10}{'max dev':>9}   axis/normal")
    print(f"{'':<30}{'':>8}{'':>7}   {'':<10}{'(edge lengths)':>19}")
    for arg in sys.argv[1:]:
        p = Path(arg)
        if not p.is_file():
            print(f"{p.name:<30} missing")
            continue
        best, edge = scan(p)
        rope = octrope(p)
        first = True
        for k in ("mirror", "C2", "C3", "C4"):
            mean, mx, n = best[k]
            verdict = ("BUILDS" if mx < 0.5 else
                       "symmetrize first" if mean < 1.5 else "no")
            lead = f"{p.name:<30}{rope:>8.2f}{edge:>7.3f}" if first else " " * 45
            print(f"{lead}   {k:<10}{mean:>10.2f}{mx:>9.2f}   "
                  f"{np.round(n,3)}  {verdict}")
            first = False
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
