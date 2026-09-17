#!/usr/bin/env python3
"""Increase a link's resolution by periodic spline interpolation.

The author's own workflow (ridgerunner(1), EXAMPLES) is a resolution ladder:
minimize at low resolution, refine, minimize again, repeat -- "ridgerunner is
usually called first on a relatively low resolution file to reach an
approximately tight configuration, then the resolution is increased with
splinevect --mps and the curve is run again." The paper puts numbers on it: 2
vertices per unit of ropelength, then 4, then 8, with most of the runtime in
the last stage.

splinevect is NOT in plCurve 11.2.2 -- the package builds plcurvature,
ropelength, struts, mrlocs, randompolygon and the pd tools, and nothing else.
The library still has plc_convert_to_spline / plc_convert_from_spline, and a
periodic cubic spline reproduces it. The "--mps" of the original is
minrad-preserving; plain splining can leave a vertex whose turning angle
violates MinRad, but ridgerunner's own MinRad constraint repairs that in the
first steps, so the distinction costs a little work rather than correctness.

Two choices that matter, both defaults here:

  * VERTICES ARE ALLOCATED IN PROPORTION TO COMPONENT LENGTH, not to the
    existing counts, so every component ends with the same edge length. Counts
    drift away from that during a run -- measured on the 4_LC mirror result,
    component A held 38.7% of the vertices for 34.0% of the length while C held
    10.5% for 14.7% -- and refinement is the natural moment to correct it.

  * SAMPLING IS BY ARCLENGTH, which preserves a mirror symmetry in which each
    component maps to ITSELF: the reflection fixes arclength, so a symmetric
    curve sampled uniformly in arclength stays symmetric. Pass --check-mirror to
    measure that rather than trust it.

Usage:
    python3 refine_xyz.py in.xyz out.xyz --vu 10            # target v/u
    python3 refine_xyz.py in.xyz out.xyz --vertices 1200    # or an explicit total
    python3 refine_xyz.py in.xyz --report                   # just say where it is
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked
# out -- these scripts previously carried an absolute path to one machine.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_xyz, write_vect  # noqa: E402


def octrope(path: Path) -> float:
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / "x.vect"
        write_vect(vect, read_xyz(path))
        return float(subprocess.run(["ropelength", "-q", str(vect)],
                                    capture_output=True, text=True, timeout=300).stdout)


def minrad(comps) -> float:
    worst = np.inf
    for c in comps:
        c = np.asarray(c)
        e = np.roll(c, -1, 0) - c
        ln = np.linalg.norm(e, axis=1)
        u = e / ln[:, None]
        cos = np.clip((u * np.roll(u, 1, 0)).sum(1), -1, 1)
        t = np.tan(np.maximum(np.arccos(cos), 1e-12) / 2)
        worst = min(worst, float((np.minimum(np.roll(ln, 1), ln) / (2 * t)).min()))
    return worst


def clearance(comps) -> float:
    comps = [np.asarray(c) for c in comps]
    return min(float(np.linalg.norm(comps[i][:, None, :] - comps[j][None, :, :], axis=2).min())
               for i in range(len(comps)) for j in range(i + 1, len(comps)))


def subdivide(points: np.ndarray, n: int) -> np.ndarray:
    """Place n points uniformly by arclength ON the existing polygon.

    This changes nothing geometrically -- the curve is identical, so length,
    thickness, minrad and ropelength are all preserved exactly -- and simply
    hands ridgerunner more degrees of freedom to use. A spline instead SMOOTHS,
    which sounds better but overshoots at tight turns: measured on the 4_LC
    mirror result, cubic splining to 10 v/u dropped minrad 0.500 -> 0.434, and
    since thickness is min(dcsd/2, minrad) that inflated ropelength 113.59 ->
    130.87 before a single step. The paper's splinevect is "minrad-preserving"
    precisely to avoid this; plain splining is not.
    """
    closed = np.vstack([points, points[:1]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    want = np.linspace(0, arc[-1], n, endpoint=False)
    idx = np.clip(np.searchsorted(arc, want, side="right") - 1, 0, len(points) - 1)
    t = (want - arc[idx]) / np.maximum(seg[idx], 1e-15)
    return closed[idx] + t[:, None] * (closed[idx + 1] - closed[idx])


def resample(points: np.ndarray, n: int) -> np.ndarray:
    """Periodic cubic spline through a closed polygon, sampled uniformly by arclength."""
    closed = np.vstack([points, points[:1]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    spline = CubicSpline(s, closed, bc_type="periodic")
    # one pass of arclength equalisation: the spline's parameter is the ORIGINAL
    # chord length, which is not the spline's own arclength, so a naive uniform
    # sample in s is not uniform along the curve
    dense = spline(np.linspace(0, total, max(8 * n, 4000), endpoint=False))
    d = np.linalg.norm(np.diff(np.vstack([dense, dense[:1]]), axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(d)])
    want = np.linspace(0, arc[-1], n, endpoint=False)
    idx = np.searchsorted(arc, want, side="right") - 1
    idx = np.clip(idx, 0, len(dense) - 1)
    return dense[idx]


def relax_minrad(comps, target, rounds=4000, step=0.25):
    """Round only the corners that violate MinRad, leaving the rest alone.

    MinRad is min(|e_prev|,|e_next|) / (2 tan(theta/2)), so it scales with EDGE
    LENGTH: refine a polygon k-fold at a fixed turning angle and minrad falls by
    k. That is why plain refinement breaks it, and why the paper's splinevect is
    "minrad-preserving". Measured on the 4_LC mirror result refined 721 -> 1136:
    subdivision gave minrad 0.297, cubic splining 0.434, against 0.500 before.

    The violation is local -- 9 of 721 vertices sat within 10% of the limit,
    median 1.01 -- so this applies Laplacian smoothing ONLY at offending
    vertices and their immediate neighbours, which rounds those corners into
    arcs without touching the rest of the curve.
    """
    out = [np.array(c, dtype=float) for c in comps]
    for _ in range(rounds):
        worst = np.inf
        moved = False
        for c in out:
            e = np.roll(c, -1, 0) - c
            ln = np.linalg.norm(e, axis=1)
            u = e / ln[:, None]
            cos = np.clip((u * np.roll(u, 1, 0)).sum(1), -1, 1)
            t = np.tan(np.maximum(np.arccos(cos), 1e-12) / 2)
            mr = np.minimum(np.roll(ln, 1), ln) / (2 * t)
            worst = min(worst, float(mr.min()))
            bad = np.where(mr < target)[0]
            if not len(bad):
                continue
            moved = True
            hit = np.unique(np.concatenate([(bad - 1) % len(c), bad, (bad + 1) % len(c)]))
            lap = 0.5 * (np.roll(c, 1, 0) + np.roll(c, -1, 0)) - c
            c[hit] += step * lap[hit]
        if not moved:
            break
    return out, worst


def mirror_residual(comps, axis=2) -> float:
    P = np.vstack([np.asarray(c) for c in comps])
    X = P - P.mean(0)
    Q = X.copy()
    Q[:, axis] *= -1
    d, _ = cKDTree(X).query(Q)
    return float(d.mean() / np.linalg.norm(X, axis=1).mean())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path, nargs="?")
    ap.add_argument("--vu", type=float, help="Target vertices per unit ropelength.")
    ap.add_argument("--vertices", type=int, help="Target total vertex count.")
    ap.add_argument("--report", action="store_true", help="Report resolution and exit.")
    ap.add_argument("--check-mirror", action="store_true",
                    help="Measure the z=0 mirror residual before and after.")
    ap.add_argument("--mode", choices=("subdivide", "spline"), default="subdivide",
                    help="subdivide: place points on the existing polygon -- geometry,\n"
                         "  and therefore ropelength and minrad, exactly unchanged.\n"
                         "spline: smooth through the vertices -- can violate MinRad.")
    ap.add_argument("--fix-minrad", type=float, metavar="TARGET",
                    help="Round violating corners until minrad reaches TARGET.\n"
                         "  Use the input's own minrad (0.5 at the MinRad limit).")
    ap.add_argument("--resymmetrize", action="store_true",
                    help="Restore an exact z=0 mirror afterwards. Arclength sampling only\n"
                         "  preserves it when sampling starts at a fixed point of the\n"
                         "  reflection, which is not generally true.")
    args = ap.parse_args()

    comps = [np.array(c) for c in read_xyz(args.input)]
    rope = octrope(args.input)
    lengths = [float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in comps]
    total = sum(lengths)
    diameter = 2 * total / rope
    n0 = sum(len(c) for c in comps)
    edges = np.concatenate([np.linalg.norm(np.roll(c, -1, 0) - c, axis=1) for c in comps])

    print(f"input  {args.input.name}")
    print(f"  ropelength {rope:.4f}   vertices {n0}   v/u {n0/rope:.2f}")
    print(f"  mean d/D {edges.mean()/diameter:.4f}   max d/D {edges.max()/diameter:.4f}"
          f"   edge max/min {edges.max()/edges.min():.2f}")
    print(f"  per component {'/'.join(str(len(c)) for c in comps)}"
          f"   length shares {'/'.join(f'{100*x/total:.1f}' for x in lengths)}")
    if args.report:
        return 0

    if args.vertices:
        target = args.vertices
    elif args.vu:
        target = int(round(args.vu * rope))
    else:
        print("give --vu or --vertices", file=sys.stderr)
        return 1

    # allocate in proportion to LENGTH so every component gets the same edge length
    want = [max(12, int(round(target * x / total))) for x in lengths]
    before_mirror = mirror_residual(comps) if args.check_mirror else None
    fn = subdivide if args.mode == "subdivide" else resample
    out = [fn(c, n) for c, n in zip(comps, want)]
    if args.fix_minrad:
        out, got = relax_minrad(out, args.fix_minrad)
        print(f"  minrad repair -> {got:.4f} (target {args.fix_minrad})")
    if args.resymmetrize:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "sm", "/Users/diliu/Documents/Claude_tmp/symmetrize_mirror.py")
        sm = importlib.util.module_from_spec(spec)
        sys.argv = ["sm"]
        spec.loader.exec_module(sm)
        centre = np.vstack(out).mean(0)
        shifted = [c - centre for c in out]
        mirrored = [c @ sm.MIRROR.T for c in shifted]
        fixed = []
        for i, m in enumerate(mirrored):
            _, flip, off = sm.best_alignment(shifted[i], m)
            src = m[::-1] if flip else m
            src = np.roll(src, -off, axis=0)
            fixed.append(0.5 * (shifted[i] + src) + centre)
        out = fixed
    write_xyz(args.output, [c.tolist() for c in out], 9)

    rope2 = octrope(args.output)
    e2 = np.concatenate([np.linalg.norm(np.roll(c, -1, 0) - c, axis=1) for c in out])
    l2 = [float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in out]
    d2 = 2 * sum(l2) / rope2
    n2 = sum(len(c) for c in out)
    print(f"\noutput {args.output.name}")
    print(f"  ropelength {rope2:.4f}   vertices {n2}   v/u {n2/rope2:.2f}")
    print(f"  mean d/D {e2.mean()/d2:.4f}   max d/D {e2.max()/d2:.4f}"
          f"   edge max/min {e2.max()/e2.min():.2f}")
    print(f"  per component {'/'.join(str(len(c)) for c in out)}")
    print(f"  minrad {minrad(out):.4f} (was {minrad(comps):.4f})   "
          f"clearance {clearance(out)/d2:.4f} D (was {clearance(comps)/diameter:.4f} D)")
    if args.check_mirror:
        print(f"  mirror residual {mirror_residual(out):.2e} (was {before_mirror:.2e})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
