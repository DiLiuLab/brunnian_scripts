#!/usr/bin/env python3
"""Squeeze the empty axial gaps out of a link, then let RidgeRunner repair it.

A tightened link can end up as two locally-ideal clusters held apart by a long
bridge component. RidgeRunner will not close that gap: shortening the bridge
first requires its wrap regions to give, which is locally uphill, so descent
sits there forever. Measured on 4_LC: 70000 steps moved the gap from 6.32 to
6.50 tube diameters while ropelength fell 52 units, all of it from elsewhere.

The map here is t -> f(t) along one axis, f piecewise linear and strictly
increasing: identity inside each cluster's axial span, slope `factor` across
the empty gaps between them. Because f is monotone the map is a homeomorphism
of R^3, so the LINK TYPE IS PRESERVED EXACTLY -- no HOMFLY check needed on the
output of this script (still check after re-tightening). What it does break is
thickness: strands end up closer than a tube diameter, which is the repair
RidgeRunner is for.

Only gaps *between* clusters are squeezed. The regions beyond the outermost
cluster are the bridge's own end caps and are left alone.

Usage:
    python3 compress_gap.py in.xyz out.xyz --factor 0.3 [--bridge 0]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked
# out -- these scripts previously carried an absolute path to one machine.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_xyz, write_vect  # noqa: E402


def _octrope(path: Path) -> float:
    """Ropelength, used only to recover the tube diameter D = 2*Len/Rop."""
    import subprocess
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / "x.vect"
        write_vect(vect, read_xyz(path))
        return float(subprocess.run(["ropelength", "-q", str(vect)],
                                    capture_output=True, text=True, timeout=300).stdout)


def principal_axis(points: np.ndarray) -> np.ndarray:
    _, _, vt = np.linalg.svd(points - points.mean(0), full_matrices=False)
    return vt[0]


def cluster_axis(comps, bridge, diameter, join=1.5):
    """Direction between the two clusters of non-bridge components.

    The principal axis is only ever a proxy for this, and it fails outright
    once a structure folds up. Measured on a 4_LC configuration at 187.54: the
    real separation was 3.92 tube diameters along the cluster-to-cluster
    direction, which sat 88.4 degrees from PC1 -- so a PC1 scan reported NO GAP
    on a structure that was 76% spacer.

    Components are grouped by proximity (join < 1.5 D counts as one cluster),
    which is what "cluster" means physically here: B and D touching at 1.00 D
    are one object that the bridge grips twice, not two separate ends.

    Returns None unless exactly two clusters are found -- three or more is a
    different problem and a single axis cannot express it.
    """
    others = [i for i in range(len(comps)) if i != bridge]

    def near(i, j):
        d = np.linalg.norm(comps[i][:, None, :] - comps[j][None, :, :], axis=2).min()
        return d / diameter < join

    groups = [[i] for i in others]
    merged = True
    while merged:
        merged = False
        for a in range(len(groups)):
            for b in range(a + 1, len(groups)):
                if any(near(i, j) for i in groups[a] for j in groups[b]):
                    groups[a] += groups[b]
                    groups.pop(b)
                    merged = True
                    break
            if merged:
                break
    if len(groups) != 2:
        return None, groups
    c1 = np.vstack([comps[i] for i in groups[0]]).mean(0)
    c2 = np.vstack([comps[i] for i in groups[1]]).mean(0)
    axis = c2 - c1
    return axis / np.linalg.norm(axis), groups


def cluster_spans(comps, bridge, axis, centre):
    """Axial [lo, hi] of every component except the bridge, merged where they overlap."""
    spans = []
    for i, c in enumerate(comps):
        if i == bridge:
            continue
        t = (c - centre) @ axis
        spans.append([t.min(), t.max()])
    spans.sort()
    merged = [spans[0]]
    for lo, hi in spans[1:]:
        if lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return merged


def build_map(merged, factor, blend):
    """Monotone t -> t'. Gaps shrink by `factor`; transitions are smoothed.

    A piecewise-linear squeeze has a corner where each gap meets a cluster, and
    compressing hard turns that corner into a kink: thickness collapses onto
    minrad and RidgeRunner spends the run undoing it rather than tightening.
    Measured at factor 0.15, thickness and minrad were both 0.533 -- entirely
    curvature-limited. So the derivative is ramped with a smoothstep over
    `blend` instead of switched, which keeps f in C^1.

    f is built by integrating a strictly positive derivative, so it stays a
    homeomorphism of R^3 and the link type is still preserved exactly.
    """
    lo = merged[0][0]
    hi = merged[-1][1]
    pad = max(blend, 1.0)
    grid = np.linspace(lo - pad, hi + pad, 20001)

    def smoothstep(x):                      # 0 -> 1 with zero slope at both ends
        x = np.clip(x, 0.0, 1.0)
        return x * x * (3.0 - 2.0 * x)

    # weight 1 = full scale (inside a cluster), 0 = squeezed (gap core)
    w = np.ones_like(grid)
    for (a_lo, a_hi), (b_lo, _) in zip(merged, merged[1:]):
        if blend <= 0:
            w = np.where((grid > a_hi) & (grid < b_lo), 0.0, w)
            continue
        # ramp down after a_hi, back up before b_lo; half the gap each at most
        half = 0.5 * (b_lo - a_hi)
        width = min(blend, half)
        down = 1.0 - smoothstep((grid - a_hi) / width)
        up = smoothstep((grid - (b_lo - width)) / width)
        inside = (grid > a_hi) & (grid < b_lo)
        w = np.where(inside, np.maximum(np.minimum(down + up, 1.0), 0.0), w)

    deriv = factor + (1.0 - factor) * w
    tp = np.concatenate([[0.0], np.cumsum(0.5 * (deriv[1:] + deriv[:-1]) * np.diff(grid))])
    tp += grid[0] - tp[0]                   # anchor the low end

    def f(t):
        out = np.interp(t, grid, tp)
        out = np.where(t < grid[0], t + (tp[0] - grid[0]), out)
        out = np.where(t > grid[-1], t + (tp[-1] - grid[-1]), out)
        return out
    return f


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--factor", type=float, required=True,
                    help="Gap scale: 1.0 leaves it, 0.3 keeps 30%% of the empty span.")
    ap.add_argument("--bridge", type=int, default=None,
                    help="Index of the bridge component (default: the longest).")
    ap.add_argument("--axis", choices=("pc1", "cluster"), default="cluster",
                    help="Direction to squeeze along. 'cluster' uses the vector between\n                         the two clusters of non-bridge components, which is what the\n                         gap actually lies along; 'pc1' uses the principal axis, which\n                         is only a proxy for it and fails on a folded structure.")
    ap.add_argument("--blend", type=float, default=1.0,
                    help="Smoothing width for the squeeze ramp, in the link's own\n                         units (tube diameter is 1.0). 0 reproduces the old\n                         piecewise-linear map, which kinks under hard squeezes.")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    comps = [np.array(c) for c in read_xyz(args.input)]
    lengths = [float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in comps]
    bridge = args.bridge if args.bridge is not None else int(np.argmax(lengths))

    allpts = np.vstack(comps)
    centre = allpts.mean(0)
    total = sum(lengths)
    # tube diameter, needed to decide which components count as one cluster
    diameter = 2 * total / _octrope(args.input) if args.axis == "cluster" else 1.0
    axis, groups = (None, None)
    if args.axis == "cluster":
        axis, groups = cluster_axis(comps, bridge, diameter)
        if axis is None:
            print(f"error: found {len(groups)} clusters, not 2 -- a single axis cannot "
                  "express that. Squeeze one pair at a time, or pass --axis pc1.",
                  file=sys.stderr)
            return 1
    else:
        axis = principal_axis(allpts)
    merged = cluster_spans(comps, bridge, axis, centre)
    if len(merged) < 2:
        print("error: clusters do not separate along the principal axis; "
              "nothing to squeeze", file=sys.stderr)
        return 1

    f = build_map(merged, args.factor, args.blend)
    out = []
    for c in comps:
        t = (c - centre) @ axis
        out.append(c + np.outer(f(t) - t, axis))

    new_lengths = [float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in out]
    if not args.quiet:
        gaps = [f"{b[0] - a[1]:.2f}" for a, b in zip(merged, merged[1:])]
        how = "cluster-to-cluster" if args.axis == "cluster" else "principal"
        print(f"axis {np.round(axis, 3)} ({how})   bridge component {bridge}")
        if groups:
            print("clusters " + " | ".join("".join("ABCD"[i] for i in g) for g in groups))
        print(f"cluster spans " + "  ".join(f"[{lo:.2f},{hi:.2f}]" for lo, hi in merged))
        print(f"gaps {', '.join(gaps)} -> scaled by {args.factor} (blend {args.blend})")
        print(f"length {sum(lengths):.2f} -> {sum(new_lengths):.2f} "
              f"({100 * (sum(new_lengths) / sum(lengths) - 1):+.1f}%)")
        # worst clearance after the squeeze, so the caller knows how hard the repair is
        worst = min(
            float(np.linalg.norm(out[i][:, None, :] - out[j][None, :, :], axis=2).min())
            for i in range(len(out)) for j in range(i + 1, len(out)))
        print(f"closest approach between components: {worst:.3f} "
              f"(tube diameter is 1.0)")

    write_xyz(args.output, [c.tolist() for c in out], 9)
    if not args.quiet:
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
