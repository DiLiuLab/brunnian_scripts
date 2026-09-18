#!/usr/bin/env python3
"""Contract a link's contact clusters toward their centroid, for ANY arrangement.

The axial squeeze in compress_gap.py handles two clusters with a gap between
them by compressing one coordinate. That is a special case. Once three or more
clusters sit at the corners of some shape there is no single direction of
slack, and the general move is:

    1. find every contact site between components, and cluster the SITES
    2. translate each cluster rigidly toward the clusters' common centroid
    3. carry each connecting strand along, interpolating between the
       translations of the clusters it runs between

Clustering SITES rather than components matters. A component that spans the
structure touches at two corners, so clustering components merges those corners
into one. Measured on a 4_LC configuration at 187.54: clustering components
found two clusters and reported a single 3.9 D gap, while the truth was three
clusters in a near-equilateral triangle of side 6.6-6.8 D, with two edges made
of two strands of the bridge and the third made of four strands of two OTHER
components that span between corners.

Step 3 is where the earlier attempts failed, and the interpolant is the reason:

  * By spatial distance: the displacement varies ACROSS each strand, so every
    strand crossing the transition bends. minrad went to 0.000 for every blend
    width under 3 D, at every factor.
  * By arclength, Gaussian weights: better -- clearance recovered to ~1.00 --
    but normalised Gaussians make the displacement nearly CONSTANT through the
    middle of a long strand and steep near the grips, so the shortening piles
    up at the ends and buckles them. minrad still 0.001 below factor 0.75.
  * By arclength, LINEAR between consecutive grips: the displacement gradient
    is constant along each strand, which is exactly uniform shortening -- the
    same thing the axial squeeze does, expressed intrinsically. This is what is
    implemented here.

Topology is NOT guaranteed the way a monotone axial map guarantees it, so the
result must be checked with verify_topology.py rather than assumed.

Usage:
    python3 contract_clusters.py in.xyz --scan
    python3 contract_clusters.py in.xyz out.xyz --factor 0.6
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked
# out -- these scripts previously carried an absolute path to one machine.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_xyz, write_vect


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


def contact_bands(comps, diameter, contact=1.05, gap=3):
    """Per component, the runs of vertices that touch another component."""
    out = {}
    for i, c in enumerate(comps):
        touching = np.zeros(len(c), bool)
        for j, o in enumerate(comps):
            if i == j:
                continue
            d = np.linalg.norm(c[:, None, :] - o[None, :, :], axis=2).min(1) / diameter
            touching |= d < contact
        idx = np.where(touching)[0]
        bands = []
        if len(idx):
            # append a COPY and rebind -- appending `cur` then clearing it makes
            # every band the same list object, which silently collapses all of a
            # component's bands onto one point
            cur = [int(idx[0])]
            for a, b in zip(idx, idx[1:]):
                if b - a <= gap:
                    cur.append(int(b))
                else:
                    bands.append(cur)
                    cur = [int(b)]
            bands.append(cur)
            if len(bands) > 1 and idx[0] == 0 and idx[-1] == len(c) - 1:
                bands[0] = bands[-1] + bands[0]
                bands.pop()
        out[i] = bands
    return out


def cluster_sites(comps, bands, diameter, join=3.5):
    """Group contact bands into spatial clusters. Returns (centres, band->cluster).

    COMPLETE linkage, not single linkage. Two components can press together
    along an edge rather than only at its ends -- measured here, B and D are in
    contact over 41% of B's length, in runs of 1.4, 4.2 and 2.8. Single linkage
    chains straight through those mid-edge contacts and merges the corners they
    connect, which is exactly what went wrong: it reported clusters {A,C}, {B},
    {D} instead of the three corners. Complete linkage requires every pair in a
    group to be within `join`, so a chain of grazing contacts cannot bridge two
    corners 7 D apart.
    """
    pts, keys = [], []
    for i, bs in bands.items():
        for n, b in enumerate(bs):
            pts.append(comps[i][b].mean(0))
            keys.append((i, n))
    if not pts:
        return [], {}
    pts = np.array(pts)
    label = [-1] * len(pts)
    groups = []
    for k, p in enumerate(pts):
        for g, members in enumerate(groups):
            if max(np.linalg.norm(p - pts[m]) for m in members) / diameter < join:
                members.append(k)
                label[k] = g
                break
        else:
            groups.append([k])
            label[k] = len(groups) - 1
    centres = [pts[m].mean(0) for m in groups]
    return centres, {keys[k]: label[k] for k in range(len(pts))}


def target_positions(centres, targets):
    """Move cluster centres so their pairwise distances match `targets`.

    A homothety -- every cluster moved toward the centroid by the same fraction
    of its own radius -- scales all sides by one factor and so preserves the
    triangle's SHAPE exactly. Measured over four rounds on 4_LC, the
    longest/shortest side ratio stayed at 1.10, 1.11, 1.11, 1.11: the 11%
    asymmetry present at the start survived every contraction untouched. That
    wastes slack, because the factor is limited by whichever side reaches its
    floor first while the longest side keeps room that can never be reached.

    This solves for positions directly instead, by least squares on the
    pairwise distances, starting from the current centres and keeping the
    result near them (the residual translation/rotation is removed by
    re-centring on the original centroid).
    """
    from scipy.optimize import least_squares
    n = len(centres)
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    c0 = np.array(centres, dtype=float)

    # 3 distance constraints against 9 coordinates is underdetermined -- the
    # missing 6 are the rigid motions. Pin them by penalising displacement from
    # the current centres, weighted well below the distance terms so the
    # targets are met and the solution is simply the nearest one that meets them.
    reg = 0.05

    def resid(x):
        p = x.reshape(n, 3)
        d = [np.linalg.norm(p[a] - p[b]) - targets[k] for k, (a, b) in enumerate(pairs)]
        return np.concatenate([d, reg * (p - c0).ravel()])

    sol = least_squares(resid, c0.ravel(), method="trf")
    p = sol.x.reshape(n, 3)
    return p - p.mean(0) + c0.mean(0)


def contract(comps, bands, assign, centres, factor, targets=None):
    """Rigid cluster translations, linearly interpolated along arclength.

    With `targets`, the clusters are moved to match those pairwise distances
    (shape-aware). Without, every cluster moves toward the centroid by the same
    fraction of its own radius (a homothety, which preserves shape).
    """
    centre = np.mean(centres, axis=0)
    if targets is not None:
        t = target_positions(centres, targets) - np.array(centres)
    else:
        t = np.array([(centre - c) * (1.0 - factor) for c in centres])
    out = []
    for i, c in enumerate(comps):
        n = len(c)
        seg = np.linalg.norm(np.roll(c, -1, 0) - c, axis=1)
        arc = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
        total = float(seg.sum())
        bs = bands.get(i, [])
        if not bs:
            out.append(c.copy())
            continue
        # Every vertex in a band takes that cluster's translation, and the free
        # run between two bands ramps from one to the other. Anchor at the band
        # ENDS, not their centres: anchoring at centres leaves a discontinuity
        # at each band edge, because the first free vertex is interpolated from
        # a point half a band away. Measured, that stretched a single edge from
        # 0.12 to 0.80 and drove minrad to 0.003.
        u = np.full((n, 3), np.nan)
        for k, b in enumerate(bs):
            tk = t[assign[(i, k)]]
            for v in b:
                u[v] = tk
        # walk the free runs between consecutive bands, cyclically
        free = np.where(np.isnan(u[:, 0]))[0]
        if len(free):
            runs, cur = [], [int(free[0])]
            for a, b in zip(free, free[1:]):
                if b - a == 1:
                    cur.append(int(b))
                else:
                    runs.append(cur)
                    cur = [int(b)]
            runs.append(cur)
            if len(runs) > 1 and free[0] == 0 and free[-1] == n - 1:
                runs[0] = runs[-1] + runs[0]
                runs.pop()
            for run in runs:
                before = u[(run[0] - 1) % n]
                after = u[(run[-1] + 1) % n]
                # arclength along the run, handling the cyclic wrap
                d = np.empty(len(run) + 1)
                d[0] = 0.0
                for m in range(len(run)):
                    d[m + 1] = d[m] + seg[(run[m]) % n]
                span = d[-1] if d[-1] > 0 else 1.0
                for m, v in enumerate(run):
                    w = d[m] / span
                    u[v] = (1 - w) * before + w * after
        out.append(c + u)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path, nargs="?")
    ap.add_argument("--factor", type=float)
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--join", type=float, default=3.5,
                    help="Contact bands closer than this (in D) are one cluster.")
    args = ap.parse_args()

    comps = [np.array(c) for c in read_xyz(args.input)]
    rope0 = octrope(args.input)
    total = sum(float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in comps)
    diameter = 2 * total / rope0
    mr0, cl0 = minrad(comps), clearance(comps)

    bands = contact_bands(comps, diameter)
    centres, assign = cluster_sites(comps, bands, diameter, args.join)
    cen = np.mean(centres, axis=0)
    print(f"{len(centres)} clusters of contact sites")
    for k, c in enumerate(centres):
        members = sorted({"ABCDEFGH"[i] for (i, _), g in assign.items() if g == k})
        print(f"   cluster {k+1}: {''.join(members)}   "
              f"radius {np.linalg.norm(c - cen)/diameter:.2f} D from the centroid")
    for a in range(len(centres)):
        for b in range(a + 1, len(centres)):
            print(f"   {a+1}-{b+1}: {np.linalg.norm(centres[a]-centres[b])/diameter:.2f} D")
    print(f"input: ropelength {rope0:.2f}  length {total:.2f}  "
          f"minrad {mr0:.3f}  clearance {cl0:.3f}\n")

    factors = np.arange(0.30, 1.00, 0.05) if args.scan else [args.factor]
    if factors[0] is None:
        print("give --factor or --scan", file=sys.stderr)
        return 1
    best = None
    if args.scan:
        print(f"  {'factor':>7}{'length':>9}{'vs now':>9}{'clear':>8}{'minrad':>8}   verdict")
    for f in factors:
        moved = contract(comps, bands, assign, centres, float(f))
        L = sum(float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in moved)
        cl, mr = clearance(moved), minrad(moved)
        ok = mr >= mr0 * 0.95 and cl >= diameter * 0.95
        if args.scan:
            print(f"  {f:>7.2f}{L:>9.2f}{100*(L/total-1):>8.1f}%{cl:>8.3f}{mr:>8.3f}"
                  f"   {'usable' if ok else 'damaged'}")
        if ok and (best is None or L < best[0]):
            best = (L, float(f))
    if args.output:
        f = args.factor if args.factor else (best[1] if best else None)
        if f is None:
            print("error: no usable factor", file=sys.stderr)
            return 1
        write_xyz(args.output, [c.tolist() for c in
                                contract(comps, bands, assign, centres, float(f))], 9)
        print(f"\nwrote {args.output} at factor {f:.2f}")
    elif best:
        print(f"\n  best usable: factor {best[1]:.2f}, length {best[0]:.2f} "
              f"({100*(best[0]/total-1):+.1f}%)  -> ropelength ~{rope0*best[0]/total:.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
