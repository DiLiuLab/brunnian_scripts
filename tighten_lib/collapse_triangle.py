#!/usr/bin/env python3
"""Contract a three-cluster (triangular) link radially about its centroid.

The axial squeeze in compress_gap.py assumes TWO clusters with a gap between
them, and expresses the move as a monotone map along the line joining them.
That model breaks on a structure whose clusters form a triangle: squeeze along
any one edge and the third corner projects to the middle of the compressed
slab, so it gets crushed along with the gap. There is no direction along which
two clusters separate and the third stays clear.

Measured on a 4_LC configuration at 187.54: three clusters of grip sites with
sides 6.64, 6.73, 6.76 tube diameters, corner radii 3.85/3.86/3.91 from the
centroid -- near-equilateral. Roughly 58% of the link's length sits in the three
edges. Note the edges are NOT all made of the bridge: two are two strands of
the bridge each, and the third is four strands, two each of two OTHER
components that span between corners. Clustering components by proximity misses
this entirely, because those two components touch at both corners and so look
like a single cluster; the thing to cluster is the GRIP SITES.

The move here is radial about the centroid:

    r -> f(r),   f(r) = factor * r              for r <= R1
                 f(r) = factor * R1 + (r - R1)  for r >  R1

with R1 the inner edge of the cluster shell. f is strictly increasing, so the
map is a homeomorphism of R^3 and THE LINK TYPE IS PRESERVED EXACTLY, as with
the axial squeeze. Beyond R1 the map is a radial translation, which preserves a
cluster's RADIAL extent but compresses it TANGENTIALLY by f(R)/R. That
tangential distortion is the price the axial version never had to pay, and it
is what the --scan option measures: it reports minrad and clearance so a factor
can be chosen that RidgeRunner can actually repair.

Usage:
    python3 collapse_triangle.py in.xyz --scan
    python3 collapse_triangle.py in.xyz out.xyz --factor 0.6
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
from tighten_link_xyz import read_xyz, write_xyz, write_vect  # noqa: E402


def octrope(path: Path) -> float:
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / "x.vect"
        write_vect(vect, read_xyz(path))
        return float(subprocess.run(["ropelength", "-q", str(vect)],
                                    capture_output=True, text=True, timeout=300).stdout)


def minrad(comps) -> float:
    """octrope's MinRad: min(|e_prev|,|e_next|) / (2 tan(theta/2)) over all vertices."""
    worst = np.inf
    for c in comps:
        c = np.asarray(c)
        e = np.roll(c, -1, 0) - c
        ln = np.linalg.norm(e, axis=1)
        u = e / ln[:, None]
        cos = np.clip((u * np.roll(u, 1, 0)).sum(1), -1, 1)
        t = np.tan(np.maximum(np.arccos(cos), 1e-12) / 2)
        worst = min(worst, (np.minimum(np.roll(ln, 1), ln) / (2 * t)).min())
    return float(worst)


def clearance(comps) -> float:
    comps = [np.asarray(c) for c in comps]
    return min(float(np.linalg.norm(comps[i][:, None, :] - comps[j][None, :, :], axis=2).min())
               for i in range(len(comps)) for j in range(i + 1, len(comps)))


def grip_clusters(comps, diameter, contact=1.05, join=3.5):
    """Cluster the CONTACT SITES between components, not the components.

    A component that spans the structure touches at two corners; clustering
    components would merge those corners into one. Grip sites do not have that
    problem -- each is local by construction.
    """
    sites = []
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            d = np.linalg.norm(comps[i][:, None, :] - comps[j][None, :, :], axis=2)
            close = np.argwhere(d / diameter < contact)
            for a, b in close:
                sites.append(0.5 * (comps[i][a] + comps[j][b]))
    if not sites:
        return []
    sites = np.array(sites)
    groups = [[sites[0]]]
    for p in sites[1:]:
        for g in groups:
            if min(np.linalg.norm(p - q) for q in g) / diameter < join:
                g.append(p)
                break
        else:
            groups.append([p])
    return [np.array(g).mean(0) for g in groups]


def build(comps, centre, factor, r1):
    """Radial map. Kept for comparison; it crushes the tube -- see build_rigid."""
    out = []
    for c in comps:
        v = c - centre
        r = np.linalg.norm(v, axis=1)
        fr = np.where(r <= r1, factor * r, factor * r1 + (r - r1))
        scale = np.divide(fr, r, out=np.ones_like(r), where=r > 1e-12)
        out.append(centre + v * scale[:, None])
    return out


def build_rigid(comps, corners, centre, factor, sigma):
    """Translate each cluster rigidly toward the centroid, blending in between.

    The radial map scales every separation inside the cluster shell, so the
    four parallel strands along an edge -- about 1 D apart -- are squeezed into
    each other. Measured: clearance tracked the factor almost exactly (0.50 ->
    0.569, 0.70 -> 0.744), which is uniform scaling, not gap removal.

    Here each cluster gets a single rigid translation t_k toward the centroid,
    and every curve point is displaced by a weighted average of the three,
    weights falling off as exp(-(d/sigma)^2) from each cluster centre. Inside a
    cluster one weight dominates, so the cluster translates without distortion;
    along an edge the two nearby clusters share the weight, so the edge is
    carried inward and shortens because its ENDS converge. Two strands running
    side by side on the same edge sit at almost the same place, so they receive
    almost the same displacement and their separation survives.

    Unlike the axial squeeze, this is NOT guaranteed to be a homeomorphism --
    a large enough translation can fold the displacement field. The link type
    must be checked with HOMFLY afterwards rather than assumed.
    """
    t = np.array([(centre - c) * (1.0 - factor) for c in corners])
    out = []
    for c in comps:
        d = np.stack([np.linalg.norm(c - k, axis=1) for k in corners])   # (3, n)
        w = np.exp(-(d / sigma) ** 2)
        w /= w.sum(0, keepdims=True)
        out.append(c + (w.T @ t))
    return out


def build_arclength(comps, corners, centre, factor, sigma_s, diameter, contact=1.05):
    """Rigid cluster translations, blended along each curve's own ARCLENGTH.

    Blending by spatial distance bends every strand that crosses the transition
    region: the displacement field varies across the strand, so the strand
    acquires curvature. Measured on the 187.54 configuration, that drove minrad
    to 0.000 for every blend width under 3 D, whatever the factor.

    Blending along arclength instead makes the displacement vary ALONG each
    strand rather than across it. A strand running from cluster i to cluster j
    has its displacement interpolate from t_i to t_j over its own length, which
    shortens it smoothly -- a reparametrisation, not a bend. Two strands lying
    side by side on the same edge are displaced independently but identically,
    so their separation survives as well.

    Weights come from arclength distance to the nearest grip of each cluster,
    measured cyclically along the component. A component with no grip in a
    given cluster gets no weight there.
    """
    # which grips belong to which cluster, per component
    sites = {i: [] for i in range(len(comps))}
    for i in range(len(comps)):
        for j in range(i + 1, len(comps)):
            d = np.linalg.norm(comps[i][:, None, :] - comps[j][None, :, :], axis=2)
            for a, b in np.argwhere(d / diameter < contact):
                mid = 0.5 * (comps[i][a] + comps[j][b])
                k = int(np.argmin([np.linalg.norm(mid - c) for c in corners]))
                sites[i].append((int(a), k))
                sites[j].append((int(b), k))

    t = np.array([(centre - c) * (1.0 - factor) for c in corners])
    out = []
    for i, c in enumerate(comps):
        n = len(c)
        seg = np.linalg.norm(np.roll(c, -1, 0) - c, axis=1)
        arc = np.concatenate([[0.0], np.cumsum(seg)[:-1]])
        total = seg.sum()
        w = np.zeros((len(corners), n))
        for k in range(len(corners)):
            idx = [a for a, kk in sites[i] if kk == k]
            if not idx:
                continue
            # cyclic arclength distance to the nearest grip of this cluster
            dist = np.full(n, np.inf)
            for a in idx:
                d = np.abs(arc - arc[a])
                dist = np.minimum(dist, np.minimum(d, total - d))
            w[k] = np.exp(-(dist / sigma_s) ** 2)
        col = w.sum(0)
        if (col <= 1e-12).any():          # a component touching nothing: leave it
            out.append(c.copy())
            continue
        w /= col
        out.append(c + (w.T @ t))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path, nargs="?")
    ap.add_argument("--factor", type=float, help="Interior scale; 1 leaves it alone.")
    ap.add_argument("--scan", action="store_true", help="Report the cost of a range of factors.")
    ap.add_argument("--mode", choices=("arc", "rigid", "radial"), default="arc",
                    help="arc: rigid cluster translations blended along arclength (default).\n"
                         "rigid: the same, blended by spatial distance -- bends the strands.\n"
                         "radial: scale the interior about the centroid -- crushes the tube.")
    ap.add_argument("--sigma", type=float, default=None,
                    help="Blend width for rigid mode, in tube diameters. "
                         "Default: 0.55 x the mean side.")
    args = ap.parse_args()

    comps = [np.array(c) for c in read_xyz(args.input)]
    rope0 = octrope(args.input)
    total = sum(float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in comps)
    diameter = 2 * total / rope0
    mr0, cl0 = minrad(comps), clearance(comps)

    corners = grip_clusters(comps, diameter)
    if len(corners) < 3:
        print(f"error: found {len(corners)} clusters of grip sites, not 3. "
              "Use compress_gap.py for the two-cluster case.", file=sys.stderr)
        return 1
    corners = np.array(corners)
    centre = corners.mean(0)
    radii = np.linalg.norm(corners - centre, axis=1)
    # inner edge of the cluster shell: how close any curve point in a cluster gets
    r1 = float(min(radii) - 0.5 * diameter)

    print(f"{len(corners)} clusters of grip sites, "
          f"radii {'/'.join(f'{r/diameter:.2f}' for r in radii)} D from the centroid")
    for i in range(len(corners)):
        for j in range(i + 1, len(corners)):
            print(f"   side {i+1}-{j+1}: {np.linalg.norm(corners[i]-corners[j])/diameter:.2f} D")
    print(f"input: ropelength {rope0:.2f}  length {total:.2f}  "
          f"minrad {mr0:.3f}  clearance {cl0:.3f}  D {diameter:.3f}\n")

    sides = [np.linalg.norm(corners[i] - corners[j]) / diameter
             for i in range(len(corners)) for j in range(i + 1, len(corners))]
    sigma = (args.sigma if args.sigma else 0.55 * float(np.mean(sides))) * diameter

    def make(f):
        if args.mode == "arc":
            return build_arclength(comps, corners, centre, float(f), sigma, diameter)
        if args.mode == "rigid":
            return build_rigid(comps, corners, centre, float(f), sigma)
        return build(comps, centre, float(f), r1)

    factors = (np.arange(0.30, 1.00, 0.05) if args.scan
               else [args.factor] if args.factor else [])
    if not factors.__len__():
        print("give --factor or --scan", file=sys.stderr)
        return 1

    if args.scan:
        print(f"  mode {args.mode}, sigma {sigma/diameter:.2f} D")
        print(f"  {'factor':>7}{'length':>9}{'vs now':>9}{'side':>8}{'clear':>8}{'minrad':>8}   verdict")
    best = None
    for f in factors:
        moved = make(f)
        L = sum(float(np.linalg.norm(np.roll(c, -1, 0) - c, axis=1).sum()) for c in moved)
        if args.mode == "rigid":
            side = np.linalg.norm((corners[0] + (centre-corners[0])*(1-f))
                                  - (corners[1] + (centre-corners[1])*(1-f))) / diameter
        else:
            side = np.linalg.norm(
                (centre + (corners[0]-centre) * ((f*r1 + (radii[0]-r1))/radii[0]))
                - (centre + (corners[1]-centre) * ((f*r1 + (radii[1]-r1))/radii[1]))) / diameter
        cl, mr = clearance(moved), minrad(moved)
        ok = mr >= mr0 * 0.95 and cl >= diameter * 0.90
        if args.scan:
            print(f"  {f:>7.2f}{L:>9.2f}{100*(L/total-1):>8.1f}%{side:>8.2f}{cl:>8.3f}"
                  f"{mr:>8.3f}   {'usable' if ok else 'damaged'}")
        if ok and (best is None or L < best[0]):
            best = (L, float(f), moved)

    if args.output:
        f = args.factor if args.factor else (best[1] if best else None)
        if f is None:
            print("error: no usable factor found", file=sys.stderr)
            return 1
        moved = make(f)
        write_xyz(args.output, [c.tolist() for c in moved], 9)
        print(f"\nwrote {args.output} at factor {f:.2f}")
    elif args.scan and best:
        print(f"\n  best usable: factor {best[1]:.2f}, length {best[0]:.2f} "
              f"({100*(best[0]/total-1):+.1f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
