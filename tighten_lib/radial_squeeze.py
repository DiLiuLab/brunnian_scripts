#!/usr/bin/env python3
"""Radial gap squeeze: close the hollow centre of an annular link.

The axial squeeze in compress_gap.py is topology-safe because it is
    Phi(t, y) = (f(t), y),   f strictly increasing
a homeomorphism of R^3 that leaves the transverse projection fixed. Cluster
contraction has no such guarantee -- it moves different pieces in different
directions -- which is why it needs a HOMFLY check every round, and why its
damage lands on minRad (buckled connecting strands) rather than on clearance.

For a link whose slack is the HOLE IN THE MIDDLE of a ring, the same trick
works in cylindrical coordinates about the symmetry axis:

    Phi(r, theta, z) = (g(r), theta, z),    g strictly increasing, g(0) = 0

g strictly increasing and g(0)=0 makes Phi a homeomorphism of R^3, and
g_s = (1-s)*id + s*g is a strictly increasing interpolation for every s, so
Phi is isotopic to the identity. THE LINK TYPE IS THEREFORE PRESERVED EXACTLY,
with no topology check needed on the output -- the same standing the axial
squeeze has, and which cluster contraction does not.

Two further properties for free:

  * It is theta-independent, so it commutes with every rotation about the axis.
    A Cn-symmetric link stays EXACTLY Cn -- no re-symmetrization step.
  * The damage is to CLEARANCE, not curvature. Choosing slope 1 across the
    link's own radial band makes the map a pure inward translation in r there,
    so a circumferential arc's curvature changes only from 1/r to 1/(r-delta)
    -- smooth and mild. Contrast cluster contraction, whose rigid translations
    plus arclength interpolation buckle the joining strands.

g is built by integrating a smoothstep-ramped derivative, so g is C^1 and the
join at the link's inner radius is not a kink (the lesson of compress_gap.py's
--blend: curvature, not clearance, was the binding constraint).

Usage:
    python3 radial_squeeze.py in.xyz out.xyz --factor 0.6 [--axis 0,0,1]
    python3 radial_squeeze.py in.xyz --scan
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked out.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_xyz, write_vect  # noqa: E402


def smoothstep(x):
    x = np.clip(x, 0.0, 1.0)
    return x * x * (3 - 2 * x)


def build_g(r0, factor, blend, rmax):
    """g' ramps from `factor` inside the hole to 1 across the link's band."""
    grid = np.linspace(0.0, rmax * 1.5 + 1.0, 20000)
    if blend <= 0:
        dg = np.where(grid < r0, factor, 1.0)
    else:
        dg = factor + (1.0 - factor) * smoothstep((grid - (r0 - blend)) / (2 * blend))
    g = np.concatenate([[0.0], np.cumsum(0.5 * (dg[1:] + dg[:-1]) * np.diff(grid))])
    return grid, g


def apply_squeeze(comps, axis, factor, blend):
    a = np.asarray(axis, float)
    a /= np.linalg.norm(a)
    allv = np.vstack(comps)
    centre = allv.mean(axis=0)

    # orthonormal frame with `a` as the axis
    tmp = np.array([1.0, 0.0, 0.0])
    if abs(tmp @ a) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = tmp - (tmp @ a) * a
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)

    def cyl(P):
        q = P - centre
        return q @ e1, q @ e2, q @ a

    xs, ys, _ = cyl(allv)
    rr = np.hypot(xs, ys)
    r0, rmax = rr.min(), rr.max()
    grid, g = build_g(r0, factor, blend, rmax)

    out = []
    for c in comps:
        x, y, z = cyl(c)
        r = np.hypot(x, y)
        th = np.arctan2(y, x)
        rn = np.interp(r, grid, g)
        out.append(centre + rn[:, None] * (np.cos(th)[:, None] * e1
                                           + np.sin(th)[:, None] * e2)
                   + z[:, None] * a)
    return out, r0, rmax


def apply_squeeze_about(comps, point, axis, factor, blend):
    """The squeeze about an EXPLICIT line (point + direction), for irregular
    structures where the axis is found rather than given. Same map, same
    guarantee: monotone in r about a line the link avoids -> homeomorphism."""
    ax = np.asarray(axis, float); ax = ax / np.linalg.norm(ax)
    point = np.asarray(point, float)
    allv = np.vstack([np.asarray(c, float) for c in comps]) - point
    rr = np.linalg.norm(np.cross(allv, ax), axis=1)
    grid, g = build_g(rr.min(), factor, blend, rr.max())
    out = []
    for c in comps:
        q = np.asarray(c, float) - point
        t = q @ ax
        rad = q - np.outer(t, ax)
        r = np.linalg.norm(rad, axis=1)
        rn = np.interp(r, grid, g)
        scale = np.where(r > 1e-12, rn / np.maximum(r, 1e-12), 1.0)
        out.append(point + np.outer(t, ax) + rad * scale[:, None])
    return out


def octrope(comps, tmp: Path):
    import subprocess, re
    write_vect(tmp, [np.asarray(c, float).tolist() for c in comps])
    t = subprocess.run(["ropelength", str(tmp)], capture_output=True, text=True)
    txt = t.stdout + t.stderr
    grab = lambda k: float(__import__("re").search(
        rf"(?<![A-Za-z] ){k}:\s*([0-9.eE+-]+)", txt).group(1))
    return grab("Ropelength"), grab("Thickness"), grab("minRad"), grab("minStrut")


def length_of(comps):
    return sum(float(np.linalg.norm(
        np.diff(np.vstack([np.asarray(c, float), np.asarray(c, float)[:1]]), axis=0),
        axis=1).sum()) for c in comps)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path, nargs="?")
    ap.add_argument("--factor", type=float, default=None,
                    help="slope of g inside the hole; 1.0 leaves it, 0.3 keeps 30%%")
    ap.add_argument("--blend", type=float, default=1.0,
                    help="smoothing width in the link's own units (D=1 after "
                         "normalisation); 0 gives the kinked piecewise-linear map")
    ap.add_argument("--axis", default="0,0,1")
    ap.add_argument("--scan", action="store_true")
    a = ap.parse_args()

    axis = [float(v) for v in a.axis.split(",")]
    comps = [np.asarray(c, float) for c in read_xyz(a.input)]
    tmp = a.input.with_suffix(".rs.vect")
    try:
        rop0, tau0, mr0, ms0 = octrope(comps, tmp)
        L0 = length_of(comps)
        _, r0, rmax = apply_squeeze(comps, axis, 1.0, a.blend)
        print(f"input: ropelength {rop0:.4f}  length {L0:.4f}  minRad {mr0:.4f}  "
              f"minStrut {ms0:.4f}")
        print(f"       inner radius {r0:.4f}  outer {rmax:.4f}  "
              f"hole = {r0 / (2 * tau0):.2f} D across the empty core\n")

        if a.scan:
            print(f"{'factor':>7}{'length':>10}{'vs now':>9}{'minRad':>9}"
                  f"{'minStrut':>10}{'Rop_est':>10}  verdict")
            for f in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
                sq, _, _ = apply_squeeze(comps, axis, f, a.blend)
                rop, tau, mr, ms = octrope(sq, tmp)
                L = length_of(sq)
                t = min(ms / 2, mr)
                print(f"{f:7.2f}{L:10.3f}{(L - L0) / L0 * 100:+8.1f}%{mr:9.4f}"
                      f"{ms:10.4f}{L / t:10.3f}  "
                      f"{'strut-limited' if mr > ms / 2 else 'CURVATURE-limited'}")
            return 0

        if a.factor is None:
            ap.error("give --factor or --scan")
        sq, _, _ = apply_squeeze(comps, axis, a.factor, a.blend)
        rop, tau, mr, ms = octrope(sq, tmp)
        write_xyz(a.output, [c.tolist() for c in sq], 12)
        print(f"wrote {a.output}")
        print(f"  ropelength {rop:.4f}  length {length_of(sq):.4f}  "
              f"minRad {mr:.4f}  minStrut {ms:.4f}")
        print("  topology preserved by construction (monotone radial map); no "
              "HOMFLY check required on THIS step")
    finally:
        tmp.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
