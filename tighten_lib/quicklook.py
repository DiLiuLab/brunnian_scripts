#!/usr/bin/env python3
"""Quicklook renders for a link .xyz: a 4-view PNG and a tube .glb.

Naming follows the existing assets_tightening/quicklooks convention:
    <tag>_<ropelength>[_<label>].png / .glb
Tubes are drawn at the file's own octrope thickness by default, so a contact
really does look like tube touching tube.  --rod-diameter overrides that with
an explicit tube diameter in the file's own units: a value below 2*tau opens
visible gaps between strands, which reads better for showing the weave.
"""
import subprocess, sys, re, argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked out.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_vect, write_xyz

COLORS = ["#2b6cb0", "#c05621", "#2f855a", "#b83280", "#6b46c1",
          "#0987a0", "#975a16", "#4a5568"]

def thickness_of(comps, tmp):
    write_vect(tmp, [c.tolist() for c in comps])
    out = subprocess.run(["ropelength", str(tmp)], capture_output=True, text=True)
    txt = out.stdout + out.stderr
    tau = re.search(r"(?<![A-Za-z] )Thickness:\s*([0-9.eE+-]+)", txt)
    rop = re.search(r"(?<![A-Za-z] )Ropelength:\s*([0-9.eE+-]+)", txt)
    return (float(tau.group(1)) if tau else None,
            float(rop.group(1)) if rop else None)

def tube_mesh(curve, radius, nseg=14):
    """Closed tube around a closed polyline. Returns (vertices, faces)."""
    P = np.asarray(curve, float)
    N = len(P)
    T = np.roll(P, -1, axis=0) - np.roll(P, 1, axis=0)
    T /= np.linalg.norm(T, axis=1)[:, None]
    # parallel transport a normal around the loop
    ref = np.array([0.0, 0.0, 1.0])
    if abs(T[0] @ ref) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    n0 = ref - (ref @ T[0]) * T[0]
    n0 /= np.linalg.norm(n0)
    Ns = np.empty_like(P); Ns[0] = n0
    for i in range(1, N):
        n = Ns[i-1] - (Ns[i-1] @ T[i]) * T[i]
        nn = np.linalg.norm(n)
        Ns[i] = n / nn if nn > 1e-12 else np.cross(T[i], ref)
    # Parallel transport does not close up: carried around a closed loop the
    # normal comes back rotated by a holonomy angle (up to ~56 deg on these
    # links). Left alone, the ring at vertex N-1 joins the ring at vertex 0
    # twisted by that angle and the wrap-around quads shear into a visible
    # seam, even though the mesh is watertight. Spread the correction evenly
    # over the loop so the frame closes; the cost is a uniform torsion of
    # theta/N per edge, the minimal-twist closure.
    #
    # One pass leaves a small residual because rotating about the tangent does
    # not commute exactly with the projection step when T turns, so iterate.
    # It converges geometrically -- three passes reach machine precision on
    # these links.
    base, phi = Ns.copy(), np.zeros(N)
    for _ in range(12):
        Ns = (base * np.cos(phi)[:, None]
              + np.cross(T, base) * np.sin(phi)[:, None])
        Ns /= np.linalg.norm(Ns, axis=1)[:, None]
        n_end = Ns[N-1] - (Ns[N-1] @ T[0]) * T[0]
        n_end /= np.linalg.norm(n_end)
        theta = np.arctan2(n_end @ np.cross(T[0], Ns[0]), n_end @ Ns[0])
        if abs(theta) < 1e-12:
            break
        phi = phi - theta * np.arange(N) / N

    B = np.cross(T, Ns)
    ang = np.linspace(0, 2*np.pi, nseg, endpoint=False)
    ring = (np.cos(ang)[None, :, None] * Ns[:, None, :] +
            np.sin(ang)[None, :, None] * B[:, None, :])
    V = (P[:, None, :] + radius * ring[0]).reshape(-1, 3) if ring.shape[0] == 1 \
        else (P[:, None, :] + radius * (np.cos(ang)[None, :, None]*Ns[:, None, :]
                                        + np.sin(ang)[None, :, None]*B[:, None, :])).reshape(-1, 3)
    F = []
    for i in range(N):
        j = (i + 1) % N
        for k in range(nseg):
            k2 = (k + 1) % nseg
            a, b = i*nseg + k, i*nseg + k2
            c, d = j*nseg + k, j*nseg + k2
            F.append([a, c, d]); F.append([a, d, b])
    return V, np.asarray(F)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--tag", default=None, help="name prefix, e.g. 5BL")
    ap.add_argument("--label", default="", help="suffix after the ropelength")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--glb", action="store_true")
    ap.add_argument("--rod-diameter", type=float, default=None,
                    help="tube DIAMETER in the file's own units; default is the "
                         "measured 2*tau, i.e. rope at true thickness")
    ap.add_argument("--glb-only", action="store_true",
                    help="write only the .glb, no PNG and no .xyz copy")
    ap.add_argument("--stem", default=None,
                    help="exact output stem, overriding the <tag>_<ropelength> convention")
    a = ap.parse_args()

    p = Path(a.input)
    comps = [np.asarray(c, float) for c in read_xyz(p)]
    tau, rop = thickness_of(comps, p.with_suffix(".ql.vect"))
    tag = a.tag or p.stem
    stem = a.stem or (f"{tag}_{rop:.2f}" + (f"_{a.label}" if a.label else ""))
    radius = a.rod_diameter / 2 if a.rod_diameter is not None else tau
    out = Path(a.outdir)

    if not a.glb_only:
        allv = np.vstack(comps)
        lo, hi = allv.min(axis=0), allv.max(axis=0)
        ctr = (lo + hi) / 2
        ext = (hi - lo) * 1.06
        ext = np.maximum(ext, ext.max() * 0.06)   # keep a flat axis from collapsing

        views = [(90, -90, "down the symmetry axis  (+z)"), (0, -90, "side  (-y)"),
                 (0, 0, "side  (+x)"), (24, -54, "oblique")]
        fig = plt.figure(figsize=(13, 9), dpi=130)
        for i, (elev, azim, title) in enumerate(views, 1):
            ax = fig.add_subplot(2, 2, i, projection="3d")
            for k, c in enumerate(comps):
                L = np.vstack([c, c[:1]])
                ax.plot(L[:, 0], L[:, 1], L[:, 2], color=COLORS[k % len(COLORS)],
                        lw=2.6, solid_capstyle="round")
            ax.set_xlim(ctr[0]-ext[0]/2, ctr[0]+ext[0]/2)
            ax.set_ylim(ctr[1]-ext[1]/2, ctr[1]+ext[1]/2)
            ax.set_zlim(ctr[2]-ext[2]/2, ctr[2]+ext[2]/2)
            # box aspect proportional to the true extents => equal scale on every
            # axis, no distortion, and a flat structure does not waste the panel
            ax.set_box_aspect(tuple(ext / ext.max()))
            ax.view_init(elev=elev, azim=azim)
            ax.set_title(title, fontsize=10, color="#444")
            ax.set_axis_off()
        fig.suptitle(f"{tag}   ropelength {rop:.3f}   tau {tau:.4f}   "
                     f"{len(comps)} components, {sum(len(c) for c in comps)} vertices"
                     + (f"\n{a.label}" if a.label else ""),
                     fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        png = out / f"{stem}.png"
        fig.savefig(png, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"  wrote {png}")

        # Keep the geometry beside its renders, under the SAME stem. The filename
        # carries the ropelength, which is how these get identified later, and the
        # png is only a preview -- the xyz is what feeds the next step.
        xyz = out / f"{stem}.xyz"
        write_xyz(xyz, [c.tolist() for c in comps], 12)
        print(f"  wrote {xyz}")

    if a.glb:
        import trimesh
        scene = trimesh.Scene()
        for k, c in enumerate(comps):
            V, F = tube_mesh(c, radius)
            m = trimesh.Trimesh(vertices=V, faces=F, process=False)
            m.fix_normals()                      # consistent outward winding
            m.vertex_normals                     # smooth normals, so the tube
                                                 # shades as a tube, not facets
            col = COLORS[k % len(COLORS)].lstrip("#")
            m.visual.face_colors = [int(col[0:2],16), int(col[2:4],16),
                                    int(col[4:6],16), 255]
            scene.add_geometry(m, node_name=f"component_{k}")
        glb = out / f"{stem}.glb"
        scene.export(glb)
        print(f"  wrote {glb}  (rod diameter {2*radius:.4f}"
              + (f", true thickness)" if a.rod_diameter is None
                 else f", true thickness is {2*tau:.4f})"))
    Path(p.with_suffix(".ql.vect")).unlink(missing_ok=True)

if __name__ == "__main__":
    main()
