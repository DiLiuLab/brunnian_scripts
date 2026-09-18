#!/usr/bin/env python3
"""Quicklook renders for a link .xyz: a 4-view PNG and a tube .glb.

Naming follows the existing assets_tightening/quicklooks convention:
    <tag>_<ropelength>[_<label>].png / .glb
Tubes are drawn at the file's own octrope thickness, so a contact really does
look like tube touching tube.
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
    a = ap.parse_args()

    p = Path(a.input)
    comps = [np.asarray(c, float) for c in read_xyz(p)]
    tau, rop = thickness_of(comps, p.with_suffix(".ql.vect"))
    tag = a.tag or p.stem
    stem = f"{tag}_{rop:.2f}" + (f"_{a.label}" if a.label else "")
    out = Path(a.outdir)

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
            V, F = tube_mesh(c, tau)
            m = trimesh.Trimesh(vertices=V, faces=F, process=False)
            col = COLORS[k % len(COLORS)].lstrip("#")
            m.visual.face_colors = [int(col[0:2],16), int(col[2:4],16),
                                    int(col[4:6],16), 255]
            scene.add_geometry(m, node_name=f"component_{k}")
        glb = out / f"{stem}.glb"
        scene.export(glb)
        print(f"  wrote {glb}")
    Path(p.with_suffix(".ql.vect")).unlink(missing_ok=True)

if __name__ == "__main__":
    main()
