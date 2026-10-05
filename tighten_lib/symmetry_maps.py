#!/usr/bin/env python3
"""Exact vertex maps and the Reynolds operator for a finite point group acting on a link.

    python3 tighten_lib/symmetry_maps.py LINK.xyz D2d          # report the maps, or why they fail
    python3 tighten_lib/symmetry_maps.py LINK.xyz C5 --axis 0,0,1 --ref 1,0,0
    python3 tighten_lib/symmetry_maps.py LINK.xyz Cn --order 7
    python3 tighten_lib/symmetry_maps.py LINK.xyz detect       # every group that holds EXACTLY

A library with no side effects; this __main__ is only a diagnostic. It is what
`slp_tighten.py` uses to move a link inside its exact symmetry subspace, and it
is meant to be the one place a vertex correspondence is built.

GROUPS, from named generators (canonical frame: principal axis z, reference x;
--axis/--ref conjugate every generator by the rotation taking z->axis, x->ref):

    C1, none, E        identity
    Cs, mirror         reflection in z            (ridgerunner's "D2" is THIS group)
    Ci                 inversion -I
    C<n>, Z/<n>Z       rotation 2pi/n about z      ("Cn"/"Z/nZ" with --order)
    C<n>v, Cpv         + mirror with normal x      (ridgerunner Cpv)
    C<n>h              + mirror with normal z
    RD<n>, D<n> (n>=3) + half turn about x         (ridgerunner RDp: the real dihedral group)
    D<n>d[:xy|:diag]   S_2n about z + mirror normal x (xy) or rotated by pi/(2n) (diag);
                       D2d:xy is exactly d2d_symmetrize's group
    S<2n>              rotation pi/n about z times the z mirror

A bare D2 is refused: ridgerunner's D2 is a single mirror (Cs); write Cs, or
RD2 for the order-4 dihedral group.

CONVENTION, shared with symmetrize_link_xyz.py, d2d_symmetrize.py and the slp
probes. X = vstack(components) in read_xyz order.

    perms[g][i] = j   <=>   M_g @ X[i] == X[j]
    composition       pi_{gh} = pi_g[pi_h]
    Reynolds          (R V)[i] = 1/|G| sum_g V[pi_g[i]] @ M_g      (row vectors)
    mats[0] = I, perms[0] = arange(N); rely on no other ordering.

REFUSES RATHER THAN GUESSES. A map is built only if every vertex lands within
`tol` mean edges (default 1e-6) of a vertex, the nearest-vertex map is a
bijection, it sends whole components to whole components, every component's
index map is a shift or a reversal (so edges go to edges), and the composed maps
satisfy every relation of the group -- a check symmetrize_link_xyz.py does not
make on its live path. A failure says which test failed and what to do: project
first, centre the link, or resample an orbit to one count (equivariant_resample).

The link is NOT centred here. Callers centre it (X -= X.mean(0)) first; the
diagnostic below does so unless --no-centre is given.
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz  # noqa: E402
from symmetrize_link_xyz import align_to, reflection_matrix, rotation_matrix  # noqa: E402

D2_MESSAGE = ("ridgerunner's D2 is a single mirror (Cs); write Cs, or RD2 for the "
              "order-4 dihedral group")


class SymmetryMapError(ValueError):
    """The requested group does not act exactly on this link's vertices (or the spec is bad)."""


# --------------------------------------------------------------------------- #
# group specifications
# --------------------------------------------------------------------------- #

@dataclass(frozen=True)
class GroupSpec:
    kind: str          # C1 Cs Ci Cn Cnv Cnh RDn Dnd S2n detect
    n: int = 1
    frame: str = "xy"  # Dnd only

    @property
    def label(self) -> str:
        k, n = self.kind, self.n
        if k in ("C1", "Cs", "Ci", "detect"):
            return k
        if k == "Cn":
            return f"C{n}"
        if k == "Cnv":
            return f"C{n}v"
        if k == "Cnh":
            return f"C{n}h"
        if k == "RDn":
            return f"RD{n}"
        if k == "S2n":
            return f"S{2 * n}"
        if k == "Dnd":
            return f"D{n}d" + (":diag" if self.frame == "diag" else "")
        return k

    @property
    def order(self) -> int:
        k, n = self.kind, self.n
        return {"C1": 1, "Cs": 2, "Ci": 2, "Cn": n, "Cnv": 2 * n, "Cnh": 2 * n,
                "RDn": 2 * n, "S2n": 2 * n, "Dnd": 4 * n}.get(k, 0)


def _count(token: str, order: int | None, spec: str) -> int:
    if token in ("n", "p"):
        if order is None:
            raise SymmetryMapError(f"{spec!r} needs an order (e.g. order=5, or write C5 / C5v)")
        return int(order)
    n = int(token)
    if order is not None and int(order) != n:
        raise SymmetryMapError(f"{spec!r} says n={n} but order={order} was also given")
    return n


def parse_spec(spec, order: int | None = None) -> GroupSpec:
    """Parse a group name into GroupSpec(kind, n, frame). See the module docstring."""
    if isinstance(spec, GroupSpec):
        return spec
    s = str(spec).strip()
    low = s.lower()
    if low in ("c1", "none", "e", "identity"):
        return GroupSpec("C1")
    if low in ("cs", "mirror"):
        return GroupSpec("Cs")
    if low == "ci":
        return GroupSpec("Ci")
    if low == "detect":
        return GroupSpec("detect")
    if low == "d2":
        raise SymmetryMapError(D2_MESSAGE)
    if low == "cplanes":
        raise SymmetryMapError("ridgerunner's cplanes is not supported here; name the point "
                               "group it should be (e.g. C2v, Cs, Ci)")
    m = re.fullmatch(r"z/(\d+|[np])z", low)
    if m:
        n = _count(m.group(1), order, s)
        if n < 1:
            raise SymmetryMapError(f"{s!r}: the order must be at least 1")
        return GroupSpec("C1") if n == 1 else GroupSpec("Cn", n)
    m = re.fullmatch(r"(rd|c|d|s)(\d+|[np])(v|h|d)?(?::(xy|diag))?", low)
    if not m:
        raise SymmetryMapError(
            f"unknown group {s!r}: use C1, Cs, Ci, C<n>, Z/<n>Z, C<n>v, C<n>h, RD<n>, D<n> (n>=3), "
            "D<n>d[:xy|:diag], S<2n> or detect")
    head, num, tail, frame = m.groups()
    n = _count(num, order, s)
    if n < 1:
        raise SymmetryMapError(f"{s!r}: n must be at least 1")
    if frame is not None and not (head == "d" and tail == "d"):
        raise SymmetryMapError(f"{s!r}: the :xy/:diag frame suffix applies to D<n>d only")
    if head == "c":
        if tail is None:
            return GroupSpec("C1") if n == 1 else GroupSpec("Cn", n)
        if tail == "v":
            return GroupSpec("Cnv", n)
        if tail == "h":
            return GroupSpec("Cnh", n)
        raise SymmetryMapError(f"{s!r}: there is no C<n>d; did you mean D<n>d?")
    if head == "rd":
        if tail is not None:
            raise SymmetryMapError(f"{s!r}: RD<n> takes no suffix")
        return GroupSpec("RDn", n)
    if head == "d":
        if tail is None:
            if n == 2:
                raise SymmetryMapError(D2_MESSAGE)
            if n < 3:
                raise SymmetryMapError(f"{s!r}: write RD{n} for the dihedral group of order {2 * n}")
            return GroupSpec("RDn", n)
        if tail == "d":
            if n < 2:
                raise SymmetryMapError(f"{s!r}: D<n>d needs n >= 2")
            return GroupSpec("Dnd", n, frame or "xy")
        raise SymmetryMapError(f"{s!r}: D<n>{tail} is not supported (supported: D<n>d, RD<n>)")
    # head == "s"
    if tail is not None:
        raise SymmetryMapError(f"{s!r}: S<2n> takes no suffix")
    if n % 2:
        raise SymmetryMapError(f"{s!r}: only even improper rotations S<2n> are supported "
                               f"(S{n} with n odd is C{n}h)")
    return GroupSpec("S2n", n // 2)


# --------------------------------------------------------------------------- #
# generators and closure
# --------------------------------------------------------------------------- #

_SNAP = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])


def _snap(M: np.ndarray) -> np.ndarray:
    """Make entries within 1e-12 of 0, +-1/2, +-1 exact (rot(z, pi/2) has cos = 6e-17)."""
    M = np.array(M, float)
    d = np.abs(M[..., None] - _SNAP)
    k = d.argmin(-1)
    hit = d.min(-1) < 1e-12
    M[hit] = _SNAP[k[hit]]
    M[M == 0] = 0.0          # no negative zeros
    return M


def frame_rotation(axis=(0, 0, 1), ref=(1, 0, 0), strict: bool = True) -> np.ndarray:
    """Rotation A with A z = axis and A x = ref projected perpendicular to the axis.

    With strict=False a ref parallel to the axis is replaced by x, or y if x is parallel too
    (for groups whose generators do not depend on the reference direction)."""
    a = np.asarray(axis, float)
    if np.linalg.norm(a) < 1e-12:
        raise SymmetryMapError("the symmetry axis is the zero vector")
    a = a / np.linalg.norm(a)
    r = np.asarray(ref, float)
    rp = r - (r @ a) * a
    if np.linalg.norm(rp) < 1e-9 and not strict:
        for alt in ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0)):
            rp = np.asarray(alt) - (np.asarray(alt) @ a) * a
            if np.linalg.norm(rp) > 1e-3:
                break
    if np.linalg.norm(rp) < 1e-9:
        raise SymmetryMapError("the reference direction is parallel to the axis")
    rp = rp / np.linalg.norm(rp)
    z, x = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])
    if np.allclose(a, z, atol=1e-15) and np.allclose(rp, x, atol=1e-15):
        return np.eye(3)
    A1 = align_to(z, a)
    x1 = A1 @ x
    ang = math.atan2(float(a @ np.cross(x1, rp)), float(x1 @ rp))
    return rotation_matrix(a, ang) @ A1


def _canonical_generators(g: GroupSpec) -> list[np.ndarray]:
    z, x = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0])
    k, n = g.kind, g.n
    rot = [rotation_matrix(z, 2 * math.pi / n)] if n > 1 else []
    if k == "C1":
        return []
    if k == "Cs":
        return [reflection_matrix(z)]
    if k == "Ci":
        return [-np.eye(3)]
    if k == "Cn":
        return rot
    if k == "Cnv":
        return rot + [reflection_matrix(x)]
    if k == "Cnh":
        return rot + [reflection_matrix(z)]
    if k == "RDn":
        return rot + [rotation_matrix(x, math.pi)]
    if k == "S2n":
        return [rotation_matrix(z, math.pi / n) @ reflection_matrix(z)]
    if k == "Dnd":
        if n == 2:      # exactly d2d_symmetrize.d2d_elements
            S4 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]], float)
            sig = (np.diag([-1.0, 1.0, 1.0]) if g.frame == "xy"
                   else np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], float))
            return [S4, sig]
        S = rotation_matrix(z, math.pi / n) @ reflection_matrix(z)
        if g.frame == "xy":
            nrm = x
        else:
            nrm = np.array([math.cos(math.pi / (2 * n)), math.sin(math.pi / (2 * n)), 0.0])
        return [S, reflection_matrix(nrm)]
    raise SymmetryMapError(f"no generators for {g.label}")


def generators(spec, axis=(0, 0, 1), ref=(1, 0, 0), order: int | None = None) -> list[np.ndarray]:
    """Generator matrices of the group, conjugated into the frame (axis, ref)."""
    g = parse_spec(spec, order)
    if g.kind == "detect":
        raise SymmetryMapError("'detect' is not a group; call detect()")
    gens = [_snap(M) for M in _canonical_generators(g)]
    uses_ref = g.kind in ("Cnv", "RDn", "Dnd")
    A = frame_rotation(axis if axis is not None else (0, 0, 1), ref if ref is not None else (1, 0, 0),
                       strict=uses_ref and ref is not None)
    if not np.array_equal(A, np.eye(3)):
        gens = [_snap(A @ M @ A.T) for M in gens]
    return gens


def close_group(gens, max_order: int = 240, atol: float = 1e-9) -> list[np.ndarray]:
    """All products of the generators (breadth first), deduplicated at atol; mats[0] = I."""
    mats = [np.eye(3)]
    frontier = [0]
    while frontier:
        new = []
        for h in frontier:
            for G in gens:
                Mn = _snap(np.asarray(G, float) @ mats[h])
                if not any(np.allclose(Mq, Mn, atol=atol) for Mq in mats):
                    mats.append(Mn)
                    new.append(len(mats) - 1)
                    if len(mats) > max_order:
                        raise SymmetryMapError(f"the group did not close within order {max_order}")
        frontier = new
    return mats


# --------------------------------------------------------------------------- #
# vertex maps
# --------------------------------------------------------------------------- #

def _layout(comps):
    counts = np.array([len(c) for c in comps])
    off = np.r_[0, np.cumsum(counts)]
    owner = np.repeat(np.arange(len(comps)), counts)
    return counts, off, owner


def mean_edge(comps) -> float:
    return float(np.mean(np.concatenate(
        [np.linalg.norm(np.roll(c, -1, 0) - c, axis=1) for c in comps])))


# An image vertex within this many mean edges of the target component's polygon counts as ON
# that curve, for the unequal-counts diagnosis only: two samplings of one smooth curve differ by
# about edge^2 x curvature / 8 (about 0.02 mean edges on a tight link), far below it.
UNEQUAL_CURVE_TOL = 0.25


def _point_polyline_dist(Q, P):
    """Distance from each point of Q to the closed polygon P (via the edges at P's nearest vertex)."""
    n = len(P)
    _, v = cKDTree(P).query(Q, k=1)
    best = np.full(len(Q), np.inf)
    for a_, b_ in ((v - 1) % n, v), (v, (v + 1) % n):
        A, B = P[a_], P[b_]
        E = B - A
        f = np.clip(((Q - A) * E).sum(1) / np.maximum((E * E).sum(1), 1e-300), 0.0, 1.0)
        best = np.minimum(best, np.linalg.norm(Q - A - f[:, None] * E, axis=1))
    return best


def _unequal_pairs(j, M, comps, counts, off, owner):
    """Components the generator maps onto a component with a different vertex count.

    Each component's image is pre-assigned to a target by a nearest-vertex majority vote. Returns
    (pairs, on_curve): pairs lists (c, t) with counts[c] != counts[t] when the vote is a
    permutation (else []), and on_curve says whether every such image lies on its target's
    polygon within UNEQUAL_CURVE_TOL mean edges, i.e. the link is symmetric as a curve and only
    the sampling differs."""
    nc = len(comps)
    sig = [int(np.bincount(owner[j[off[c]:off[c + 1]]], minlength=nc).argmax()) for c in range(nc)]
    if sorted(sig) != list(range(nc)):
        return [], False
    pairs = [(c, sig[c]) for c in range(nc) if counts[c] != counts[sig[c]]]
    if not pairs:
        return [], False
    edge = mean_edge(comps)
    on_curve = all(_point_polyline_dist(comps[c] @ M.T, comps[t]).max() <= UNEQUAL_CURVE_TOL * edge
                   for c, t in pairs)
    return pairs, on_curve


def vertex_maps(comps, spec, *, axis=None, ref=None, tol: float = 1e-6, return_info: bool = False,
                order: int | None = None, max_order: int = 240):
    """(mats, perms[, info]) for the group acting on the vertices of `comps`, or SymmetryMapError.

    info: sigma[g] (component map), orient[(g, c)] (+1 shift / -1 reversal), shift[(g, c)]
    (image index of vertex 0), orbits (component orbits), edge (mean edge), dim (equivariant
    dimension = trace R), group (label), order.
    """
    g = parse_spec(spec, order)
    comps = [np.asarray(c, float) for c in comps]
    counts, off, owner = _layout(comps)
    nc = len(comps)
    X = np.vstack(comps)
    N = len(X)
    edge = mean_edge(comps)
    gens = generators(g, axis, ref)
    lim = tol * edge
    tree = cKDTree(X)
    cen = X.mean(0)
    cen_off = max((float(np.linalg.norm(G @ cen - cen)) for G in gens), default=0.0)

    def centre_hint(unequal_counts):
        # Unequal counts within an orbit move the vertex mean off the symmetry centre on their
        # own, so the hint would be false there.
        if cen_off <= lim or unequal_counts:
            return ""
        return (f"; the vertex mean {np.array2string(cen, precision=4)} is {cen_off / edge:.3g} "
                "mean edges from the symmetry centre -- centre the link first (X -= X.mean(0))")

    def perm_of(M, label):
        d, j = tree.query(X @ M.T, k=1)
        # Pre-assign each component's image to a target component (nearest-vertex majority
        # vote) before the vertex-exact test, so an orbit whose members have unequal vertex
        # counts is reported as that (the advice is a resample, not a projection or centring).
        pairs, on_curve = _unequal_pairs(j, M, comps, counts, off, owner)
        if pairs and on_curve:
            c, t = pairs[0]
            raise SymmetryMapError(
                f"{label}: component {c} ({counts[c]} vertices) maps onto component {t} "
                f"({counts[t]} vertices) as a curve; unequal counts in one orbit -- resample the "
                "orbit to one count (symmetrize_link_xyz --beads, or "
                "symmetry_maps.equivariant_resample)")
        worst = float(d.max())
        if worst > lim:
            raise SymmetryMapError(
                f"{label}: worst vertex {worst / edge:.3g} mean edges from its image (tol {tol}) "
                f"-- not vertex-exact; project first{centre_hint(bool(pairs))}")
        if len(np.unique(j)) != N:
            raise SymmetryMapError(f"{label}: the nearest-vertex map is not a bijection")
        sigma = owner[j[off[:-1]]]
        if not np.array_equal(owner[j], sigma[owner]):
            raise SymmetryMapError(f"{label}: the vertex map splits a component across two")
        if sorted(sigma.tolist()) != list(range(nc)):
            raise SymmetryMapError(f"{label}: the component map {sigma.tolist()} is not a permutation")
        # A bijection that maps every component into one component forces equal counts around
        # each cycle of sigma, so unequal counts cannot reach this point (caught above).
        for c in range(nc):
            t, n = int(sigma[c]), int(counts[c])
            k = j[off[c]:off[c + 1]] - off[t]
            if n < 3:
                continue
            dlt = int((k[1] - k[0]) % n)
            if dlt not in (1, n - 1):
                raise SymmetryMapError(f"{label}: component {c}'s index map is not a shift or a reversal")
            o = 1 if dlt == 1 else -1
            if not np.array_equal(k, (k[0] + o * np.arange(n)) % n):
                raise SymmetryMapError(f"{label}: component {c}'s index map is not k -> a + o*k "
                                       "(edges are not mapped to edges)")
        return j

    gp = []
    for i, M in enumerate(gens):
        lab = (f"{g.label} generator {i} (det {np.linalg.det(M):+.0f}, "
               f"trace {np.trace(M):+.3f})")
        gp.append((M, perm_of(M, lab)))

    mats, perms = [np.eye(3)], [np.arange(N)]
    frontier = [0]
    while frontier:
        new = []
        for h in frontier:
            for G, pG in gp:
                Mn = _snap(G @ mats[h])
                pn = pG[perms[h]]
                hit = next((q for q, Mq in enumerate(mats) if np.allclose(Mq, Mn, atol=1e-9)), None)
                if hit is None:
                    mats.append(Mn)
                    perms.append(pn)
                    new.append(len(mats) - 1)
                    if len(mats) > max_order:
                        raise SymmetryMapError(f"the group did not close within order {max_order}")
                elif not np.array_equal(perms[hit], pn):
                    raise SymmetryMapError("generator vertex maps violate a group relation "
                                           "(the maps are not a group action)")
        frontier = new

    info = {"sigma": {}, "orient": {}, "shift": {}, "edge": edge, "group": g.label,
            "order": len(mats)}
    for gi, (M, p) in enumerate(zip(mats, perms)):
        r = float(np.abs(X @ M.T - X[p]).max())
        if r > lim:
            raise SymmetryMapError(f"element {gi}: composed map residual {r / edge:.3g} mean edges "
                                   f"(tol {tol})")
        sig = owner[p[off[:-1]]]
        info["sigma"][gi] = sig
        for c in range(nc):
            n = int(counts[c])
            k = p[off[c]:off[c + 1]] - off[sig[c]]
            info["orient"][(gi, c)] = 1 if n < 2 or (k[1] - k[0]) % n == 1 else -1
            info["shift"][(gi, c)] = int(k[0])
    orbits, seen = [], set()
    for c in range(nc):
        if c in seen:
            continue
        orb = sorted({int(info["sigma"][gi][c]) for gi in range(len(mats))})
        seen |= set(orb)
        orbits.append(orb)
    info["orbits"] = orbits
    info["dim"] = float(sum(np.trace(M) * np.count_nonzero(p == np.arange(N))
                            for M, p in zip(mats, perms)) / len(mats))
    return (mats, perms, info) if return_info else (mats, perms)


# --------------------------------------------------------------------------- #
# Reynolds operator and the equivariant basis
# --------------------------------------------------------------------------- #

def reynolds(V, mats, perms) -> np.ndarray:
    """(R V)[i] = 1/|G| sum_g V[pi_g[i]] @ M_g  (V is N x 3)."""
    V = np.asarray(V, float)
    if len(mats) == 1:
        return V.copy()
    out = np.zeros_like(V)
    for M, p in zip(mats, perms):
        out += V[p] @ np.asarray(M, float)
    return out / len(mats)


def reynolds_matrix(N: int, mats, perms) -> sparse.csr_matrix:
    """R as a sparse 3N x 3N matrix acting on V.ravel()."""
    G = len(mats)
    rows, cols, vals = [], [], []
    ar = np.arange(N)
    for M, pi in zip(mats, perms):
        Mt = np.asarray(M, float).T / G
        for a in range(3):
            for b in range(3):
                if abs(Mt[a, b]) > 1e-15:
                    rows.append(3 * ar + a)
                    cols.append(3 * np.asarray(pi) + b)
                    vals.append(np.full(N, Mt[a, b]))
    return sparse.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(3 * N, 3 * N))


def equivariant_basis(R) -> sparse.csr_matrix:
    """Orthonormal basis (3N x m, sparse) of range(R), built block by block over vertex orbits."""
    R = sparse.csr_matrix(R)
    N3 = R.shape[0]
    if R.nnz == N3 and np.allclose(R.diagonal(), 1.0):
        return sparse.identity(N3, format="csr")
    nc, lab = sparse.csgraph.connected_components(R + R.T, directed=False)
    order = np.argsort(lab, kind="stable")
    starts = np.r_[0, np.cumsum(np.bincount(lab, minlength=nc))]
    ri, ci, vv = [], [], []
    m = 0
    for c in range(nc):
        idx = order[starts[c]:starts[c + 1]]
        blk = R[idx][:, idx].toarray()
        w, V = np.linalg.eigh(0.5 * (blk + blk.T))
        for k in np.where(w > 0.5)[0]:
            nz = np.abs(V[:, k]) > 0
            ri.append(idx[nz])
            ci.append(np.full(int(nz.sum()), m))
            vv.append(V[nz, k])
            m += 1
    if m == 0:
        return sparse.csr_matrix((N3, 0))
    return sparse.csr_matrix((np.concatenate(vv), (np.concatenate(ri), np.concatenate(ci))),
                             shape=(N3, m))


def symmetry_error(X, mats, perms) -> float:
    """max_g |X M_g^T - X[pi_g]|_inf (absolute units)."""
    X = np.asarray(X, float)
    return max(float(np.abs(X @ np.asarray(M, float).T - X[p]).max()) for M, p in zip(mats, perms))


# --------------------------------------------------------------------------- #
# equivariant resample
# --------------------------------------------------------------------------- #

def _cumarc(P):
    e = np.linalg.norm(np.roll(P, -1, 0) - P, axis=1)
    return np.r_[0.0, np.cumsum(e)[:-1]], float(e.sum())


def _sample(P, t, interp):
    s, L = _cumarc(P)
    if interp == "linear":
        Q = np.vstack([P, P[:1]])
        ss = np.r_[s, L]
        return np.column_stack([np.interp(np.mod(t, L), ss, Q[:, k]) for k in range(3)])
    cs = CubicSpline(np.r_[s, L], np.vstack([P, P[:1]]), bc_type="periodic")
    return cs(np.mod(t, L))


def equivariant_resample(comps, spec, total: int, interp: str = "spline", min_count: int = 8, *,
                         axis=None, ref=None, order: int | None = None, tol: float = 1e-6):
    """Resample a vertex-exactly symmetric link to about `total` vertices, keeping the symmetry exact.

    Per component orbit: the representative r (lowest index) is resampled by arclength with a
    count that is a multiple of the number k of rotations in its stabiliser, on a grid that
    starts at a reflection's fixed point when the stabiliser holds one; every other member is
    the image of r. One strict Reynolds pass over freshly built maps removes float residue.
    Returns (new_comps, report) with report rows (orbit, n_old, n_new, k, n_reflections).
    """
    comps = [np.asarray(c, float) for c in comps]
    mats, perms, info = vertex_maps(comps, spec, axis=axis, ref=ref, tol=tol, return_info=True,
                                    order=order)
    nc, G = len(comps), len(mats)
    sig = [info["sigma"][gi] for gi in range(G)]
    Ls = [_cumarc(c)[1] for c in comps]
    Ltot = sum(Ls)
    newc = [None] * nc
    report = []
    for orb in info["orbits"]:
        r = orb[0]
        n_old = len(comps[r])
        s_r, L = _cumarc(comps[r])
        stab = [gi for gi in range(G) if sig[gi][r] == r]
        shifts = {int(info["shift"][(gi, r)]) for gi in stab if info["orient"][(gi, r)] == 1}
        refl = [gi for gi in stab if info["orient"][(gi, r)] == -1]
        k = len(shifts)
        n = max(float(min_count), total * L / Ltot)
        n = int(k * max(1, round(n / k)))
        b = 0.5 * s_r[int(info["shift"][(refl[0], r)])] if refl else 0.0
        Yr = _sample(comps[r], b + np.arange(n) * L / n, interp)
        for c in orb:
            gi = next(q for q in range(G) if sig[q][r] == c)
            o = info["orient"][(gi, r)]
            newc[c] = (Yr @ np.asarray(mats[gi], float).T)[(o * np.arange(n)) % n]
        report.append((list(orb), n_old, n, k, len(refl)))
    mats2, perms2 = vertex_maps(newc, spec, axis=axis, ref=ref, tol=tol, order=order)
    Xs = reynolds(np.vstack(newc), mats2, perms2)
    off = np.r_[0, np.cumsum([len(c) for c in newc])]
    return [Xs[off[i]:off[i + 1]] for i in range(nc)], report


# --------------------------------------------------------------------------- #
# detection
# --------------------------------------------------------------------------- #

DETECT_CANDIDATES = (["C1", "Cs", "Ci"] + [f"C{n}" for n in range(2, 13)]
                     + [f"C{n}v" for n in range(2, 9)] + [f"C{n}h" for n in range(2, 7)]
                     + [f"RD{n}" for n in range(2, 7)]
                     + ["D2d", "D2d:diag", "D3d", "D4d", "S4", "S6", "S8"])


def detect(comps, candidates=None, tol: float = 1e-6, axis=None, ref=None):
    """Every candidate group that acts vertex-exactly, as (label, order, dim), sorted by order.

    Never chooses one: the caller decides. Canonical frame unless axis/ref are given.
    """
    out = []
    for i, name in enumerate(candidates or DETECT_CANDIDATES):
        try:
            mats, perms, info = vertex_maps(comps, name, axis=axis, ref=ref, tol=tol, return_info=True)
        except SymmetryMapError:
            continue
        out.append((len(mats), i, parse_spec(name).label, info["dim"]))
    out.sort()
    return [(lab, n, dim) for n, _, lab, dim in out]


# --------------------------------------------------------------------------- #
# diagnostic main
# --------------------------------------------------------------------------- #

def _vec(text):
    v = [float(t) for t in text.split(",")]
    if len(v) != 3:
        raise argparse.ArgumentTypeError("expected X,Y,Z")
    return v


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("group", help="group name (see above), or 'detect'")
    ap.add_argument("--order", type=int, default=None, help="n for Cn / Cnv / Cpv / RDp / Z/nZ")
    ap.add_argument("--axis", type=_vec, default=None, help="principal axis X,Y,Z (default 0,0,1)")
    ap.add_argument("--ref", type=_vec, default=None, help="reference direction X,Y,Z (default 1,0,0)")
    ap.add_argument("--tol", type=float, default=1e-6, help="vertex tolerance in mean edges (default 1e-6)")
    ap.add_argument("--no-centre", action="store_true", help="do not subtract the vertex mean first")
    a = ap.parse_args(argv)
    comps = [np.asarray(c, float) for c in read_xyz(a.input)]
    X = np.vstack(comps)
    cen = X.mean(0)
    if not a.no_centre:
        comps = [c - cen for c in comps]
    print(f"{a.input.name}: {len(comps)} components, counts {[len(c) for c in comps]}, "
          f"vertex mean {np.array2string(cen, precision=6)}{' (subtracted)' if not a.no_centre else ''}")
    try:
        spec = parse_spec(a.group, a.order)
    except SymmetryMapError as e:
        print(f"  error: {e}", file=sys.stderr)
        return 2
    if spec.kind == "detect":
        hits = detect(comps, tol=a.tol, axis=a.axis, ref=a.ref)
        print(f"  exact groups (tol {a.tol} mean edges), by order:")
        for lab, n, dim in hits:
            print(f"    {lab:9s} |G| {n:3d}   equivariant dim {dim:.0f}")
        return 0
    try:
        mats, perms, info = vertex_maps(comps, spec, axis=a.axis, ref=a.ref, tol=a.tol, return_info=True)
    except SymmetryMapError as e:
        print(f"  REFUSED: {e}")
        return 3
    N = len(np.vstack(comps))
    print(f"  {info['group']}: |G| = {len(mats)}, mean edge {info['edge']:.4g}, orbits {info['orbits']}")
    for gi, M in enumerate(mats):
        print(f"   g{gi:<3d} det {np.linalg.det(M):+.0f} tr {np.trace(M):+.3f}  "
              f"sigma {info['sigma'][gi].tolist()}  "
              f"orient {[info['orient'][(gi, c)] for c in range(len(comps))]}")
    V = np.random.default_rng(0).standard_normal((N, 3))
    W = np.random.default_rng(1).standard_normal((N, 3))
    RV = reynolds(V, mats, perms)
    print(f"  |R R V - R V| = {np.abs(reynolds(RV, mats, perms) - RV).max():.2e}   "
          f"<RV,W> - <V,RW> = {abs((RV * W).sum() - (V * reynolds(W, mats, perms)).sum()):.2e}")
    print(f"  symmetry error {symmetry_error(np.vstack(comps), mats, perms):.2e}   "
          f"equivariant dim {info['dim']:.0f} of {3 * N}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
