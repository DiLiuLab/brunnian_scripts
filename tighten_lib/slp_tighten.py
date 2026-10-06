#!/usr/bin/env python3
"""Symmetry-equivariant sequential-LP (SLP) polish for a tightened link.

    python3 tighten_lib/slp_tighten.py BEST.xyz -o slp/slp.xyz --group C7
    python3 tighten_lib/slp_tighten.py BEST.xyz -o slp/slp.xyz --group C7 --gate-ref REF.xyz
    python3 tighten_lib/slp_tighten.py BEST.xyz --group C7 --diagnose   # first-order feasible fraction
    python3 tighten_lib/slp_tighten.py BEST.xyz --group detect          # which candidate groups hold exactly
    python3 tighten_lib/slp_tighten.py BEST.xyz --group D2d --dry-run   # group, dimension, plan

  The 8BL recipe (298.242 -> 291.745 in 618 s, residual 0.307 -> 0.005, 2486 struts, D2d exact,
  HOMFLY SAME; an earlier build's 291.745 read residual 0.063. The original hand-run slp_B chain
  reached 291.748 at residual 0.010):

    python3 tighten_lib/slp_tighten.py 8BL_D2d_final_298.242.xyz -o slp/8BL.xyz \\
        --group D2d --passes 2 --resample-vu 6 --gate-ref 8BL_D2d_final_298.242.xyz

WHAT IT DOES. A sequential linear program over per-vertex displacements,
restricted to the exact symmetry subspace: the range of the Reynolds operator
over the vertex maps that tighten_lib/symmetry_maps.py builds (C1, Cs, Ci, Cn,
Cnv, Cnh, RDn, Dnd, S2n; C1 means no symmetry and all of R^3N). Each LP
minimises the linearised length subject to
  - every pair of edges nearer than 2tau + 3delta + 0.04tau not coming closer
    than 2tau (linearised, closest points held fixed);
  - both branches of Rawdon's polygonal minRad staying >= tau at every vertex
    within 0.12tau of it;
  - every edge staying within [0.8, 1.25] x its component's mean edge (an edge
    already outside may not move further out);
  - a trust region |dx_k| <= delta.
The tau of the first two rows starts as the pass's starting tau, so a tau that
an accepted step dipped (by up to --tau-slack) is pulled back up. If that LP is
infeasible at an iterate whose own octrope tau is lower, the rows are rebuilt
at min(target, the iterate's tau), where dx = 0 is feasible, and that becomes
the target for the rest of the pass (floor_rebuilds counts these). A pass then
reaches --dmin because no step gains, not because the rows asked for more than
the trust region could restore. The acceptance floor stays fixed at the pass's
starting tau x (1 - --tau-slack).
The step is projected back onto the exact symmetry, refined by a second-order
correction LP, and ACCEPTED only if octrope's ropelength falls, octrope's
thickness stays within --tau-slack of the pass's starting thickness, AND a PL
isotopy certificate proves that the straight-line homotopy from the old polygon
to the new one passes through embedded polygons only. delta grows x1.5 on an
accept and halves on a reject. All scale constants are in units of tau.

WHY. Under --Symmetry the patched ridgerunner is forced onto the Animation
stepper, which "won't converge to low residual". A symmetric run CAN park well
above residual 0.1 (8BL RD2 legs 0.31-0.43, 7BL C7 0.322, lark C5 0.263), and
there the floor was the stepper's, not the configuration's; it is not universal
(4BL_wider under C2v stopped on residual 0.0146). On 8BL D2d, after every
geometric move had returned nothing, the hand-run SLP chain took
298.242 -> 291.748 with the residual 0.307 -> 0.010 (this tool's recipe: 291.745). One
pass of this tool, the default, took 7BL C7 274.365 -> 272.816, lark C5
210.041 -> 202.282 and 4_LC Cs 113.590 -> 113.049, all link-gated. Below a
residual of about 0.1 it buys little; use the default-stepper polish there.

EDGE BOUNDS ARE NOT OPTIONAL. Without them an end-ring edge collapsed (edge
ratio 125, ropelength 297.035) and HOMFLY still said SAME. Turning them off
requires --allow-edge-collapse.

IT DOES NOT CHECK TOPOLOGY. The certificate proves each accepted step is an
isotopy of the polygon, but gate the output anyway: pass --gate-ref REF (exit 4
unless HOMFLY SAME in the cycle's generic frame), or let tighten_cycle.py --slp
gate it. The gate is link-gated (HOMFLY); the diagram may change -- check with
a diagram-isomorphism tool before quoting a per-diagram minimum.

INPUT. A tightened file (longest edge <= 0.5 D and >= 50 struts; otherwise
descend it first) that is vertex-exactly symmetric under --group. A file that is
only nearly symmetric (within --project-tol mean edges) is projected once,
certified, and refused (exit 3) if the projection costs more than
--max-projection-cost of ropelength. The link is centred for the run and the
centre is restored (its symmetric part, when the input was symmetric about the
origin) on output.

PASSES. --passes P splits --max-minutes over P passes (unused time rolls over).
Between passes --resample-vu V resamples equivariantly to V vertices per unit
ropelength (certified as an isotopy through a common refinement); each pass
re-reads tau, the edge means and the maps from its own input. The result is
the lowest-ropelength pass best. --robustness-vu 4.5,6,8 re-reads the output
resampled to each v/u (reported, not required).

OUTPUT. OUT.xyz (15 decimals, rescaled about the centre to tau = --final-tau,
default 0.49997, the normalised window of tighten_cycle.py), OUT_trace.csv (one
row per accepted step), OUT_pass{k}.xyz (each pass's best, unscaled),
OUT_pass{k}_input.xyz (the polygon pass k started from, after any projection
or resample: the floor of min_edge_frac, re-read by tighten_cycle.py) and
OUT_cur.xyz (every --save-every accepted steps). -o may not be the --gate-ref
file (exit 2); the reference is read before anything is written. The last line is

    SLP RESULT status=ok group=D2d order=8 rop_in=... rop_out=... gain=... tau_out=...
    minrad_out=... N_in=... N_out=... passes=... best_pass=... accepted=... trials=...
    cert_rejects=... tau_rejects=... lp_fail=... floor_rebuilds=... soc=... on_bounds=...
    edge_ratio=... min_edge_frac=... min_edge_frac_in=... sym_err=... projection_cost=...
    cert_margin=... max_disp_over_tau=... tau_min_ratio=... gate=... [robust=...]
    wall_s=... note="link-gated (HOMFLY); ..."

(one line; robust= only with --robustness-vu, gate=not-run without --gate-ref).

min_edge_frac is the shortest edge of each component over min(elo x mean edge,
shortest edge) of the pass input the LP ran on (OUT_pass{best_pass}_input.xyz):
>= 1 with the bounds on, and the collapse guard to gate on. min_edge_frac_in is
the same against the tool's input (rescaled to tau_out); it falls below 1 after
a resample to a higher v/u by design.

EXIT CODES. 0 ok (including no gain); 2 usage, or an input that is not
tightened; 3 symmetry refused; 4 gate not SAME; 5 octrope parse failure.

Temporary files (the scratch .vect, the generic-frame gate copies, anything the
octrope binaries drop in their cwd) go in --work-dir or a temporary directory,
never in the current directory; `residual` is not run here. Needs
numpy, scipy >= 1.6 (HiGHS) and octrope's `ropelength` and `struts` on PATH.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import linprog, nnls
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tighten_link_xyz import TightenError, read_xyz, write_xyz, write_vect  # noqa: E402
import symmetry_maps as sm  # noqa: E402

NOTE = ('note="link-gated (HOMFLY); the diagram may change -- check with a diagram-isomorphism '
        'tool before quoting a per-diagram minimum"')

EXIT_OK, EXIT_USAGE, EXIT_SYMMETRY, EXIT_GATE, EXIT_OCTROPE = 0, 2, 3, 4, 5


class Refused(Exception):
    def __init__(self, code: int, reason: str, message: str):
        super().__init__(message)
        self.code, self.reason = code, reason


# --------------------------------------------------------------------------- #
# polygon topology and geometry
# --------------------------------------------------------------------------- #

class Topo:
    """Closed components laid end to end: off, comp (owner), nxt, prv."""

    def __init__(self, counts):
        self.counts = [int(c) for c in counts]
        n = self.counts
        self.off = np.r_[0, np.cumsum(n)].astype(int)
        self.N = int(self.off[-1])
        self.comp = np.repeat(np.arange(len(n)), n)
        self.nxt = np.concatenate([o + (np.arange(m) + 1) % m for o, m in zip(self.off[:-1], n)])
        self.prv = np.concatenate([o + (np.arange(m) - 1) % m for o, m in zip(self.off[:-1], n)])

    def split(self, X):
        return [np.asarray(X[self.off[k]:self.off[k + 1]]) for k in range(len(self.counts))]


def edge_lengths(X, topo):
    return np.linalg.norm(X[topo.nxt] - X, axis=1)


def comp_mean_per_vertex(X, topo):
    el = edge_lengths(X, topo)
    means = np.array([el[topo.comp == k].mean() for k in range(len(topo.counts))])
    return means[topo.comp], means, np.array([el[topo.comp == k].min() for k in range(len(topo.counts))])


def seg_seg(a0, a1, b0, b1):
    """Closest points of segment pairs: s, t, distance, unit normal from the b point to the a point."""
    u = a1 - a0
    v = b1 - b0
    w = a0 - b0
    A = (u * u).sum(1)
    B = (u * v).sum(1)
    Cc = (v * v).sum(1)
    Dd = (u * w).sum(1)
    E = (v * w).sum(1)
    den = A * Cc - B * B
    ok = den > 1e-14 * np.maximum(A * Cc, 1e-300)
    s = np.where(ok, (B * E - Cc * Dd) / np.where(ok, den, 1.0), 0.0)
    s = np.clip(s, 0.0, 1.0)
    t = (B * s + E) / np.maximum(Cc, 1e-300)
    t0 = t < 0
    t1 = t > 1
    t = np.clip(t, 0.0, 1.0)
    s = np.where(t0, np.clip(-Dd / np.maximum(A, 1e-300), 0, 1), s)
    s = np.where(t1, np.clip((B - Dd) / np.maximum(A, 1e-300), 0, 1), s)
    dv = (a0 + s[:, None] * u) - (b0 + t[:, None] * v)
    d = np.linalg.norm(dv, axis=1)
    return s, t, d, dv / np.maximum(d, 1e-15)[:, None]


def arc_positions(X, topo):
    """Arclength of each edge midpoint along its component, and each component's length."""
    el = edge_lengths(X, topo)
    arc = np.zeros(topo.N)
    Lc = np.zeros(len(topo.counts))
    for k in range(len(topo.counts)):
        e = el[topo.off[k]:topo.off[k + 1]]
        arc[topo.off[k]:topo.off[k + 1]] = np.r_[0.0, np.cumsum(e)[:-1]] + 0.5 * e
        Lc[k] = e.sum()
    return arc, Lc, el


def contact_pairs(X, topo, cutoff, arc_excl):
    """Edge pairs (edge i = X[i] -> X[nxt i]) nearer than cutoff; same-component pairs closer
    than arc_excl along the curve are excluded. Returns i, j, s, t, d, unit normal."""
    arc, Lc, el = arc_positions(X, topo)
    nxt, comp = topo.nxt, topo.comp
    tree = cKDTree(0.5 * (X + X[nxt]))
    pr = tree.query_pairs(cutoff + el.max() + 1e-12, output_type="ndarray")
    if len(pr) == 0:
        z = np.zeros(0)
        return z.astype(int), z.astype(int), z, z, z, np.zeros((0, 3))
    i, j = pr[:, 0], pr[:, 1]
    same = comp[i] == comp[j]
    da = np.abs(arc[i] - arc[j])
    da = np.minimum(da, Lc[comp[i]] - da)
    keep = ~(same & (da < arc_excl))
    i, j = i[keep], j[keep]
    s, t, d, nv = seg_seg(X[i], X[nxt[i]], X[j], X[nxt[j]])
    k = d < cutoff
    return i[k], j[k], s[k], t[k], d[k], nv[k]


def contact_rows(N, nxt, i, j, s, t, nv):
    """Sparse rows a_p with d(distance_p) = a_p . dx (Danskin: s, t held fixed)."""
    P = len(i)
    rows = np.repeat(np.arange(P), 12)
    verts = np.stack([i, nxt[i], j, nxt[j]], 1)
    coef = np.stack([1 - s, s, -(1 - t), -t], 1)
    cols = (verts[:, :, None] * 3 + np.arange(3)[None, None, :]).reshape(P, 12)
    vals = (coef[:, :, None] * nv[:, None, :]).reshape(P, 12)
    return sparse.csr_matrix((vals.ravel(), (rows, cols.ravel())), shape=(P, 3 * N))


def _mr(T):
    a = T[:, 1] - T[:, 0]
    b = T[:, 2] - T[:, 1]
    la = np.linalg.norm(a, axis=1)
    lb = np.linalg.norm(b, axis=1)
    ua = a / la[:, None]
    ub = b / lb[:, None]
    th = np.linalg.norm(np.cross(ua, ub), axis=1) / (1 + (ua * ub).sum(1))
    with np.errstate(divide="ignore", invalid="ignore"):
        return la / (2 * th), lb / (2 * th)


def minrad_parts(X, topo):
    """Rawdon's polygonal minRad at each vertex, both branches: |e_prev|, |e_next| / (2 tan(theta/2))."""
    return _mr(np.stack([X[topo.prv], X, X[topo.nxt]], 1))


def minrad_rows(X, topo, v, h):
    """Central-difference gradients of both minRad branches at vertices v (9 coordinates each)."""
    N = len(X)
    p, q = topo.prv[v], topo.nxt[v]
    trip = np.stack([X[p], X[v], X[q]], 1)
    G0 = np.zeros((len(v), 9))
    G1 = np.zeros((len(v), 9))
    for k in range(9):
        T1 = trip.copy()
        T2 = trip.copy()
        T1[:, k // 3, k % 3] += h
        T2[:, k // 3, k % 3] -= h
        r1a, r1b = _mr(T1)
        r2a, r2b = _mr(T2)
        G0[:, k] = (r1a - r2a) / (2 * h)
        G1[:, k] = (r1b - r2b) / (2 * h)
    vv = np.stack([p, v, q], 1)
    cols = (vv[:, :, None] * 3 + np.arange(3)[None, None, :]).reshape(len(v), 9)
    rows = np.repeat(np.arange(len(v)), 9)
    return (sparse.csr_matrix((G0.ravel(), (rows, cols.ravel())), shape=(len(v), 3 * N)),
            sparse.csr_matrix((G1.ravel(), (rows, cols.ravel())), shape=(len(v), 3 * N)))


def edge_rows(X, topo):
    """Rows u_e . (dx_nxt - dx_e): the linearised change of each edge length."""
    N = topo.N
    e = X[topo.nxt] - X
    l = np.linalg.norm(e, axis=1)
    u = e / l[:, None]
    rr = np.repeat(np.arange(N), 6)
    cc = np.stack([3 * topo.nxt[:, None] + np.arange(3), 3 * np.arange(N)[:, None] + np.arange(3)],
                  1).reshape(N, 6)
    return sparse.csr_matrix((np.concatenate([u, -u], 1).ravel(), (rr, cc.ravel())), shape=(N, 3 * N)), l


def length_grad(X, topo):
    a = X - X[topo.prv]
    b = X[topo.nxt] - X
    return a / np.linalg.norm(a, axis=1)[:, None] - b / np.linalg.norm(b, axis=1)[:, None]


# --------------------------------------------------------------------------- #
# the PL isotopy certificate
# --------------------------------------------------------------------------- #

def _chain_monotone(Xa, Xb, topo, i, j, max_edges):
    """True if the shorter run of edges from edge i to edge j (same component) has every edge
    direction, at both ends of the homotopy, in one open hemisphere: then every intermediate
    chain is strictly monotone along that direction, hence an embedded arc, so i and j never meet."""
    c = topo.comp[i]
    n = topo.counts[c]
    o = topo.off[c]
    a, b = i - o, j - o
    fw = (b - a) % n
    if fw > n - fw:
        a, b, fw = b, a, n - fw
    if fw + 1 > max_edges:
        return False
    idx = o + (a + np.arange(fw + 1)) % n
    nx = topo.nxt[idx]
    Ea = Xa[nx] - Xa[idx]
    Eb = Xb[nx] - Xb[idx]
    Ua = Ea / np.linalg.norm(Ea, axis=1)[:, None]
    Ub = Eb / np.linalg.norm(Eb, axis=1)[:, None]
    v = Ua.sum(0) + Ub.sum(0)
    nv = np.linalg.norm(v)
    if nv < 1e-12:
        return False
    v = v / nv
    return bool(min((Ua @ v).min(), (Ub @ v).min()) > 1e-9)


def isotopy_certificate(X, Y, topo, chain_max: int = 64):
    """Sufficient condition that the straight-line homotopy H_t = (1-t)X + tY, t in [0, 1], is
    through embedded polygons. Returns (ok, worst_margin, why, n_chain).

    With E_e = X[nxt e] - X[e], dE_e its change and m_e = max(|dX_e|, |dX_nxt e|):
      (a) |dE_e| < |E_e| for every edge, so no edge degenerates;
      (b) segdist_X(e, f) > m_e + m_f for every non-adjacent pair whose midpoints are within
          2 max|dX| + max|E| (a sound superset: farther pairs cannot close the gap), since each
          point of edge e moves at most m_e; a same-component pair that fails (b) may instead be
          cleared by a chain test: every edge between them, at both ends, lies in one open
          hemisphere, so that run stays a monotone, embedded arc;
      (c) theta_v + phi_prev(v) + phi_v < pi with phi_e = arcsin(|dE_e|/|E_e|), so adjacent edges
          share only their common vertex.
    """
    X = np.asarray(X, float)
    Y = np.asarray(Y, float)
    nxt, prv = topo.nxt, topo.prv
    dX = np.linalg.norm(Y - X, axis=1)
    E = X[nxt] - X
    lE = np.linalg.norm(E, axis=1)
    dE = np.linalg.norm((Y[nxt] - Y) - E, axis=1)
    if np.any(dE >= lE):
        return False, float((lE - dE).min()), "an edge could degenerate", 0
    m = np.maximum(dX, dX[nxt])
    worst = math.inf
    n_chain = 0
    tree = cKDTree(0.5 * (X + X[nxt]))
    pr = tree.query_pairs(2 * float(m.max()) + float(lE.max()) + 1e-12, output_type="ndarray")
    if len(pr):
        i, j = pr[:, 0], pr[:, 1]
        adj = (nxt[i] == j) | (nxt[j] == i)
        i, j = i[~adj], j[~adj]
        if len(i):
            _, _, d, _ = seg_seg(X[i], X[nxt[i]], X[j], X[nxt[j]])
            marg = d - m[i] - m[j]
            bad = np.where(marg <= 0)[0]
            if len(bad) > 20000:
                return False, float(marg.min()), "too many near pairs for the chain test", 0
            for k in bad:
                if topo.comp[i[k]] != topo.comp[j[k]] or not _chain_monotone(X, Y, topo, int(i[k]),
                                                                            int(j[k]), chain_max):
                    return False, float(marg[k]), "non-adjacent edges could meet", n_chain
                n_chain += 1
            good = np.ones(len(i), bool)
            good[bad] = False
            if good.any():
                worst = float(marg[good].min())
    ua = E[prv] / lE[prv][:, None]
    ub = E / lE[:, None]
    th0 = np.arccos(np.clip((ua * ub).sum(1), -1.0, 1.0))
    phi = np.arcsin(np.clip(dE / lE, 0.0, 1.0))
    if np.any(th0 + phi[prv] + phi >= np.pi):
        return False, worst, "a turning angle could reach pi", n_chain
    return True, worst, "", n_chain


def _closest_params(P, Q):
    """Arclength parameter (along closed polygon P) of the closest point of P to each point of Q."""
    n = len(P)
    nx = np.r_[np.arange(1, n), 0]
    E = P[nx] - P
    lE = np.linalg.norm(E, axis=1)
    s0 = np.r_[0.0, np.cumsum(lE)[:-1]]
    k = min(16, n)
    _, cand = cKDTree(0.5 * (P + P[nx])).query(Q, k=k)
    cand = np.atleast_2d(cand).reshape(len(Q), k)
    best_d = np.full(len(Q), np.inf)
    best_u = np.zeros(len(Q))
    for col in range(k):
        e = cand[:, col]
        w = Q - P[e]
        f = np.clip((w * E[e]).sum(1) / (lE[e] ** 2), 0.0, 1.0)
        d = np.linalg.norm(w - f[:, None] * E[e], axis=1)
        better = d < best_d
        best_d[better] = d[better]
        best_u[better] = s0[e[better]] + f[better] * lE[e[better]]
    return best_u, float(lE.sum()), s0


def _pl_eval(knots, vals, t, L):
    """Closed piecewise-linear curve through vals at increasing parameters knots (period L)."""
    kk = np.r_[knots[-1] - L, knots, knots[0] + L]
    vv = np.vstack([vals[-1:], vals, vals[:1]])
    tt = (t - knots[0]) % L + knots[0]
    return np.column_stack([np.interp(tt, kk, vv[:, c]) for c in range(3)])


def certify_resample(old_comps, new_comps, merge: float = 1e-9):
    """Certify that a resample is a PL isotopy, through a common refinement.

    Each new vertex gets the arclength parameter of its closest point on the old component;
    these must be strictly cyclically increasing (one turn, same orientation). Both polygons are
    then evaluated at the union of the old vertex parameters and these (a new parameter within
    merge*L of an old vertex is moved onto it), which leaves both curves unchanged and gives them
    one combinatorics, and isotopy_certificate is applied to the pair.
    Returns (ok, worst_margin, why).
    """
    A_all, B_all, counts = [], [], []
    for c, (P, Q) in enumerate(zip(old_comps, new_comps)):
        P = np.asarray(P, float)
        Q = np.asarray(Q, float)
        u, L, s0 = _closest_params(P, Q)
        u = u % L
        # a new parameter within merge*L of an old vertex (cyclically: L is 0) is moved onto it
        grid = np.r_[s0, L]
        pos = np.clip(np.searchsorted(grid, u), 1, len(grid) - 1)
        lo, hi = grid[pos - 1], grid[pos]
        near = np.where(u - lo <= hi - u, lo, hi)
        u = np.where(np.abs(u - near) < merge * L, near % L, u)
        r = int(np.argmin(u))
        u = np.roll(u, -r)
        Qr = np.roll(Q, -r, axis=0)
        gaps = np.diff(np.r_[u, u[0] + L])
        if np.any(gaps <= 0):
            return False, -math.inf, (f"component {c}: the new vertices are not in strictly increasing "
                                      "cyclic order along the old curve")
        T = np.unique(np.r_[s0, u])
        A_all.append(_pl_eval(s0, P, T, L))
        B_all.append(_pl_eval(u, Qr, T, L))
        counts.append(len(T))
    topo = Topo(counts)
    ok, worst, why, _ = isotopy_certificate(np.vstack(A_all), np.vstack(B_all), topo)
    return ok, worst, why


# --------------------------------------------------------------------------- #
# octrope
# --------------------------------------------------------------------------- #

def _grab(text: str, key: str):
    m = re.search(rf"(?<![A-Za-z] ){key}:\s*([0-9.eE+-]+)", text)
    return float(m.group(1)) if m else None


class Octrope:
    """octrope's binaries on a scratch .vect, always with cwd in the scratch directory."""

    def __init__(self, tmp: Path):
        self.tmp = Path(tmp).resolve()
        self.calls = 0
        self.failures = 0

    def _vect(self, X, topo, name="slp.vect"):
        p = self.tmp / name
        write_vect(p, [c.tolist() for c in topo.split(X)])
        return p

    def measure(self, X, topo):
        self.calls += 1
        p = self._vect(X, topo)
        try:
            out = subprocess.run(["ropelength", str(p)], capture_output=True, text=True,
                                 cwd=str(self.tmp)).stdout
        except OSError:
            self.failures += 1
            return None
        vals = {k: _grab(out, key) for k, key in (("rop", "Ropelength"), ("tau", "Thickness"),
                                                 ("minRad", "minRad"), ("minStrut", "minStrut"))}
        if None in vals.values() or not all(math.isfinite(v) for v in vals.values()):
            self.failures += 1
            return None
        return vals

    def struts(self, X, topo, radius):
        p = self._vect(X, topo)
        try:
            out = subprocess.run(["struts", "-s", "-n", "-r", f"{radius:.9f}", str(p)], capture_output=True,
                                 text=True, cwd=str(self.tmp))
        except OSError:
            return None
        txt = out.stdout + out.stderr
        m = re.search(r"Sorting (\d+) struts", txt) or re.search(r"(\d+) struts", txt)
        return int(m.group(1)) if m else None


def strut_radius(tau, minstrut):
    """tighten_cycle.strut_radius: nominal 0.5 (+1e-4) for a normalised file, else minStrut/2."""
    try:
        from tighten_cycle import strut_radius as sr
        return sr(tau, minstrut)
    except Exception:  # pragma: no cover - same rule, kept local if the driver cannot import
        if abs(tau - 0.5) <= 1e-3:
            return max(0.5, minstrut / 2) * (1 + 1e-4)
        return (minstrut / 2) * (1 + 1e-4)


# --------------------------------------------------------------------------- #
# the symmetric problem
# --------------------------------------------------------------------------- #

# octrope prints tau with 6 decimals (+-5e-7), so 2 tau_eff is uncertain by 1e-6: rows within
# OCT_RESOLUTION x tau_t (2e-6 at tau 0.5) below their floor are at it, not violating it.
OCT_RESOLUTION = 4e-6


def _at_floor(rhs, tol):
    """Right-hand sides a hair below 0 (within octrope's print resolution) become 0; tol 0 is a no-op."""
    rhs = np.asarray(rhs, float).copy()
    if tol > 0:
        rhs[(rhs < 0) & (rhs > -tol)] = 0.0
    return rhs

class Problem:
    """Everything that is fixed during one pass: maps, basis, targets, edge means."""

    def __init__(self, X, topo, spec, a, oc, moving=None):
        self.topo = topo
        self.a = a
        self.spec = spec
        comps = topo.split(X)
        self.mats, self.perms, self.info = sm.vertex_maps(comps, spec, axis=a.axis, ref=a.ref,
                                                          tol=a.tol, return_info=True)
        N = topo.N
        self.trivial = len(self.mats) == 1
        if self.trivial:
            B = sparse.identity(3 * N, format="csr")
        else:
            B = sm.equivariant_basis(sm.reynolds_matrix(N, self.mats, self.perms))
        self.box = self.trivial          # B is a column selection of I: trust region as bounds
        if moving is not None:
            keep = np.repeat(np.isin(topo.comp, moving), 3)
            Bc = B.tocsc()
            cols = [k for k in range(Bc.shape[1])
                    if keep[Bc.indices[Bc.indptr[k]:Bc.indptr[k + 1]]].all()]
            B = Bc[:, cols].tocsr()
        self.B = B
        self.Bt = B.T.tocsr()
        self.m = B.shape[1]
        o = oc.measure(X, topo)
        if o is None:
            raise Refused(EXIT_OCTROPE, "octrope", "octrope could not measure the pass input")
        self.o0 = o
        self.tau = o["tau"]
        self.D = 2 * self.tau
        self.emean, self.comp_mean, self.comp_min = comp_mean_per_vertex(X, topo)

    def project(self, Y):
        return Y if self.trivial else sm.reynolds(Y, self.mats, self.perms)

    def tau_eff(self, tau_x=None):
        """The floor the LP rows are built on: the pass's starting tau_t, or the iterate's own
        octrope tau when that is lower (acceptance allows it to dip to tau_t(1 - tau_slack))."""
        return self.tau if tau_x is None else min(self.tau, float(tau_x))

    def system(self, X, delta, tau_x=None):
        """Linearised constraints at X for trust radius delta (rows reused when delta shrinks).

        The contact rows ask for 2 tau_eff and the curvature rows for tau_eff, tau_eff =
        min(tau_t, tau_x); tau_x is octrope's tau of X, None meaning tau_t itself.
          - run_pass builds at its target, tau_t at first (tau_x None). After an accepted step
            dipped tau below it those rows ask the LP to restore it, which keeps the iterate
            off the acceptance floor tau_t (1 - tau_slack).
          - If that LP is infeasible, run_pass rebuilds with tau_x = tau(X), where dx = 0 is
            feasible, and lowers its target to that tau_eff for the rest of the pass. Before
            this, every failure halved delta and the pass stopped at dmin with gain left (7BL
            C7 273.945, lark C5 207.19).
          - Building at the iterate's tau ALWAYS was measured and lost: with nothing pulling tau
            back up it sank to the floor within a few dozen steps, tau rejects halved delta, and
            the passes stopped at dmin in 10-40 s (8BL recipe 298.133, 4_LC Cs 113.540, lark
            209.952, against 291.745, 113.045, 207.188 with rows at tau_t).
          - Rebuilding per iterate but returning to tau_t afterwards was also measured: the
            restore demand grew with the drift, pinned delta near it, and 7BL stopped at dmin at
            273.676 (189 s), where a restart still gained 0.80 in 2 minutes; with the target
            lowered it went on to 272.816, and lark C5 to 202.282 (against 202.383). A target
            of min(tau_t, tau(X) + 0.25 delta) lost to it on lark (203.247) and 7BL (272.865).
        octrope prints tau to 6 decimals, so in the rebuild a row that sits below its floor by less
        than that resolution is read as AT its floor (may not get closer) rather than as violated."""
        a, t, B = self.a, self.topo, self.B
        tau_t = self.tau
        tau = self.tau_eff(tau_x)
        D = 2 * tau
        # only the rebuild promises dx = 0 feasible; the tau_t rows keep their exact right-hand
        # sides (a hair-negative row there is a real, if tiny, ask to restore tau)
        at_floor = OCT_RESOLUTION * tau_t if tau_x is not None else 0.0
        g = length_grad(X, t).ravel()
        buf = a.buf * tau_t
        i, j, s_, t_, d, nv = contact_pairs(X, t, D + 3 * delta + a.contact_slack * tau_t, math.pi * tau_t)
        A = contact_rows(t.N, t.nxt, i, j, s_, t_, nv)
        rows = [-(A @ B)]
        rhs = [_at_floor(d - D - buf, at_floor)]
        ra, rb = minrad_parts(X, t)
        cv = np.where(np.minimum(ra, rb) < tau + a.curv_window * tau_t)[0]
        if len(cv):
            A0, A1 = minrad_rows(X, t, cv, a.fd_step * tau_t)
            rows += [-(A0 @ B), -(A1 @ B)]
            rhs += [_at_floor(ra[cv] - tau - buf / 2, at_floor), _at_floor(rb[cv] - tau - buf / 2, at_floor)]
        if not a.allow_edge_collapse:
            Ae, l = edge_rows(X, t)
            Ae = Ae @ B
            rows += [Ae, -Ae]
            rhs += [np.maximum(a.ehi * self.emean, l) - l, l - np.minimum(a.elo * self.emean, l)]
        n_fixed = sum(r.shape[0] for r in rows)
        if not self.box:
            rows += [B, -B]
        Aub = sparse.vstack(rows).tocsr()
        return dict(g=g, c=self.Bt @ g, A=Aub, rhs=np.concatenate(rhs), n_fixed=n_fixed,
                    pairs=len(i), curv=len(cv), tau_eff=tau)

    def solve(self, S, delta):
        if self.box:
            return linprog(S["c"], A_ub=S["A"], b_ub=S["rhs"], bounds=(-delta, delta), method="highs")
        b = np.concatenate([S["rhs"], np.full(S["A"].shape[0] - S["n_fixed"], delta)])
        return linprog(S["c"], A_ub=S["A"], b_ub=b, bounds=(None, None), method="highs")

    def step(self, X, y):
        return self.project(X + (self.B @ y).reshape(self.topo.N, 3))


# --------------------------------------------------------------------------- #
# one pass
# --------------------------------------------------------------------------- #

def run_pass(k, X, P: Problem, oc, a, deadline, trace, save_cur, t0):
    """Iterate SLP steps from X until delta < dmin, the deadline, --iters or a stall."""
    tau = P.tau
    best = P.o0["rop"]
    best_o = dict(P.o0)
    delta = a.delta0 * tau
    dmax, dmin = a.dmax * tau, a.dmin * tau
    st = dict(accepted=0, trials=0, cert_rejects=0, tau_rejects=0, lp_fail=0, soc=0,
              floor_rebuilds=0, cert_margin=math.inf, max_disp=0.0, tau_min=tau, chain=0, stop="")
    hist = [best]
    tau_rows = None          # the rows' target tau: tau_t until an LP at it is infeasible
    L = float(edge_lengths(X, P.topo).sum())
    print(f"pass {k}: start rop {best:.6f} tau {tau:.6f} L {L:.6f} |G| {len(P.mats)} "
          f"eq-dim {P.m} of {3 * P.topo.N} N {P.topo.N}", flush=True)
    it = 0
    while True:
        if it >= a.iters:
            st["stop"] = "iters"
            break
        if time.time() >= deadline:
            st["stop"] = "budget"
            break
        # Rows at the current target (tau_t at first): when an accepted step has dipped tau below
        # it they ask the LP to restore it, which keeps the iterate off the acceptance floor.
        S = P.system(X, delta, tau_rows)
        rebuilt = False
        accepted = False
        while True:
            if time.time() >= deadline:
                st["stop"] = "budget"
                break
            st["trials"] += 1
            lp = P.solve(S, delta)
            if lp.status != 0 and not rebuilt and best_o["tau"] < P.tau_eff(tau_rows):
                # Infeasible because this iterate's tau sits below the target by more than the
                # trust region can restore: rebuild the rows at tau_eff = min(tau_t, tau(X)),
                # where dx = 0 is feasible, instead of halving delta down to dmin with gain left;
                # and make tau_eff the target from here on, so later iterates restore toward a
                # level they can reach (the acceptance floor stays tau_t (1 - tau_slack)).
                S = P.system(X, delta, best_o["tau"])
                tau_rows = P.tau_eff(best_o["tau"])
                rebuilt = True
                st["floor_rebuilds"] += 1
                lp = P.solve(S, delta)
            if lp.status != 0:
                st["lp_fail"] += 1
                delta *= 0.5
                if delta < dmin:
                    break
                continue
            Y = P.step(X, lp.x)
            oy = oc.measure(Y, P.topo)
            tag = ""
            if oy is not None and not a.no_soc:
                S2 = P.system(Y, delta, tau_rows)
                lp2 = P.solve(S2, 0.5 * delta)
                if lp2.status != 0 and rebuilt and oy["tau"] < P.tau_eff(tau_rows):
                    S2 = P.system(Y, delta, oy["tau"])      # this iterate already runs at tau_eff
                    lp2 = P.solve(S2, 0.5 * delta)
                if lp2.status == 0:
                    Y2 = P.step(Y, lp2.x)
                    o2 = oc.measure(Y2, P.topo)
                    if o2 is not None and o2["rop"] < oy["rop"]:
                        Y, oy, tag = Y2, o2, " soc"
            if oy is not None and oy["rop"] < best - 1e-7:
                if oy["tau"] < tau * (1 - a.tau_slack):
                    st["tau_rejects"] += 1
                else:
                    ok, margin, why, nch = isotopy_certificate(X, Y, P.topo)
                    if not ok:
                        st["cert_rejects"] += 1
                        print(f"  it {it}: certificate refused the step ({why}, margin {margin:.3g})",
                              flush=True)
                    else:
                        disp = float(np.linalg.norm(Y - X, axis=1).max()) / tau
                        Ly = float(edge_lengths(Y, P.topo).sum())
                        st["accepted"] += 1
                        st["soc"] += bool(tag)
                        st["chain"] += nch
                        st["cert_margin"] = min(st["cert_margin"], margin)
                        st["max_disp"] = max(st["max_disp"], disp)
                        st["tau_min"] = min(st["tau_min"], oy["tau"])
                        n_acc = st["accepted"]
                        if n_acc <= 5 or n_acc % 100 == 0:
                            print(f"it {it:5d} delta {delta:.2e} lin dL {lp.fun:+.5f} real dL {Ly - L:+.5f} "
                                  f"tau {oy['tau']:.6f} minRad {oy['minRad']:.4f} rop {oy['rop']:.6f} "
                                  f"({oy['rop'] - best:+.5f}) pairs {S['pairs']}{tag} cert {margin:.3g} "
                                  f"disp/tau {disp:.3g} [{time.time() - t0:.0f}s]", flush=True)
                        trace.writerow([k, it, f"{delta:.6e}", f"{oy['rop']:.9f}", f"{oy['tau']:.9f}",
                                        f"{oy['minRad']:.9f}", f"{Ly:.9f}", S["pairs"], S["curv"],
                                        int(bool(tag)), f"{margin:.6g}", f"{disp:.6g}",
                                        f"{time.time() - t0:.2f}"])
                        X, best, best_o, L = Y, oy["rop"], oy, Ly
                        hist.append(best)
                        if a.save_every and n_acc % a.save_every == 0:
                            save_cur(X)
                        delta = min(1.5 * delta, dmax)
                        accepted = True
                        break
            delta *= 0.5
            if delta < dmin:
                break
        it += 1
        if st["stop"]:
            break
        if not accepted and delta < dmin:
            st["stop"] = "dmin"
            break
        w = a.stall_window
        if w and len(hist) > w and hist[-w - 1] - hist[-1] < a.stall_eps:
            st["stop"] = "stall"
            break
    print(f"pass {k}: stop ({st['stop']}) rop {best:.6f} tau {best_o['tau']:.6f} "
          f"minRad {best_o['minRad']:.6f} accepted {st['accepted']} trials {st['trials']} "
          f"cert_rejects {st['cert_rejects']} tau_rejects {st['tau_rejects']} lp_fail {st['lp_fail']} "
          f"floor_rebuilds {st['floor_rebuilds']} [{time.time() - t0:.0f}s]", flush=True)
    return X, best_o, st


# --------------------------------------------------------------------------- #
# diagnose: the first-order feasible cone (fo_cone, in sparse form)
# --------------------------------------------------------------------------- #

def diagnose(X, P: Problem, margin_tau=0.02):
    t, tau, B = P.topo, P.tau, P.B
    g = length_grad(X, t).ravel()
    i, j, s, tt, d, nv = contact_pairs(X, t, P.D + margin_tau * tau, math.pi * tau)
    A = contact_rows(t.N, t.nxt, i, j, s, tt, nv)
    ra, rb = minrad_parts(X, t)
    r = np.minimum(ra, rb)
    act = np.where(r < tau * (1 + 1e-4))[0]
    rows = [A]
    if len(act):
        A0, A1 = minrad_rows(X, t, act, P.a.fd_step * tau)
        rows += [A0, A1]
    Acone = sparse.vstack(rows).tocsr()
    AB = (Acone @ B).tocsr()
    gn = float(np.linalg.norm(g))
    free = float(np.abs(g).sum())
    print(f"diagnose: contact pairs within D + {margin_tau}tau: {len(i)}; curvature-active vertices "
          f"{len(act)}; min Rawdon radius / tau {r.min() / tau:.6f}; eq-dim {P.m}")
    if P.box:
        lp = linprog(P.Bt @ g, A_ub=-AB, b_ub=np.zeros(AB.shape[0]), bounds=(-1, 1), method="highs")
    else:
        Aub = sparse.vstack([-AB, B, -B]).tocsr()
        lp = linprog(P.Bt @ g, A_ub=Aub, b_ub=np.r_[np.zeros(AB.shape[0]), np.ones(2 * B.shape[0])],
                     bounds=(None, None), method="highs")
    box = -lp.fun / free if lp.status == 0 else float("nan")
    print(f"  LP box (|dx_k| <= 1): status {lp.status}; max -g.dx {(-lp.fun if lp.status == 0 else float('nan')):.5f}"
          f" of unconstrained {free:.3f}: ratio {box:.4f}")
    frac = None
    if AB.shape[0] * P.m <= 2e7:
        M = AB.toarray()
        gv = P.Bt @ g
        lam, res = nnls(M.T, gv, maxiter=50 * max(M.shape[0], 1))
        frac = res / gn
        print(f"  equivariant NNLS: |g| {gn:.5f}, unresolved |g - A^T lam| {res:.5f}: feasible first-order "
              f"fraction {100 * frac:.1f}% (about 0 means critical; 13.7% at the 8BL 298.242 plateau)")
    else:
        print(f"  equivariant NNLS skipped: {AB.shape[0]} rows x {P.m} columns > 2e7")
    return box, frac


# --------------------------------------------------------------------------- #
# driver
# --------------------------------------------------------------------------- #

@contextmanager
def _workdir(path):
    if path:
        p = Path(path).resolve()
        p.mkdir(parents=True, exist_ok=True)
        yield p
    else:
        with tempfile.TemporaryDirectory(prefix="slp_tighten_") as td:
            yield Path(td)


@contextmanager
def _chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def _same_file(p, q) -> bool:
    """p and q name one file (same resolved path, or the same inode when both exist)."""
    if p is None or q is None:
        return False
    p, q = Path(p), Path(q)
    if p.resolve() == q.resolve():
        return True
    try:
        return p.exists() and q.exists() and os.path.samefile(p, q)
    except OSError:
        return False


def _vec(text):
    try:
        v = [float(x) for x in text.split(",")]
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected X,Y,Z, got {text!r}")
    if len(v) != 3:
        raise argparse.ArgumentTypeError(f"expected X,Y,Z, got {text!r}")
    return v


def _floats(text):
    try:
        return [float(x) for x in text.split(",") if x.strip()]
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected a comma list of numbers, got {text!r}")


def _userpath(s: str) -> Path:
    """argparse type: a path with ~ expanded, so a typed ~/... means the home folder, not ./~/..."""
    return Path(s).expanduser()


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=_userpath, help="tightened .xyz link")
    ap.add_argument("-o", "--output", type=_userpath, help="output .xyz (required unless --diagnose, --dry-run "
                    "or --group detect)")
    ap.add_argument("--gui", action="store_true",
                    help="open the parameter window (also what running with no arguments does)")
    g = ap.add_argument_group("symmetry")
    g.add_argument("--group", required=True, help="C1|none, Cs, Ci, C<n>, Z/<n>Z, C<n>v, C<n>h, RD<n>, "
                   "D<n> (n>=3), D<n>d[:xy|:diag], S<2n>, or detect (list which of 38 candidate "
                   "groups hold exactly, canonical frame unless --axis/--ref, and stop). "
                   "No default: C1 only when asked")
    g.add_argument("--order", type=int, default=None, help="n for a generic spelling (Cn, Cnv, Cpv, RDp, Z/nZ)")
    g.add_argument("--axis", type=_vec, default=None, help="principal axis X,Y,Z (default 0,0,1)")
    g.add_argument("--ref", type=_vec, default=None, help="reference direction X,Y,Z (default 1,0,0)")
    g.add_argument("--tol", type=float, default=1e-6, help="strict vertex-map tolerance, mean edges (1e-6)")
    g.add_argument("--project-tol", type=float, default=0.05,
                   help="a file within this many mean edges of exact is projected once (default 0.05)")
    g.add_argument("--max-projection-cost", type=float, default=0.05,
                   help="refuse (exit 3) if that projection costs more ropelength than this (0.05)")
    g.add_argument("--components", default=None,
                   help="comma list of components allowed to move; must be a union of orbits")
    r = ap.add_argument_group("run")
    r.add_argument("--max-minutes", type=float, default=15.0, help="wall budget for all passes (15)")
    r.add_argument("--passes", type=int, default=1, help="passes (1); the 8BL recipe uses 2")
    r.add_argument("--resample-vu", type=float, default=0.0,
                   help="equivariant resample to this v/u between passes (0 = off; the 8BL recipe uses 6)")
    r.add_argument("--robustness-vu", type=_floats, default=None,
                   help="comma list of v/u: resample the output to each and report octrope's ropelength")
    r.add_argument("--iters", type=int, default=10 ** 7, help="cap on outer iterations per pass")
    r.add_argument("--no-soc", action="store_true", help="turn off the second-order-correction LP")
    r.add_argument("--stall-window", type=int, default=200, help="accepted steps in the stall window (200)")
    r.add_argument("--stall-eps", type=float, default=5e-4,
                   help="stop a pass when the window gained less ropelength than this (5e-4)")
    r.add_argument("--save-every", type=int, default=100, help="write OUT_cur.xyz every N accepts (100)")
    r.add_argument("--final-tau", type=float, default=0.49997, help="rescale the output to this tau")
    r.add_argument("--work-dir", type=_userpath, default=None, help="scratch directory (default: a temp dir)")
    r.add_argument("--gate-ref", type=_userpath, default=None,
                   help="HOMFLY-gate the output against this file in the cycle's generic frame (exit 4 "
                        "unless SAME)")
    r.add_argument("--diagnose", action="store_true", help="print the first-order feasible cone and exit")
    r.add_argument("--dry-run", action="store_true", help="print the group, dimension and plan and exit")
    c = ap.add_argument_group("step constants (units of the pass's starting tau)")
    c.add_argument("--delta0", type=float, default=0.02, help="initial trust radius (0.02)")
    c.add_argument("--dmax", type=float, default=0.06, help="largest trust radius (0.06)")
    c.add_argument("--dmin", type=float, default=2e-7, help="stop below this trust radius (2e-7)")
    c.add_argument("--buf", type=float, default=0.0, help="extra clearance demanded of contacts (0)")
    c.add_argument("--contact-slack", type=float, default=0.04, help="contact rows within 2tau+3delta+this (0.04)")
    c.add_argument("--curv-window", type=float, default=0.12, help="curvature rows within tau+this (0.12)")
    c.add_argument("--fd-step", type=float, default=2e-7, help="finite-difference step for minRad (2e-7)")
    c.add_argument("--tau-slack", type=float, default=1e-3,
                   help="accept only if tau >= (1 - this) x the pass's starting tau (1e-3)")
    c.add_argument("--elo", type=float, default=0.8, help="edge lower bound, x component mean (0.8)")
    c.add_argument("--ehi", type=float, default=1.25, help="edge upper bound, x component mean (1.25)")
    c.add_argument("--allow-edge-collapse", action="store_true",
                   help="drop the edge bounds. Without them an edge collapsed with HOMFLY still SAME")
    return ap


def _result(fields: dict, note=True) -> str:
    return "SLP RESULT " + " ".join(f"{k}={v}" for k, v in fields.items()) + (" " + NOTE if note else "")


GUI_BLURB = ("The SLP polish finishes a link that RidgeRunner has already tightened and on which "
             "it has stopped improving (a plateau). Pick the tightened .xyz, choose the symmetry "
             "group the run used (or first type detect in the group field and press 'Run now "
             "(short jobs)' to list the exact groups), name the output, and set --gate-ref to the "
             "cycle's reference (results/round0_input_sym.xyz; for a file not made by a cycle, the "
             "input itself) so the topology is checked. Every field has a  ?  that explains it. "
             "'Run now (short jobs)' runs only the quick modes (detect, --dry-run, --diagnose; "
             "--diagnose can take up to a minute, and the window waits) and refuses a real run; "
             "'Run in background' is for a real run (minutes), logged next to the output with "
             ".xyz replaced by .log. Paths may start with ~. Do not hand the output back to "
             "RidgeRunner.")

GUI_EXAMPLES = {
    "input": "~/tightening/my_link/cycle/BEST.xyz\n\nA TIGHTENED file at RidgeRunner scale "
             "(thickness about 0.5) that is exactly symmetric under the group you choose. "
             "A raw drawing is refused: tighten it first.",
    "output": "~/tightening/my_link/slp/polished.xyz\n\nThe run log goes next to it, as "
              "polished.log. Not needed for --group detect, --dry-run or --diagnose.",
    "group": "detect   (press 'Run now' to list the groups that hold exactly)\n\nthen e.g. "
             "C5 for a 5-loop link run under Z/5Z, D2d for the folded 8BL band, Ci, Cs, "
             "C2v, or C1 for no symmetry.",
    "gate_ref": "~/tightening/my_link/cycle/results/round0_input_sym.xyz\n\nThe cycle's own "
                "reference (round0_input.xyz if the cycle did not symmetrize); for a file not made "
                "by a cycle, the input itself. The output is HOMFLY-checked against it and the run "
                "exits 4 if the link changed. It must not be the output file.",
    "max_minutes": "15\n\nMost 15-minute one-pass runs on the project's links stopped on this "
                   "budget while still descending; continue such a run from its output.",
    "passes": "2 with --resample-vu 6: the 8BL D2d recipe, which reached 291.745 from 298.242.",
    "resample_vu": "6\n\nResample between passes to this many vertices per unit ropelength "
                   "(symmetric and certified). 0 = never.",
    "robustness_vu": "4.5,6,8\n\nAfter the run, resample the output to each value and print "
                     "the ropelength, as a check that the gain is not a polygon artefact. "
                     "Compare with the INPUT resampled the same way.",
}


def _gui_log(app):
    """Run in background: log next to the output file."""
    out = app.value("output")
    if not out:
        return None, ("Set -o/--output for a real run: the log is written next to it. "
                      "For --group detect, --dry-run or --diagnose use 'Run now' instead.")
    out = Path(out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    return out.with_suffix(".log"), None


def _gui_quick(app):
    """Run now is only for the modes that finish in seconds."""
    if app.value("group").strip().lower() == "detect":
        return True, None
    flags = [r.dest for r in app.rows if r.is_flag and r.var.get()]
    if "dry_run" in flags or "diagnose" in flags:
        return True, None
    return False, ("Run now is only for --group detect, --dry-run and --diagnose, which finish in "
                   "seconds. A real run takes minutes: use Run in background, which logs next to "
                   "the output and keeps going if you close the window.")


def gui(parser) -> int:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    try:
        from tighten_gui import launch
    except ImportError as exc:          # a Python built without Tk
        print(f"error: no tkinter in this Python ({exc}); run with --help for the command line.",
              file=sys.stderr)
        return EXIT_USAGE
    return launch(parser, Path(__file__).resolve(), blurb=GUI_BLURB, examples=GUI_EXAMPLES,
                  browse={"input": "open", "output": "save", "gate_ref": "open", "work_dir": "dir"},
                  log_for=_gui_log, quick_run=_gui_quick)


def main(argv=None) -> int:
    ap = build_parser()
    raw = sys.argv[1:] if argv is None else list(argv)
    # A launcher, not a run: no arguments ON THE COMMAND LINE, or --gui. main([]) called from code
    # keeps argparse's usage error rather than opening a window.
    if (argv is None and not raw) or "--gui" in raw:
        return gui(ap)
    a = ap.parse_args(argv)
    if a.gui:                           # an abbreviation such as --gu
        return gui(ap)
    t0 = time.time()
    try:
        spec = sm.parse_spec(a.group, a.order)
    except sm.SymmetryMapError as e:
        print(f"error: {e}", file=sys.stderr)
        return EXIT_USAGE
    if spec.kind != "detect" and not (a.diagnose or a.dry_run) and a.output is None:
        ap.error("-o/--output is required (except with --diagnose, --dry-run or --group detect)")
    if not a.allow_edge_collapse and not (0 < a.elo < 1 < a.ehi < math.inf):
        ap.error("edge bounds need 0 < --elo < 1 < --ehi; dropping them requires --allow-edge-collapse")
    if a.passes < 1 or a.max_minutes <= 0:
        ap.error("--passes must be >= 1 and --max-minutes > 0")
    if not a.input.is_file():
        print(f"error: there is no file called {a.input}", file=sys.stderr)
        return EXIT_USAGE
    gating = a.gate_ref is not None and spec.kind != "detect" and not (a.diagnose or a.dry_run)
    if gating:
        if not a.gate_ref.is_file():
            print(f"error: there is no --gate-ref file called {a.gate_ref}", file=sys.stderr)
            return EXIT_USAGE
        if _same_file(a.output, a.gate_ref):
            print(f"error: -o and --gate-ref are the same file ({a.gate_ref}): the output would replace "
                  "the reference and the gate would compare it with itself. Write the output elsewhere",
                  file=sys.stderr)
            return EXIT_USAGE
    if a.output is not None and spec.kind != "detect" and not (a.diagnose or a.dry_run) \
            and _same_file(a.output, a.input):
        print(f"warning: -o is the input file ({a.input}); it is read first and will be overwritten",
              file=sys.stderr, flush=True)
    with _workdir(a.work_dir) as tmp:
        try:
            return run(a, spec, tmp, t0)
        except Refused as e:
            print(f"refused: {e}", file=sys.stderr)
            print(_result({"status": "refused", "reason": e.reason, "group": spec.label,
                           "wall_s": f"{time.time() - t0:.0f}"}, note=False), flush=True)
            return e.code


def run(a, spec, tmp: Path, t0: float) -> int:
    try:
        comps = [np.asarray(c, float) for c in read_xyz(a.input)]
    except TightenError as e:
        print(f"error: {e}", file=sys.stderr)
        return EXIT_USAGE
    topo = Topo([len(c) for c in comps])
    X_in = np.vstack(comps)
    c0 = X_in.mean(0)
    X = X_in - c0
    oc = Octrope(tmp)

    # The gate's reference is read, and its generic-frame copy made, BEFORE anything is written:
    # an output that aliases the reference (refused in main) or any of OUT_*.xyz can then never
    # turn the gate into a comparison of the output with itself.
    ref_gen = None
    if a.gate_ref is not None and spec.kind != "detect" and not (a.diagnose or a.dry_run):
        from tighten_cycle import to_generic_frame
        gdir = tmp / "gate_reference"
        gdir.mkdir(parents=True, exist_ok=True)
        try:
            ref_gen = to_generic_frame(a.gate_ref.resolve(), gdir / "reference_generic.xyz")
        except (TightenError, OSError, ValueError) as e:
            print(f"error: --gate-ref {a.gate_ref}: {e}", file=sys.stderr)
            return EXIT_USAGE

    if spec.kind == "detect":
        hits = sm.detect(topo.split(X), tol=a.tol, axis=a.axis, ref=a.ref)
        print(f"{a.input.name}: exact groups (tol {a.tol} mean edges, centred, canonical frame unless "
              "--axis/--ref), by order:")
        for lab, n, dim in hits:
            print(f"  {lab:9s} |G| {n:3d}   equivariant dim {dim:.0f} of {3 * topo.N}")
        print("pass the one the run used with --group; nothing is chosen for you")
        return EXIT_OK

    o_in = oc.measure(X, topo)
    if o_in is None:
        raise Refused(EXIT_OCTROPE, "octrope", f"octrope could not measure {a.input}")
    rop_in = o_in["rop"]
    el = edge_lengths(X, topo)
    max_dD = float(el.max()) / (2 * o_in["tau"])
    n_struts = oc.struts(X, topo, strut_radius(o_in["tau"], o_in["minStrut"]))
    if n_struts is None:
        raise Refused(EXIT_OCTROPE, "octrope", "could not parse octrope's struts output")
    print(f"input {a.input.name}: rop {rop_in:.6f} tau {o_in['tau']:.6f} minRad {o_in['minRad']:.6f} "
          f"minStrut {o_in['minStrut']:.6f} struts {n_struts} longest edge/D {max_dD:.3f} N {topo.N} "
          f"counts {topo.counts}", flush=True)
    if max_dD > 0.5 or n_struts < 50:
        raise Refused(EXIT_USAGE, "not_tightened",
                      f"{a.input.name} is not tightened (longest edge {max_dD:.3f} D, {n_struts} struts; "
                      "needs <= 0.5 D and >= 50). Descend it first (tighten_cycle.py --initial-steps "
                      "30000, or ridgerunner); the SLP polishes a plateau, it does not untangle a layout")

    # --- exact symmetry, or one certified projection -------------------------
    projection_cost = 0.0
    try:
        sm.vertex_maps(topo.split(X), spec, axis=a.axis, ref=a.ref, tol=a.tol)
    except sm.SymmetryMapError as strict:
        if a.project_tol <= a.tol:
            raise Refused(EXIT_SYMMETRY, "symmetry", f"{spec.label}: {strict}")
        try:
            mats_p, perms_p = sm.vertex_maps(topo.split(X), spec, axis=a.axis, ref=a.ref, tol=a.project_tol)
        except sm.SymmetryMapError as loose:
            raise Refused(EXIT_SYMMETRY, "symmetry",
                          f"{spec.label}: {strict}; and not within --project-tol {a.project_tol}: {loose}")
        Xp = sm.reynolds(X, mats_p, perms_p)
        ok, marg, why, _ = isotopy_certificate(X, Xp, topo)
        if not ok:
            raise Refused(EXIT_SYMMETRY, "symmetry", f"{spec.label}: the projection is not certified ({why})")
        op = oc.measure(Xp, topo)
        if op is None:
            raise Refused(EXIT_OCTROPE, "octrope", "octrope could not measure the projected input")
        projection_cost = op["rop"] - rop_in
        print(f"projected onto exact {spec.label} (move {np.abs(Xp - X).max():.2e}): rop {op['rop']:.6f} "
              f"(cost {projection_cost:+.6f})", flush=True)
        if projection_cost > a.max_projection_cost:
            raise Refused(EXIT_SYMMETRY, "projection_cost",
                          f"projecting onto {spec.label} costs {projection_cost:.4f} > "
                          f"--max-projection-cost {a.max_projection_cost}")
        X = Xp
        try:
            sm.vertex_maps(topo.split(X), spec, axis=a.axis, ref=a.ref, tol=a.tol)
        except sm.SymmetryMapError as e:
            raise Refused(EXIT_SYMMETRY, "symmetry", f"{spec.label} after projection: {e}")

    moving = None
    try:
        P = Problem(X, topo, spec, a, oc)
    except sm.SymmetryMapError as e:
        raise Refused(EXIT_SYMMETRY, "symmetry", f"{spec.label}: {e}")
    if a.components:
        try:
            moving = sorted({int(x) for x in a.components.split(",") if x.strip()})
        except ValueError:
            print(f"error: --components wants a comma list of integers, got {a.components!r}", file=sys.stderr)
            return EXIT_USAGE
        if any(c < 0 or c >= len(topo.counts) for c in moving):
            print(f"error: --components {moving} out of range 0..{len(topo.counts) - 1}", file=sys.stderr)
            return EXIT_USAGE
        for orb in P.info["orbits"]:
            if 0 < len(set(orb) & set(moving)) < len(orb):
                print(f"error: --components must be a union of orbits; orbit {orb} is split", file=sys.stderr)
                return EXIT_USAGE
        P = Problem(X, topo, spec, a, oc, moving=moving)

    G = len(P.mats)
    print(f"group {spec.label}: |G| {G}, equivariant dim {P.m} of {3 * topo.N} (trace R {P.info['dim']:.0f}), "
          f"orbits {P.info['orbits']}" + (f", moving {moving}" if moving else ""), flush=True)

    if a.dry_run:
        print(f"plan: {a.passes} pass(es) in {a.max_minutes:g} min (budget split evenly, unused time rolls "
              f"over); start tau {P.tau:.6f}; delta0 {a.delta0 * P.tau:.3e} dmax {a.dmax * P.tau:.3e} "
              f"dmin {a.dmin * P.tau:.3e}; edge bounds "
              + ("OFF (--allow-edge-collapse)" if a.allow_edge_collapse else f"[{a.elo}, {a.ehi}] x component mean")
              + f"; soc {'off' if a.no_soc else 'on'}")
        if a.passes > 1:
            print("  between passes: " + (f"equivariant resample to v/u {a.resample_vu:g} "
                                          f"(about {round(a.resample_vu * rop_in)} vertices now)"
                                          if a.resample_vu > 0 else "restart from the pass best (no resample)"))
        print(f"  gate: {a.gate_ref if a.gate_ref else 'none (gate the output yourself)'}")
        return EXIT_OK

    if a.diagnose:
        diagnose(X, P)
        return EXIT_OK

    # --- passes ---------------------------------------------------------------
    out = a.output
    out.parent.mkdir(parents=True, exist_ok=True)
    stem = out.with_suffix("")
    # Restore the centre. If the input was symmetric about the origin (its centre's
    # non-fixed part is inside tolerance) restore only the part fixed by the group,
    # so the output is exactly symmetric about the origin, not about c0.
    c_fix = sm.reynolds(c0[None, :], P.mats, [np.zeros(1, int)] * G)[0]
    centre = c_fix if np.linalg.norm(c0 - c_fix) <= a.tol * P.info["edge"] else c0

    def write(path, Y, tp):
        write_xyz(path, [(c + centre).tolist() for c in tp.split(Y)], 15)

    deadline = t0 + 60.0 * a.max_minutes
    tf = open(Path(str(stem) + "_trace.csv"), "w", newline="")
    trace = csv.writer(tf)
    trace.writerow(["pass", "iteration", "delta", "rop", "tau", "minrad", "length", "pairs", "curv_rows",
                    "soc", "cert_margin", "max_disp_over_tau", "t"])
    totals = dict(accepted=0, trials=0, cert_rejects=0, tau_rejects=0, lp_fail=0, soc=0, chain=0,
                  floor_rebuilds=0)
    cert_margin, max_disp, tau_min = math.inf, 0.0, 1.0
    bests = []
    cur, cur_topo = X, topo
    for k in range(1, a.passes + 1):
        now = time.time()
        if now >= deadline:
            print(f"pass {k}: no time left", flush=True)
            break
        if k > 1:
            try:
                P = Problem(cur, cur_topo, spec, a, oc, moving=moving)
            except sm.SymmetryMapError as e:
                print(f"pass {k}: maps fail on the pass input ({e}); stopping passes", flush=True)
                break
        # The input this pass's LP runs on (after any projection or resample): its edges set the
        # floor that min_edge_frac is measured against, and tighten_cycle.py re-reads it to make its
        # own collapse-guard reading when the vertex count has changed.
        write(Path(f"{stem}_pass{k}_input.xyz"), cur, cur_topo)
        budget = (deadline - now) / (a.passes - k + 1)
        Xk, ok_, st = run_pass(k, cur, P, oc, a, now + budget, trace,
                               lambda Y, tp=cur_topo: write(Path(str(stem) + "_cur.xyz"), Y, tp), t0)
        tf.flush()
        for key in totals:
            totals[key] += st[key]
        cert_margin = min(cert_margin, st["cert_margin"])
        max_disp = max(max_disp, st["max_disp"])
        tau_min = min(tau_min, st["tau_min"] / P.tau)
        write(Path(f"{stem}_pass{k}.xyz"), Xk, cur_topo)
        bests.append(dict(rop=ok_["rop"], X=Xk, topo=cur_topo, P=P, k=k))
        if k == a.passes:
            break
        if a.resample_vu > 0:
            total = int(round(a.resample_vu * ok_["rop"]))
            try:
                new, rep = sm.equivariant_resample(cur_topo.split(Xk), spec, total, axis=a.axis, ref=a.ref,
                                                   tol=a.tol)
            except sm.SymmetryMapError as e:
                print(f"resample refused ({e}); stopping passes", flush=True)
                break
            okc, marg, why = certify_resample(cur_topo.split(Xk), new)
            ntop = Topo([len(c) for c in new])
            orr = oc.measure(np.vstack(new), ntop)
            print(f"resample to v/u {a.resample_vu:g}: N {cur_topo.N} -> {ntop.N} (target {total}), rop "
                  f"{ok_['rop']:.6f} -> {orr['rop'] if orr else float('nan'):.6f}; certificate "
                  f"{'ok' if okc else 'FAILED: ' + why} (margin {marg:.3g}); orbits "
                  + ", ".join(f"{o}:{n0}->{n1}" for o, n0, n1, _, _ in rep), flush=True)
            if not okc or orr is None:
                print("resample not certified; stopping passes", flush=True)
                break
            cur, cur_topo = np.vstack(new), ntop
        else:
            cur = Xk
    tf.close()
    if not bests:
        raise Refused(EXIT_USAGE, "no_time", "no pass ran (--max-minutes too small?)")

    # --- result ---------------------------------------------------------------
    b = min(bests, key=lambda r: r["rop"])
    Xb, tb, Pb = b["X"], b["topo"], b["P"]
    ob = oc.measure(Xb, tb)
    if ob is None:
        raise Refused(EXIT_OCTROPE, "octrope", "octrope could not measure the best pass")
    scale = a.final_tau / ob["tau"]
    Xs = Xb * scale
    os_ = oc.measure(Xs, tb)
    if os_ is None:
        raise Refused(EXIT_OCTROPE, "octrope", "octrope could not measure the rescaled output")
    try:
        sm.vertex_maps(tb.split(Xs), spec, axis=a.axis, ref=a.ref, tol=a.tol)
    except sm.SymmetryMapError as e:
        raise Refused(EXIT_SYMMETRY, "symmetry", f"the output lost exact {spec.label}: {e}")
    write(out, Xs, tb)
    # what was written, re-read
    wc = [np.asarray(c, float) for c in read_xyz(out)]
    Xw = np.vstack(wc)
    try:
        mw, pw = sm.vertex_maps(wc, spec, axis=a.axis, ref=a.ref, tol=a.tol)
        sym_err = sm.symmetry_error(Xw, mw, pw)
        sym_frame = "file"
    except sm.SymmetryMapError:
        cw = [c - Xw.mean(0) for c in wc]
        mw, pw = sm.vertex_maps(cw, spec, axis=a.axis, ref=a.ref, tol=a.tol)
        sym_err = sm.symmetry_error(np.vstack(cw), mw, pw)
        sym_frame = "centred"
    elw = edge_lengths(Xw, tb)
    edge_ratio = max(float(elw[tb.comp == c].max() / elw[tb.comp == c].min()) for c in range(len(tb.counts)))
    elb = edge_lengths(Xb, tb)
    on_bounds = 0 if a.allow_edge_collapse else int(((elb >= a.ehi * Pb.emean * (1 - 1e-6))
                                                     | (elb <= a.elo * Pb.emean * (1 + 1e-6))).sum())
    # collapse guard: shortest edge per component over min(elo x mean, shortest) of the pass input the
    # LP actually ran on (>= 1 by construction with the bounds on), and, for information, of the tool's
    # input (below 1 after a resample to a higher v/u, by design)
    floor_pass = np.minimum(a.elo * Pb.comp_mean, Pb.comp_min)
    min_edge_frac = min(float(elb[tb.comp == c].min() / floor_pass[c]) for c in range(len(tb.counts)))
    _, in_mean, in_min = comp_mean_per_vertex(X_in, topo)
    floor_in = np.minimum(a.elo * in_mean, in_min) * (a.final_tau / o_in["tau"])
    min_edge_frac_in = min(float(elw[tb.comp == c].min() / floor_in[c]) for c in range(len(tb.counts)))
    print(f"output {out.name}: rop {os_['rop']:.6f} ({os_['rop'] - rop_in:+.6f}) tau {os_['tau']:.6f} "
          f"minRad {os_['minRad']:.6f} N {tb.N} (pass {b['k']}); symmetry error {sym_err:.1e} ({sym_frame} "
          f"frame); edge ratio {edge_ratio:.3f}; {on_bounds} edges on the bounds", flush=True)
    if edge_ratio > 1.15:
        print("  edge ratio > 1.15: redistribute before the next ridgerunner leg", flush=True)

    robust = []
    for v in (a.robustness_vu or []):
        nv = 0
        try:
            rs, _ = sm.equivariant_resample(tb.split(Xs), spec, int(round(v * os_["rop"])), axis=a.axis,
                                            ref=a.ref, tol=a.tol)
            nv = sum(len(c) for c in rs)
            ro = oc.measure(np.vstack(rs), Topo([len(c) for c in rs]))
        except sm.SymmetryMapError:
            ro = None
        val = ro["rop"] if ro else float("nan")
        robust.append(f"{v:g}:{val:.4f}")
        print(f"  robustness: resampled to v/u {v:g} (N {nv}) rop {val:.6f} ({val - os_['rop']:+.4f})",
              flush=True)

    verdict = "not-run"
    if ref_gen is not None:
        from tighten_cycle import homfly_verdict
        gout = tmp / "gate_candidate"          # its own directory: no name can collide with the reference
        gout.mkdir(parents=True, exist_ok=True)
        out_abs = out.resolve()                # before the chdir: -o may be relative
        with _chdir(tmp):
            verdict = homfly_verdict(ref_gen, out_abs, gout)
        print(f"gate: HOMFLY {verdict} against {a.gate_ref.name} (generic frame, reference read before "
              "any output was written)", flush=True)

    fields = {
        "status": "ok" if verdict in ("SAME", "not-run") else "gate_not_same",
        "group": spec.label, "order": G,
        "rop_in": f"{rop_in:.6f}", "rop_out": f"{os_['rop']:.6f}", "gain": f"{rop_in - os_['rop']:.6f}",
        "tau_out": f"{os_['tau']:.6f}", "minrad_out": f"{os_['minRad']:.6f}",
        "N_in": topo.N, "N_out": tb.N, "passes": len(bests), "best_pass": b["k"],
        "accepted": totals["accepted"], "trials": totals["trials"], "cert_rejects": totals["cert_rejects"],
        "tau_rejects": totals["tau_rejects"], "lp_fail": totals["lp_fail"],
        "floor_rebuilds": totals["floor_rebuilds"], "soc": totals["soc"],
        "on_bounds": on_bounds, "edge_ratio": f"{edge_ratio:.3f}", "min_edge_frac": f"{min_edge_frac:.4f}",
        "min_edge_frac_in": f"{min_edge_frac_in:.4f}",
        "sym_err": f"{sym_err:.1e}", "projection_cost": f"{projection_cost:.6f}",
        "cert_margin": f"{cert_margin:.3g}", "max_disp_over_tau": f"{max_disp:.3g}",
        "tau_min_ratio": f"{tau_min:.6f}", "gate": verdict,
    }
    if robust:
        fields["robust"] = ",".join(robust)
    fields["wall_s"] = f"{time.time() - t0:.0f}"
    print(_result(fields), flush=True)
    return EXIT_OK if verdict in ("SAME", "not-run") else EXIT_GATE


if __name__ == "__main__":
    raise SystemExit(main())
