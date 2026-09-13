#!/usr/bin/env python3
"""Tighten a knot or link stored as an .xyz file with SONO.

SONO (Shrink-On-No-Overlaps, Pieranski 1998) is a different algorithm from
RidgeRunner's constrained gradient descent, which is exactly what makes it
worth having. Agreement between the two is real evidence that a
configuration is tight; disagreement usually means the RidgeRunner run
stalled at a non-critical point rather than that SONO found something
better.

The model is a chain of hard spheres. Every component is a closed loop of
beads of diameter D, with consecutive centres held at a common segment
length d. The tube of radius D/2 around the chain is the rope, so

    ropelength = length / (D / 2) = 2 * N * d / D

for N beads in total. One iteration restores three constraints

    segment length   |p_i - p_i+1| = d
    curvature        turn angle <= 2 asin(d / D)
    no overlaps      |p_i - p_j| >= D for non-local pairs

and whenever all three hold, d shrinks a little. That is the whole method:
shrink on no overlaps.

THE LENGTH DISTRIBUTION IS FROZEN. Every component shares one d, so
component i has length n_i * d and the ratios L_i : L_j are fixed by the
bead counts chosen at setup. That is a property of the algorithm, not a
shortcoming of this script -- SONO has no mechanism to move length between
components. RidgeRunner minimises total ropelength with per-component
lengths as free variables, so it redistributes freely. Which behaviour you
want depends on the question:

    fixed      the right model when every component must end up the same
               length, e.g. four equal-length strands for synthesis
    free       the right model when you want the true ideal ropelength of
               an asymmetric link, where the tight shape has unequal
               components

Use --beads-per-component to set the frozen ratio explicitly, or
--reallocate to hill-climb over allocations and let the distribution find
its own optimum.

Internal ropelength is a bead-model number. The polygon that comes out is
also measured with octrope's `ropelength` when that binary is available,
which is the figure comparable to RidgeRunner and to the literature.
"""

from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from tighten_link_xyz import TightenError, read_xyz, write_vect, write_xyz

# Bead diameter. Everything is measured in these units, and ropelength is
# scale invariant, so there is no reason to expose this as an option.
BEAD_DIAMETER = 1.0


class Chain:
    """Every component concatenated into one (N, 3) array, plus bookkeeping.

    Keeping all beads in a single array is what makes the overlap query one
    KD-tree build instead of one per component pair, which matters because
    that query dominates the per-iteration cost.
    """

    def __init__(self, components: list[list[tuple[float, float, float]]]):
        self.sizes = [len(c) for c in components]
        self.P = np.array([p for c in components for p in c], dtype=float)

        bounds = np.cumsum([0] + self.sizes)
        self.starts, self.ends = bounds[:-1], bounds[1:]

        # Component id per bead, and the bead's index within its own
        # component. Both are needed to tell a genuine contact from two
        # beads that are merely close because they are neighbours along the
        # chain.
        self.cid = np.repeat(np.arange(len(self.sizes)), self.sizes)
        self.pos = np.concatenate([np.arange(n) for n in self.sizes])
        self.n_of = np.array(self.sizes)[self.cid]

        # Cyclic successor and predecessor, since every component is closed.
        self.nxt = np.empty(len(self.P), dtype=int)
        self.prv = np.empty(len(self.P), dtype=int)
        for s, e in zip(self.starts, self.ends):
            idx = np.arange(s, e)
            self.nxt[idx] = np.roll(idx, -1)
            self.prv[idx] = np.roll(idx, 1)

    def components(self) -> list[list[tuple[float, float, float]]]:
        return [
            [tuple(p) for p in self.P[s:e]] for s, e in zip(self.starts, self.ends)
        ]

    def rebuild(self, components: list[list[tuple[float, float, float]]]) -> None:
        """Adopt a new discretisation of the same curve, in place.

        Callers hold a reference to this object, so a re-discretisation has
        to mutate rather than return a replacement.
        """
        self.__init__(components)

    def component_lengths(self) -> list[float]:
        delta = self.P[self.nxt] - self.P
        seg = np.linalg.norm(delta, axis=1)
        return [float(seg[s:e].sum()) for s, e in zip(self.starts, self.ends)]


def clamped(d: float, diameter: float) -> bool:
    """True when d >= D, where the bead model stops describing a tube.

    turn_max clamps d / diameter at 1 to keep asin finite. Once that clamp
    is active the admissible turn angle is pi, which is to say no curvature
    constraint at all, and beads spaced further apart than they are wide no
    longer form a continuous tube. Runs in this regime look self-consistent
    and quietly pass strands through each other.
    """
    return d >= diameter


def turn_max(d: float, diameter: float) -> float:
    """Largest turn angle a chain of chord d can take inside a tube of radius D/2.

    A circle of radius D/2 subtends 2 asin(d / D) per chord of length d, and
    anything sharper would pinch the tube shut locally.
    """
    return 2.0 * math.asin(min(1.0, d / diameter))


def local_skip(d: float, diameter: float) -> int:
    """Index gap below which two beads are allowed to sit closer than D.

    Beads that are near each other only because they are close along the
    chain must not be pushed apart, or the curvature operator and the
    overlap operator fight each other forever. On the most sharply curved
    admissible arc the chord only reaches D after pi*D/2 of arclength, so
    that is the honest threshold rather than the more obvious D.
    """
    return max(1, int(math.ceil(math.pi * diameter / (2.0 * d))))


def enforce_segments(P: np.ndarray, nxt: np.ndarray, d: float, iters: int = 2) -> None:
    """Pull every consecutive pair back to separation exactly d."""
    for _ in range(iters):
        delta = P[nxt] - P
        dist = np.maximum(np.linalg.norm(delta, axis=1), 1e-12)
        corr = (0.5 * (dist - d) / dist)[:, None] * delta
        # Jacobi style: both endpoints move by half the error, using
        # positions from the start of the sweep.
        np.add.at(P, np.arange(len(P)), corr)
        np.add.at(P, nxt, -corr)


def turn_angles(P: np.ndarray, prv: np.ndarray, nxt: np.ndarray) -> np.ndarray:
    a, b = P[prv] - P, P[nxt] - P
    na = np.maximum(np.linalg.norm(a, axis=1), 1e-12)
    nb = np.maximum(np.linalg.norm(b, axis=1), 1e-12)
    cos_interior = np.clip(np.einsum("ij,ij->i", a, b) / (na * nb), -1.0, 1.0)
    return np.pi - np.arccos(cos_interior)


def enforce_curvature(
    P: np.ndarray, prv: np.ndarray, nxt: np.ndarray, d: float,
    diameter: float, relax: float = 0.5, margin: float = 0.95,
) -> float:
    """Straighten vertices that turn too sharply. Returns the worst violation.

    The correction is the exact one. For a vertex whose two edges both have
    length d, the distance to the midpoint M of its neighbours is
    d * sin(turn / 2), and sliding the vertex a fraction t toward M scales
    that distance by (1 - t) -- so the t that lands exactly on a target
    angle is 1 - sin(target / 2) / sin(turn / 2).

    Driving the vertex by a gain proportional to the raw excess instead
    looks reasonable and is a trap: near the limit the excess is tiny, so a
    mild residual kink heals asymptotically slowly while still counting as
    infeasible, which starves the shrink step and stalls the whole run far
    from the tight configuration.
    """
    limit = turn_max(d, diameter)
    turn = turn_angles(P, prv, nxt)
    bad = turn > limit
    if not np.any(bad):
        return 0.0

    sin_now = np.sin(np.clip(turn, 1e-12, np.pi) / 2.0)
    step = np.zeros(len(P))
    step[bad] = relax * (1.0 - math.sin(margin * limit / 2.0) / sin_now[bad])
    np.clip(step, 0.0, 1.0, out=step)

    mid = 0.5 * (P[prv] + P[nxt])
    P += step[:, None] * (mid - P)
    return float(max(0.0, (turn - limit).max()))


def overlap_pairs(chain: Chain, P: np.ndarray, diameter: float, skip: int):
    """Non-local bead pairs closer than D, as (i, j, distance)."""
    tree = cKDTree(P)
    pairs = tree.query_pairs(r=diameter, output_type="ndarray")
    if len(pairs) == 0:
        return np.empty(0, int), np.empty(0, int), np.empty(0, float)

    i, j = pairs[:, 0], pairs[:, 1]
    same = chain.cid[i] == chain.cid[j]
    gap = np.abs(chain.pos[i] - chain.pos[j])
    gap = np.minimum(gap, chain.n_of[i] - gap)  # closed loops, so gaps wrap
    keep = ~same | (gap > skip)

    i, j = i[keep], j[keep]
    if len(i) == 0:
        return i, j, np.empty(0, float)
    dist = np.linalg.norm(P[j] - P[i], axis=1)
    return i, j, dist


def remove_overlaps(
    chain: Chain, P: np.ndarray, diameter: float, skip: int, relax: float = 1.0
) -> float:
    """Push overlapping non-local pairs apart to D. Returns the worst overlap."""
    i, j, dist = overlap_pairs(chain, P, diameter, skip)
    if len(i) == 0:
        return 0.0
    safe = np.maximum(dist, 1e-12)
    delta = P[j] - P[i]
    push = (relax * 0.5 * (diameter - safe) / safe)[:, None] * delta
    np.add.at(P, i, -push)
    np.add.at(P, j, push)
    return float((diameter - dist).max())


def violations(chain: Chain, P: np.ndarray, d: float, diameter: float) -> tuple[float, float, float]:
    """Worst segment, curvature and overlap residual, all relative, without touching P.

    Relative matters. The admissible turn angle 2 asin(d / D) shrinks along
    with d, so an absolute radian bound silently tightens as the run
    proceeds until nothing can satisfy it and the shrink stops dead.
    """
    seg = np.linalg.norm(P[chain.nxt] - P, axis=1)
    seg_rel = float(np.abs(seg - d).max()) / d

    limit = turn_max(d, diameter)
    turn = turn_angles(P, chain.prv, chain.nxt)
    curv_rel = float(max(0.0, (turn - limit).max())) / limit

    _, _, dist = overlap_pairs(chain, P, diameter, local_skip(d, diameter))
    lap_rel = float(max(0.0, (diameter - dist).max())) / diameter if len(dist) else 0.0
    return seg_rel, curv_rel, lap_rel


def resample_closed(points: list[tuple[float, float, float]], n: int) -> list[tuple[float, float, float]]:
    """Re-place n equally spaced beads along a closed polyline."""
    P = np.asarray(points, dtype=float)
    loop = np.vstack([P, P[:1]])
    seg = np.linalg.norm(np.diff(loop, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = s[-1]
    if total <= 0:
        raise TightenError("a component has zero length")
    want = np.linspace(0.0, total, n, endpoint=False)
    out = np.empty((n, 3))
    for axis in range(3):
        out[:, axis] = np.interp(want, s, loop[:, axis])
    return [tuple(p) for p in out]


def allocate_beads(components, total: int, minimum: int = 12) -> list[int]:
    """Split a bead budget across components in proportion to their length.

    This allocation *is* the frozen length ratio, so it is the single most
    consequential default in the script. Proportional-to-input keeps the
    shape you handed in; pass --beads-per-component to impose your own.
    """
    lengths = []
    for comp in components:
        P = np.asarray(comp, dtype=float)
        loop = np.vstack([P, P[:1]])
        lengths.append(float(np.linalg.norm(np.diff(loop, axis=0), axis=1).sum()))

    scale = total / sum(lengths)
    counts = [max(minimum, int(round(L * scale))) for L in lengths]
    return counts


def autoscale(chain: Chain, diameter: float, d: float, cap: float = 100.0) -> float:
    """Scale the configuration up until no non-local beads overlap.

    Only overlaps are worth scaling for. A merely kinked curve is something
    the curvature operator smooths out by itself at fixed d, whereas scaling
    to clear the single sharpest kink lets a little input noise blow the
    configuration up by orders of magnitude -- and since SONO only ever
    shrinks, every bit of that has to be walked back one step at a time.
    """
    P = chain.P
    _, _, dist = overlap_pairs(chain, P, diameter, local_skip(d, diameter))
    if len(dist) == 0:
        return d

    closest = float(dist.min())
    if closest <= 0:
        raise TightenError("input has coincident beads on non-adjacent strands")
    factor = min(cap, 1.05 * diameter / closest)
    if factor <= 1.0:
        return d
    P *= factor
    return d * factor


def cap_displacement(P: np.ndarray, before: np.ndarray, max_move: float) -> None:
    """Clamp how far any bead travelled, in place.

    This is what keeps the topology. Hard-sphere repulsion is resolved only
    at discrete times, so a bead that jumps further than the tube diameter
    can land on the far side of another strand; the overlap operator then
    cheerfully pushes it apart in its new, wrong position and the link comes
    undone. Nothing downstream notices, because the resulting configuration
    is perfectly self-consistent -- it is just a different link, with a
    correspondingly meaningless ropelength. Capping every step well below D
    removes the tunnelling route entirely.
    """
    delta = P - before
    dist = np.linalg.norm(delta, axis=1)
    hot = dist > max_move
    if np.any(hot):
        P[hot] = before[hot] + delta[hot] * (max_move / dist[hot])[:, None]


def relax_once(
    chain: Chain, d: float, diameter: float, skip: int, inner: int, max_move: float
) -> None:
    """One round of constraint restoration: leashes, curvature, overlaps."""
    P = chain.P
    for _ in range(inner):
        before = P.copy()
        enforce_segments(P, chain.nxt, d)
        enforce_curvature(P, chain.prv, chain.nxt, d, diameter)
        remove_overlaps(chain, P, diameter, skip)
        cap_displacement(P, before, max_move)


def sono(
    chain: Chain,
    steps: int,
    shrink: float,
    inner: int,
    tol: float,
    patience: int,
    min_shrink: float,
    seg_slack: float,
    curv_tol: float,
    settle: int,
    presmooth: float,
    max_move: float,
    diameter: float = BEAD_DIAMETER,
    coarsen_below: float = 0.0,
    coarsen_floor: int = 40,
    coarsen_factor: float = 0.5,
    check_topology: bool = False,
    report=None,
) -> tuple[float, int]:
    """Run the shrink-on-no-overlaps loop. Returns (final d, steps used).

    The shrink lowers the *target* segment length and leaves the leash
    operator to haul the curve in. That is not an implementation detail: the
    leash sweep moves every vertex toward its two neighbours, which is a
    discrete curve-shortening flow, and curve shortening is the only thing
    here that removes wiggle.

    Contracting the coordinates uniformly instead is much faster per step
    and quietly wrong. Uniform scaling preserves angles, so it cannot take a
    wiggle out, while the admissible turn angle 2 asin(d / D) keeps
    tightening as d falls -- the curve crinkles up and locks solid with
    curvature saturated everywhere. On a round unknot that shows up as the
    total turning climbing from 2 pi to many multiples of it.

    The shrink gate is likewise deliberate. Only genuine constraint
    violations -- curvature and overlap -- may block it. The segment
    residual is the shortening flow doing its job, so it gates only through
    a loose slack that keeps the leash from falling too far behind.
    """
    P = chain.P
    seg = np.linalg.norm(P[chain.nxt] - P, axis=1)
    # Start from the tightest component's spacing, not the mean over all
    # beads. A shared d means any component whose beads are currently spaced
    # wider than d has to contract and any spaced narrower has to grow --
    # and growth is precisely what the leash operator is worst at, since it
    # has to push a component outward against everything it is linked
    # through. Taking the minimum leaves every component contracting, which
    # is the motion SONO is built around. It matters most when the bead
    # allocation is forced far from the input's own proportions.
    d = min(
        float(seg[s:e].mean()) for s, e in zip(chain.starts, chain.ends)
    )
    d = autoscale(chain, diameter, d)

    # A run that starts with d >= D has no curvature constraint at all (see
    # clamped()), so it will happily drive strands through each other while
    # reporting a converged-looking ropelength. Say so rather than let the
    # numbers look real.
    if clamped(d, diameter):
        print(
            f"  WARNING: d/D = {d / diameter:.2f} >= 1 after autoscale. The curvature "
            "limit is clamped\n           at pi, so nothing prevents strands passing "
            "through each other and the\n           ropelength below is not "
            "trustworthy. Use more beads.",
            flush=True,
        )

    # Reference for the topology check. Taken before any shrinking, so a
    # later mismatch means this run broke the link rather than inheriting a
    # broken input.
    reference = homfly_signature(chain.components()) if check_topology else None
    if check_topology and reference is None:
        print("  note: could not read a HOMFLY reference, so coarsening will not be "
              "topology-checked", flush=True)

    def admissible() -> tuple[bool, tuple[float, float, float]]:
        seg_err, curv_err, lap_err = violations(chain, P, d, diameter)
        ok = curv_err < curv_tol and lap_err < tol and seg_err < seg_slack
        return ok, (seg_err, curv_err, lap_err)

    # Settle before shrinking. Inputs are rarely admissible -- resampling
    # noise and coarse sampling both leave kinks past the curvature limit --
    # and no amount of shrink-rate tuning fixes that, it simply needs
    # relaxation at fixed d before there is anything to shrink.
    #
    # The curvature operator alone is not enough here, because it only moves
    # vertices that are already over the limit. Coordinate noise comparable
    # to the bead spacing leaves a scatter of kinks that trade places faster
    # than they clear, and the residual then oscillates forever just above
    # tolerance. A plain Laplacian pass over every vertex kills that in a few
    # sweeps; it tapers to nothing so it cannot smooth away real features.
    for k in range(settle):
        if admissible()[0]:
            break
        if presmooth > 0.0:
            before = P.copy()
            weight = presmooth * (1.0 - k / settle)
            mid = 0.5 * (P[chain.prv] + P[chain.nxt])
            P += weight * (mid - P)
            # Smoothing is blind to contacts, so it needs the same cap.
            cap_displacement(P, before, max_move * diameter)
        relax_once(
            chain, d, diameter, local_skip(d, diameter), inner, max_move * diameter
        )

    eps, stall, shrank = shrink, 0, False
    for step in range(1, steps + 1):
        relax_once(
            chain, d, diameter, local_skip(d, diameter), inner, max_move * diameter
        )
        feasible, (seg_err, curv_err, lap_err) = admissible()

        if feasible:
            d *= 1.0 - eps
            stall, shrank = 0, True
            # Let the step recover after a backoff, or one hard patch early
            # on would throttle the whole remaining run.
            eps = min(eps * 1.05, shrink)
        else:
            stall += 1
            # Only a shrink that was actually taken and then failed is
            # evidence the step is too big. Backing off because the curve has
            # not finished settling would stop the run before it starts.
            if stall >= patience and shrank:
                eps *= 0.5
                stall, shrank = 0, False
                if eps < min_shrink:
                    return d, step

        # Adaptive coarsening. d falls as the run proceeds, so a chain that
        # started at a safe d/D becomes progressively over-resolved: every
        # step costs more while the shrink per step is fractional in d and
        # independent of the bead count. Dropping to a fraction of the beads
        # raises d by the reciprocal, which restores d/D and cuts the cost,
        # and it leaves ropelength exactly where it was because
        # 2*N*d/D = 2L/D is invariant under re-discretisation.
        if coarsen_below and d / diameter < coarsen_below:
            sizes = chain.sizes
            halved = [
                max(coarsen_floor, int(round(n * coarsen_factor))) for n in sizes
            ]
            if halved != sizes:
                lengths_before = chain.component_lengths()
                snapshot = chain.components()
                chain.rebuild([
                    resample_closed(c, n) for c, n in zip(snapshot, halved)
                ])
                P = chain.P
                # Same curve, fewer beads: d follows the new spacing so that
                # N*d, and therefore the ropelength, is unchanged.
                d_new = sum(lengths_before) / len(P)

                # Check before continuing. Dropping resolution is exactly when
                # a strand can be cut through, and carrying on from a broken
                # link wastes the rest of the run and hides which step did it.
                if reference is not None:
                    now = homfly_signature(chain.components())
                    if now != reference:
                        chain.rebuild(snapshot)
                        P = chain.P
                        coarsen_below = 0.0  # no further attempts
                        if report:
                            reason = "split into separate components" if now is None \
                                else "changed link type"
                            print(
                                f"  coarsening {sum(sizes)} -> {sum(halved)} beads "
                                f"{reason}; reverted and stopped coarsening",
                                flush=True,
                            )
                        continue

                d = d_new
                if report:
                    print(
                        f"  coarsened {sum(sizes)} -> {len(P)} beads, "
                        f"d/D now {d / diameter:.3f}",
                        flush=True,
                    )

        if report and step % report == 0:
            rope = 2.0 * len(P) * d / diameter
            print(
                f"  step {step:6d}  d {d:.6f}  ropelength {rope:9.4f}  "
                f"shrink {eps:.2e}  seg {seg_err:.1e} curv {curv_err:.1e} lap {lap_err:.1e}",
                flush=True,
            )

    return d, steps


def linking_numbers(components) -> dict[tuple[int, int], int]:
    """Gauss linking number for each pair of components, rounded to an integer.

    Worth computing on every run. A shrink that outruns the contact
    resolution pushes one strand through another, and the configuration that
    results is entirely self-consistent -- correct segment lengths, no
    overlaps, a converged ropelength. It is simply a different link. On the
    Hopf link the giveaway is Lk falling from -1 to 0 while the reported
    ropelength drops below anything the real link can achieve, because what
    is actually being measured is two separate circles.

    This is a necessary check, not a sufficient one: linking numbers are
    blind to Brunnian changes, where every pairwise Lk is zero to begin
    with. Re-identify anything important with the SnapPy scripts.
    """
    out: dict[tuple[int, int], int] = {}
    for a in range(len(components)):
        for b in range(a + 1, len(components)):
            A = np.asarray(components[a], dtype=float)
            B = np.asarray(components[b], dtype=float)
            mid_a, step_a = 0.5 * (A + np.roll(A, -1, axis=0)), np.roll(A, -1, axis=0) - A
            mid_b, step_b = 0.5 * (B + np.roll(B, -1, axis=0)), np.roll(B, -1, axis=0) - B
            r = mid_a[:, None, :] - mid_b[None, :, :]
            norm = np.maximum(np.linalg.norm(r, axis=2), 1e-12)
            cross = np.cross(step_a[:, None, :], step_b[None, :, :])
            total = np.einsum("ijk,ijk->ij", r, cross) / norm**3
            out[(a, b)] = int(round(float(total.sum()) / (4.0 * np.pi)))
    return out


def homfly_signature(components, timeout: int = 120) -> str | None:
    """HOMFLYPT polynomial of a configuration, via plCurve's knottype.

    Linking numbers cannot police these links: every pair is unlinked before
    and after a break, so a change that keeps them unlinked is invisible.
    HOMFLY sees it. Returns None when knottype is unavailable or the link has
    come apart into split components, which knottype refuses outright -- and
    that refusal is itself a break, so callers should treat None as failure
    once a reference has been established.
    """
    binary = shutil.which("knottype") or str(Path("~/.local/bin/knottype").expanduser())
    if not Path(binary).is_file():
        return None
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / "check.vect"
        write_vect(vect, components)
        try:
            out = subprocess.run(
                [binary, "-h", "-t", str(timeout), "--seed", "1", str(vect)],
                capture_output=True, text=True, timeout=timeout + 60,
            )
        except (subprocess.TimeoutExpired, OSError):
            return None
    match = re.search(r"Homfly polynomial:\((.*)\)", out.stdout)
    return match.group(1).strip() if match else None


def describe_linking(lk: dict[tuple[int, int], int]) -> str:
    return "  ".join(f"Lk({a},{b})={v:+d}" for (a, b), v in sorted(lk.items())) or "single component"


def octrope_ropelength(components) -> float | None:
    """Measure the finished polygon with octrope, the independent checker.

    The bead model's own number depends on its discretisation; this one is
    what compares to RidgeRunner and to published values.
    """
    binary = shutil.which("ropelength") or shutil.which(
        str(Path("~/.local/bin/ropelength").expanduser())
    )
    if not binary:
        return None
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / "sono.vect"
        write_vect(vect, components)
        try:
            out = subprocess.run(
                [binary, "-q", str(vect)], capture_output=True, text=True, timeout=300
            )
        except (OSError, subprocess.TimeoutExpired):
            return None
    for token in out.stdout.split():
        try:
            return float(token)
        except ValueError:
            continue
    return None


def report_lengths(chain: Chain, label: str) -> None:
    lengths = chain.component_lengths()
    total = sum(lengths)
    print(f"{label}: {len(lengths)} components, total length {total:.4f}")
    for i, (n, L) in enumerate(zip(chain.sizes, lengths)):
        print(f"  comp {i}: {n:5d} beads  len {L:9.4f}   {100 * L / total:6.2f}% of total")


def solve(components, counts, args, report=None) -> tuple[Chain, float, float]:
    """Resample to the given bead allocation and run SONO. Returns (chain, d, rope)."""
    resampled = [resample_closed(c, n) for c, n in zip(components, counts)]
    chain = Chain(resampled)
    d, _ = sono(
        chain,
        steps=args.steps,
        shrink=args.shrink,
        inner=args.inner,
        tol=args.tol,
        patience=args.patience,
        min_shrink=args.min_shrink,
        seg_slack=args.seg_slack,
        curv_tol=args.curv_tol,
        settle=args.settle,
        presmooth=args.presmooth,
        max_move=args.max_move,
        coarsen_below=getattr(args, "coarsen_below", 0.0),
        coarsen_floor=getattr(args, "coarsen_floor", 40),
        coarsen_factor=getattr(args, "coarsen_factor", 0.5),
        check_topology=getattr(args, "check_topology", False),
        report=report,
    )
    rope = 2.0 * len(chain.P) * d / BEAD_DIAMETER
    return chain, d, rope


def allocation_grid(n_components: int, increment: float, total: int) -> list[list[int]]:
    """Every distribution of length across components on a fixed grid.

    These are ordered compositions: component identity matters, so
    [1, 2, 3, 4] and [4, 3, 2, 1] are different links rather than one
    relabelled. For n components at increment 1/k there are C(k-1, n-1) of
    them -- 84 for four components at 0.1, 969 at 0.05.

    Every component gets at least one increment, so the largest share any
    one can take is 1 - (n-1) * increment: 0.7 for four components at 0.1.
    """
    units = int(round(1.0 / increment))
    if units < n_components:
        raise TightenError(
            f"increment {increment} leaves only {units} units for "
            f"{n_components} components"
        )

    grid: list[list[int]] = []

    def walk(prefix: list[int], remaining: int, left: int) -> None:
        if left == 1:
            grid.append(prefix + [remaining])
            return
        for take in range(1, remaining - (left - 1) + 1):
            walk(prefix + [take], remaining - take, left - 1)

    walk([], units, n_components)

    # Convert unit shares to bead counts, fixing rounding drift on the last
    # component so every allocation spends exactly the same bead budget.
    out = []
    for share in grid:
        counts = [max(1, int(round(total * u / units))) for u in share]
        counts[-1] += total - sum(counts)
        if min(counts) >= 1:
            out.append(counts)
    return out


def sweep_allocations(components, args) -> list[dict]:
    """Run one short SONO solve per grid point and keep every structure.

    This is a seed generator, not a search, and the distinction matters. A
    trial rarely reaches a feasible configuration, and that is fine: its
    value is as a starting point for a proper minimiser, not as an answer.
    The bead-model ropelength of an infeasible run is just the initial d and
    means nothing -- one such run reported 66.68 where octrope measured
    382.25 on the polygon it produced -- so nothing here is ranked or
    discarded on that number. Every structure is written out, and octrope,
    which measures the polygon that actually came out, is recorded alongside.

    An earlier hill-climb over the same moves was actively misleading: it
    accepted whichever move lowered the bead number, which meant climbing
    towards more infeasible allocations, and it kept only its final answer.
    The seed that turned out to be worth having was a discarded waypoint
    part way along that march.
    """
    trial = argparse.Namespace(**vars(args))
    trial.steps = args.realloc_steps

    total = sum(allocate_beads(components, args.beads))
    grid = allocation_grid(len(components), args.sweep_increment, total)

    seed_dir = args.seed_dir or args.input.with_name(args.input.stem + "_seeds")
    seed_dir.mkdir(parents=True, exist_ok=True)

    print(f"Sweeping {len(grid)} allocations at increment {args.sweep_increment} "
          f"({args.realloc_steps} steps each), writing seeds to {seed_dir}")

    results: list[dict] = []
    for index, counts in enumerate(grid, start=1):
        chain, d, rope = solve(components, counts, trial)
        produced = chain.components()
        path = seed_dir / ("seed_" + "-".join(str(c) for c in counts) + ".xyz")
        write_xyz(path, produced, args.decimals)

        measured = octrope_ropelength(produced)
        lk = linking_numbers(produced)
        row = {
            "counts": counts,
            "path": path,
            "bead_rope": rope,
            "octrope": measured,
            "linking": lk,
        }
        results.append(row)
        shown = f"{measured:9.3f}" if measured is not None else "        -"
        print(f"  [{index:>3}/{len(grid)}] {str(counts):<26} "
              f"octrope {shown}   bead {rope:9.3f}", flush=True)

    ranked = [r for r in results if r["octrope"] is not None]
    ranked.sort(key=lambda r: r["octrope"])
    if ranked:
        print()
        print("Best seeds by octrope (the figure that survives an infeasible bead run):")
        for row in ranked[:5]:
            print(f"  {str(row['counts']):<26} octrope {row['octrope']:9.3f}   {row['path'].name}")
    print()
    print("These are seeds, not answers. Refine the promising ones with "
          "tighten_link_xyz.py,")
    print("and check any you act on with verify_topology.py.")
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Tighten a knot or link in an .xyz file with SONO.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", type=Path, help="Input .xyz path.")
    parser.add_argument("-o", "--output", type=Path, help="Output .xyz path (default: <input>_sono.xyz).")

    beads = parser.add_argument_group("bead allocation (this sets the frozen length ratio)")
    beads.add_argument("--beads", type=int, default=600, help="Total bead budget, split proportionally to input component lengths.")
    beads.add_argument("--beads-per-component", help="Explicit per-component counts, comma separated, e.g. 150,150,150,150 to force equal lengths.")
    beads.add_argument("--reallocate", action="store_true", help="Sweep every length-ratio distribution on a grid, writing one seed structure per allocation. This is a generator, not a search: nothing is discarded, and nothing is ranked on the bead-model number.")
    beads.add_argument("--sweep-increment", type=float, default=0.1, help="Grid step for each component's share of the total length. 0.1 gives 84 allocations for four components, 0.05 gives 969. Every component gets at least one increment, so the largest share possible is 1 - (n-1) * increment.")
    beads.add_argument("--seed-dir", type=Path, help="Where to write the swept structures (default: <input>_seeds).")
    beads.add_argument("--realloc-steps", type=int, default=8000, help="SONO steps per swept allocation. These are seeds rather than answers, so a short budget is the point: it buys more of the grid. Raising it does not rescue an allocation that never takes a feasible step.")

    solver = parser.add_argument_group("solver")
    solver.add_argument("-s", "--steps", type=int, default=40000, help="Maximum shrink steps.")
    solver.add_argument("--shrink", type=float, default=1e-3, help="Initial fractional shrink of d per feasible step.")
    solver.add_argument("--min-shrink", type=float, default=1e-7, help="Stop once the adaptive shrink falls below this.")
    solver.add_argument("--inner", type=int, default=3, help="Constraint-restoration sweeps per step.")
    solver.add_argument("--tol", type=float, default=2e-3, help="Relative overlap tolerance; this one sets the thickness, so keep it tight.")
    solver.add_argument("--curv-tol", type=float, default=2e-2, help="Relative turn-angle overshoot tolerated. Looser than --tol on purpose: the curvature sweep updates every vertex at once, so neighbours keep knocking each other back over a tight bound.")
    solver.add_argument("--coarsen-below", type=float, default=0.0, help="Halve the bead count whenever d/D falls below this, resampling the same curve. d only falls as a run proceeds, so a chain that started safe becomes over-resolved and slow; halving restores d/D and halves the cost per step, and leaves ropelength unchanged because 2*N*d/D = 2L/D. 0 disables. Try 0.3.")
    solver.add_argument("--check-topology", action="store_true", help="Verify the link type with HOMFLY after every coarsening, reverting the step and stopping further coarsening if it broke. Costs about a second per check and saves the rest of a run that would otherwise continue from a broken link.")
    solver.add_argument("--coarsen-factor", type=float, default=0.5, help="Fraction of the beads kept at each coarsening. 0.5 halves; 0.7 is gentler and gives more, smaller steps before d/D reaches the degenerate regime. Ratios between components are preserved.")
    solver.add_argument("--coarsen-floor", type=int, default=40, help="Never take a component below this many beads when coarsening.")
    solver.add_argument("--max-move", type=float, default=0.05, help="Largest distance a bead may move per sweep, as a fraction of the tube diameter. This is the guard against strands passing through each other; raise it only if you verify the topology afterwards.")
    solver.add_argument("--presmooth", type=float, default=0.15, help="Laplacian smoothing applied during settle, tapering to zero. Set 0 for input you know is already clean.")
    solver.add_argument("--settle", type=int, default=5000, help="Relaxation rounds at fixed d before shrinking starts.")
    solver.add_argument("--seg-slack", type=float, default=0.02, help="How far the leash operator may lag before the shrink pauses.")
    solver.add_argument("--patience", type=int, default=12, help="Infeasible steps tolerated before backing off the shrink.")
    solver.add_argument("--report", type=int, default=1000, help="Print progress every N steps (0 to silence).")
    solver.add_argument("--decimals", type=int, default=9, help="Decimal places in the output xyz.")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        components = read_xyz(args.input)
        print(f"Read {args.input}: {len(components)} component(s), "
              f"{sum(len(c) for c in components)} vertices")

        if args.beads_per_component:
            counts = [int(v) for v in args.beads_per_component.split(",")]
            if len(counts) != len(components):
                raise TightenError(
                    f"--beads-per-component has {len(counts)} entries "
                    f"but the file has {len(components)} components"
                )
        else:
            counts = allocate_beads(components, args.beads)

        before_lk = linking_numbers(components)

        if args.reallocate:
            # A sweep produces one structure per allocation rather than a
            # single answer, so it finishes here.
            sweep_allocations(components, args)
            return 0

        print(f"Running SONO with bead counts {counts}:")
        chain, d, rope = solve(components, counts, args, report=args.report or None)

        report_lengths(chain, "Tightened")
        print(f"Bead-model ropelength: {rope:.4f}  (segment length d = {d:.6f})")

        after = linking_numbers(chain.components())
        if before_lk != after:
            print("TOPOLOGY CHANGED -- this result is not usable.")
            print(f"  before: {describe_linking(before_lk)}")
            print(f"  after:  {describe_linking(after)}")
            print("  A strand passed through another during the run. Re-run with a")
            print("  smaller --max-move and/or a smaller --shrink.")
        elif after:
            print(f"Topology preserved: {describe_linking(after)}")

        measured = octrope_ropelength(chain.components())
        if measured is None:
            print("octrope `ropelength` not found, so no independent check was made.")
        else:
            print(f"octrope ropelength:    {measured:.4f}  <- the comparable figure")

        out = args.output or args.input.with_name(args.input.stem + "_sono.xyz")
        write_xyz(out, chain.components(), args.decimals)
        print(f"Wrote {out}")
    except TightenError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
