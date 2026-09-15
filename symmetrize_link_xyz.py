#!/usr/bin/env python3
"""Find a link's symmetry elements, project onto them, and report the damage.

Given an .xyz link that is approximately symmetric, this finds the symmetry
elements themselves, checks each one against the components, averages the
structure over the group, and writes it back out in a canonical frame.

    centre -> find the elements -> validate each -> match components and
    vertices -> average over the group -> align to z -> check the topology

Groups: Cs (one mirror), Cn (n-fold rotation), Cnv (n-fold rotation plus n
mirrors through the axis; C2v is n=2). Elements are found automatically, or
pinned with --axis and --normal.

Three things decide whether a projection is sound, and all three are easy to
get wrong in ways that look like success.

Correspondence. Before two vertices can be averaged the script has to know
which maps to which -- which component, what offset around the loop, and
which way round it runs. Get it wrong and the output is still exactly
symmetric, just no longer the same link. The arclength report is the check
that catches it: a sound projection barely changes it.

Group structure. The vertex permutations have to compose the way the group
does, so they are built from generators and composed rather than matched
element by element. Matching independently yields maps that are individually
plausible and jointly not a group, and the average then comes back invariant
under one element and stubbornly off under the rest.

Representability. A rotation that carries a component onto itself does so by
a 1/n turn of its index circle, which exists only when the vertex count
divides by n -- a 307-vertex component has no half turn at all -- so
components are resampled to suit before anything is matched.

Topology last. Averaging is a geometric projection with no idea the strands
are solid, so linking numbers are checked before and after. For Brunnian
links every pairwise linking number is zero and the check is vacuous; it
says so, and you should re-identify with the SnapPy scripts.
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

from sono_link_xyz import describe_linking, linking_numbers, resample_closed
from tighten_link_xyz import TightenError, read_xyz, write_xyz


def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues rotation about a unit axis."""
    axis = axis / max(np.linalg.norm(axis), 1e-15)
    x, y, z = axis
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + math.sin(angle) * K + (1.0 - math.cos(angle)) * (K @ K)


def reflection_matrix(normal: np.ndarray) -> np.ndarray:
    """Householder reflection in the plane through the origin with this normal."""
    n = normal / max(np.linalg.norm(normal), 1e-15)
    return np.eye(3) - 2.0 * np.outer(n, n)


def align_to(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotation carrying unit vector `source` onto unit vector `target`."""
    a = source / max(np.linalg.norm(source), 1e-15)
    b = target / max(np.linalg.norm(target), 1e-15)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    if s < 1e-12:
        # Parallel, or antiparallel and needing a half turn about any
        # perpendicular axis.
        if float(a @ b) > 0:
            return np.eye(3)
        perp = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        return rotation_matrix(np.cross(a, perp), math.pi)
    return rotation_matrix(v / s, math.atan2(s, float(a @ b)))


def sphere_directions(count: int) -> np.ndarray:
    """Roughly even directions on the sphere, for multi-start.

    Only half the sphere is needed: n and -n give rotations by the same
    angle in opposite senses, and a structure with one is symmetric under
    the other.
    """
    i = np.arange(count) + 0.5
    phi = np.arccos(1.0 - i / count)  # upper hemisphere
    theta = np.pi * (1.0 + 5.0**0.5) * i
    return np.column_stack(
        [np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)]
    )


def symmetry_error(P: np.ndarray, tree: cKDTree, M: np.ndarray, scale: float) -> float:
    """RMS distance from the transformed cloud back to the original, relative to scale.

    A nearest-neighbour (Chamfer) distance is what makes the axis search
    tractable: it compares the link to its image as point sets, so it needs
    no knowledge of which component maps to which or where each one's
    parametrisation starts. Working that out is a separate problem, and a
    much easier one once the axis is already known.
    """
    moved = P @ M.T
    dist, _ = tree.query(moved, k=1)
    return float(np.sqrt(np.mean(dist**2)) / scale)


def find_axis(P: np.ndarray, order: int, starts: int, refine: int) -> tuple[np.ndarray, float]:
    """Search for the best p-fold rotation axis through the origin."""
    tree = cKDTree(P)
    scale = float(np.sqrt(np.mean(np.sum(P**2, axis=1)))) or 1.0
    angle = 2.0 * math.pi / order

    def score(direction: np.ndarray) -> float:
        return symmetry_error(P, tree, rotation_matrix(direction, angle), scale)

    # Seed with the principal axes as well as the grid; a symmetry axis of
    # a link is very often one of them, and starting there costs nothing.
    _, _, vh = np.linalg.svd(P - P.mean(axis=0), full_matrices=False)
    candidates = np.vstack([vh, sphere_directions(starts)])
    scored = sorted(((score(d), d) for d in candidates), key=lambda t: t[0])

    best_err, best_dir = scored[0]
    for _, direction in scored[:refine]:
        # Two angles are enough to parametrise a direction, and Nelder-Mead
        # suits an objective built on nearest-neighbour distances, which is
        # continuous but not smooth.
        phi0 = math.acos(max(-1.0, min(1.0, direction[2])))
        theta0 = math.atan2(direction[1], direction[0])

        def objective(angles: np.ndarray) -> float:
            phi, theta = angles
            return score(
                np.array(
                    [
                        math.sin(phi) * math.cos(theta),
                        math.sin(phi) * math.sin(theta),
                        math.cos(phi),
                    ]
                )
            )

        out = minimize(objective, [phi0, theta0], method="Nelder-Mead",
                       options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 2000})
        if out.fun < best_err:
            phi, theta = out.x
            best_err = float(out.fun)
            best_dir = np.array(
                [
                    math.sin(phi) * math.cos(theta),
                    math.sin(phi) * math.sin(theta),
                    math.cos(phi),
                ]
            )
    return best_dir / np.linalg.norm(best_dir), best_err


def find_plane(P: np.ndarray, starts: int, refine: int) -> tuple[np.ndarray, float]:
    """Search for the best mirror plane through the origin, returning its normal."""
    tree = cKDTree(P)
    scale = float(np.sqrt(np.mean(np.sum(P**2, axis=1)))) or 1.0

    def score(normal: np.ndarray) -> float:
        return symmetry_error(P, tree, reflection_matrix(normal), scale)

    _, _, vh = np.linalg.svd(P - P.mean(axis=0), full_matrices=False)
    candidates = np.vstack([vh, sphere_directions(starts)])
    scored = sorted(((score(d), d) for d in candidates), key=lambda t: t[0])

    best_err, best_dir = scored[0]
    for _, direction in scored[:refine]:
        phi0 = math.acos(max(-1.0, min(1.0, direction[2])))
        theta0 = math.atan2(direction[1], direction[0])

        def objective(angles: np.ndarray) -> float:
            phi, theta = angles
            return score(
                np.array(
                    [
                        math.sin(phi) * math.cos(theta),
                        math.sin(phi) * math.sin(theta),
                        math.cos(phi),
                    ]
                )
            )

        out = minimize(objective, [phi0, theta0], method="Nelder-Mead",
                       options={"xatol": 1e-10, "fatol": 1e-14, "maxiter": 2000})
        if out.fun < best_err:
            phi, theta = out.x
            best_err = float(out.fun)
            best_dir = np.array(
                [
                    math.sin(phi) * math.cos(theta),
                    math.sin(phi) * math.sin(theta),
                    math.cos(phi),
                ]
            )
    return best_dir / np.linalg.norm(best_dir), best_err


def match_components(components, M: np.ndarray) -> list[int]:
    """Which component each one lands on under M, by nearest-neighbour distance."""
    trees = [cKDTree(np.asarray(c, dtype=float)) for c in components]
    sigma = []
    for comp in components:
        moved = np.asarray(comp, dtype=float) @ M.T
        errs = [float(np.sqrt(np.mean(tree.query(moved, k=1)[0] ** 2))) for tree in trees]
        sigma.append(int(np.argmin(errs)))

    if sorted(sigma) != list(range(len(components))):
        raise TightenError(
            f"the group does not permute the components (got {sigma}). The axis is "
            "probably wrong, or the structure is too far from symmetric to match."
        )
    return sigma


def orbits_of(sigma: list[int]) -> list[list[int]]:
    """Cycles of the component permutation."""
    seen, cycles = set(), []
    for start in range(len(sigma)):
        if start in seen:
            continue
        cycle, node = [], start
        while node not in seen:
            seen.add(node)
            cycle.append(node)
            node = sigma[node]
        cycles.append(cycle)
    return cycles


def best_alignment(A: np.ndarray, B: np.ndarray, involutive: bool = False) -> tuple[np.ndarray, float]:
    """Index map carrying A onto B, with its residual cost.

    A rotation or reflection maps a closed component onto another one as a
    curve, but says nothing about where the stored parametrisation happens to
    start or which way round it runs. Both have to be recovered before the
    two can be averaged vertex by vertex, and getting either wrong silently
    averages unrelated points -- which looks like success, because the result
    is still exactly symmetric, just no longer the same link.

    The map itself is returned rather than an (offset, orientation) pair.
    Rebuilding it from those two downstream is how this went wrong before:
    the search scored candidates built by reversing and rolling the array,
    the caller reconstructed them as (orient * i + k) % n, and for reflected
    components the two differ by one, so every such component was averaged
    against the wrong partners and collapsed.
    """
    n = len(A)
    idx = np.arange(n)
    candidates = []
    if involutive:
        # An involution of the group must get an involution of the indices,
        # or the permutations do not compose into a group action and the
        # average stops being a projection -- it comes out exactly invariant
        # under one element and stuck well off the rest. The reversing maps
        # i -> (k - i) are involutions for every k; of the shifts, only k = 0
        # and the half turn are.
        candidates = [(k - idx) % n for k in range(n)]
        candidates.append(idx.copy())
        if n % 2 == 0:
            candidates.append((idx + n // 2) % n)
    else:
        candidates = [(idx + k) % n for k in range(n)]
        candidates += [(k - idx) % n for k in range(n)]

    best_map, best_cost = idx, float("inf")
    for mapped in candidates:
        cost = float(np.sum((A - B[mapped]) ** 2))
        if cost < best_cost:
            best_map, best_cost = mapped, cost
    return best_map, best_cost


def build_permutation(components, sigma: list[int], M: np.ndarray, involutive: bool = False) -> np.ndarray:
    """Vertex-level permutation induced by the group element."""
    sizes = [len(c) for c in components]
    starts = np.cumsum([0] + sizes)[:-1]
    perm = np.empty(sum(sizes), dtype=int)

    for c, comp in enumerate(components):
        target = sigma[c]
        n = sizes[c]
        if sizes[target] != n:
            raise TightenError(
                f"component {c} maps to {target} but they have {n} and "
                f"{sizes[target]} vertices; resample the orbit to a common count"
            )
        A = np.asarray(comp, dtype=float) @ M.T
        B = np.asarray(components[target], dtype=float)
        mapped, _ = best_alignment(A, B, involutive and sigma[c] == c)
        perm[starts[c] + np.arange(n)] = starts[target] + mapped
    return perm


def symmetrize(P: np.ndarray, perm: np.ndarray, M: np.ndarray, order: int) -> np.ndarray:
    """Average every vertex over its group orbit.

    With M v_i = v_(perm i) for an exactly symmetric link, each vertex can be
    written from any point of its orbit by rotating back, so the mean of
    those p estimates is the projection onto the symmetric subspace.
    """
    total = np.zeros_like(P)
    inverse = M.T  # orthogonal, so the transpose is the inverse
    current = np.arange(len(P))
    for m in range(order):
        total += P[current] @ np.linalg.matrix_power(inverse, m).T
        current = perm[current]

    if not np.array_equal(current, np.arange(len(P))):
        raise TightenError(
            f"the vertex permutation does not have order {order}; the matched "
            "correspondence is inconsistent, so averaging would mix unrelated points"
        )
    return total / order


def _loop_length(points) -> float:
    """Arclength of a closed polyline; the sanity check on a projection."""
    P = np.asarray(points, dtype=float)
    return float(np.linalg.norm(np.diff(np.vstack([P, P[:1]]), axis=0), axis=1).sum())


def orthonormal_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two unit vectors spanning the plane perpendicular to `axis`."""
    seed = np.array([1.0, 0.0, 0.0])
    if abs(float(axis @ seed)) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    u = seed - float(seed @ axis) * axis
    u /= np.linalg.norm(u)
    return u, np.cross(axis, u)


def find_vertical_mirror(P: np.ndarray, axis: np.ndarray, samples: int = 720) -> tuple[np.ndarray, float]:
    """Best mirror whose plane contains `axis`, i.e. whose normal is perpendicular to it.

    Searching only the circle of admissible normals rather than the whole
    sphere is what makes the group structure hold by construction: a mirror
    found freely will sit at some arbitrary angle to the rotation axis and
    the two will not close into a group.
    """
    tree = cKDTree(P)
    scale = float(np.sqrt(np.mean(np.sum(P**2, axis=1)))) or 1.0
    u, v = orthonormal_basis(axis)

    def normal_at(t: float) -> np.ndarray:
        return math.cos(t) * u + math.sin(t) * v

    grid = np.linspace(0.0, math.pi, samples, endpoint=False)
    scores = [(symmetry_error(P, tree, reflection_matrix(normal_at(t)), scale), t) for t in grid]
    best_err, best_t = min(scores)

    out = minimize(
        lambda a: symmetry_error(P, tree, reflection_matrix(normal_at(a[0])), scale),
        [best_t], method="Nelder-Mead",
        options={"xatol": 1e-12, "fatol": 1e-15, "maxiter": 2000},
    )
    if out.fun < best_err:
        best_err, best_t = float(out.fun), float(out.x[0])
    return normal_at(best_t), best_err


def group_elements(kind: str, order: int, axis: np.ndarray, normal: np.ndarray):
    """Every matrix of the requested point group, as (label, matrix)."""
    elements = [("E", np.eye(3))]
    if kind == "Ci":
        # Inversion, r -> -r. Unlike every other element here it has no free
        # parameters: no axis to search for, no plane to orient, only a centre,
        # and for a symmetric structure that has to be the centroid.
        elements.append(("i", -np.eye(3)))
    if kind in ("Cn", "Cnv"):
        for k in range(1, order):
            elements.append((f"C{order}^{k}", rotation_matrix(axis, 2.0 * math.pi * k / order)))
    if kind == "Cs":
        elements.append(("sigma", reflection_matrix(normal)))
    if kind == "Cnv":
        # Mirror planes in Cnv sit 180/n degrees apart, so their normals do too.
        for k in range(order):
            nk = rotation_matrix(axis, math.pi * k / order) @ normal
            elements.append((f"sigma_v{k}", reflection_matrix(nk)))
    return elements


def component_map(components, M: np.ndarray, no_permute: bool) -> list[int]:
    """Which component each one lands on, or the identity when pinned.

    Pinning is the right call whenever the symmetry is known to carry each
    component onto itself. Free matching then has nothing to gain and a real
    way to lose: a drifted component can score marginally better against a
    neighbour than against itself, and the resulting map is not a
    permutation, so a group that is genuinely present gets rejected.
    """
    if no_permute:
        return list(range(len(components)))
    return match_components(components, M)


def inspect_element(components, M: np.ndarray, scale: float, no_permute: bool = False):
    """Symmetry error, the component map, and how clean that match is.

    The component map is the real test, and a plain symmetry error will not
    tell you when it fails. A transform can look passable as a point cloud
    and still send one component onto a different one -- at which point it
    is not a symmetry of the link at all, whatever the error says.
    """
    trees = [cKDTree(np.asarray(c, dtype=float)) for c in components]
    mapping, worst, margin = [], 0.0, float("inf")
    for i, comp in enumerate(components):
        moved = np.asarray(comp, dtype=float) @ M.T
        errs = [float(np.sqrt(np.mean(t.query(moved, k=1)[0] ** 2))) / scale for t in trees]
        if no_permute:
            mapping.append(i)
            worst = max(worst, errs[i])
            others = [e for j, e in enumerate(errs) if j != i]
            if others:
                margin = min(margin, min(others) - errs[i])
        else:
            order = np.argsort(errs)
            mapping.append(int(order[0]))
            worst = max(worst, errs[order[0]])
            if len(errs) > 1:
                margin = min(margin, errs[order[1]] - errs[order[0]])
    valid = sorted(mapping) == list(range(len(components)))
    return mapping, worst, margin, valid


def build_parser() -> argparse.ArgumentParser:
    """Every option, in one place, so the GUI and the CLI cannot drift apart."""
    parser = argparse.ArgumentParser(
        description="Find a link's symmetry elements, project onto them, and report the damage.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input", type=Path, nargs="?", help="Input .xyz path.")
    parser.add_argument("--gui", action="store_true", help="Open the Tkinter GUI. This is also the default when no arguments are given.")
    parser.add_argument("-o", "--output", type=Path, help="Output .xyz path (default: <input>_sym.xyz).")
    parser.add_argument("--group", default="Cs", help="Point group to look for: Cs (one mirror), Ci (inversion, r -> -r), Cn (n-fold rotation), or Cnv (n-fold rotation plus n mirrors through the axis; C2v is n=2). Aliases: Zp=Cn, mirror=Cs.")
    parser.add_argument("-p", "--order", type=int, default=2, help="n, for Cn and Cnv.")
    parser.add_argument("--axis", help='Rotation axis "x,y,z" instead of searching. A symmetric structure often has several valid axes, and the search returns whichever scores best, not the one you had in mind.')
    parser.add_argument("--normal", help='Mirror normal "x,y,z" instead of searching (for Cs; for Cnv it is projected perpendicular to the axis).')
    parser.add_argument("--tolerance", type=float, default=0.10, help="Reject an element whose component match is worse than this relative error, or whose map is not a permutation.")
    parser.add_argument("--no-permute", action="store_true", help="Require every symmetry element to carry each component onto itself, instead of matching components freely. Use when you know the symmetry passes through the components.")
    parser.add_argument("--force", action="store_true", help="Symmetrize over the requested group even if elements fail validation. Expect a large correction and check the topology.")
    parser.add_argument("--beads", type=int, help="Resample every component to this many vertices. Applied automatically when an orbit's components differ.")
    parser.add_argument("--starts", type=int, default=200, help="Directions tried in the coarse search.")
    parser.add_argument("--refine", type=int, default=5, help="Best coarse directions handed to local refinement.")
    parser.add_argument("--decimals", type=int, default=9, help="Decimal places in the output xyz.")
    return parser


def gui_defaults() -> argparse.Namespace:
    """Parse a placeholder CLI so the GUI starts from the documented defaults."""
    return build_parser().parse_args(["__gui_placeholder__"])


def symmetrize_run(args: argparse.Namespace, say=print) -> int:
    """Do the work. `say` receives every line, so the GUI can capture it."""

    kind = {"Zp": "Cn", "mirror": "Cs"}.get(args.group, args.group)
    if kind not in ("Cs", "Ci", "Cn", "Cnv"):
        say(f"error: unknown group {args.group!r}; use Cs, Ci, Cn or Cnv")
        return 1

    try:
        components = read_xyz(args.input)
        say(f"Read {args.input}: {len(components)} component(s), "
              f"{sum(len(c) for c in components)} vertices")

        before_lk = linking_numbers(components)
        P = np.array([p for c in components for p in c], dtype=float)
        centre = P.mean(axis=0)
        components = [[tuple(np.asarray(p) - centre) for p in c] for c in components]
        P = P - centre
        scale = float(np.sqrt(np.mean(np.sum(P**2, axis=1)))) or 1.0

        def parse(text, label):
            vec = np.array([float(v) for v in text.split(",")], dtype=float)
            if vec.shape != (3,) or np.linalg.norm(vec) == 0:
                raise TightenError(f"--{label} needs three numbers and must be non-zero")
            return vec / np.linalg.norm(vec)

        axis = normal = None
        if kind in ("Cn", "Cnv"):
            if args.order < 2:
                raise TightenError("--order must be at least 2")
            if args.axis:
                axis = parse(args.axis, "axis")
            else:
                axis, err = find_axis(P, args.order, args.starts, args.refine)
                say(f"Found C{args.order} axis: [{axis[0]:+.6f} {axis[1]:+.6f} {axis[2]:+.6f}]  err {err:.3e}")
        if kind == "Cnv":
            if args.normal:
                normal = parse(args.normal, "normal")
                normal = normal - float(normal @ axis) * axis
                normal /= np.linalg.norm(normal)
            else:
                normal, err = find_vertical_mirror(P, axis)
                say(f"Found mirror normal: [{normal[0]:+.6f} {normal[1]:+.6f} {normal[2]:+.6f}]  err {err:.3e}")
        if kind == "Cs":
            if args.normal:
                normal = parse(args.normal, "normal")
            else:
                normal, err = find_plane(P, args.starts, args.refine)
                say(f"Found mirror normal: [{normal[0]:+.6f} {normal[1]:+.6f} {normal[2]:+.6f}]  err {err:.3e}")

        if kind == "Ci":
            say("Inversion has no axis or plane to find; the centre is the centroid.")

        elements = group_elements(kind, args.order, axis, normal)
        say(f"\n{kind} has {len(elements)} elements. Validating each against the components:")
        bad = []
        for label, M in elements:
            mapping, worst, margin, valid = inspect_element(components, M, scale, args.no_permute)
            verdict = "ok" if (valid and worst <= args.tolerance) else "REJECTED"
            if verdict == "REJECTED" and label != "E":
                bad.append(label)
            note = "" if valid else "  <- not a permutation"
            say(f"  {label:<10} map {mapping}  worst err {worst:.4f}  margin {margin:.4f}  [{verdict}]{note}")

        if bad and not args.force:
            say(f"\n{kind} is not supported by this structure: {', '.join(bad)} failed.")
            say("A rejected element either maps a component onto a different one, or")
            say("misses by more than --tolerance. Averaging over it would blend")
            say("unrelated strands together rather than remove drift.")
            say("Re-run with a smaller group, or --force if you are sure.")
            return 1

        sigma_by_element = []
        for label, M in elements:
            if label == "E":
                continue
            sigma_by_element.append((label, M, component_map(components, M, args.no_permute)))

        # Components sharing an orbit must share a vertex count.
        sizes = [len(c) for c in components]
        counts = list(sizes)
        reach = {i: {i} for i in range(len(components))}
        for _, _, sig in sigma_by_element:
            for i, j in enumerate(sig):
                reach[i] |= {j}
        changed = True
        while changed:
            changed = False
            for i in range(len(components)):
                merged = set(reach[i])
                for j in list(reach[i]):
                    merged |= reach[j]
                if merged != reach[i]:
                    reach[i] = merged
                    changed = True
        for i, orbit in reach.items():
            target = args.beads or max(sizes[c] for c in orbit)
            for c in orbit:
                counts[c] = target
        # A rotation-like element that fixes a component maps it onto itself by
        # a 1/n turn of the index circle, which only exists when the count
        # divides by n. With 307 vertices and a C2 there is no half turn at
        # all, so the matcher is forced onto some unrelated map and the
        # component collapses. This binds ONLY on self-mapped components: one
        # swapped with a partner is matched against a different index circle
        # and any offset will do.
        turn_order = 2 if kind == "Ci" else args.order
        if kind in ("Cn", "Cnv", "Ci"):
            for _, _, sig in sigma_by_element:
                for c, target in enumerate(sig):
                    if target == c and counts[c] % turn_order:
                        counts[c] += turn_order - (counts[c] % turn_order)
            # Orbit members still have to agree with each other.
            for i, orbit in reach.items():
                target = max(counts[c] for c in orbit)
                for c in orbit:
                    counts[c] = target
        if args.beads or counts != sizes:
            say(f"\nResampling components {sizes} -> {counts}")
            components = [resample_closed(c, n) for c, n in zip(components, counts)]

        P = np.array([p for c in components for p in c], dtype=float)

        # Build the permutations from GENERATORS and compose them. Matching
        # every element independently looks equivalent and is not: the
        # resulting maps need not satisfy pi_gh = pi_g . pi_h, and averaging
        # over an inconsistent action projects onto nothing in particular --
        # the output comes back exactly invariant under a single element and
        # stubbornly off under the rest, however many times it is re-run.
        identity = np.arange(len(P))
        def perm_for(M):
            return build_permutation(
                components, component_map(components, M, args.no_permute), M,
                # Only a true reflection (det -1) needs an involutive index
                # map. A C2 rotation also squares to the identity but carries
                # a component onto itself the other way round -- preserving
                # the direction of travel, by a half-turn of the index circle
                # -- so forcing a reversing map on it mangles the component.
                involutive=bool(
                    np.allclose(M @ M, np.eye(3), atol=1e-9)
                    and np.linalg.det(M) < 0
                ),
            )

        with_perms = []
        if kind in ("Cs", "Ci"):
            S = reflection_matrix(normal) if kind == "Cs" else -np.eye(3)
            with_perms = [(np.eye(3), identity), (S, perm_for(S))]
        else:
            R = rotation_matrix(axis, 2.0 * math.pi / args.order)
            pi_R = perm_for(R)
            powers, cur, M_cur = [], identity, np.eye(3)
            for _ in range(args.order):
                powers.append((M_cur, cur))
                cur, M_cur = pi_R[cur], R @ M_cur
            with_perms = list(powers)
            if kind == "Cnv":
                S = reflection_matrix(normal)
                pi_S = perm_for(S)
                for M_k, pi_k in powers:
                    with_perms.append((S @ M_k, pi_S[pi_k]))

        # A reflection is its own inverse, so its vertex map has to be too. If
        # it is not, the matcher paired a component with itself in an
        # orientation-preserving way and the averaging will not be a
        # projection.
        for M, perm in with_perms:
            if (np.allclose(M @ M, np.eye(3), atol=1e-9) and np.linalg.det(M) < 0
                    and not np.array_equal(perm[perm], identity)):
                say("  warning: a reflection's vertex map is not an involution; "
                      "the projection may be imperfect")

        # q_i = mean over the group of g^-1 applied to the vertex g sends i to.
        total = np.zeros_like(P)
        for M, perm in with_perms:
            total += P[perm] @ M          # rows, so v @ M applies M^-1 = M.T
        Q = total / len(with_perms)

        moved = float(np.max(np.linalg.norm(Q - P, axis=1)))
        say(f"\nLargest vertex move: {moved:.6f}  ({100 * moved / scale:.2f}% of scale)")
        lengths_in = [_loop_length(c) for c in components]
        lengths_out = [
            _loop_length(Q[s:e])
            for s, e in zip(np.cumsum([0] + [len(c) for c in components])[:-1],
                            np.cumsum([0] + [len(c) for c in components])[1:])
        ]
        worst = max(abs(b - a) / a for a, b in zip(lengths_in, lengths_out))
        say("Component arclength, input -> output:")
        for i, (la, lb) in enumerate(zip(lengths_in, lengths_out)):
            say(f"  comp {i}: {la:9.4f} -> {lb:9.4f}  ({100 * (lb - la) / la:+.1f}%)")
        if worst > 0.05:
            say(f"WARNING: arclength changed by up to {100 * worst:.1f}%. A sound projection")
            say("  barely changes it; this much means vertices were averaged against the")
            say("  wrong partners, or the group is not a symmetry of this structure.")
        for label, M in elements:
            if label == "E":
                continue
            say(f"  residual error under {label}: "
                  f"{symmetry_error(Q, cKDTree(Q), M, scale):.3e}")

        # Put the symmetry on z so the output lands in a canonical frame: a
        # rotation axis becomes z, and a lone mirror's normal becomes z so its
        # plane is xy. Nothing downstream requires this, it just makes results
        # from different runs directly comparable.
        if kind == "Ci":
            # Nothing to point anywhere: inversion fixes only the centre, which
            # is already the origin, and every orientation of the structure is
            # equally canonical.
            pass
        else:
            A = align_to(axis if kind in ("Cn", "Cnv") else normal, np.array([0.0, 0.0, 1.0]))

            if kind == "Cnv":
                # Also bring the first mirror normal onto x. ridgerunner takes
                # the azimuth of the first mirror from --SymmetryRef, which
                # defaults to 1,0,0, so a Cnv structure whose mirrors sit at
                # some arbitrary azimuth gets rejected for a reason nothing in
                # the output would explain. Fixing it here lets the two tools
                # compose on their defaults.
                image = A @ normal
                azimuth = math.atan2(image[1], image[0])
                A = rotation_matrix(np.array([0.0, 0.0, 1.0]), -azimuth) @ A

            Q = Q @ A.T

        bounds = np.cumsum([0] + [len(c) for c in components])
        out_components = [[tuple(p) for p in Q[s:e]] for s, e in zip(bounds[:-1], bounds[1:])]

        after_lk = linking_numbers(out_components)
        if before_lk != after_lk:
            say("TOPOLOGY CHANGED -- this result is not usable.")
            say(f"  before: {describe_linking(before_lk)}")
            say(f"  after:  {describe_linking(after_lk)}")
        elif after_lk:
            say(f"Topology preserved: {describe_linking(after_lk)}")
            if all(v == 0 for v in after_lk.values()):
                say("  NOTE: every pairwise linking number is zero, so this check is")
                say("  vacuous -- it cannot see a Brunnian change. That is the normal")
                say("  case for the links this repo screens. Re-identify the result")
                say("  with determine_brunnian_borromean.py before trusting it.")

        out = args.output or args.input.with_name(args.input.stem + "_sym.xyz")
        write_xyz(out, out_components, args.decimals)
        say(f"Wrote {out}")
        if kind == "Ci":
            say("Nothing to align: inversion fixes only the centre, now the origin.")
        else:
            axis_name = "rotation axis" if kind in ("Cn", "Cnv") else "mirror normal"
            say(f"The {axis_name} is aligned to z in the output.")
            if kind == "Cnv":
                say("The first mirror normal is aligned to x.")
    except TightenError as exc:
        say(f"error: {exc}")
        return 1
    return 0


# Help text behind the "?" buttons in the GUI. Kept at module level so it can
# be checked against the parser without starting Tk.
HELP = {
    "input": ("Input .xyz", """One "x y z" per line, with a blank line between components.
Every component is treated as a CLOSED loop, so do not repeat the
first vertex at the end.

A leading label column ("C 1.0 2.0 3.0") and the two-line atom-count
header of a molecular .xyz are both tolerated and skipped."""),

    "output": ("Output .xyz", """Where the symmetrized link is written. Leave blank for
<input>_sym.xyz.

Picking an input auto-fills this with the group appended, e.g.
link_sphere_4BL_wider_noeq_from_best_C2v.xyz, which keeps a
directory of results readable at a glance."""),

    "group": ("Point group", """Cs   one mirror plane.
Ci   inversion, r -> -r about the centroid.
Cn   n-fold rotation about an axis.
Cnv  n-fold rotation PLUS n mirrors that contain the axis.

Ci has no axis or plane to find -- inversion is a fixed matrix with
no free parameters, so the search fields are ignored. It requires an
AMPHICHIRAL link: inversion reverses orientation, so a centrally
symmetric embedding is isotopic to its own mirror image. A chiral
link cannot have it in any embedding.

C2v is Cnv with n = 2: a C2 axis and two perpendicular mirrors
through it. The C2 is the product of the two mirrors, so the three
are not independent -- if the search finds two perpendicular mirrors
whose cross product matches the best C2 axis, the group is real
rather than three coincidences."""),

    "order": ("Order n", """The n in Cn and Cnv. Ignored for Cs.

  n = 2  ->  C2v (two mirrors, one C2)
  n = 3  ->  C3v, and Cn with n = 3 is the Borromean case

Components are resampled so each vertex count divides by n, because
a rotation carrying a component onto itself does so by a 1/n turn of
its index circle. A 307-vertex component has no half turn at all, so
it is resampled to 308."""),

    "axis": ("Rotation axis", """Blank to search for it. Otherwise "x,y,z", e.g.

    1,1,1        the body diagonal
    0,0,1        the z axis

Worth pinning when several valid axes exist. Borromean rings have a
C3 about EVERY body diagonal, so the search returns whichever scores
best, which need not be the one you meant."""),

    "normal": ("Mirror normal", """Blank to search for it. Otherwise "x,y,z" -- the normal to
the plane, not a direction lying in it.

For Cnv it is projected perpendicular to the rotation axis, since a
vertical mirror has to contain that axis for the group to close."""),

    "starts": ("Coarse directions", """How many directions on the sphere the coarse search tries
before local refinement. 200 is ample for a clean structure; raise it
if a symmetry you believe in is being missed.

Ignored when you pin the axis or normal."""),

    "refine": ("Refined candidates", """How many of the best coarse directions get handed to the
local Nelder-Mead refinement. Raising it costs little and helps when
the coarse grid lands near several competing minima."""),

    "no_permute": ("Map each component onto itself", """Require every symmetry element to carry each
component onto itself, instead of matching components freely.

Turn this ON when the symmetry passes through the components. It is
forced whenever the component arclengths all differ, since no
isometry can exchange curves of different length.

Turn it OFF when the symmetry EXCHANGES components. Inversion often
does: on one 4-component link it swaps components 0 and 2 and fixes
1 and 3, and pinning drives the match from 0.003 to 0.644, i.e. from
exact to absent. Two components with near-equal arclength are the
tell that a swap is coming.

Free matching has a real failure mode: a drifted component can score
marginally better against a NEIGHBOUR than against itself, the map
then is not a permutation, and a group that is genuinely present gets
rejected."""),

    "force": ("Symmetrize anyway", """Average over the requested group even when elements were
rejected.

This does not find a hidden symmetry -- it DEFORMS the link until the
missing elements hold. Real cost, measured: forcing C2v on a
structure that only had one mirror took 4.6% off one component's
arclength and moved vertices by 39% of the structure's scale.

Use it when the symmetry is a design target you are imposing. Do not
use it when you are reporting what a structure is."""),

    "tolerance": ("Tolerance", """An element is rejected when its component match is worse
than this, relative to the structure's own scale, or when its
component map is not a permutation.

For calibration, on real structures: a genuine element scores around
0.01 to 0.05, and an absent one around 0.5. The gap is wide, so 0.10
separates them comfortably."""),

    "beads": ("Resample to N beads", """Force every component to this many vertices. Blank
resamples only where it is needed:

  - components sharing an orbit must share a count
  - a component fixed by an n-fold rotation needs a count divisible
    by n

Raising it costs accuracy nothing; it mainly changes output size."""),

    "decimals": ("Output decimals", """Decimal places in the written .xyz. The default of 9 is
far beyond what the geometry justifies and costs only file size."""),
}


def run_gui() -> None:
    """Launch a Tkinter GUI for setting the symmetrization parameters."""
    try:
        import queue
        import threading
        import tkinter as tk
        from tkinter import filedialog
        from tkinter.scrolledtext import ScrolledText
    except Exception as exc:  # pragma: no cover - only when Tkinter is absent.
        raise RuntimeError(f"Tkinter GUI is not available in this Python environment: {exc}")

    defaults = gui_defaults()

    root = tk.Tk()
    root.title("symmetrize_link_xyz: project a link onto its symmetry")
    root.geometry("1020x800")

    fields = {
        "input": tk.StringVar(),
        "output": tk.StringVar(),
        "group": tk.StringVar(value="Cnv"),
        "order": tk.StringVar(value="2"),
        "axis": tk.StringVar(),
        "normal": tk.StringVar(),
        "tolerance": tk.StringVar(value=str(defaults.tolerance)),
        "beads": tk.StringVar(),
        "starts": tk.StringVar(value=str(defaults.starts)),
        "refine": tk.StringVar(value=str(defaults.refine)),
        "decimals": tk.StringVar(value=str(defaults.decimals)),
    }
    flags = {
        "no_permute": tk.BooleanVar(value=True),
        "force": tk.BooleanVar(value=False),
    }

    def show_help(key: str) -> None:
        title, body = HELP[key]
        win = tk.Toplevel(root)
        win.title(f"Help: {title}")
        win.transient(root)
        win.geometry("620x360")
        tk.Label(win, text=title, anchor="w", font=("Helvetica", 14, "bold")).pack(
            fill="x", padx=12, pady=(12, 4)
        )
        text = ScrolledText(win, wrap="word", font=("Courier New", 11),
                            relief="flat", background="#f4f8fb")
        text.pack(fill="both", expand=True, padx=12, pady=(0, 8))
        text.insert("1.0", body)
        text.config(state="disabled")
        tk.Button(win, text="Close", command=win.destroy).pack(pady=(0, 12))
        win.bind("<Escape>", lambda _e: win.destroy())

    def help_dot(parent, key: str):
        """A '?' that opens the explanation for one argument.

        Deliberately a Label and not a Button: on macOS the Aqua theme
        ignores a Button's background, so a light blue one renders grey.
        """
        dot = tk.Label(parent, text=" ? ", fg="#0b5c8a", bg="#bfe3f5",
                       font=("Helvetica", 11, "bold"), cursor="hand2",
                       relief="raised", borderwidth=1)
        dot.bind("<Button-1>", lambda _e, k=key: show_help(k))
        dot.bind("<Enter>", lambda _e: dot.config(bg="#8fd0ee"))
        dot.bind("<Leave>", lambda _e: dot.config(bg="#bfe3f5"))
        return dot

    def add_row(parent, row, label, var, hint="", browse=None, width=None, help_key=None):
        tk.Label(parent, text=label, anchor="w").grid(
            row=row, column=0, sticky="w", padx=(8, 6), pady=3
        )
        entry = tk.Entry(parent, textvariable=var, width=width or (52 if browse else 18))
        entry.grid(row=row, column=1, sticky="w", padx=(0, 6), pady=3)
        if browse is not None:
            tk.Button(parent, text="Browse...", command=browse).grid(
                row=row, column=2, sticky="w", padx=(0, 6), pady=3
            )
        if help_key:
            help_dot(parent, help_key).grid(row=row, column=3, sticky="w", padx=(4, 4), pady=3)
        if hint:
            tk.Label(parent, text=hint, anchor="w", fg="#555555").grid(
                row=row, column=4, sticky="w", padx=(4, 8), pady=3
            )
        return entry

    def pick_input():
        chosen = filedialog.askopenfilename(
            title="Select an .xyz link file",
            filetypes=[("xyz files", "*.xyz"), ("all files", "*.*")],
        )
        if chosen:
            fields["input"].set(chosen)
            if not fields["output"].get().strip():
                source = Path(chosen)
                group = fields["group"].get()
                tag = group if group != "Cnv" else f"C{fields['order'].get()}v"
                fields["output"].set(str(source.with_name(f"{source.stem}_{tag}.xyz")))

    def pick_output():
        chosen = filedialog.asksaveasfilename(
            title="Symmetrized .xyz output", defaultextension=".xyz",
            filetypes=[("xyz files", "*.xyz"), ("all files", "*.*")],
        )
        if chosen:
            fields["output"].set(chosen)

    files = tk.LabelFrame(root, text="Files")
    files.pack(fill="x", padx=10, pady=(10, 4))
    add_row(files, 0, "Input .xyz", fields["input"], browse=pick_input, help_key="input")
    add_row(files, 1, "Output .xyz", fields["output"], browse=pick_output,
            hint="blank = <input>_sym.xyz", help_key="output")

    grp = tk.LabelFrame(root, text="Point group")
    grp.pack(fill="x", padx=10, pady=4)
    tk.Label(grp, text="Group", anchor="w").grid(row=0, column=0, sticky="w", padx=(8, 6), pady=3)
    row = tk.Frame(grp)
    row.grid(row=0, column=1, columnspan=3, sticky="w")
    for value, text in (("Cs", "Cs - one mirror"),
                        ("Ci", "Ci - inversion"),
                        ("Cn", "Cn - n-fold rotation"),
                        ("Cnv", "Cnv - rotation + n mirrors (C2v is n=2)")):
        tk.Radiobutton(row, text=text, variable=fields["group"], value=value).pack(side="left", padx=(0, 12))
    help_dot(row, "group").pack(side="left")
    add_row(grp, 1, "Order n", fields["order"], hint="for Cn and Cnv; ignored for Cs", help_key="order")

    elems = tk.LabelFrame(root, text="Elements (blank = find them automatically)")
    elems.pack(fill="x", padx=10, pady=4)
    add_row(elems, 0, "Rotation axis", fields["axis"], width=28,
            hint='"x,y,z"; pin it when several valid axes exist', help_key="axis")
    add_row(elems, 1, "Mirror normal", fields["normal"], width=28,
            hint='"x,y,z"; for Cnv it is projected perpendicular to the axis', help_key="normal")
    add_row(elems, 2, "Coarse directions", fields["starts"], hint="--starts", help_key="starts")
    add_row(elems, 3, "Refined candidates", fields["refine"], hint="--refine", help_key="refine")

    match = tk.LabelFrame(root, text="Matching")
    match.pack(fill="x", padx=10, pady=4)
    permute_row = tk.Frame(match)
    permute_row.grid(row=0, column=0, columnspan=5, sticky="w", padx=8, pady=2)
    tk.Checkbutton(
        permute_row, anchor="w",
        text="Every element maps each component onto itself (--no-permute)",
        variable=flags["no_permute"],
    ).pack(side="left")
    help_dot(permute_row, "no_permute").pack(side="left", padx=(6, 0))

    force_row = tk.Frame(match)
    force_row.grid(row=1, column=0, columnspan=5, sticky="w", padx=8, pady=2)
    tk.Checkbutton(
        force_row, anchor="w",
        text="Symmetrize even if elements are rejected (--force) - deforms the link to "
             "manufacture a symmetry it does not have",
        variable=flags["force"],
    ).pack(side="left")
    help_dot(force_row, "force").pack(side="left", padx=(6, 0))
    add_row(match, 2, "Tolerance", fields["tolerance"],
            hint="reject an element whose component match is worse than this", help_key="tolerance")
    add_row(match, 3, "Resample to N beads", fields["beads"],
            hint="blank = only as needed; orbits and rotations force some resampling", help_key="beads")
    add_row(match, 4, "Output decimals", fields["decimals"], help_key="decimals")

    buttons = tk.Frame(root)
    buttons.pack(fill="x", padx=10, pady=(6, 2))
    run_button = tk.Button(buttons, text="Symmetrize")
    run_button.pack(side="left")
    tk.Button(buttons, text="Clear log", command=lambda: log_widget.delete("1.0", "end")).pack(side="left", padx=(8, 0))
    tk.Button(buttons, text="Quit", command=root.destroy).pack(side="right")

    status_var = tk.StringVar(value="")
    tk.Label(root, textvariable=status_var, anchor="w").pack(fill="x", padx=10, pady=(2, 2))

    log_widget = ScrolledText(root, height=20, wrap="none", font=("Courier New", 10))
    log_widget.pack(fill="both", expand=True, padx=10, pady=(0, 10))

    log_queue: queue.Queue = queue.Queue()

    def append_log(text: str) -> None:
        log_widget.insert("end", text + "\n")
        log_widget.see("end")

    def optional_number(key, label, cast):
        text = fields[key].get().strip()
        if not text:
            return None
        try:
            return cast(text)
        except ValueError:
            raise TightenError(f"{label}: expected a number, got {text!r}")

    def build_namespace() -> argparse.Namespace:
        source = fields["input"].get().strip()
        if not source:
            raise TightenError("Choose an input .xyz file first.")
        output = fields["output"].get().strip()

        # Start from the parser's own defaults and override only what the GUI
        # exposes, so an option added to the CLI later cannot go missing here.
        settings = gui_defaults()
        settings.input = Path(source)
        settings.output = Path(output) if output else None
        settings.group = fields["group"].get()
        settings.order = optional_number("order", "Order n", int) or 2
        settings.axis = fields["axis"].get().strip() or None
        settings.normal = fields["normal"].get().strip() or None
        settings.tolerance = optional_number("tolerance", "Tolerance", float) or defaults.tolerance
        settings.beads = optional_number("beads", "Resample to N beads", int)
        settings.starts = optional_number("starts", "Coarse directions", int) or defaults.starts
        settings.refine = optional_number("refine", "Refined candidates", int) or defaults.refine
        settings.decimals = optional_number("decimals", "Output decimals", int) or defaults.decimals
        settings.no_permute = flags["no_permute"].get()
        settings.force = flags["force"].get()
        return settings

    def drain() -> None:
        try:
            while True:
                item = log_queue.get_nowait()
                if item is None:
                    run_button.config(state="normal")
                    status_var.set("Done.")
                    return
                append_log(item)
        except queue.Empty:
            pass
        root.after(80, drain)

    def start() -> None:
        try:
            settings = build_namespace()
        except TightenError as exc:
            append_log(f"error: {exc}")
            return
        run_button.config(state="disabled")
        status_var.set("Searching for symmetry elements; the axis search takes a few seconds...")

        def work() -> None:
            try:
                symmetrize_run(settings, say=log_queue.put)
            except Exception as exc:
                log_queue.put(f"error: {exc}")
            finally:
                log_queue.put(None)

        threading.Thread(target=work, daemon=True).start()
        root.after(80, drain)

    run_button.config(command=start)
    append_log("Pick an .xyz file, choose a group, and press Symmetrize.")
    append_log("")
    append_log("Read three things in the output before trusting a result:")
    append_log("  every element marked [ok], arclength change near 0.0%, residuals ~1e-16.")
    append_log("A rejected element means that group is not in your structure.")
    root.mainloop()



def main(argv: list[str] | None = None) -> int:
    """Entry point for CLI and GUI modes."""
    args_in = sys.argv[1:] if argv is None else list(argv)
    args = build_parser().parse_args(args_in)

    if not args_in or args.gui:
        try:
            run_gui()
            return 0
        except Exception as exc:
            print(f"GUI could not be started: {exc}", file=sys.stderr)
            print("Use CLI mode, for example:", file=sys.stderr)
            print("  python3 symmetrize_link_xyz.py link.xyz --group Cnv -p 2 --no-permute",
                  file=sys.stderr)
            return 1

    if args.input is None:
        print("error: CLI mode needs an input .xyz file, unless --gui is used.", file=sys.stderr)
        return 1
    return symmetrize_run(args)


if __name__ == "__main__":
    raise SystemExit(main())
