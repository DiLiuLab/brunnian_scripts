#!/usr/bin/env python3
"""Scale one component of a link, to test whether it is jammed at the wrong length.

RidgeRunner's --MangleMode perturbs a link with a global torus rotation, which
never specifically trades length between components. This does that directly:
it scales one component about its own centroid, so that component gets longer
(or shorter) while the others are left alone. Re-tightening afterwards shows
whether the tightening had pushed that component to the wrong size.

Growing a component is the safe direction -- it gets looser relative to its
neighbours. Either way the scaling is swept continuously and the minimum
distance between components is checked along the way: if that distance never
approaches zero, no strand passed through another and the link type is intact.

    python3 perturb_component.py in.xyz out.xyz --component 2 --scale 2.5
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# tighten_link_xyz lives in the repository root, one level up from this
# package. Resolved from __file__ so it works wherever the repo is checked out.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import tighten_link_xyz as t  # noqa: E402


def component_lengths(components) -> list[float]:
    out = []
    for comp in components:
        a = np.asarray(comp)
        out.append(float(np.linalg.norm(a - np.roll(a, -1, axis=0), axis=1).sum()))
    return out


def min_inter_component_distance(components, only: int | None = None) -> float:
    """Closest approach between two different components.

    With `only` set, just the pairs involving that component. The global
    minimum is useless for this check: in a tight link some other pair is
    already at contact distance, so it never moves when one component is
    scaled, and would mask a strand passing through.
    """
    arrays = [np.asarray(c) for c in components]
    best = float("inf")
    for i in range(len(arrays)):
        for j in range(i + 1, len(arrays)):
            if only is not None and only not in (i, j):
                continue
            d = np.linalg.norm(arrays[i][:, None, :] - arrays[j][None, :, :], axis=2)
            best = min(best, float(d.min()))
    return best


def scaled(components, index: int, factor: float):
    """Scale one component about its own centroid, leaving the others alone."""
    out = [list(c) for c in components]
    a = np.asarray(components[index])
    centre = a.mean(axis=0)
    out[index] = [tuple(p) for p in (centre + (a - centre) * factor)]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("input", type=Path)
    ap.add_argument("output", type=Path)
    ap.add_argument("--component", type=int, required=True, help="0-based component index.")
    ap.add_argument("--scale", type=float, required=True, help="Length factor, e.g. 2.5.")
    ap.add_argument("--samples", type=int, default=40, help="Sweep steps checked for pass-through.")
    ap.add_argument("--decimals", type=int, default=9)
    args = ap.parse_args()

    components = t.read_xyz(args.input)
    if not 0 <= args.component < len(components):
        print(f"error: component {args.component} out of range (0-{len(components)-1})", file=sys.stderr)
        return 1

    before = component_lengths(components)
    total = sum(before)
    start_gap = min_inter_component_distance(components, only=args.component)

    # Sweep the scaling and watch the closest approach between components. A
    # strand passing through another would drive this to zero on the way.
    worst = start_gap
    for k in range(1, args.samples + 1):
        factor = 1.0 + (args.scale - 1.0) * k / args.samples
        worst = min(
            worst,
            min_inter_component_distance(
                scaled(components, args.component, factor), only=args.component
            ),
        )

    result = scaled(components, args.component, args.scale)
    after = component_lengths(result)
    total_after = sum(after)

    print(f"component {args.component}: length {before[args.component]:.3f} -> "
          f"{after[args.component]:.3f}  ({100*before[args.component]/total:.2f}% -> "
          f"{100*after[args.component]/total_after:.2f}% of total)")
    print("  share now: " + "  ".join(f"{100*v/total_after:5.2f}%" for v in after))
    print(f"  clearance of component {args.component} to the others: {start_gap:.4f} "
          f"at the start, {worst:.4f} at its worst during the sweep")
    if worst <= 1e-6:
        print("  REFUSING: components touch during the sweep, so the link type would change.",
              file=sys.stderr)
        return 1
    if worst < 0.25 * start_gap:
        print("  WARNING: the gap closed a lot; treat the resulting link type as unverified.")

    t.write_xyz(args.output, result, args.decimals)
    print(f"  wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
