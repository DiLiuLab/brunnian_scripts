#!/usr/bin/env python3
"""Print everything you need to know about a link file, in plain language.

    python3 tighten_lib/measure_link.py mylink.xyz
    python3 tighten_lib/measure_link.py mylink.xyz --symmetry 5

Every number is measured from the file itself with octrope, never read from a
run log. Each line ends with a plain-language verdict so you can tell at a
glance whether a configuration is healthy, and the last block tells you what
the next step should be.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from tighten_cycle import (measure, symmetry_deviation, symmetrize_is_safe,  # noqa: E402
                           hole_wall_radius, marginal_sensitivity)
from tighten_link_xyz import read_xyz  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", type=Path)
    ap.add_argument("--symmetry", type=int, default=None, metavar="N",
                    help="check for N-fold rotational symmetry about z (e.g. 5)")
    a = ap.parse_args()

    if not a.input.is_file():
        print(f"There is no file called {a.input}", file=sys.stderr)
        return 1

    comps = [np.asarray(c, float) for c in read_xyz(a.input)]
    m = measure(a.input)

    print(f"FILE      {a.input.name}")
    print(f"          {len(comps)} components, {m.vertices} points "
          f"({'/'.join(str(len(c)) for c in comps)})")
    print()
    print(f"ROPELENGTH {m.rop:10.3f}   <-- the number we are trying to make smaller")
    print(f"  length   {m.length:10.3f}   total length of all the rope")
    print(f"  thickness{m.tau:10.4f}   radius of the tube; ropelength = length / thickness")
    print()

    # --- health checks, each with a verdict -------------------------------
    print("HEALTH CHECKS")
    lim = "contact (strands touching)" if m.strut_limited else "CURVATURE (a sharp corner)"
    print(f"  what sets the thickness   {lim}")
    if not m.strut_limited:
        print("      ^ a corner is limiting the tube, not strands touching. On a freshly drawn")
        print("        layout this is normal and a descent fixes it. On a file that came out of")
        print("        a squeeze or contraction it means the move was too aggressive.")
    print(f"  minRad / thickness        {m.minrad_over_tau:6.3f}   "
          + ("plenty of room in the corners" if m.minrad_over_tau > 1.2 else
             "getting tight in the corners" if m.minrad_over_tau > 1.05 else
             "almost no room left in the corners"))
    print(f"  contacts (struts)         {m.struts if m.struts is not None else '-':>6}   "
          + ("healthy" if (m.struts or 0) > 200 else "very few -- is this run degenerate?"))
    if m.residual is None:
        print("  residual                       -   not measurable on this file (its thickness")
        print("                                     is not at the standard 0.5), which is normal")
        print("                                     right after a squeeze or contraction")
    else:
        print(f"  residual                  {m.residual:6.3f}   "
              + ("very close to ideal" if m.residual < 0.05 else
                 "near critical" if m.residual < 0.1 else
                 "still far from ideal -- there is room to descend"))
    print(f"  points per unit length    {m.vu:6.2f}   "
          + ("fine" if m.vu >= 4 else
             "coarse -- consider refining" if m.vu >= 2 else
             "VERY coarse; descend first, it improves on its own"))
    print(f"  longest edge / diameter   {m.max_dD:6.3f}   "
          + ("fine" if m.max_dD <= 0.5 else
             "the polygon is cutting corners across its own tube"))
    print()

    # --- symmetry ---------------------------------------------------------
    if a.symmetry:
        dev = symmetry_deviation(a.input, a.symmetry)
        if dev is None:
            print(f"SYMMETRY  this file does not have {a.symmetry} components, so a "
                  f"{a.symmetry}-fold check does not apply")
        else:
            safe, why = symmetrize_is_safe(dev)
            print(f"SYMMETRY  C{a.symmetry} about the z-axis")
            print(f"  worst point is {dev[1]:.3f} half-edges from where symmetry would put it")
            print(f"  {dev[2] * 100:.1f}% of points are further than half an edge")
            print(f"  -> {'SAFE to make it exactly symmetric now' if safe else 'DO NOT symmetrise yet -- descend first'}")
            if not safe:
                print("     (symmetrising a shape that is too far from symmetric can destroy")
                print("      the link. Descend it first, then check again.)")
        print()

    # --- what to do next --------------------------------------------------
    print("WHAT TO DO NEXT")
    hole = hole_wall_radius(a.input, m.tau)
    if m.max_dD > 0.5 or (m.struts or 0) < 50:
        print("  This link is not tightened yet. Run the cycle with --initial-steps 30000")
        print("  so RidgeRunner descends it before any geometric move is attempted.")
    elif not m.strut_limited:
        print("  A corner is limiting the thickness. Let RidgeRunner descend this file to")
        print("  repair it before trying another geometric move.")
    else:
        print("  This looks like a tightened configuration. Empty tunnel through the middle:")
        if hole <= 0.25:
            print(f"  wall radius {hole:.2f} D -- the hole is closed; the squeeze has "
                  f"nothing to take.")
        else:
            # An open hole is NOT a budget. What decides it is whether the hole is
            # held open by slack or by the packing of the strands around it, and
            # only the marginal cost of an infinitesimal squeeze can tell those
            # apart. Reporting the wall radius alone told a reader to squeeze a
            # 7-component structure whose strands were already in mutual contact
            # at a median 1.001 D -- the same structure the cycle driver refuses.
            import tempfile
            with tempfile.TemporaryDirectory() as td:
                sens = marginal_sensitivity(a.input, Path(td), 0.99, 0.5)
            if sens is None:
                print(f"  wall radius {hole:.2f} D, but a probe squeeze removes nothing "
                      f"measurable -- treat the hole as closed.")
            elif sens > 2.30:
                print(f"  wall radius {hole:.2f} D, but squeezing it costs {sens:.2f} of")
                print("  thickness per unit of length gained. The strands around the hole are")
                print("  already packed against each other, so there is nowhere for them to")
                print("  go. NO SQUEEZE -- an open hole on its own is not a budget.")
                print("  The cycle will refuse it too. Ask me before spending time here.")
            else:
                print(f"  wall radius {hole:.2f} D, and squeezing it costs only {sens:.2f} of")
                print("  thickness per unit of length gained -- there is real slack here.")
                print("  Run the cycle with --choose-move and it will pick the factor.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
