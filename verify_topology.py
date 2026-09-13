#!/usr/bin/env python3
"""Check that tightened .xyz files are still the same link as their original.

Tightening and perturbation can push one strand through another. The result
is entirely self-consistent -- correct thickness, converged ropelength, a
plausible-looking curve -- but it is a different link, and its ropelength is
meaningless as an answer about the original one. The giveaway is usually a
ropelength far below anything the real link can achieve.

Pairwise linking numbers do not catch this for the links in this project:
every pair here has Lk = 0 already, so a change that keeps them unlinked is
invisible. This uses the HOMFLYPT polynomial instead, via plCurve's
`knottype`, which does see such changes.

    python3 verify_topology.py original.xyz candidate1.xyz candidate2.xyz

HOMFLY is a necessary check, not a sufficient one: distinct links can share a
polynomial, so a match is strong evidence rather than proof. A mismatch is
conclusive -- the link changed.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from tighten_link_xyz import TightenError, read_xyz, write_vect

HOMFLY_LINE = re.compile(r"Homfly polynomial:\((.*)\)")


def find_knottype(explicit: str | None) -> str:
    for candidate in (explicit, "knottype"):
        if candidate and shutil.which(candidate):
            return shutil.which(candidate)
    fallback = Path("~/.local/bin/knottype").expanduser()
    if fallback.is_file():
        return str(fallback)
    raise TightenError(
        "knottype not found. It ships with libplCurve; put it on PATH or pass --knottype."
    )


def homfly(path: Path, binary: str, timeout: int, seed: int) -> str:
    """HOMFLYPT polynomial of the link in an .xyz file, as a string."""
    components = read_xyz(path)
    with tempfile.TemporaryDirectory() as tmp:
        vect = Path(tmp) / f"{path.stem}.vect"
        write_vect(vect, components)
        try:
            result = subprocess.run(
                [binary, "-h", "-t", str(timeout), "--seed", str(seed), str(vect)],
                capture_output=True,
                text=True,
                timeout=timeout + 60,
            )
        except subprocess.TimeoutExpired:
            raise TightenError(f"{path.name}: knottype timed out")

    match = HOMFLY_LINE.search(result.stdout)
    if not match:
        tail = (result.stdout + result.stderr).strip().splitlines()[-3:]
        raise TightenError(f"{path.name}: no HOMFLY in knottype output: {tail}")
    return match.group(1).strip()


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("reference", type=Path, help="The link the others should still be.")
    ap.add_argument("candidates", type=Path, nargs="+")
    ap.add_argument("--timeout", type=int, default=300, help="Seconds allowed per HOMFLY (default 300).")
    ap.add_argument("--seed", type=int, default=1, help="knottype random seed; HOMFLY is an invariant, so this only fixes the projection.")
    ap.add_argument("--knottype", help="Path to the knottype binary.")
    args = ap.parse_args(argv)

    try:
        binary = find_knottype(args.knottype)
        reference = homfly(args.reference, binary, args.timeout, args.seed)
    except TightenError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    print(f"reference {args.reference.name}")
    print(f"  HOMFLY {reference}\n")

    changed = 0
    for candidate in args.candidates:
        try:
            poly = homfly(candidate, binary, args.timeout, args.seed)
        except TightenError as exc:
            print(f"  ERROR  {candidate.name}: {exc}")
            changed += 1
            continue
        if poly == reference:
            print(f"  SAME   {candidate.name}")
        else:
            changed += 1
            print(f"  CHANGED {candidate.name}")
            print(f"          HOMFLY {poly}")

    print()
    if changed:
        print(f"{changed} of {len(args.candidates)} differ from the reference and must not be "
              "compared against it on ropelength.")
    else:
        print("All candidates share the reference polynomial; no link change detected.")
    return 1 if changed else 0


if __name__ == "__main__":
    raise SystemExit(main())
