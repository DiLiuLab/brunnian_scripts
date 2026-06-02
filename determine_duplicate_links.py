#!/usr/bin/env python3
"""Determine duplicate links by comparing SnapPy link exteriors.

This master script replaces determine_duplicate_links*.py. It keeps the
essential method from the original scripts based on Prof. Dunfield's
duplicate_links.txt: compare non-geometric exteriors, randomize triangulations,
and check whether an isomorphism extends to the link.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from pathlib import Path

import snappy


def isosig(manifold):
    return manifold.triangulation_isosig(
        ignore_cusp_ordering=True,
        ignore_curve_orientations=True,
        ignore_orientation=True,
    )


def isomorphic_triangulations_preserving_links(A, B) -> bool:
    isomorphisms = A.isomorphisms_to(B)
    return any(isom.extends_to_link for isom in isomorphisms)


def same_link_exterior_no_geometry(A, B, max_tries: int = 100000) -> bool:
    T = A.exterior(with_hyperbolic_structure=False)
    S = B.exterior(with_hyperbolic_structure=False)
    target = isosig(T)
    for _ in range(max_tries):
        if isosig(S) == target:
            break
        S.randomize()
    return isomorphic_triangulations_preserving_links(T, S)


def read_link_groups(path: Path) -> list[list[str]]:
    groups: list[list[str]] = []
    current: list[str] = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            if current:
                groups.append(current)
                current = []
            continue
        if stripped.startswith("#"):
            continue
        current.append(stripped)
    if current:
        groups.append(current)
    return groups


def compare_group(group_id: int, link_names: list[str], max_tries: int) -> list[dict]:
    rows = []
    checked = set()
    for left, right in itertools.combinations(link_names, 2):
        pair = tuple(sorted((left, right)))
        if pair in checked:
            continue
        row = {"group": group_id, "link_a": pair[0], "link_b": pair[1], "duplicate": "", "error": ""}
        try:
            A = snappy.Link(pair[0])
            B = snappy.Link(pair[1])
            row["duplicate"] = same_link_exterior_no_geometry(A, B, max_tries=max_tries)
        except Exception as exc:
            row["error"] = str(exc)
        rows.append(row)
        checked.add(pair)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["group", "link_a", "link_b", "duplicate", "error"]
    with path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def print_rows(rows: list[dict]) -> None:
    for row in rows:
        if row["error"]:
            print(f"[group {row['group']}] {row['link_a']} vs {row['link_b']}: ERROR: {row['error']}")
        elif row["duplicate"]:
            print(f"[group {row['group']}] {row['link_a']} vs {row['link_b']}: DUPLICATE FOUND")
        else:
            print(f"[group {row['group']}] {row['link_a']} vs {row['link_b']}: not equivalent")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Compare candidate duplicate links.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "--links",
        nargs="+",
        help="One candidate group, with each SnapPy link string quoted if it contains spaces.",
    )
    mode.add_argument(
        "--input-file",
        type=Path,
        help="Text file of candidate groups. Use blank lines to separate groups.",
    )
    parser.add_argument("--max-tries", type=int, default=100000, help="Triangulation randomization tries.")
    parser.add_argument("--output", type=Path, help="Optional CSV output path.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    groups = [args.links] if args.links else read_link_groups(args.input_file)

    rows = []
    for group_id, link_names in enumerate(groups, start=1):
        rows.extend(compare_group(group_id, link_names, max_tries=args.max_tries))

    if args.output:
        write_csv(args.output, rows)
    else:
        print_rows(rows)


if __name__ == "__main__":
    main()
