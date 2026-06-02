#!/usr/bin/env python3
"""Determine Brunnian and Borromean links with SnapPy.

This master script replaces the older determine_brunnian*.py,
determine_borromean_nR.py, and is_brunnian*.py scripts.

Link strings are passed directly to snappy.Link(...), so SnapPy names such as
L14n63195 and DT codes such as "DT: [(-8,-12,16),...]" are both valid when
quoted on the command line or listed in an input text file.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import sys
from pathlib import Path
from typing import Iterable

import snappy


METHOD_NR = "nr"
METHOD_SIMPLIFY = "simplify"
CSV_FIELDNAMES = [
    "source",
    "query",
    "link_name",
    "method",
    "components",
    "crossings",
    "volume",
    "solution_type",
    "mark",
    "is_brunnian",
    "is_borromean",
    "error",
]


def is_unlink_exterior(link) -> bool:
    """Use the num_relators test recommended in the original nR scripts."""
    return link.exterior().fundamental_group().num_relators() == 0


def is_unlink_simplify(link) -> bool:
    """Use diagram simplification, the method from the original non-nR scripts."""
    link_copy = link.copy()
    link_copy.simplify("global")
    return len(link_copy.crossings) == 0


def is_unlink(link, method: str) -> bool:
    if method == METHOD_NR:
        return is_unlink_exterior(link)
    if method == METHOD_SIMPLIFY:
        return is_unlink_simplify(link)
    raise ValueError(f"Unknown unlink method: {method}")


def is_brunnian(link, method: str) -> bool:
    n = len(link.link_components)
    # Must have at least 3 components, and the full link must be nontrivial.
    if n < 3:
        return False
    if is_unlink(link, method):
        return False

    # Test all (n-1)-component sublinks. A Brunnian link becomes an unlink
    # whenever any single component is removed.
    for i in range(n):
        sub = link.sublink([j for j in range(n) if j != i])
        if not is_unlink(sub, method):
            return False
    return True


def is_borromean(link, method: str) -> bool:
    n = len(link.link_components)
    # In these scripts, "Borromean" means the full link is nontrivial and every
    # 2-component sublink is an unlink.
    if n < 3:
        return False
    if is_unlink(link, method):
        return False

    for i, j in itertools.combinations(range(n), 2):
        sub = link.sublink([i, j])
        if not is_unlink(sub, method):
            return False
    return True


def safe_value(callable_obj, default=""):
    try:
        return callable_obj()
    except Exception:
        return default


def link_metadata(link, link_string: str, source: str, method: str) -> dict:
    exterior = safe_value(link.exterior, default=None)
    return {
        "source": source,
        "query": link_string,
        "link_name": safe_value(link.name, default=link_string),
        "method": method,
        "components": len(link.link_components),
        "crossings": len(link.crossings),
        "volume": safe_value(lambda: round(exterior.volume(), 6), default="") if exterior else "",
        "solution_type": safe_value(exterior.solution_type, default="") if exterior else "",
        "mark": "",
    }


def classify_link(link_string: str, method: str, tests: Iterable[str], source: str) -> dict:
    link = snappy.Link(link_string)
    row = link_metadata(link, link_string, source, method)

    if "brunnian" in tests:
        row["is_brunnian"] = is_brunnian(link, method)
    if "borromean" in tests:
        row["is_borromean"] = is_borromean(link, method)
    return row


def classify_ht_exterior(M, crossings: int, method: str, tests: Iterable[str]) -> dict:
    link_name = M.name()
    link = snappy.Link(link_name)
    row = {
        "source": "ht_table",
        "query": link_name,
        "link_name": link_name,
        "method": method,
        "components": M.num_cusps(),
        "crossings": crossings,
        "volume": safe_value(lambda: round(M.volume(), 6), default=""),
        # For hyperbolic links, the original scripts noted that a typical
        # solution_type is "all tetrahedra positively oriented".
        "solution_type": safe_value(M.solution_type, default=""),
        "mark": "",
    }

    if "brunnian" in tests:
        row["is_brunnian"] = is_brunnian(link, method)
    if "borromean" in tests:
        row["is_borromean"] = is_borromean(link, method)
    return row


def read_link_strings(path: Path) -> list[str]:
    links = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            links.append(stripped)
    return links


def requested_tests(value: str) -> list[str]:
    if value == "both":
        return ["brunnian", "borromean"]
    return [value]


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    assign_potential_duplicate_marks(rows)
    with path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(csv_row(row))


def csv_row(row: dict) -> dict:
    return {name: row.get(name, "") for name in CSV_FIELDNAMES}


def row_matches(row: dict) -> bool:
    return row.get("error") or row.get("is_brunnian") is True or row.get("is_borromean") is True


def duplicate_group_key(row: dict):
    if row.get("error"):
        return None
    try:
        volume = round(float(row["volume"]), 6)
    except (KeyError, TypeError, ValueError):
        return None
    return (str(row.get("components", "")), str(row.get("crossings", "")), volume)


def assign_potential_duplicate_marks(rows: list[dict]) -> None:
    """Mark quick candidate duplicate groups by component count, crossings, and volume.

    This preserves the original scripts' rough grouping idea: rows with the same
    number of components, same crossing number, and equal rounded volume are
    marked as potentially identical links. Use determine_duplicate_links.py for
    the stronger exterior-isomorphism check.
    """
    groups: dict[tuple[str, str, float], list[dict]] = {}
    for row in rows:
        row["mark"] = ""
        key = duplicate_group_key(row)
        if key is not None:
            groups.setdefault(key, []).append(row)

    mark_count = 0
    for group_rows in groups.values():
        if len(group_rows) < 2:
            continue
        mark_count += 1
        mark = "*" * mark_count
        for row in group_rows:
            row["mark"] = mark


def annotate_csv_with_marks(path: Path) -> None:
    with path.open(newline="") as csv_file:
        rows = list(csv.DictReader(csv_file))
    assign_potential_duplicate_marks(rows)
    with path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(csv_row(row))


def print_rows(rows: list[dict]) -> None:
    for row in rows:
        if row.get("error"):
            print(f"{row.get('query', '')}: ERROR: {row['error']}")
            continue
        labels = []
        if "is_brunnian" in row:
            labels.append(f"Brunnian={row['is_brunnian']}")
        if "is_borromean" in row:
            labels.append(f"Borromean={row['is_borromean']}")
        print(f"{row['query']}: {', '.join(labels)}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Determine whether links are Brunnian and/or Borromean."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--link-string", help="Single SnapPy link string, including quoted DT codes.")
    mode.add_argument("--input-file", type=Path, help="Text file with one SnapPy link string per line.")
    mode.add_argument("--screen-ht", action="store_true", help="Screen SnapPy HTLinkExteriors.")

    parser.add_argument(
        "--property",
        choices=["brunnian", "borromean", "both"],
        default="both",
        help="Which property to test. Default: both.",
    )
    parser.add_argument(
        "--method",
        choices=[METHOD_NR, METHOD_SIMPLIFY],
        default=METHOD_NR,
        help="Unlink test: 'nr' uses fundamental_group().num_relators(); "
        "'simplify' uses simplify('global') and crossing count. Default: nr.",
    )
    parser.add_argument("--crossings-start", type=int, default=3, help="First HT crossing number.")
    parser.add_argument("--crossings-end", type=int, default=14, help="Last HT crossing number, inclusive.")
    parser.add_argument("--min-components", type=int, default=3, help="Minimum components for HT screening.")
    parser.add_argument("--output", type=Path, help="CSV output path. Required for --screen-ht.")
    parser.add_argument(
        "--only-matches",
        action="store_true",
        help="For CSV output, keep only rows matching at least one requested property.",
    )
    parser.add_argument("--verbose", action="store_true", help="Print HT screening progress to stderr.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    tests = requested_tests(args.property)
    rows: list[dict] = []

    if args.screen_ht and not args.output:
        parser.error("--output is required for --screen-ht")

    if args.link_string:
        try:
            rows.append(classify_link(args.link_string, args.method, tests, "argument"))
        except Exception as exc:
            rows.append({"source": "argument", "query": args.link_string, "method": args.method, "error": str(exc)})

    if args.input_file:
        for link_string in read_link_strings(args.input_file):
            try:
                rows.append(classify_link(link_string, args.method, tests, str(args.input_file)))
            except Exception as exc:
                rows.append({"source": str(args.input_file), "query": link_string, "method": args.method, "error": str(exc)})

    if args.screen_ht:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("w", newline="") as csv_file:
            writer = csv.DictWriter(csv_file, fieldnames=CSV_FIELDNAMES)
            writer.writeheader()
            for crossings in range(args.crossings_start, args.crossings_end + 1):
                if args.verbose:
                    print(f"Screening HT links with {crossings} crossings...", file=sys.stderr, flush=True)
                for M in snappy.HTLinkExteriors(crossings=crossings):
                    try:
                        if M.num_cusps() >= args.min_components:
                            row = classify_ht_exterior(M, crossings, args.method, tests)
                        else:
                            continue
                    except Exception as exc:
                        row = {
                            "source": "ht_table",
                            "query": safe_value(M.name, default=""),
                            "method": args.method,
                            "crossings": crossings,
                            "error": str(exc),
                        }
                    if not args.only_matches or row_matches(row):
                        writer.writerow(csv_row(row))
                        csv_file.flush()
        annotate_csv_with_marks(args.output)
        return

    if args.only_matches:
        rows = [row for row in rows if row_matches(row)]

    if args.output:
        write_csv(args.output, rows)
    else:
        print_rows(rows)


if __name__ == "__main__":
    main()
