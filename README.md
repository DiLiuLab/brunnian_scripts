# Brunnian and Borromean Link Screening

This repository contains cleaned SnapPy scripts for determining Brunnian and
Borromean links and for checking candidate duplicate links.

## Master Scripts

- `determine_brunnian_borromean.py`
  - tests Brunnian links, Borromean links, or both
  - accepts a single `--link-string`, an `--input-file`, or `--screen-ht`
  - accepts SnapPy link names and quoted DT-code strings as link strings
  - supports `--method simplify` and `--method nr`
- `determine_duplicate_links.py`
  - checks candidate duplicate links by comparing non-geometric exteriors
  - accepts one group via `--links` or groups from `--input-file`

## Examples

See `Examples/README.md` for runnable examples and input files.

Single link:

```bash
python3 determine_brunnian_borromean.py --link-string L14n63195 --property both --method nr
```

Single DT code:

```bash
python3 determine_brunnian_borromean.py --link-string "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" --property brunnian --method nr
```

HT screening:

```bash
python3 determine_brunnian_borromean.py --screen-ht --property both --method nr --only-matches --output Run_results/ht_screening_brunnian_borromean_nr.csv
```

## Screening Results

`Run_results/` includes HT-table screening outputs for crossings 3 through 14:

- `ht_screening_brunnian_simplify.csv`: 53 Brunnian matches
- `ht_screening_borromean_simplify.csv`: 65 Borromean matches
- `ht_screening_brunnian_borromean_simplify.csv`: combined simplify output
- `ht_screening_brunnian_nr.csv`: 53 Brunnian matches
- `ht_screening_borromean_nr.csv`: 65 Borromean matches
- `ht_screening_brunnian_borromean_nr.csv`: combined nR output

For this range, the `simplify` and `nr` methods produced identical link-name
sets for both Brunnian and Borromean screening.

The `mark` column preserves the original quick duplicate-candidate grouping:
matching stars (`*`, `**`, etc.) indicate rows with the same component count,
crossing number, and rounded volume. Use `determine_duplicate_links.py` for the
stronger exterior-isomorphism check.

## Archived Originals

The original scripts were moved to `previous/`, which is ignored by git. That
folder also contains notes describing how the master scripts replaced the older
files.
