# Examples

Run a single link name:

```bash
python3 determine_brunnian_borromean.py --link-string L14n63195 --property both --method nr
```

Run a single DT code. Quote the whole DT string so the shell passes it as one
argument:

```bash
python3 determine_brunnian_borromean.py --link-string "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" --property brunnian --method nr
```

Run a text file of link strings:

```bash
python3 determine_brunnian_borromean.py --input-file Examples/brunnian_borromean_links.txt --property both --method nr
```

Screen SnapPy's HT table and write Brunnian matches:

```bash
python3 determine_brunnian_borromean.py --screen-ht --property brunnian --method simplify --only-matches --output Run_results/ht_screening_brunnian_simplify.csv
```

Screen SnapPy's HT table and write Borromean matches:

```bash
python3 determine_brunnian_borromean.py --screen-ht --property borromean --method simplify --only-matches --output Run_results/ht_screening_borromean_simplify.csv
```

Screen both properties in one pass:

```bash
python3 determine_brunnian_borromean.py --screen-ht --property both --method simplify --only-matches --output Run_results/ht_screening_brunnian_borromean_simplify.csv
```

Run the same combined screen with the `nR` method for comparison:

```bash
python3 determine_brunnian_borromean.py --screen-ht --property both --method nr --only-matches --output Run_results/ht_screening_brunnian_borromean_nr.csv
```

Screening CSVs include a `mark` column. Matching stars (`*`, `**`, etc.)
indicate quick candidate duplicate groups with the same component count,
crossing number, and rounded volume; use `determine_duplicate_links.py` to
verify those candidates by exterior isomorphism.

Check duplicate-link candidate groups:

```bash
python3 determine_duplicate_links.py --input-file Examples/duplicate_candidate_groups.txt --output Run_results/duplicate_candidate_groups.csv
```
