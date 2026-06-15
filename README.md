# Brunnian and Borromean Link Screening

This repository contains cleaned SnapPy scripts for determining Brunnian and
Borromean links and for checking candidate duplicate links.

## Get The Code

Clone the repository the first time you want a local copy:

```bash
git clone https://github.com/DiLiuLab/brunnian_scripts.git
cd brunnian_scripts
```

`git clone` downloads the repository and creates a local `brunnian_scripts`
folder. The `cd` command moves your terminal into that folder so you can run
the scripts.

If you already cloned the repository, update your local copy with:

```bash
git pull
```

`git pull` fetches the newest commits from GitHub and merges them into your
current local branch. Run it from inside the `brunnian_scripts` folder before
using the scripts if you want the latest README, examples, and results.

## Properties

A link is **Brunnian** if the full link is nontrivial, but removing any one
component gives an unlink. Equivalently, every proper sublink is trivial. The
classical example is the Borromean rings.

This repository uses the common convention that Brunnian links have at least 3
components. Some formal knot-theory definitions allow the n=2 case; under that
broader convention, a 2-component link is Brunnian if the full link is
nontrivial and each individual component is an unknot. For example, the Hopf
link can be Brunnian under the broad n=2-allowed convention, but not under the
n>=3 convention used here.

For the n>=3 convention, the code checks that the full link is nontrivial and
that every delete-one-component sublink is an unlink. Checking only delete-one
sublinks is sufficient: any smaller proper sublink is contained in one of those
delete-one sublinks, and a sublink of an unlink is also an unlink.

A link is treated as **Borromean** here if the full link is nontrivial and every
2-component sublink is an unlink. This includes the classical Borromean rings
and extends the same pairwise-unlinked condition to links with more than three
components. For n>3, this generalized Borromean condition is weaker than the
Brunnian condition, because a 3-component or larger proper sublink may still be
nontrivial even when every pair of components is unlinked.

The n=2 Brunnian edge case is not the same thing as being a 2-component prime
link. Under the broad n=2 convention, Brunnian means "nontrivial link with both
components unknotted." Prime means the link cannot be decomposed as a
nontrivial connected sum. These are different properties; for instance, Knot
Atlas describes `L10a108` as two interlinked trefoil knots, so its components
are not unknots.

## HT Table

The HT table is SnapPy's
[`HTLinkExteriors`](https://www.math.uic.edu/t3m/SnapPy/censuses.html#snappy.HTLinkExteriors)
census from the Hoste-Thistlethwaite link tables. It contains link exteriors
for links up to 14 crossings, together with data accessible through SnapPy such
as the link name, number of cusps/components, volume, triangulation
information, DT code, and solution type. The screening mode in this repository
iterates through this table by crossing number and filters to links with at
least 3 components.

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

Duplicate-link check from command-line link strings:

```bash
python3 determine_duplicate_links.py --links L14n49573 L14n50824
```

Duplicate-link check from command-line DT-code strings. Quote each DT code so
the shell passes it as one argument:

```bash
python3 determine_duplicate_links.py --links "DT: [(-8,-12,16),(-24,-22,-28,-26),(-10,-14,-2),(-20,-6,-18,-4)]" "DT: [(-6,16),(26,-2,34,32,-40,-30,4,38,-28,-36),(-8,-22,18,12),(10,-24,20,-14)]"
```

Duplicate-link check from an input file. In the file, use one link string per
line and blank lines to separate independent candidate groups:

```bash
python3 determine_duplicate_links.py --input-file Examples/duplicate_candidate_groups.txt
```

Duplicate-link check with CSV output:

```bash
python3 determine_duplicate_links.py --input-file Examples/duplicate_candidate_groups.txt --output Run_results/duplicate_candidate_groups.csv
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

## References

- H. N. Howards, [Brunnian Spheres](https://users.wfu.edu/howards/papers/brunnianspheres.pdf).
- B. S. Mangum and T. Stanford, [Brunnian links are determined by their complements](https://eudml.org/doc/121302), Algebraic & Geometric Topology 1, 143-152, 2001.
- C. Liang and K. Mislow, [On Borromean links](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/liangmislow.pdf), Journal of Mathematical Chemistry 16, 27-35, 1994.
- Knot Atlas, [L10a108 Quick Notes](https://katlas.org/wiki/L10a108_Quick_Notes).
