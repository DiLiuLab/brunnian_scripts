# tighten_lib

The geometric and diagnostic tools for tightening links. `tighten_link_xyz.py`
stays in the repository root: it is the driver everything else supports, and it
is what a user runs first.

Each script here imports the driver with

```python
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from tighten_link_xyz import read_xyz, write_xyz, write_vect
```

resolved from `__file__`, so the repository can be checked out anywhere. They
previously carried an absolute path to one machine, which is why they only ever
ran from the scratch tree they were written in.

Run each with `--help`; `blocks2k.py`, `extract_best.py` and
`symmetrize_mirror.py` print their docstring instead, as they take positional
arguments only.

## Moving a trapped configuration

| script | what it does |
|---|---|
| `compress_gap.py` | the axial gap squeeze. A monotone coordinate map, therefore a homeomorphism, therefore **provably** cannot change the link — the only geometric move here with that guarantee. `--blend` ramps the derivative so the join is not a kink |
| `contract_clusters.py` | the general case: N clusters at the corners of any shape, moved rigidly toward their centroid with the strands rebuilt by arclength interpolation between band ends. Band assignment uses **complete** linkage — single linkage chains through mid-edge contacts and collapses everything into one cluster. **Not monotone, so HOMFLY-check every round** |
| `collapse_triangle.py` | the single-purpose precursor of the above, kept for the record |
| `perturb_component.py` | displaces one component, for testing whether a configuration is a genuine minimum or a saddle |

## Symmetry

| script | what it does |
|---|---|
| `symmetrize_mirror.py` | makes a link exactly mirror-symmetric by matching components explicitly and searching the cyclic offset. Necessary because `plc_build_symmetry` has no distance tolerance: it maps each vertex with the bare matrix, takes the NEAREST vertex, and fails only on a collision — so a curve well away from symmetric yields a perfectly injective but entirely wrong map, reporting "Error before symmetrizing 0, after 0" while it does so |
| `detect_symmetry.py` | scans mirror / C2 / C3 / C4 over a Fibonacci sphere, reporting deviation in **edge lengths**, the unit `plc_build_symmetry` actually cares about |

See also `symmetrize_link_xyz.py` in the root, which projects onto Cs / Ci / Cn / Cnv and writes the canonical frame the patched ridgerunner's `--SymmetryAxis` and `--SymmetryRef` defaults expect.

## Resolution

| script | what it does |
|---|---|
| `refine_xyz.py` | changes v/u. `--fix-minrad` is what makes refinement viable: MinRad = min(\|e_prev\|,\|e_next\|)/(2·tan(θ/2)) scales with edge length, so refining k-fold at a fixed turning angle divides minrad by k. Plain subdivision took a 113.59 structure to 190.98. The violation is local — 9 of 721 vertices — so rounding only those corners preserves minrad. `--vertices` hits an exact total, which matters when normalising several structures to a common count |

## Diagnostics

| script | what it does |
|---|---|
| `strut_free.py` | measures the longest strut-free arc, which decides `--Timewarp`. The flag rescales the gradient on free sections, and the man page adds "there is a performance cost in turning this on if there are no such segments" — worth it below roughly 60% contact, pure overhead above 90% |
| `extract_best.py` | pulls out **both** the lowest-ropelength and the most nearly critical snapshot, because the two diverge and the dominance rule cannot choose between them. Scores residual by windowed median, never by single step — a single-step residual minimum is a spike from a strut rebuild |
| `track_best.py` | reports how much a run's true optimum cost to snapshot granularity, so a missed optimum is a known quantity rather than a silent one |
| `blocks2k.py` | per-N-step ropelength record across a set of arms |
| `plot_rr_progress.py` | ropelength-versus-step plots from run logs |

## Two cautions that cost real time to learn

`strutcount.dat` has **three** columns — step, struts, minrad-locations — and the
second is the strut count.

A trace value and a file measurement **disagree**, and the file is
authoritative: a trace reports one step, and the step a file was written from is
not necessarily the last logged one.
