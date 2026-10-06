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

## Polishing a plateau

| script | what it does |
|---|---|
| `slp_tighten.py` | run with no arguments (or `--gui`) for a parameter window with file pickers, *Run now* for `--group detect`/`--dry-run`/`--diagnose`, and *Run in background* logging to `<output>.log`. The **equivariant sequential LP**: per-vertex moves in the exact symmetry subspace, linearised contact, minRad and edge-length constraints, a trust region and an exact projection; a step is accepted only on an octrope decrease plus a PL isotopy certificate that the straight-line step stays embedded. For a symmetric run parked well above residual 0.1 by the forced Animation stepper (8BL RD2 legs 0.31-0.43, 7BL C7 0.322, lark C5 0.263; not every symmetric run parks: 4BL_wider under C2v stopped on residual 0.0146). On 8BL D2d the hand-run slp_B chain took 298.242 -> 291.748 (residual 0.307 -> 0.010). This tool's recipe `--group D2d --passes 2 --resample-vu 6` reaches 291.745 in 618 s (residual 0.005), and one pass, the default, took 7BL (C7) 274.365 -> 272.816, lark (C5) 210.041 -> 202.282 and 4_LC (Cs) 113.590 -> 113.049. The designer's prototype it was built from had reached 273.910, 204.141 and 113.036 on those three. `--diagnose` reports the first-order feasible fraction (13.7% at 298.242; ~0 means critical), and `--group detect` which of its 38 candidate groups the file holds exactly, in the canonical frame (axis z, reference x) unless `--axis`/`--ref` are given. Edge bounds 0.8/1.25 are on by default because without them an edge collapsed with HOMFLY still SAME. **Gate its output** (`--gate-ref`, or `tighten_cycle.py --slp`), and note the gate is on the link, not the diagram: the diagram may change, so check it with a diagram-isomorphism tool before quoting a per-diagram minimum |

## Symmetry

| script | what it does |
|---|---|
| `symmetrize_mirror.py` | makes a link exactly mirror-symmetric by matching components explicitly and searching the cyclic offset. Necessary because `plc_build_symmetry` has no distance tolerance: it maps each vertex with the bare matrix, takes the NEAREST vertex, and fails only on a collision — so a curve well away from symmetric yields a perfectly injective but entirely wrong map, reporting "Error before symmetrizing 0, after 0" while it does so |
| `detect_symmetry.py` | scans mirror / C2 / C3 / C4 over a Fibonacci sphere, reporting deviation in **edge lengths**, the unit `plc_build_symmetry` actually cares about |
| `symmetry_maps.py` | exact vertex maps and the Reynolds operator for the point groups C1, Cs, Ci, Cn, Cnv, Cnh, RDn, Dnd (D2d) and S2n, built from named generators (Dnh, D2h and the polyhedral groups are not built). **Refuses rather than guesses**: tolerance 1e-6 mean edges, per-component maps must be shifts or reversals, and the maps must satisfy every group relation, a check `symmetrize_link_xyz.py` does not make on its live path. A bare `D2` is refused, because ridgerunner's D2 is a single mirror (write `Cs`, or `RD2` for the order-4 dihedral group). Also an equivariant resample. A library with a diagnostic `__main__`: `symmetry_maps.py FILE detect` lists which of its 38 candidate groups (C1, Cs, Ci, C2-C12, C2v-C8v, C2h-C6h, RD2-RD6, D2d in both frames, D3d, D4d, S4, S6, S8) hold exactly, in the canonical frame (axis z, reference x) unless `--axis`/`--ref` are given |

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
| `radial_squeeze.py` | the **radial hole squeeze**: closes the void *enclosed* by a ring rather than a gap between clusters. `Φ(r,θ,z) = (g(r),θ,z)` with `g` strictly increasing and `g(0)=0` is a homeomorphism of R³ isotopic to the identity, so like `compress_gap.py` it **provably cannot change the link** — and being θ-independent it preserves `Cn` and the z-mirror exactly. `apply_squeeze_about()` does the same about an arbitrary line. Its damage is azimuthal *proximity*, the kind RidgeRunner repairs, not curvature |
| `find_axis.py` | finds the squeeze axis when there is no symmetry to hand you one. Sweeps directions on a Fibonacci sphere, solves the largest-empty-circle problem on each projection (Voronoi vertices inside the Delaunay hull), then ranks the widest tunnels by **the length a probe squeeze actually removes** — clearance alone is not leverage, since an axis the curve merely parallels shortens nothing. **Not for symmetric links.** It does not reliably recover the symmetry axis: on an exactly `C7` structure it returned an axis 6.4 degrees off `z`, and squeezing about that broke the symmetry and aborted the round. Use the symmetry axis whenever there is one |
| `measure_link.py` | one command for "what is this file and what should I do with it": every quantity with a plain-language verdict, plus a `WHAT TO DO NEXT` block. Written for newcomers; also the fastest way to check a result. A file with minRad = minStrut/2 = tau (each within 1e-4) and more than 200 struts now reads "contact AND curvature both active" at any residual. Below residual 0.1 it is called converged and not to be re-descended. At 0.1 or above it reads NOT converged, its residual line says "room left for another SLP pass", and it is told to run another `slp_tighten.py` pass rather than a RidgeRunner descent, which pays the restart cost below; when the residual cannot be measured it says so. Before 2026-10-05 every such file was told to re-descend, the hand-run 291.748 included. A contact-limited SLP output gets the old wording: lark's 202.282 (minRad/tau 1.034, residual 0.136) still reads "there is room to descend", and the restart cost applies to it all the same |
| `quicklook.py` | renders any file to `.png` (four views), `.glb` (a tube at the file's own thickness, so a contact really looks like tube touching tube) and a `.xyz` copy, all under one ropelength-named stem |
| `blocks2k.py` | per-N-step ropelength record across a set of arms |
| `plot_rr_progress.py` | ropelength-versus-step plots from run logs |

## A note on `contract_clusters.py` and folds

A free run between two clusters can be asked to move sideways **further than it
is long**, and linear interpolation then folds the strand back on itself.
Measured at `f = 0.90` on one link: ten of twenty runs were three vertices
(0.298 D) joining clusters whose translations differed by 0.32–0.42 D, giving a
required displacement gradient of 1.09–1.40. Turning angles reached 161.7° and
minRad fell from 0.533 to 0.0017, taking ropelength from 215 to 63 557 — while
HOMFLY returned **SAME**, so no topology check catches it.

`--max-gradient` (default 0.5) now widens the ramp into the adjacent bands until
the gradient is safe, with a smoothstep profile. `contract_clusters.py.bak` is
the pre-fix version, kept so the failure stays reproducible. Note the guard costs
a little of the rigidity the method relied on — up to 0.16 D of intra-cluster
distortion when it fires, against exactly zero before — and it does not address
the *other* way minRad dies, which is edge compression with no angular change at
all.

## A note on what `slp_tighten.py` hands back

Its result is **link-gated, not diagram-checked**. HOMFLY SAME says the link
survived; it cannot say the diagram did, and only the hand-run 291.748 has had a
diagram check. Check the diagram with a diagram-isomorphism tool before quoting
any other output as a per-diagram minimum: lark's 202.282, for one, sits 0.12
above 5BL_square's 202.158.

The output is often **curvature-active**: minRad = minStrut/2 = tau. On 8BL the
hand-run 291.748 measures 2416 struts at residual 0.010, and this tool's 291.745
measures 2486 at 0.005; 7BL (C7, 272.816) and 4_LC (Cs, 113.049, whose input
already was) end the same way. Not always: lark (C5, 202.282) ends
contact-limited at minRad/tau 1.034, and 7BL after the cycle's 2-minute stage at
1.108. Do not re-descend either kind. The reason is the restart cost, not the
curvature: on 8BL every ridgerunner restart cost +0.36 to +0.38 at once, and a
continuation then needed 13.7k-16.3k steps to get back to its input value. A
file still above residual 0.1 (7BL 0.147, lark 0.136, 4_LC 0.141) wants another
SLP pass, not a descent.

Its edges are uneven by design. The bounds allow 0.8-1.25 of each component's
mean edge, and the 8BL recipe ends at edge ratio 1.561 with 24 edges on the
bounds. Resample to even edges before any ridgerunner leg that follows; the
8BL campaign redistributed whenever a component's edge ratio passed 1.15.

Its number depends on resolution, so quote the range with it.
`--robustness-vu 4.5,6,8` re-reads the output resampled to each v/u: the
recipe's 291.745 reads 292.230 / 291.957 / 292.185. Compare each reading with
the **input resampled the same way**, never with the output's own ropelength:
the equivariant spline resample alone costs some inputs a lot (lark +18.0 at
v/u 8, 4_LC Cs +8.7 and +10.6 at v/u 6 and 8), so a raw "+5" can be the
resample's doing, not fragility. Input -> output, each resampled to v/u 4.5 / 6 / 8:

| link | v/u 4.5 | v/u 6 | v/u 8 |
|---|--:|--:|--:|
| 8BL D2d, 298.242 -> 291.745 | 298.544 -> 292.230 | 298.482 -> 291.957 | 298.531 -> 292.185 |
| 7BL C7, 274.365 -> 272.816 | 274.715 -> 273.246 | 274.576 -> 273.006 | 274.659 -> 273.075 |
| lark C5, 210.041 -> 202.282 | 210.249 -> 202.590 | 210.196 -> 202.429 | 228.043 -> 202.956 |
| 4_LC Cs, 113.590 -> 113.049 | 113.779 -> 113.249 | 122.325 -> 117.958 | 124.202 -> 118.676 |
| 2_BC Ci, 116.569 -> 115.997 | 116.539 -> 116.229 | 116.403 -> 116.607 | 117.163 -> 117.744 |

Every gain survives except 2_BC's (an earlier build of the tool), which reads
worse than its input at v/u 6 and 8. The hand-run SLP results read
292.06-292.55 over the same range with an equal-edge resample, against
298.48-298.59 for the 298.242 input.

Gate the collapse guard on `min_edge_frac` in the `SLP RESULT` line (>= 1 while
the bounds are on), not on `min_edge_frac_in`. The latter is measured against
the tool's input, and after the recipe's v/u 6 resample it reads 0.7431 by
design: the resample shortens the mean edge by about 1350/1748 = 0.77.

## Two cautions that cost real time to learn

`strutcount.dat` has **three** columns — step, struts, minrad-locations — and the
second is the strut count.

A trace value and a file measurement **disagree**, and the file is
authoritative: a trace reports one step, and the step a file was written from is
not necessarily the last logged one.
