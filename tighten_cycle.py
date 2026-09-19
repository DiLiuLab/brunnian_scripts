#!/usr/bin/env python3
"""Automate contract-then-recondition cycles for tightening a link.

Straightforward RidgeRunner descent takes a link down to a trapped
configuration and stops. On a structure whose slack is *space between tight
clusters* rather than length in a bridge, the remaining progress comes from a
geometric move -- cluster contraction -- and not from more steps. Measured on a
5-component C5 link over five rounds:

    316.83  raw diagram
    251.05  descent, 8956 steps         -65.78   <- descent does the bulk once
    238.69  contraction f=0.87          -12.35
    237.00  descent, 5927 steps          -1.69
    230.46  contraction f=0.92           -6.54
    228.43  descent, 30000 steps         -2.04
    222.42  contraction f=0.92           -6.01
    222.42  descent, 15038 steps         +0.01   <- bought no length at all
    219.66  contraction f=0.96           -2.77

After the first descent, essentially all of the ropelength came from the
contractions. That does NOT make the descents optional, and this is the point
the script is built around: the descent that bought +0.01 took minRad/tau from
1.072 back to 1.341 and rebuilt the contact set from 785 struts to 1260. The
NEXT contraction was only feasible because of that headroom. So the cycle is

    contract  ->  (re-symmetrize)  ->  HOMFLY gate  ->  descend  ->  measure

where the descent is a RECONDITIONING step, not a polish, and the loop carries
the *reconditioned* file forward while tracking the best ropelength separately.
Feeding the raw contraction into the next contraction stacks two unrelaxed
geometric moves and stalls.

Four rules learned the hard way, all encoded here:

1.  FEASIBILITY IS minRad > minStrut/2, not minRad >= the input's minRad.
    `contract_clusters.py --scan` calls a factor "damaged" whenever minRad
    falls below the input's, which is far more conservative than the actual
    constraint. Thickness is min(minStrut/2, minRad); a contraction is fine as
    long as it stays *strut*-limited. Scanning on the real criterion beat the
    tool's own pick by 2.8 units in one round and 2.2 in another. Below the
    boundary, thickness collapses to minRad and ropelength jumps (one measured
    case: f=0.955 -> minRad 0.476 -> ropelength 230.2 instead of 219.7).

2.  MEASURE STRUTS AND RESIDUAL AT THE NOMINAL RADIUS 0.5 -- BUT ONLY ON A
    NORMALISED FILE. RidgeRunner normalises to a nominal 0.5 and counts contacts
    within its overstep tolerance of *that*; the measured thickness drifts just
    below it, and both binaries want struts strictly shorter than 2r. Key them
    to the measured tau and a healthy 1200-strut configuration reads as
    "residual 0.703, 210 struts" -- indistinguishable from a degenerate run.
    The converse trap is worse: a file that is NOT normalised (a freshly
    contracted one, tau 0.4926) accepts the nominal radius happily and returns
    "residual 0.0000" from an over-determined system, which reads as perfectly
    critical. Residual is therefore reported as UNAVAILABLE off-nominal, and
    --auto-flags treats unavailable as far-from-critical.

3.  ROTATE OUT OF THE CANONICAL FRAME BEFORE ANY HOMFLY CHECK. The symmetrizer
    writes the rotation axis on z, knottype's projection looks down it, and the
    structure doubles onto its own image. Reference and candidate both go
    through one fixed generic rotation here so the comparison is like-for-like.

4.  CHECK TOPOLOGY EVERY ROUND. Cluster contraction moves different pieces in
    different directions and is a monotone map of nothing, so unlike the axial
    squeeze it carries no topological guarantee. A failed gate aborts the round
    and keeps the previous configuration.

Usage:

    python3 tighten_cycle.py input.xyz --work-dir /path/outside/dropbox \
        --rounds 6 --steps 30000 --symmetry Z/5Z --sym-group Cn --sym-order 5

    python3 tighten_cycle.py input.xyz --work-dir /tmp/run --dry-run

    # continue a cycle whose driver died (crash, reboot, kill). The per-round
    # .xyz files are the durable record; the round index, the running best and
    # the ledger live only in the driver's memory until the cycle ends.
    python3 tighten_cycle.py --resume --work-dir /path/outside/dropbox \
        --rounds 8 --symmetry Z/5Z --sym-group Cn --sym-order 5
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict, field, fields
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent
LIB = REPO / "tighten_lib"
sys.path.insert(0, str(REPO))
from tighten_link_xyz import read_xyz, write_xyz, write_vect  # noqa: E402

# One fixed rotation, applied to the reference and to every candidate, so no
# symmetry axis of a Cn/Cs/Ci structure lands on z. See rule 3.
_A, _B, _C = 0.7234561, 0.4312789, 1.1123457


def _rot(a: float, axis: int) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    if axis == 0:
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == 1:
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


GENERIC = _rot(_A, 2) @ _rot(_B, 0) @ _rot(_C, 1)


class CycleError(RuntimeError):
    pass


# --------------------------------------------------------------------------- #
# measurement
# --------------------------------------------------------------------------- #

@dataclass
class Measurement:
    path: str
    rop: float
    tau: float
    minrad: float
    minstrut: float
    length: float
    vertices: int
    struts: int | None
    residual: float | None
    mean_dD: float
    max_dD: float
    vu: float
    minrad_over_tau: float
    strut_radius: float

    @property
    def strut_limited(self) -> bool:
        """True when thickness is set by contact, not by curvature (rule 1)."""
        return self.minrad > self.minstrut / 2


def _run(cmd: list[str], cwd: Path | None = None) -> str:
    p = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, cwd=cwd)
    return p.stdout + p.stderr


def _grab(text: str, key: str) -> float | None:
    # A word boundary so "Diameter Ropelength" does not match "Ropelength".
    m = re.search(rf"(?<![A-Za-z] ){key}:\s*([0-9.eE+-]+)", text)
    return float(m.group(1)) if m else None


NOMINAL_TAU = 0.5
# RidgeRunner's printed overstep tolerance is 1e-4 of thickness. A file that sits
# further than this from the nominal 0.5 is NOT a normalised ridgerunner output,
# and keying the strut radius to the nominal over-counts its contact set.
NOMINAL_WINDOW = 1e-3


def is_normalised(tau: float) -> bool:
    return abs(tau - NOMINAL_TAU) <= NOMINAL_WINDOW


def strut_radius_for_minstrut(minstrut: float) -> float:
    """Radius that reliably includes the binding strut: (minStrut/2)*(1+1e-4).

    THIS IS THE ROBUST RULE. Both `struts` and `residual` look for struts
    STRICTLY shorter than 2r, so any radius keyed to a nominal value fails
    whenever minStrut lands exactly on 2r. Measured three times today on the
    same underlying mistake:

      tau 0.499953, minStrut 0.999906 -> r=tau gave 0 struts / residual 1
      tau 0.492554, minStrut 0.985100 -> r=0.5  gave 1484 struts / residual 0
      tau 0.500000, minStrut 1.000000 -> r=0.5  gave 1 strut  / residual 1

    That last reading looks exactly like a collapsed run on a configuration
    carrying ~1000 contacts. Keying to the file's own minStrut removes the
    special-casing: it reproduces 0.5000 on the files where the nominal happened
    to work, and 0.50005 where it did not.

    Caveat worth knowing: on a configuration with a dense band of near-contacts
    just above the minimum, the count stays radius-sensitive (1038 at +1e-4,
    2208 at +1e-3 on one file). Quote a range, and treat strutcount.dat column 2
    as the authority for a run.
    """
    return (minstrut / 2) * 1.0001


def strut_radius(tau: float, minstrut: float) -> float:
    """The radius to hand `struts` and `residual`.

    Two requirements, and only this rule meets both:

      NON-EMPTY  -- must exceed minStrut/2, since both binaries want struts
                    STRICTLY shorter than 2r. A radius sitting exactly on
                    minStrut/2 returns an empty constraint set, which reads as a
                    collapsed run.
      COMPARABLE -- must be the SAME threshold across a lineage, or strut counts
                    cannot be compared between rounds. Keying purely to each
                    file's own minStrut moves the threshold per file: one
                    configuration dropped from 1160 struts to 695 purely because
                    its minStrut sat a fraction below the nominal.

    So: a normalised file gets the nominal 0.5 nudged up by the overstep
    tolerance (0.50005) -- fixed for every such file, and always above
    minStrut/2 in practice. An off-nominal file (raw layout, freshly contracted)
    has no shared threshold to use, so it falls back to its own minStrut.
    """
    if is_normalised(tau):
        return max(NOMINAL_TAU, minstrut / 2) * (1 + 1e-4)
    return (minstrut / 2) * (1 + 1e-4)


def strut_radius_for(tau: float) -> float:
    """Rule 2: nominal 0.5 for a ridgerunner-normalised file, else the file's own.

    The window used to be `0.49 < tau < 0.51`, which was far too loose. A freshly
    contracted file measured tau = 0.492554 with minStrut 0.9851; that passed the
    old window, so the radius became 0.5, counting every contact out to 1.0 D
    (1484 of them, left over from the PRE-contraction spacing) and handing the
    residual solver an over-determined system that returned exactly 0.0. A false
    "residual 0.0000" then told --auto-flags the configuration was critical and
    selected --no-eq for a file that had just had strain injected into it.
    """
    return NOMINAL_TAU if is_normalised(tau) else tau * 1.0001


def measure(path: Path, want_struts: bool = True) -> Measurement:
    comps = [np.asarray(c, float) for c in read_xyz(path)]
    vect = path.with_suffix(".cyc.vect")
    write_vect(vect, [c.tolist() for c in comps])
    try:
        out = _run(["ropelength", str(vect)])
        rop = _grab(out, "Ropelength")
        tau = _grab(out, "Thickness")
        minrad = _grab(out, "minRad")
        minstrut = _grab(out, "minStrut")
        if None in (rop, tau, minrad, minstrut):
            raise CycleError(f"could not parse octrope output for {path.name}:\n{out}")

        r = strut_radius(tau, minstrut)
        struts = residual = None
        if want_struts:
            st = _run(["struts", "-s", "-n", "-r", f"{r:.9f}", str(vect)])
            m = re.search(r"Sorting (\d+) struts", st) or re.search(r"(\d+) struts", st)
            struts = int(m.group(1)) if m else None
            # Only a normalised file has a meaningful ridgerunner residual. On a
            # freshly contracted file the contact set straddles two different
            # spacings (the new closest approach and the pre-contraction 1.0 D),
            # so every radius is wrong: the file's own tau finds 7 struts at
            # residual 0.98, the nominal finds 1484 at residual 0.0. Report it as
            # unavailable rather than pick one and be believed.
            if is_normalised(tau):
                residual = _grab(_run(["residual", "-r", f"{r:.9f}", str(vect)]),
                                 "Residual")

        edges = np.concatenate(
            [np.linalg.norm(np.diff(np.vstack([c, c[:1]]), axis=0), axis=1) for c in comps]
        )
        L = float(edges.sum())
        N = int(sum(len(c) for c in comps))
        D = 2 * tau
        return Measurement(
            path=str(path), rop=rop, tau=tau, minrad=minrad, minstrut=minstrut,
            length=L, vertices=N, struts=struts, residual=residual,
            mean_dD=float(edges.mean() / D), max_dD=float(edges.max() / D),
            vu=N / rop, minrad_over_tau=minrad / tau, strut_radius=r,
        )
    finally:
        vect.unlink(missing_ok=True)


def to_generic_frame(src: Path, dst: Path) -> Path:
    comps = [np.asarray(c, float) for c in read_xyz(src)]
    write_xyz(dst, [(c @ GENERIC.T).tolist() for c in comps], 12)
    return dst


def symmetry_deviation(path: Path, order: int) -> float | None:
    """Deviation of the link from exact Cn about z, in units of a HALF MEAN EDGE.

    Returns (median, max, fraction_over_half_edge), or None only when the
    component count is not `order`.

    PARAMETRISATION-FREE, and that matters. An earlier version compared vertex
    lists under every component match, cyclic offset and orientation, and gave up
    (returning None) whenever the components had different vertex counts. On a
    7-component layout drawn with C7 intent but sampled at 83-86 vertices per
    component, that None was read as "no symmetry" -- when the structure's
    component centroids in fact form a regular heptagon to 1.13% in radius and
    1.96 degrees in angle. A declined measurement is not a negative result.

    The measure here rotates the link and asks, for each vertex, how far it is
    from the curve (not from another vertex): edges are densified first, because
    the vertex spacing is far larger than any symmetry error worth detecting.

    The unit is the half mean edge because that is the threshold that decides
    whether `plc_build_symmetry` will work: it has NO distance tolerance, takes
    the nearest vertex to each image, and silently bakes in a wrong
    correspondence when the curve is further off than about half an edge.
    Measured reference points: an exactly-symmetric layout came in at 0.003
    half-edges and symmetrizing it was free; a raw diagram at 0.61 median /
    4.94 max half-edges with 33.6% of vertices over the threshold is the regime
    where symmetrizing the raw 4_LC diagram produced the 4-component unlink.
    """
    from scipy.spatial import cKDTree

    comps = [np.asarray(c, float) for c in read_xyz(path)]
    if len(comps) != order:
        return None
    origin = np.array([c.mean(axis=0) for c in comps]).mean(axis=0)
    verts = np.vstack(comps) - origin
    edges = np.concatenate(
        [np.linalg.norm(np.diff(np.vstack([c, c[:1]]), axis=0), axis=1) for c in comps])
    half_edge = float(edges.mean()) / 2

    step = half_edge / 10
    dense = []
    for c in comps:
        loop = np.vstack([c, c[:1]]) - origin
        for a, b in zip(loop[:-1], loop[1:]):
            n = max(2, int(np.linalg.norm(b - a) / step))
            dense.append(a + np.outer(np.linspace(0, 1, n, endpoint=False), b - a))
    tree = cKDTree(np.vstack(dense))

    th = 2 * np.pi / order
    c, s = np.cos(th), np.sin(th)
    d, _ = tree.query(verts @ np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]).T)
    return (float(np.median(d)) / half_edge, float(d.max()) / half_edge,
            float((d > half_edge).mean()))


def symmetrize_is_safe(dev: tuple[float, float, float] | None,
                       max_dev: float = 0.8) -> tuple[bool, str]:
    """Decide whether to symmetrize now or descend first.

    The standing rule is "tighten first, symmetrize after", but the rule is
    really about DISTANCE FROM THE SYMMETRY, not about raw-versus-tightened:
    measure the deviation and decide.

    GATE ON THE MAX, NOT THE MEDIAN. plc_build_symmetry maps each vertex with
    the bare matrix and takes the NEAREST vertex to the image; it fails only when
    two vertices claim the same target. The failure is therefore per-vertex and
    local: one region further off than about half an edge is enough to bake in a
    wrong correspondence, however good the median is. Conversely a curve whose
    WORST vertex is inside half an edge cannot pick a wrong target anywhere, even
    if its median is not tiny.

    Measured calibration, in half-mean-edge units:
        5BL_square raw      max 0.00,  0.0% over  -> safe, and it was free
        7BL_square_fixed    max 0.43,  0.0% over  -> safe (median 0.16 is
                                                      irrelevant; nothing is
                                                      near the threshold)
        raw 4_LC (guide)    max 11.3, 30.5% over  -> symmetrizing gave the UNLINK

    Default 0.8 rather than 1.0 keeps a margin below the half-edge threshold.
    """
    if dev is None:
        return False, "component count does not match the requested order"
    med, mx, frac = dev
    if mx <= max_dev and frac <= 0.0:
        return True, (f"worst vertex {mx:.3f} half-edges from its image and no "
                      f"vertex over the threshold (median {med:.3f}): the "
                      f"nearest-vertex map cannot go wrong; safe to symmetrize now")
    return False, (f"worst vertex {mx:.2f} half-edges from its image, "
                   f"{frac * 100:.1f}% of vertices over half an edge (median "
                   f"{med:.3f}): plc_build_symmetry's nearest-vertex map can pick "
                   f"the wrong target. Descend first, then symmetrize.")


# --------------------------------------------------------------------------- #
# external steps
# --------------------------------------------------------------------------- #

def homfly_reference(ref_xyz: Path, scratch: Path) -> str:
    gen = to_generic_frame(ref_xyz, scratch / "reference_generic.xyz")
    out = _run([sys.executable, REPO / "verify_topology.py", gen, gen])
    m = re.search(r"HOMFLY (.+)", out)
    if not m:
        raise CycleError(f"could not compute a reference polynomial:\n{out}")
    return m.group(1).strip()


def homfly_verdict(ref_generic: Path, candidate: Path, scratch: Path) -> str:
    """Returns 'SAME', 'DIFFERENT' or 'UNKNOWN' (uncomputable, not a change)."""
    gen = to_generic_frame(candidate, scratch / (candidate.stem + "_generic.xyz"))
    out = _run([sys.executable, REPO / "verify_topology.py", ref_generic, gen])
    for verdict in ("DIFFERENT", "UNKNOWN", "SAME"):
        if re.search(rf"^\s*{verdict}\b", out, re.M):
            return verdict
    return "UNKNOWN"


def contract(src: Path, dst: Path, factor: float, join: float | None = None) -> None:
    cmd = [sys.executable, LIB / "contract_clusters.py", src, dst, "--factor", factor]
    if join is not None:
        cmd += ["--join", join]
    _run(cmd)
    if not dst.exists():
        raise CycleError(f"contract_clusters.py wrote nothing for factor {factor}")


def symmetrize(src: Path, dst: Path, group: str, order: int, axis: str,
               no_permute: bool = False) -> str:
    cmd = [sys.executable, REPO / "symmetrize_link_xyz.py", src, "-o", dst,
           "--group", group, "-p", order, "--axis", axis]
    if no_permute:
        cmd.append("--no-permute")
    out = _run(cmd)
    if not dst.exists():
        raise CycleError(f"symmetrize_link_xyz.py failed:\n{out}")
    m = re.search(r"Largest vertex move:\s*([0-9.eE+-]+)", out)
    return m.group(1) if m else "?"


def _build_descend_cmd(src: Path, dst: Path, work_dir: Path, steps: int,
                       symmetry: str | None, extra: list[str], snapshot: int,
                       stop_res: float, wrapper: list[str] | None = None) -> list[str]:
    cmd = [sys.executable, REPO / "tighten_link_xyz.py", src, "-o", dst,
           "--work-dir", work_dir, "-s", steps, "--stop-res", stop_res,
           "--snapshot-interval", snapshot, "--select", "both", "--keep-run-dir"]
    if symmetry:
        cmd += ["--symmetry", symmetry]
    cmd += list(wrapper or [])
    for e in extra:
        cmd.append(f"--rr-arg={e}")
    return [str(c) for c in cmd]


def descend(src: Path, dst: Path, work_dir: Path, steps: int, symmetry: str | None,
            extra: list[str], snapshot: int, stop_res: float,
            wrapper: list[str] | None = None) -> Path:
    out = _run(_build_descend_cmd(src, dst, work_dir, steps, symmetry, extra,
                                  snapshot, stop_res, wrapper))
    if not dst.exists():
        raise CycleError(f"tighten_link_xyz.py produced no output:\n{out[-2000:]}")
    return dst


def _block_medians(trace: np.ndarray | None, block: int) -> dict[int, float]:
    if trace is None or trace.size == 0:
        return {}
    s, v = trace[:, 0], trace[:, 1]
    out = {}
    for lo in range(0, int(s[-1]), block):
        seg = v[(s >= lo) & (s < lo + block)]
        if len(seg):
            out[lo + block] = float(np.median(seg))
    return out


def _ridgerunner_pid(vect_stem: str) -> int | None:
    """The wrapper's ridgerunner child, found by the .vect it was handed.

    We signal the CHILD, not the wrapper: tighten_link_xyz.py catches a
    signal-killed child as RunCancelled and still runs its dominance selection,
    writing the best saved configuration. Killing the wrapper would skip that.
    """
    out = _run(["pgrep", "-f", f"ridgerunner.*{vect_stem}.vect"])
    pids = [int(t) for t in out.split() if t.isdigit()]
    return pids[-1] if pids else None


ALL_PLATEAU = frozenset({"minrad", "struts", "ropelength"})


def resolve_plateau_on(setting, move: str | None, factor: float | None,
                       hard_f: float) -> tuple[frozenset[str], str]:
    """Decide which signals gate this round's stop.

    An explicit --plateau-on wins. Otherwise "auto" keys on how badly the round's
    move broke the configuration, because that decides what the descent is FOR.

    A hard squeeze (small f -- f is the slope of the radial map, so SMALLER is
    more aggressive) collapses the thickness and hands the descent a large
    repair job whose output feeds another move. There, minRad and the strut
    count are the signals that matter: they say the structure is fit to be
    squeezed again. Waiting for ropelength to converge as well spends thousands
    of steps polishing a configuration that is about to be deliberately broken.

    A gentle squeeze barely perturbs anything, so the round is mostly a descent
    and its ropelength IS the product. Contraction likewise preserves tau, and
    an initial descent has no move at all. Those all get the full gate.

    The caveat, measured and worth stating: truncation PROPAGATES. The
    7-component link's round 1 stopped early at ropelength 276.080, so round 2
    began from a less converged input than it needed to, and a later
    continuation with no move at all recovered 0.84 units from that same file.
    Cheap intermediate rounds are a deliberate trade of result quality for
    throughput, not a free saving -- run the final round, or any round you will
    report, with the full gate.
    """
    if setting != "auto":
        return setting, "as given"
    if move == "squeeze" and factor is not None and factor <= hard_f:
        return (ALL_PLATEAU - {"ropelength"},
                f"auto: squeeze f={factor:.2f} <= {hard_f} is a repair round, "
                f"ropelength not gated")
    if move == "squeeze" and factor is not None:
        return ALL_PLATEAU, (f"auto: squeeze f={factor:.2f} > {hard_f} barely "
                             f"perturbs, so ropelength is the product")
    return ALL_PLATEAU, f"auto: {move or 'no move'} preserves tau"


def descend_auto(src: Path, dst: Path, work_dir: Path, symmetry: str | None,
                 extra: list[str], snapshot: int, stop_res: float,
                 block: int, patience: int, min_steps: int, max_steps: int,
                 minrad_eps: float, strut_eps: float, rop_eps: float = 1e-4,
                 plateau_on: frozenset[str] = frozenset(
                     {"minrad", "struts", "ropelength"}),
                 poll: float = 20.0,
                 wrapper: list[str] | None = None) -> tuple[Path, str]:
    """Run one descent and stop when reconditioning AND ropelength have plateaued.

    minRad and the strut count say when the THICKNESS REPAIR is finished, which
    is what gates the next geometric move's feasibility. They were the original
    stopping signal, calibrated on a round where minRad's block median peaked at
    step 3000 and only oscillated for 12000 more while the strut count peaked at
    8000.

    Keying the stop on those two alone was wrong, and measurably so. A repaired
    thickness does not mean a converged ropelength: the 7-component link's round
    stopped at step 12000 with minRad and struts both flat for three blocks, at
    ropelength 276.080, and restarting that exact file with no geometric move at
    all -- no squeeze, no mirror, no refinement -- recovered 0.84 units within
    8500 steps and was still falling at ~0.21 per 2000 with no sign of tapering.
    The corner margin ROSE over that stretch, 1.023 to 1.164, so the descent was
    still reconditioning too; the plateau test had simply read the oscillation of
    an already-high minRad as convergence. Every round of that link stopped the
    same way, at 7000-16000 steps against a 30000 cap, so the per-round gains
    were systematically understated and rounds were judged failures on truncated
    evidence.

    So ropelength must flatten as well. It is tracked the opposite way from the
    other two -- improvement means going DOWN -- and needs a much tighter
    threshold, because ropelength changes by far less in relative terms: a run
    actively descending moves ~4e-4 of its own value per 1000-step block, while
    a genuinely converged one moves under 1e-5. The 1e-4 default sits between
    those with an order of magnitude either side.

    Which signals GATE the stop is selectable, because the right answer depends
    on what the round is for. A round whose job is to make the next geometric
    move feasible only needs the thickness repaired -- minRad and struts -- and
    waiting for ropelength to converge as well spends thousands of steps
    polishing a configuration that is about to be deliberately broken again. A
    final round, or one being measured as a result, needs all three. All three
    are always measured and printed; `plateau_on` chooses which must be flat.

    The run is launched once and killed in place, rather than chained in blocks:
    restarting would re-autoscale and rebuild the contact set every block.
    """
    if not plateau_on:
        raise CycleError("--plateau-on needs at least one signal")
    import signal
    import time

    # Send the child's stdout to a FILE, never a pipe. Piping it and not
    # draining it deadlocks: ridgerunner prints a line per step, the OS pipe
    # buffer (~64 KB, about 1300 lines) fills, the wrapper blocks on write(),
    # and because it is then not reading ridgerunner's output ridgerunner blocks
    # too. Measured before this fix: a run showed 5h31m wall for 1m33s of CPU,
    # frozen at step 2000, with ridgerunner's own walltime.dat stopped at 84s.
    # The polling loop below reads the TRACE FILES, so it never needs the pipe.
    work_dir.mkdir(parents=True, exist_ok=True)
    child_log = work_dir / "descend_stdout.log"
    with child_log.open("w") as sink:
        proc = subprocess.Popen(
            _build_descend_cmd(src, dst, work_dir, max_steps, symmetry, extra,
                               snapshot, stop_res, wrapper),
            stdout=sink, stderr=subprocess.STDOUT, text=True)

    best_mr = best_st = -np.inf
    best_rop = np.inf                 # ropelength improves DOWNWARD
    stale_mr = stale_st = stale_rop = 0
    seen: set[int] = set()
    reason = f"reached --max-steps {max_steps}"
    try:
        while proc.poll() is None:
            time.sleep(poll)
            mr = _block_medians(read_trace(work_dir, "minrad"), block)
            st = _block_medians(read_trace(work_dir, "strutcount"), block)
            rp = _block_medians(read_trace(work_dir, "ropelength"), block) or {}
            for k in sorted(mr):
                if k in seen or k not in st:
                    continue
                seen.add(k)
                improved_mr = mr[k] > best_mr * (1 + minrad_eps)
                improved_st = st[k] > best_st * (1 + strut_eps)
                stale_mr = 0 if improved_mr else stale_mr + 1
                stale_st = 0 if improved_st else stale_st + 1
                best_mr, best_st = max(best_mr, mr[k]), max(best_st, st[k])
                if k in rp:
                    improved_rop = rp[k] < best_rop * (1 - rop_eps)
                    stale_rop = 0 if improved_rop else stale_rop + 1
                    best_rop = min(best_rop, rp[k])
                    rop_note = f"  rop {rp[k]:.4f} (best {best_rop:.4f}, stale {stale_rop})"
                else:
                    # no ropelength row for this block: do not let a missing
                    # trace silently satisfy the new condition
                    stale_rop = 0
                    rop_note = "  rop -"
                print(f"      step {k:6d}  minRad {mr[k]:.4f} (best {best_mr:.4f}, "
                      f"stale {stale_mr})  struts {st[k]:5.0f} (best {best_st:.0f}, "
                      f"stale {stale_st}){rop_note}")
                gates = {"minrad": stale_mr, "struts": stale_st,
                         "ropelength": stale_rop}
                if k >= min_steps and all(gates[g] >= patience for g in plateau_on):
                    pid = _ridgerunner_pid(Path(src).stem)
                    if pid is None:
                        reason = (f"plateau at step {k} but the ridgerunner process "
                                  f"could not be found; letting it run on")
                        print(f"      {reason}")
                        seen.update(range(k, max_steps + block, block))
                        break
                    names = ", ".join(sorted(plateau_on))
                    print(f"      plateau: {names} flat for {patience} blocks; "
                          f"stopping at step {k}")
                    reason = f"plateau ({names}) at step {k}"
                    import os
                    os.kill(pid, signal.SIGTERM)
                    break
        proc.wait()
    finally:
        if proc.poll() is None:
            proc.terminate()
    if not dst.exists():
        raise CycleError("tighten_link_xyz.py produced no output during auto descent")
    return dst, reason


def refine(src: Path, dst: Path, vu: float, fix_minrad: float,
           mode: str = "subdivide") -> str:
    """Raise resolution to a v/u target, repairing MinRad as we go.

    Refining does not move the curve, so the turning angle at a corner is
    preserved while the adjacent edges are divided: a k-fold refinement divides
    minRad by k. Measured in the guide, 721 -> 1136 vertices took minRad from
    0.500 to 0.297 by plain subdivision and 0.434 by splining, both inadmissible
    against the constraint minRad >= lambda*tau. The violation is not spread out
    (1.2% of vertices were near the limit there), so rounding only the offending
    corners restores it -- which is what --fix-minrad does.

    `subdivide` rather than `spline` by default: the guide chose it on edge
    distribution (max/min 1.01 vs 1.11), because a refined file is an INPUT to a
    run and what a run cares about is the edge distribution, not a fractionally
    lower starting ropelength.
    """
    out = _run([sys.executable, LIB / "refine_xyz.py", src, dst,
                "--vu", vu, "--mode", mode, "--fix-minrad", fix_minrad])
    if not dst.exists():
        raise CycleError(f"refine_xyz.py failed:\n{out[-1500:]}")
    return out.strip().splitlines()[-1] if out.strip() else "refined"


def next_rung(vu: float, ladder: list[float]) -> float | None:
    """Lowest ladder rung the current v/u has not reached."""
    for rung in sorted(ladder):
        if vu < rung * 0.98:      # 2% tolerance so we do not refine for nothing
            return rung
    return None


def maybe_refine(path: Path, results: Path, tag: str, args, ref_generic: Path,
                 scratch: Path) -> tuple[Path, str | None]:
    """Climb one ladder rung if this configuration is below it.

    Deliberately AFTER a descent, never before one. d/D is an output, fixed by
    Rop/(2N), so it falls for free as a run shortens the curve: a raw 7-component
    layout at v/u 0.128 needs 4653 vertices (7.8x) to reach v/u 1.0 at its own
    inflated ropelength, and gets there with no refinement at all once the first
    descent brings Rop under 596.
    """
    if not args.vu_ladder:
        return path, None
    m = measure(path, want_struts=False)
    rung = next_rung(m.vu, args.vu_ladder)
    if rung is None:
        return path, None
    dst = results / f"{tag}_refined_vu{rung:g}.xyz"
    print(f"  refining: v/u {m.vu:.2f} -> {rung:g} "
          f"({m.vertices} vertices -> ~{round(rung * m.rop)})")
    note = refine(path, dst, rung, args.refine_fix_minrad, args.refine_mode)
    print(f"    {note}")
    if args.sym_group:
        sym = dst.with_name(dst.stem + "_sym.xyz")
        moved = symmetrize(dst, sym, args.sym_group, args.sym_order,
                           args.sym_axis, args.no_permute)
        print(f"    re-symmetrized after refinement; largest vertex move {moved}")
        dst = sym
    m2 = measure(dst)
    verdict = homfly_verdict(ref_generic, dst, scratch)
    print(f"    refined: rop {m2.rop:.4f}  v/u {m2.vu:.2f}  minRad/tau "
          f"{m2.minrad_over_tau:.3f}  max d/D {m2.max_dD:.4f}  HOMFLY {verdict}")
    if verdict == "DIFFERENT":
        print("    refinement changed the link; discarding it and keeping the "
              "coarser file")
        return path, "refinement changed the link; discarded"
    if m2.minrad_over_tau < 1.0:
        print("    refinement left minRad below the constraint boundary even "
              "after --fix-minrad; discarding it")
        return path, "refinement violated minRad; discarded"
    return dst, f"refined to v/u {rung:g}"


def read_trace(run_root: Path, name: str) -> np.ndarray | None:
    """Robust trace read. A live run leaves its last line half-written, which
    makes np.loadtxt raise on the whole file."""
    hits = list(run_root.glob(f"*.rr/logfiles/{name}.dat"))
    if not hits:
        return None
    rows = []
    for line in hits[0].read_text(errors="replace").splitlines():
        f = line.split()
        if len(f) >= 2:
            try:
                rows.append([float(f[0]), float(f[1])])
            except ValueError:
                pass
    return np.array(rows) if rows else None


# --------------------------------------------------------------------------- #
# the contraction scan
# --------------------------------------------------------------------------- #

@dataclass
class ScanRow:
    factor: float
    length: float
    minrad: float
    minstrut: float
    rop_est: float
    feasible: bool


def scan(src: Path, scratch: Path, lo: float, hi: float, step: float,
         join: float | None = None) -> list[ScanRow]:
    """Sweep contraction factors and rank on the real feasibility criterion."""
    rows: list[ScanRow] = []
    f = lo
    while f <= hi + 1e-9:
        cand = scratch / f"scan_{f:.3f}.xyz"
        try:
            contract(src, cand, round(f, 4), join)
            m = measure(cand, want_struts=False)
            tau = min(m.minstrut / 2, m.minrad)
            rows.append(ScanRow(round(f, 4), m.length, m.minrad, m.minstrut,
                                m.length / tau, m.minrad > m.minstrut / 2))
        except CycleError:
            pass
        finally:
            cand.unlink(missing_ok=True)
        f += step
    return rows


def best_factor(rows: list[ScanRow]) -> ScanRow | None:
    ok = [r for r in rows if r.feasible]
    return min(ok, key=lambda r: r.rop_est) if ok else None


def scan_is_noop(rows: list[ScanRow], tol: float = 1e-6) -> bool:
    """True when every factor returned the SAME length, i.e. the contraction did
    nothing at all.

    contract_clusters.py works by finding contact bands and translating the
    clusters they form. An input that is not yet tightened has no contacts to
    find -- measured on a raw 7-component layout with 1 strut, every factor from
    0.80 to 0.95 returned length 1048.644 and minRad 0.5185, identical to the
    input. The loop would otherwise report this as "gain too small", which sends
    you looking at --min-gain instead of at the missing descent.
    """
    if len(rows) < 2:
        return False
    lens = [r.length for r in rows]
    return (max(lens) - min(lens)) <= tol * max(1.0, abs(lens[0]))


def timewarp_worthwhile(path: Path) -> tuple[bool, str]:
    """Ask tighten_lib/strut_free.py whether any strut-free stretch remains.

    --Timewarp only accelerates sections carrying no strut and no kink, and the
    man page warns of a cost when there are none. Measured across this campaign
    the verdict flipped from "WORTH IT" (3.0% of vertices in contact, 4.3 D
    longest free arc) to "pure overhead" (85.4% in contact, no free run at all)
    over the course of one link, so it is not a set-once flag.
    """
    tool = LIB / "strut_free.py"
    if not tool.is_file():
        return False, "strut_free.py not available; leaving --Timewarp off"
    out = _run([sys.executable, tool, path])
    for line in out.splitlines():
        if path.name in line:
            if "WORTH IT" in line:
                return True, line.strip()
            return False, line.strip()
    return False, "could not read a strut_free verdict; leaving --Timewarp off"


def choose_rr_flags(m: Measurement, path: Path, dd_threshold: float
                    ) -> tuple[list[str], list[str], str]:
    """Pick (wrapper flags, ridgerunner extras, rationale) from the measurement.

    This is the guide's Step 2 table:
        max d/D > 0.5        -> --eq --Timewarp   (the curve chords across its
                                                   own tube; no contact set worth
                                                   keeping)
        residual > 0.1       -> --eq --Timewarp   (far from critical; conditioning
                                                   is what limits descent)
        residual < 0.1 and
        max d/D acceptable   -> --no-eq           (the contact set is the asset)
    with the standing rule that --Timewarp is never used without --eq, and the
    added check that it is only worth anything while a free stretch exists.
    """
    coarse = m.max_dD > dd_threshold
    far = m.residual is None or m.residual > 0.1
    if not (coarse or far):
        return [], [], (f"max d/D {m.max_dD:.4f} <= {dd_threshold} and residual "
                        f"{m.residual:.4f} < 0.1: --no-eq, no Timewarp")
    if m.residual is None:
        why_res = ("residual not measurable (tau %.6f is off the nominal 0.5, so "
                   "the constraint set is ambiguous) -- treated as far from "
                   "critical" % m.tau)
    else:
        why_res = f"residual {m.residual:.4f} > 0.1"
    why = []
    if coarse:
        why.append(f"max d/D {m.max_dD:.4f} > {dd_threshold}")
    if far:
        why.append(why_res)
    tw, tw_note = timewarp_worthwhile(path)
    extras = ["--Timewarp"] if tw else []
    return (["--eq"], extras,
            f"{' and '.join(why)}: --eq" + (" --Timewarp" if tw else "")
            + f"  [{tw_note}]")


def hole_wall_radius(path: Path, tau: float) -> float:
    """Radius of the empty cylinder about z: min centreline radius minus tau.

    This is the hole-squeeze's whole budget. Measured on the structure that
    settled the question: a 5-component annulus stuck at 213.86 for three rounds
    had a wall radius of 1.44 D -- an empty cylinder 3 D across that uniform
    cluster contraction could not touch, because its binding contact formed out
    in the ring body (radial compression at r ~ 3) long before the core filled.
    """
    comps = [np.asarray(c, float) for c in read_xyz(path)]
    P = np.vstack(comps)
    P = P - P.mean(axis=0)
    return float(np.hypot(P[:, 0], P[:, 1]).min() - tau)


@dataclass
class SqueezeRow:
    factor: float
    length: float
    minrad: float
    minstrut: float
    margin: float          # minRad / (minStrut/2): >1 = strut-limited = repairable
    admissible: bool
    axis: tuple | None = None    # None = z through the centroid
    point: tuple | None = None
    # kept last so the positional axis/point call sites stay valid
    d_len: float = 0.0     # fraction of length the squeeze removes
    d_tau: float = 0.0     # fraction of thickness it destroys
    sens: float | None = None   # d_tau / d_len -- see scan_squeeze


def scan_squeeze(src: Path, scratch: Path, factors: list[float], blend: float,
                 margin: float, axis=None, point=None) -> list[SqueezeRow]:
    """Sweep hole-squeeze factors, reporting both the margin and the SENSITIVITY.

    A squeeze is judged the opposite way from a contraction: its immediate
    ropelength is meaningless (it deliberately breaks thickness for RidgeRunner
    to repair), so the admissibility test is only that minRad stays above
    minStrut/2 with some margin -- curvature must never become the binding
    constraint, because proximity damage is repairable and curvature damage is
    what stalls rounds.

    That margin is a RATIO, and a ratio cannot see a squeeze that crushes minRad
    and minStrut together: both fall, the quotient holds, and the gate waves
    through a move that has destroyed the thickness. Worse, the ratio is not
    monotone in the factor. Measured on the 7-component link at 276.080, with
    the driver's own blend:

        f     dL/L   dtau/tau   sens   margin   gate says
        0.95  1.09%     2.8%    2.59    1.01    rejected
        0.80  4.26%    12.4%    2.92    0.99    rejected
        0.40 11.96%    36.3%    3.03    0.98    rejected
        0.30 13.67%    47.5%    3.47    1.12    ADMISSIBLE
        0.20 15.27%    62.3%    4.08    1.45    ADMISSIBLE

    The margin sags in the middle and climbs again at the extreme, so the gate
    rejects every gentle factor and admits only the most violent ones -- and
    `choose_move` then takes the hardest admissible. Both f=0.20 rounds run from
    that configuration LOST (+1.12 and +1.06) while four squeezes at sens 1.75
    to 3.17 all won.

    SENSITIVITY is the quantity that separates them:

        sens = (dtau/tau) / (dL/L)

    how much thickness the move destroys per unit of length it removes. It is
    the slack-versus-packing question made measurable: if the hole is genuine
    slack the strands move inward freely and tau barely stirs; if the strands
    ringing the hole are already packed against each other, moving them inward
    forces contact and tau collapses. Over every squeeze this project has run,
    wins scored 1.75 / 2.14 / 2.71 / 3.17 and the one distinct losing
    configuration scored 4.08.

    It is REPORTED, not gated on: the ordering is trustworthy and physically
    motivated, but a threshold drawn from four wins and one distinct loss is
    not, and changing the measurement and the selection rule in the same step
    would make the next result impossible to attribute. Read the column; if it
    is above roughly 3.5, expect the round to lose.
    """
    sys.path.insert(0, str(LIB))
    from radial_squeeze import apply_squeeze, apply_squeeze_about
    comps = [np.asarray(c, float) for c in read_xyz(src)]
    m_in = measure(src, want_struts=False)
    len0, tau0 = m_in.length, m_in.tau
    rows = []
    for f in factors:
        if point is not None:
            sq = apply_squeeze_about(comps, point, axis, f, blend)
        else:
            sq, _, _ = apply_squeeze(comps, [0, 0, 1], f, blend)
        cand = scratch / f"sq_{f:.3f}.xyz"
        write_xyz(cand, [c.tolist() for c in sq], 12)
        try:
            m = measure(cand, want_struts=False)
            ratio = m.minrad / (m.minstrut / 2)
            d_len = (len0 - m.length) / len0 if len0 else 0.0
            d_tau = (tau0 - m.tau) / tau0 if tau0 else 0.0
            # a factor that removes nothing has no meaningful sensitivity
            sens = (d_tau / d_len) if d_len > 1e-9 else None
            rows.append(SqueezeRow(f, m.length, m.minrad, m.minstrut, ratio,
                                   ratio >= margin,
                                   tuple(axis) if axis is not None else None,
                                   tuple(point) if point is not None else None,
                                   d_len=d_len, d_tau=d_tau, sens=sens))
        finally:
            cand.unlink(missing_ok=True)
    return rows


def choose_move(current: Path, m_in, args, scratch: Path):
    """Pick this round's geometric move from the measured structure.

    Returns (move, payload, rationale):
      ('squeeze',  SqueezeRow, why)   -- hole-squeeze at that factor
      ('contract', ScanRow,    why)   -- cluster contraction at that factor
      (None,       None,       why)   -- neither move can help

    The rule, from the measured cases that decided it:

      1. Both moves close empty space; they differ in WHICH space. The
         hole-squeeze targets the void ENCLOSED by the ring; contraction the
         space BETWEEN clusters. Measure both budgets, in that order.
      2. Prefer the squeeze whenever it has budget: it is topology-safe by
         construction (monotone radial map -- no HOMFLY gate), preserves Cn and
         the z-mirror exactly, and its damage is proximity, the repairable kind.
         Head to head on the same 213.86 structure, contraction's best
         clearance-bound outcome was 215.7 (a certain loss) while the squeeze
         reached 207.4 within 2800 steps of repair.
      3. Contraction earns its keep where the slack is genuinely inter-cluster:
         thin necklace-type rings (its 2082 -> 360 run on a 7-component
         necklace), or when there is no enclosed void. Its immediate Rop_est
         must beat the input by --min-gain, and it pays with a HOMFLY gate.
    """
    ax = pt = None
    if args.find_axis:
        sys.path.insert(0, str(LIB))
        from find_axis import find_axis
        comps = [np.asarray(c, float) for c in read_xyz(current)]
        found = find_axis(comps, ndirs=args.axis_dirs,
                          probe=0.5, blend=args.squeeze_blend)
        if found is None:
            hole = -1.0
            parts = ["--find-axis: no empty tunnel found through the structure"]
        else:
            ax, pt, clear, rem = found
            hole = clear - m_in.tau
            parts = [f"--find-axis: axis ({ax[0]:+.3f},{ax[1]:+.3f},{ax[2]:+.3f}), "
                     f"wall radius {hole:.3f} D, probe squeeze removes "
                     f"{rem * 100:.1f}% of length"]
    else:
        hole = hole_wall_radius(current, m_in.tau)
        parts = [f"hole wall radius {hole:.3f} D about z "
                 f"(budget threshold {args.min_hole})"]

    if hole >= args.min_hole:
        srows = scan_squeeze(current, scratch, args.squeeze_factors,
                             args.squeeze_blend, args.squeeze_margin,
                             axis=ax, point=pt)
        ok = [r for r in srows if r.admissible]
        for r in srows:
            flag = ""
            if r.sens is not None and r.sens > 3.5:
                # reported, not enforced -- see scan_squeeze
                flag = "  <- sens high, expect a loss"
            parts.append(f"  squeeze f={r.factor:.2f}: len {r.length:9.3f}  "
                         f"dL {r.d_len * 100:5.2f}%  dtau {r.d_tau * 100:5.1f}%  "
                         f"sens {fmt(r.sens, '5.2f')}  "
                         f"minRad/(clear/2) {r.margin:5.2f}  "
                         f"{'admissible' if r.admissible else 'too kinked'}{flag}")
        if ok:
            pick = min(ok, key=lambda r: r.factor)   # hardest safe squeeze
            sym_note = ("axis found by search" if args.find_axis
                        else "symmetry preserved")
            parts.append(f"-> SQUEEZE f={pick.factor:.2f}: removes "
                         f"{(m_in.length - pick.length) / m_in.length * 100:.1f}% of "
                         f"length, topology guaranteed, {sym_note}")
            return "squeeze", pick, "\n  ".join(parts)
        parts.append("no admissible squeeze factor despite the open hole")
    else:
        parts.append("hole closed; squeeze has no budget")

    rows = scan(current, scratch, args.factor_lo, args.factor_hi,
                args.factor_step, args.join)
    if scan_is_noop(rows):
        return None, None, "\n  ".join(parts + [
            "contraction is a no-op (no contact clusters found); needs a descent"])
    pick = best_factor(rows)
    if pick is not None and m_in.rop - pick.rop_est >= args.min_gain:
        parts.append(f"-> CONTRACT f={pick.factor:.3f}: predicted rop "
                     f"{pick.rop_est:.3f} (gain {m_in.rop - pick.rop_est:.3f})")
        return "contract", pick, "\n  ".join(parts)
    if pick is None:
        parts.append("no feasible contraction factor")
    else:
        parts.append(f"best contraction gain {m_in.rop - pick.rop_est:.3f} < "
                     f"--min-gain {args.min_gain}")
    return None, None, "\n  ".join(parts + [
        "neither move has budget: refine (--vu-ladder) or accept this basin"])


# --------------------------------------------------------------------------- #
# the cycle
# --------------------------------------------------------------------------- #

@dataclass
class Round:
    index: int
    move: str = "contract"
    factor: float | None = None
    rop_contracted: float | None = None
    rop_descended: float | None = None
    gain: float | None = None
    homfly: str = "-"
    struts: int | None = None
    residual: float | None = None
    minrad_over_tau: float | None = None
    sym_dev: float | None = None
    note: str = ""
    kept: str = ""


def fmt(v, spec="8.3f"):
    if v is None:
        # Pad to the field width the spec asks for, so a row of missing values
        # lines up with the header instead of collapsing into "- - -".
        width = re.match(r"[+-]?(\d+)", spec)
        return "-".rjust(int(width.group(1)) if width else 1)
    return format(v, spec)


@dataclass
class Resumed:
    """What could be recovered from an interrupted cycle's work directory."""
    src: Path
    origin: Path                 # round0 input, the reference's true source
    current: Path                # what the next round starts from
    start_round: int
    best_rop: float
    best_path: Path
    rounds: list


def round_output(results: Path, i: int) -> Path | None:
    """The file round i carried forward, or None if it never finished.

    A round can end on any of three files depending on which options were in
    play, and they are produced in this order, so the most derived one that
    exists is the one that became the next round's input.
    """
    for pat in (f"round{i}_refined_vu*_sym.xyz", f"round{i}_refined_vu*.xyz",
                f"round{i}_descended_sym.xyz", f"round{i}_descended.xyz"):
        hits = sorted(results.glob(pat))
        if hits:
            return hits[-1]
    return None


def resume_state(work: Path, results: Path, args) -> Resumed:
    """Rebuild a cycle's state from what an interrupted run left on disk.

    The driver holds the round index, the running best and the ledger rows in
    memory and writes the ledger only when the whole cycle finishes, so a crash
    loses every one of them while the per-round .xyz files survive. Those files
    are the only durable record, and this reconstructs the state from them,
    using the ledger merely to enrich it when one happens to exist.

    Two things it deliberately does NOT do:

      * It does not re-derive the reference polynomial from the resumed
        configuration. Doing so would adopt whatever the interrupted run had
        drifted to as the new definition of "unchanged", which is exactly the
        check the gate exists to make. The caller rebuilds the reference from
        `origin` -- round 0's input, which is preserved.

      * It does not resume a descent that was cut off part-way. A round whose
        move was applied but whose descent never produced an output is redone
        from that round's input, and the steps already spent are reported as
        lost rather than quietly dropped.
    """
    # The reference was recorded from whatever round 0 carried forward, which is
    # the symmetrized input when the up-front symmetrize ran. Rebuilding from the
    # unsymmetrized one would fail the integrity check for no real reason.
    origin = next((results / name for name in
                   ("round0_input_sym.xyz", "round0_input.xyz")
                   if (results / name).is_file()), None)
    if origin is None:
        raise CycleError(
            f"--resume found no round0 input under {results}. Point --work-dir "
            f"at a directory an earlier cycle actually wrote, or drop --resume "
            f"to start a new cycle.")

    last = 0
    while round_output(results, last + 1) is not None:
        last += 1

    # a round whose move landed but whose descent did not
    partial = sorted(results.glob(f"round{last + 1}_*.xyz"))
    if partial:
        spent = ""
        run_dir = work / f"run_round{last + 1}"
        traces = sorted(run_dir.glob("*.rr/logfiles/ropelength.dat"))
        if traces:
            vals = []
            for ln in traces[-1].read_text().splitlines(keepends=True):
                if not ln.endswith("\n"):
                    continue          # a crash leaves a half-written final line
                r = ln.split()
                if len(r) == 2:
                    try:
                        vals.append((int(r[0]), float(r[1])))
                    except ValueError:
                        pass
            if vals:
                spent = (f"; its descent had reached step {vals[-1][0]} "
                         f"(best {min(v for _, v in vals):.4f}) -- those steps are lost, "
                         f"resume is at round granularity")
        print(f"  round {last + 1} was interrupted after its move but before it "
              f"produced an output{spent}.\n  Redoing it from round {last}'s result; "
              f"the partial files are left in place.")

    # prefer the ledger's own record of each round, fall back to the files
    rounds: list = []
    ledger = work / "ledger.csv"
    if ledger.is_file():
        with ledger.open(newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    idx = int(row["index"])
                except (KeyError, TypeError, ValueError):
                    continue
                if idx > last:
                    continue
                rd = Round(index=idx)
                for k, v in row.items():
                    if k == "index" or not hasattr(rd, k):
                        continue
                    if v in ("", "None"):
                        continue
                    f = next(fl for fl in fields(Round) if fl.name == k)
                    try:
                        setattr(rd, k, int(v) if f.type.startswith("int")
                                else float(v) if "float" in f.type else v)
                    except ValueError:
                        setattr(rd, k, v)
                rounds.append(rd)
    have = {r.index for r in rounds}
    for i in range(1, last + 1):
        if i in have:
            continue
        out = round_output(results, i)
        m = measure(out, want_struts=False)
        rounds.append(Round(index=i, rop_descended=m.rop, kept=out.name,
                            note="reconstructed on resume (no ledger row)"))
    rounds.sort(key=lambda r: r.index)

    # the running best, measured rather than trusted
    best_path, best_rop = origin, measure(origin, want_struts=False).rop
    for i in range(1, last + 1):
        out = round_output(results, i)
        rop = measure(out, want_struts=False).rop
        if rop < best_rop:
            best_rop, best_path = rop, out

    current = round_output(results, last) if last else origin
    return Resumed(src=origin, origin=origin, current=current,
                   start_round=last + 1, best_rop=best_rop,
                   best_path=best_path, rounds=rounds)


def run_cycle(args) -> int:
    work = Path(args.work_dir).expanduser().resolve()
    scratch = work / "scratch"
    results = work / "results"
    for d in (work, scratch, results):
        d.mkdir(parents=True, exist_ok=True)

    resumed = None
    if args.resume:
        print(f"work dir      {work}")
        print("resuming from what the previous run left on disk")
        resumed = resume_state(work, results, args)
        src = resumed.src
        current = resumed.current
        print(f"  round 0 input {resumed.origin.name}")
        print(f"  completed rounds 1..{resumed.start_round - 1}, "
              f"best so far {resumed.best_rop:.4f} ({resumed.best_path.name})")
        print(f"  next round    {resumed.start_round} of {args.rounds}, starting "
              f"from {Path(current).name}")
        mo = measure(resumed.origin, want_struts=False)
        print(f"  cycle started at ropelength {mo.rop:.4f}; this invocation "
              f"resumes at {measure(current, want_struts=False).rop:.4f}")
        if resumed.start_round > args.rounds:
            print(f"\nnothing to do: rounds 1..{args.rounds} are already complete. "
                  f"Raise --rounds to continue.")
            return 0
    else:
        src = Path(args.input).expanduser().resolve()
        current = results / "round0_input.xyz"
        shutil.copy(src, current)
        print(f"work dir      {work}")
        print(f"input         {src}")

    m0 = measure(current)
    label = "resume at" if resumed else "start"
    print(f"{label:<14}ropelength {m0.rop:.4f}  "
          f"minRad/tau {m0.minrad_over_tau:.3f}  struts {m0.struts}  "
          f"v/u {m0.vu:.2f}  max d/D {m0.max_dD:.4f}")

    # Symmetrize the INPUT only if it is already close enough to its symmetry.
    # The rule "tighten first, symmetrize after" exists because
    # plc_build_symmetry has no distance tolerance, so it is really a rule about
    # distance: measure, then decide. A layout drawn symmetric can be projected
    # immediately and for free; one in the raw-4_LC band must be descended first
    # or symmetrizing destroys the link.
    if args.sym_group and args.sym_order and not resumed:
        dev = symmetry_deviation(current, args.sym_order)
        safe, why = symmetrize_is_safe(dev)
        print(f"\nsymmetry      C{args.sym_order} deviation: "
              + ("unmeasurable" if dev is None else
                 f"median {dev[0]:.4f}, max {dev[1]:.3f}, {dev[2] * 100:.2f}% over "
                 f"half an edge")
              + f"\n              {why}")
        if safe and not args.dry_run:
            pre_sym = results / "round0_input_sym.xyz"
            moved = symmetrize(current, pre_sym, args.sym_group, args.sym_order,
                               args.sym_axis, args.no_permute)
            m_sym = measure(pre_sym)
            print(f"              symmetrized up front: largest vertex move {moved}, "
                  f"ropelength {m0.rop:.4f} -> {m_sym.rop:.4f}")
            current, m0 = pre_sym, m_sym

    # The reference is built from round 0's input, NOT from the configuration a
    # resumed run happens to be sitting on. Re-deriving it here would adopt any
    # drift the earlier rounds introduced as the new definition of "unchanged",
    # which is the one thing the gate exists to catch.
    ref_src = resumed.origin if resumed else current
    if resumed:
        print("\nreference polynomial: rebuilding from round 0's input")
    else:
        print("\nrecording the reference polynomial (once, before anything is modified)")
    ref_generic = to_generic_frame(ref_src, scratch / "reference_generic.xyz")
    if args.dry_run:
        print("  [dry run] skipped")
    else:
        poly = homfly_reference(ref_src, scratch)
        stored = work / "reference_homfly.txt"
        if resumed and stored.is_file():
            was = stored.read_text().strip()
            if was and was != poly:
                raise CycleError(
                    "the reference polynomial rebuilt from round 0's input does not "
                    "match the one this work directory recorded. Either --work-dir "
                    "points at a different cycle's tree, or round0_input.xyz has been "
                    "edited. Refusing to resume against a reference that is not the "
                    "one the earlier rounds were gated on.")
            print("  matches the polynomial recorded by the interrupted run")
        stored.write_text(poly + "\n")
        print(f"  {poly[:96]}{'...' if len(poly) > 96 else ''}")

    if resumed:
        best_rop, best_path = resumed.best_rop, resumed.best_path
        rounds: list[Round] = list(resumed.rounds)
        m0 = measure(resumed.origin)       # report progress against the true start
    else:
        best_rop, best_path = m0.rop, current
        rounds = []

    # ------------------------------------------------------------------ #
    # Initial descent. The cycle STARTS with a contraction, and a contraction
    # needs contact clusters to translate -- which a raw diagram does not have.
    # Measured on a raw 7-component layout (max d/D 3.99, 1 strut): every
    # factor returned the input's own length unchanged. So an untightened input
    # has to be descended once before the loop can do anything.
    # ------------------------------------------------------------------ #
    if m0.max_dD > args.dd_threshold and not args.initial_steps and not args.dry_run:
        print(f"\nWARNING: max d/D {m0.max_dD:.4f} > {args.dd_threshold}: this input "
              f"is not tightened.\n         The contraction has no contact clusters "
              f"to find and will be a no-op.\n         Pass --initial-steps to "
              f"descend first.")

    if args.initial_steps and not args.dry_run and not resumed:
        wrap, extras, why = ([], args.rr_arg, "flags as given")
        if args.auto_flags:
            wrap, extras, why = choose_rr_flags(m0, current, args.dd_threshold)
            extras = extras + args.rr_arg
        print(f"\n{'=' * 72}\ninitial descent: {args.initial_steps} steps\n{'=' * 72}")
        print(f"  flags: {why}")
        pre = results / "round0_descended.xyz"
        descend(current, pre, work / "run_initial", args.initial_steps,
                args.symmetry, extras, args.snapshot_interval, args.stop_res, wrap)
        if args.sym_group:
            symmetrize(pre, pre.with_name("round0_descended_sym.xyz"),
                       args.sym_group, args.sym_order, args.sym_axis, args.no_permute)
            pre = pre.with_name("round0_descended_sym.xyz")
        m_pre = measure(pre)
        verdict = homfly_verdict(ref_generic, pre, scratch)
        print(f"  after initial descent: rop {m_pre.rop:.4f}  minRad/tau "
              f"{m_pre.minrad_over_tau:.3f}  struts {m_pre.struts}  "
              f"max d/D {m_pre.max_dD:.4f}  HOMFLY {verdict}")
        if verdict == "DIFFERENT":
            print("  STOP: the initial descent changed the link")
            return 1
        current = pre
        if m_pre.rop < best_rop:
            best_rop, best_path = m_pre.rop, pre
        current, rnote = maybe_refine(current, results, "round0", args,
                                      ref_generic, scratch)
        if rnote:
            print(f"  {rnote}")

    for i in range(resumed.start_round if resumed else 1, args.rounds + 1):
        rd = Round(index=i)
        print(f"\n{'=' * 72}\nround {i}\n{'=' * 72}")

        m_in = measure(current)
        print(f"  input   rop {m_in.rop:.4f}  minRad/tau {m_in.minrad_over_tau:.3f}  "
              f"struts {m_in.struts}  residual {fmt(m_in.residual, '.4f')}")

        if m_in.vu < args.min_vu:
            rd.note = (f"v/u {m_in.vu:.2f} below --min-vu {args.min_vu}; refine with "
                       f"tighten_lib/refine_xyz.py --fix-minrad before continuing")
            print(f"  STOP: {rd.note}")
            rounds.append(rd)
            break

        pick = None
        did_squeeze = False
        if args.choose_move:
            move, payload, why = choose_move(current, m_in, args, scratch)
            print("  move selection:\n  " + why)
            if move is None:
                rd.note = why.splitlines()[-1].strip()
                print(f"  STOP: {rd.note}")
                rounds.append(rd)
                break
            if args.dry_run:
                rd.move, rd.factor = move, payload.factor
                rd.note = f"dry run: would {move} at f={payload.factor:.3f}"
                print(f"  [dry run] stopping before applying the {move}")
                rounds.append(rd)
                break
            if move == "squeeze":
                sys.path.insert(0, str(LIB))
                from radial_squeeze import apply_squeeze, apply_squeeze_about
                comps = [np.asarray(c, float) for c in read_xyz(current)]
                if payload.point is not None:
                    print(f"  squeeze axis ({payload.axis[0]:+.3f},"
                          f"{payload.axis[1]:+.3f},{payload.axis[2]:+.3f}) "
                          f"through ({payload.point[0]:+.2f},"
                          f"{payload.point[1]:+.2f},{payload.point[2]:+.2f})")
                    sq = apply_squeeze_about(comps, payload.point, payload.axis,
                                             payload.factor, args.squeeze_blend)
                else:
                    sq, _, _ = apply_squeeze(comps, [0, 0, 1], payload.factor,
                                             args.squeeze_blend)
                contracted = results / f"round{i}_squeezed_f{payload.factor:.3f}.xyz"
                write_xyz(contracted, [c.tolist() for c in sq], 12)
                rd.move, rd.factor = "squeeze", payload.factor
                if args.sym_group:
                    moved = symmetrize(contracted,
                                       contracted.with_name(contracted.stem + "_sym.xyz"),
                                       args.sym_group, args.sym_order,
                                       args.sym_axis, args.no_permute)
                    contracted = contracted.with_name(contracted.stem + "_sym.xyz")
                    print(f"  re-symmetrized after squeeze; largest vertex move {moved}")
                m_c = measure(contracted)
                rd.rop_contracted = m_c.rop
                # A monotone radial map is a homeomorphism: the link type
                # cannot change, so no pre-descent HOMFLY gate is needed.
                # The post-descent gate still runs as usual.
                rd.homfly = "SAFE"
                print(f"  squeezed    rop {m_c.rop:.4f} (immediate; thickness repair "
                      f"pending)  minRad/(clear/2) "
                      f"{m_c.minrad / (m_c.minstrut / 2):.2f}  topology guaranteed")
                did_squeeze = True
            else:
                pick = payload
                rd.move = "contract"

        if not did_squeeze and pick is None:
            print(f"  scanning factors {args.factor_lo}..{args.factor_hi} step {args.factor_step}")
            rows = scan(current, scratch, args.factor_lo, args.factor_hi, args.factor_step,
                        args.join)
            for r in rows:
                mark = "feasible" if r.feasible else "curvature-limited"
                print(f"    f={r.factor:.3f}  len {r.length:9.3f}  minRad {r.minrad:.4f}  "
                      f"Rop~{r.rop_est:9.3f}  {mark}")

            if scan_is_noop(rows):
                rd.note = (f"contraction is a no-op: every factor returned the input's "
                           f"own length, so no contact clusters were found. The input "
                           f"has {m_in.struts} struts and max d/D {m_in.max_dD:.4f}; it "
                           f"needs a descent (--initial-steps), not a contraction")
                print(f"  STOP: {rd.note}")
                rounds.append(rd)
                break

            pick = best_factor(rows)
            if pick is None:
                rd.note = "no feasible contraction factor: the window is closed"
                print(f"  STOP: {rd.note}")
                rounds.append(rd)
                break
            if m_in.rop - pick.rop_est < args.min_gain:
                rd.note = (f"best predicted gain {m_in.rop - pick.rop_est:.3f} < "
                           f"--min-gain {args.min_gain}")
                print(f"  STOP: {rd.note}")
                rounds.append(rd)
                break

        if not did_squeeze:
            rd.factor = pick.factor
            print(f"  chose f={pick.factor:.3f}  predicted rop {pick.rop_est:.3f} "
                  f"(gain {m_in.rop - pick.rop_est:.3f})")

            if args.dry_run:
                rd.note = "dry run: stopped before contracting"
                rounds.append(rd)
                break

            contracted = results / f"round{i}_contracted_f{pick.factor:.3f}.xyz"
            contract(current, contracted, pick.factor, args.join)

            if args.sym_group:
                moved = symmetrize(contracted, contracted.with_name(contracted.stem + "_sym.xyz"),
                                   args.sym_group, args.sym_order, args.sym_axis,
                                   args.no_permute)
                contracted = contracted.with_name(contracted.stem + "_sym.xyz")
                print(f"  re-symmetrized ({args.sym_group}, p={args.sym_order}); "
                      f"largest vertex move {moved}")

            m_c = measure(contracted)
            rd.rop_contracted = m_c.rop
            print(f"  contracted  rop {m_c.rop:.4f}  minRad/tau {m_c.minrad_over_tau:.3f}  "
                  f"struts {m_c.struts}  strut-limited {m_c.strut_limited}")

            rd.homfly = homfly_verdict(ref_generic, contracted, scratch)
            print(f"  HOMFLY gate: {rd.homfly}")
            if rd.homfly == "DIFFERENT":
                rd.note = "link changed; round abandoned, previous configuration kept"
                print(f"  STOP: {rd.note}")
                rounds.append(rd)
                break
            if rd.homfly == "UNKNOWN" and args.strict_topology:
                rd.note = "HOMFLY uncomputable and --strict-topology set; stopping"
                print(f"  STOP: {rd.note}")
                rounds.append(rd)
                break

        if m_c.rop < best_rop:
            best_rop, best_path = m_c.rop, contracted

        descended = results / f"round{i}_descended.xyz"
        run_dir = work / f"run_round{i}"
        wrap, extras = [], args.rr_arg
        if args.auto_flags:
            wrap, auto_extras, why = choose_rr_flags(m_c, contracted, args.dd_threshold)
            extras = auto_extras + args.rr_arg
            print(f"  flags: {why}")
        if args.auto_steps:
            gate, gate_why = resolve_plateau_on(
                args.plateau_on, rd.move, rd.factor, args.plateau_hard_f)
            print(f"  reconditioning: auto (cap {args.max_steps}, block "
                  f"{args.step_block}, patience {args.patience})"
                  + (f", --symmetry {args.symmetry}" if args.symmetry else ""))
            print(f"  plateau gate: {', '.join(sorted(gate))}  [{gate_why}]")
            _, why = descend_auto(
                contracted, descended, run_dir, args.symmetry, extras,
                args.snapshot_interval, args.stop_res, args.step_block,
                args.patience, args.min_steps, args.max_steps,
                args.minrad_eps, args.strut_eps, args.rop_eps,
                gate, wrapper=wrap)
            rd.note = why
            print(f"  stopped: {why}")
        else:
            print(f"  reconditioning: {args.steps} steps"
                  + (f", --symmetry {args.symmetry}" if args.symmetry else ""))
            descend(contracted, descended, run_dir, args.steps,
                    args.symmetry, extras, args.snapshot_interval, args.stop_res,
                    wrap)

        m_d = measure(descended)
        rd.rop_descended = m_d.rop
        rd.struts, rd.residual = m_d.struts, m_d.residual
        rd.minrad_over_tau = m_d.minrad_over_tau
        if args.sym_order:
            rd.sym_dev = symmetry_deviation(descended, args.sym_order)
        print(f"  descended   rop {m_d.rop:.4f}  minRad/tau {m_d.minrad_over_tau:.3f}  "
              f"struts {m_d.struts}  residual {fmt(m_d.residual, '.4f')}")
        print(f"              reconditioning: minRad/tau {m_c.minrad_over_tau:.3f} -> "
              f"{m_d.minrad_over_tau:.3f}, struts {m_c.struts} -> {m_d.struts}")

        verdict = homfly_verdict(ref_generic, descended, scratch)
        if verdict == "DIFFERENT":
            rd.note = "link changed during descent; keeping the contracted file"
            rd.kept = contracted.name
            print(f"  STOP: {rd.note}")
            rounds.append(rd)
            break

        if m_d.rop < best_rop:
            best_rop, best_path = m_d.rop, descended

        # Carry the RECONDITIONED file forward even when it is not the shortest:
        # its minRad and strut headroom are what make the next contraction
        # feasible. Best-by-ropelength is tracked separately.
        current = descended
        rd.kept = descended.name
        rd.gain = (m_in.rop - m_d.rop)
        current, rnote = maybe_refine(current, results, f"round{i}", args,
                                      ref_generic, scratch)
        if rnote:
            rd.note = (rd.note + "; " if rd.note else "") + rnote
            rd.kept = Path(current).name
        rounds.append(rd)

    # ---------------------------------------------------------------- report
    print(f"\n{'=' * 72}\nsummary\n{'=' * 72}")
    hdr = (f"{'rd':>3} {'move':>8} {'f':>6} {'contracted':>11} {'descended':>10} "
           f"{'net':>7} {'struts':>7} {'resid':>7} {'mr/tau':>7} {'HOMFLY':>8}  note")
    print(hdr)
    for r in rounds:
        print(f"{r.index:3d} {r.move[:8]:>8} {fmt(r.factor, '6.3f')} {fmt(r.rop_contracted, '11.3f')} "
              f"{fmt(r.rop_descended, '10.3f')} {fmt(r.gain, '+7.3f')} "
              f"{'-' if r.struts is None else r.struts:>7} {fmt(r.residual, '7.4f')} "
              f"{fmt(r.minrad_over_tau, '7.3f')} {r.homfly:>8}  {r.note}")
    print(f"\nstart {m0.rop:.4f}  ->  best {best_rop:.4f}  "
          f"({(best_rop - m0.rop) / m0.rop * 100:+.2f}%)")
    print(f"best configuration: {best_path}")

    final = work / "BEST.xyz"
    shutil.copy(best_path, final)
    print(f"copied to {final}")

    with (work / "ledger.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(Round(0)).keys()))
        w.writeheader()
        for r in rounds:
            w.writerow(asdict(r))
    (work / "ledger.json").write_text(json.dumps(
        {"input": str(src), "start_rop": m0.rop, "best_rop": best_rop,
         "best_path": str(best_path), "rounds": [asdict(r) for r in rounds]}, indent=2))
    print(f"ledger: {work / 'ledger.csv'}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", nargs="?",
                   help="starting .xyz; blank lines separate components. Omit it "
                        "with --resume, which takes its input from the work dir.")
    p.add_argument("--work-dir", required=True,
                   help="scratch tree. Keep it OUT of Dropbox: RidgeRunner rewrites a "
                        "multi-megabyte constraint matrix continuously.")
    p.add_argument("--gui", action="store_true",
                   help="open a window to fill in these options instead of "
                        "typing them. Also the default when the script is run "
                        "with no arguments at all. Needs tkinter, which ships "
                        "with Python; every field in the window carries a '?' "
                        "that shows this same help text.")
    p.add_argument("--resume", action="store_true",
                   help="continue an interrupted cycle in --work-dir instead of "
                        "starting a new one. The driver keeps the round index, the "
                        "running best and the ledger rows in memory and writes the "
                        "ledger only at the end, so a crash loses all of it while the "
                        "per-round .xyz files survive; this rebuilds the state from "
                        "those files. It restarts at the first round that produced no "
                        "output, rebuilds the HOMFLY reference from round 0's input "
                        "(not from the resumed configuration) and refuses to continue "
                        "if that reference disagrees with the recorded one. A descent "
                        "cut off part-way is NOT resumed: its round is redone and the "
                        "spent steps are reported.")
    p.add_argument("--rounds", type=int, default=6, help="max cycles (default 6)")
    p.add_argument("--steps", type=int, default=30000,
                   help="fixed reconditioning steps per round (default 30000). "
                        "Ignored under --auto-steps.")
    p.add_argument("--stop-res", type=float, default=1e-4,
                   help="stop a descent early if RidgeRunner's residual falls "
                        "below this (default 1e-4, the author's convergence "
                        "benchmark). Our runs rarely reach it -- typical "
                        "finished values are 0.2 to 0.4 -- so in practice the "
                        "step cap or the plateau test ends a descent first.")
    p.add_argument("--snapshot-interval", type=int, default=1000,
                   help="write a snapshot every N steps (default 1000). The "
                        "best configuration is recovered from these, so the "
                        "interval bounds how much of a mid-run optimum can be "
                        "missed; extract_best.py reports that penalty rather "
                        "than hiding it. Smaller costs disk, not compute.")

    a = p.add_argument_group(
        "automatic step count",
        "Stop each descent when RECONDITIONING plateaus instead of running a fixed "
        "budget. The descent's job in a late round is to restore the minRad and "
        "strut headroom the next contraction spends, and that finishes long before "
        "the step cap: on one measured round minRad's block median peaked at step "
        "3000 and then only oscillated for 12000 more steps, while struts peaked at "
        "8000 of 15038. The defaults are calibrated on that round and on a 400-step "
        "round which stopped INSIDE the initial contact-shedding transient and left "
        "the next round only f=0.99 feasible -- which is why --min-steps is not small.")
    a.add_argument("--auto-steps", action="store_true",
                   help="decide each round's step count from the run's own traces")
    a.add_argument("--step-block", type=int, default=1000,
                   help="trace block for the medians (default 1000). minRad jitters "
                        "hard step to step; never test it on raw values.")
    a.add_argument("--patience", type=int, default=3,
                   help="stop after this many consecutive blocks in which NEITHER "
                        "minRad nor struts beat their best so far (default 3)")
    a.add_argument("--min-steps", type=int, default=4000,
                   help="never stop before this; the first blocks shed contacts "
                        "before rebuilding them (default 4000)")
    a.add_argument("--max-steps", type=int, default=30000,
                   help="cap, and the value handed to ridgerunner (default 30000)")
    a.add_argument("--minrad-eps", type=float, default=0.01,
                   help="relative improvement in the minRad median that counts as "
                        "progress (default 0.01)")
    a.add_argument("--plateau-on", default="auto",
                   help="which signals must be flat before --auto-steps stops a "
                        "descent: 'auto' (default), or a comma-separated subset "
                        "of minrad, struts, ropelength. All three are always "
                        "measured and printed; this chooses which ones gate. "
                        "'auto' keys on the round's own move: a squeeze at "
                        "f <= --plateau-hard-f is a repair round whose output "
                        "feeds another move, so it gates on minrad and struts "
                        "only; a gentler squeeze, a contraction and an initial "
                        "descent all gate on all three. Truncation propagates, "
                        "so run any round you will REPORT with all three.")
    a.add_argument("--plateau-hard-f", type=float, default=0.5,
                   help="squeeze factors at or below this count as repair rounds "
                        "for --plateau-on auto (default 0.5). f is the slope of "
                        "the radial map, so SMALLER f is more aggressive.")
    a.add_argument("--rop-eps", type=float, default=1e-4,
                   help="relative ropelength improvement per block that still "
                        "counts as progress (default 1e-4). A descent that is "
                        "still working moves ~4e-4 of its own value per 1000 "
                        "steps; a converged one moves under 1e-5. Stopping on "
                        "minRad and struts alone left 0.84 ropelength units on "
                        "the table in a measured case, because a repaired "
                        "thickness is not a converged ropelength.")
    a.add_argument("--strut-eps", type=float, default=0.02,
                   help="same for the strut median (default 0.02)")

    g = p.add_argument_group("contraction scan")
    g.add_argument("--factor-lo", type=float, default=0.80,
                   help="hardest contraction factor to try (default 0.80). For "
                        "a CONTRACTION the factor is how far clusters move "
                        "toward their common centre, so a lower number is the "
                        "more aggressive move.")
    g.add_argument("--factor-hi", type=float, default=0.99,
                   help="gentlest contraction factor to try (default 0.99). "
                        "1.0 would be no contraction at all.")
    g.add_argument("--factor-step", type=float, default=0.01,
                   help="the window narrows as rounds proceed; 0.01 resolves it "
                        "(default 0.01)")
    g.add_argument("--join", type=float, default=None,
                   help="contract_clusters.py --join; contact bands closer than this "
                        "(in D) are one cluster")
    g.add_argument("--min-gain", type=float, default=0.5,
                   help="stop when the best predicted gain falls below this "
                        "(default 0.5)")

    s = p.add_argument_group("symmetry")
    s.add_argument("--symmetry", default=None,
                   help="ridgerunner --Symmetry for the descent, e.g. Z/5Z")
    s.add_argument("--sym-group", default=None,
                   help="symmetrize_link_xyz.py --group for re-symmetrizing after each "
                        "contraction, e.g. Cn (there is no 'C5'; use Cn with "
                        "--sym-order 5)")
    s.add_argument("--sym-order", type=int, default=None,
                   help="the n in Cn, for our own symmetry tool -- normally the "
                        "number of loops. It must agree with the n in "
                        "--symmetry Z/nZ, which is RidgeRunner's own flag.")
    s.add_argument("--sym-axis", default="0,0,1",
                   help="rotation axis as x,y,z (default 0,0,1). Files from our "
                        "drawing scripts already put the axis on z. Note this "
                        "is also the axis a hole squeeze uses: on a symmetric "
                        "link the squeeze axis MUST be the symmetry axis.")
    s.add_argument("--no-permute", action="store_true",
                   help="require every element to carry each component onto itself. "
                        "WRONG when the rotation cycles components, which it does for "
                        "n components under Cn.")

    i = p.add_argument_group(
        "initial descent and flag selection",
        "The cycle starts with a contraction, and a contraction needs contact "
        "clusters to translate. A raw diagram has none: measured on a raw "
        "7-component layout (max d/D 3.99, 1 strut) every factor from 0.80 to 0.95 "
        "returned the input's own length, unchanged. Descend once first.")
    i.add_argument("--initial-steps", type=int, default=0,
                   help="descend the input this many steps before round 1 "
                        "(default 0 = assume the input is already tightened)")
    i.add_argument("--auto-flags", action="store_true",
                   help="choose --eq / --Timewarp per round from the measured "
                        "max d/D and residual, following the guide's Step 2 table, "
                        "and consult tighten_lib/strut_free.py before enabling "
                        "Timewarp")
    i.add_argument("--vu-ladder", default="",
                   help="comma-separated v/u rungs to climb AFTER each descent, "
                        "e.g. 2,4,8 -- the ridgerunner author's own resolution "
                        "ladder (\"2 vertices per unit of ropelength, then 4, then "
                        "8, with most of the runtime in the last stage\"). Empty = "
                        "never refine. Checked after a descent and never before "
                        "one: d/D = Rop/(2N) is an output, so it improves for free "
                        "as the curve shortens -- a raw layout at v/u 0.128 needs "
                        "7.8x the vertices to reach v/u 1.0 at its own inflated "
                        "ropelength, and gets there for nothing after one descent.")
    i.add_argument("--refine-mode", choices=("subdivide", "spline"),
                   default="subdivide",
                   help="subdivide returns a near-equilateral curve (edge max/min "
                        "1.01 vs splining's 1.11); a refined file is an input to a "
                        "run, and the edge distribution is what a run cares about")
    i.add_argument("--refine-fix-minrad", type=float, default=0.5,
                   help="MinRad target passed to refine_xyz.py --fix-minrad. "
                        "Refining divides minRad by the refinement ratio, so this "
                        "is not optional (default 0.5)")
    i.add_argument("--dd-threshold", type=float, default=0.5,
                   help="max d/D above which a curve is treated as chording across "
                        "its own tube (default 0.5)")

    mv = p.add_argument_group(
        "move selection",
        "Pick each round's geometric move from the measured structure instead of "
        "always contracting. Both moves close empty space; they differ in which "
        "space: the hole-squeeze closes the void ENCLOSED by the ring (topology-"
        "safe monotone radial map, symmetry preserved, damage is repairable "
        "proximity), cluster contraction closes space BETWEEN clusters (needs a "
        "HOMFLY gate, damage lands on curvature). Measured head to head on one "
        "structure stuck at 213.86: contraction's best possible outcome was 215.7 "
        "while the squeeze reached 207.4 after 2800 steps of repair.")
    mv.add_argument("--choose-move", action="store_true",
                    help="enable the selector; without it the cycle contracts, "
                         "as before")
    mv.add_argument("--min-hole", type=float, default=0.25,
                    help="hole wall radius (in D) below which the squeeze is "
                         "considered out of budget (default 0.25)")
    mv.add_argument("--squeeze-blend", type=float, default=0.5,
                    help="smoothstep width of the radial ramp (default 0.5)")
    mv.add_argument("--squeeze-margin", type=float, default=1.10,
                    help="required minRad/(minStrut/2) after the squeeze; keeps "
                         "curvature comfortably non-binding (default 1.10)")
    mv.add_argument("--find-axis", action="store_true",
                    help="search for the best squeeze axis each round instead of "
                         "assuming z through the centroid -- for IRREGULAR "
                         "structures. Sweeps directions, solves the largest-empty-"
                         "circle problem per projection, ranks the top tunnels by "
                         "the length a probe squeeze actually removes. On a "
                         "trapped asymmetric 4-component link this found a 2.95 D "
                         "tunnel whose f=0.20 squeeze removed 30%% of the length; "
                         "on Cn-symmetric structures it recovers the symmetry "
                         "axis by itself.")
    mv.add_argument("--axis-dirs", type=int, default=160,
                    help="directions swept by --find-axis (default 160)")
    mv.add_argument("--squeeze-factors", default="0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
                    help="factors swept by the squeeze scan; the HARDEST "
                         "admissible one is taken")

    o = p.add_argument_group("other")
    o.add_argument("--rr-arg", action="append", default=[],
                   help="extra ridgerunner flag, repeatable, without the leading dashes "
                        "handled for you: --rr-arg=--Timewarp. Check "
                        "tighten_lib/strut_free.py first: Timewarp is pure overhead once "
                        "the curve has no strut-free stretch.")
    o.add_argument("--min-vu", type=float, default=0.0,
                   help="stop and recommend refinement if vertices-per-unit-ropelength "
                        "falls below this. The author's own ladder is 2/4/8; set 6-8 to "
                        "catch a resolution-limited structure.")
    o.add_argument("--strict-topology", action="store_true",
                   help="treat an uncomputable HOMFLY as a stop. Off by default: "
                        "UNKNOWN means knottype could not build a polynomial, which is "
                        "not evidence of a change.")
    o.add_argument("--dry-run", action="store_true",
                   help="measure and scan round 1, then stop without modifying anything")
    return p


def main(argv: list[str] | None = None) -> int:
    # A full cycle is an hours-long unattended job and its output is almost
    # always redirected to a file. Without this, Python block-buffers stdout and
    # the log stays empty for hours, so there is no way to watch progress.
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, OSError):
        pass
    parser = build_parser()
    raw = sys.argv[1:] if argv is None else list(argv)
    # No arguments at all, or --gui: this is a launcher, not a run. The window
    # assembles a command line and starts it detached, so it never needs the
    # validation below -- which would reject an empty argv outright.
    if not raw or "--gui" in raw:
        from tighten_gui import launch
        return launch(parser)

    args = parser.parse_args(argv)
    if not args.input and not args.resume:
        print("error: give an input .xyz, or --resume to continue a cycle "
              "already in --work-dir", file=sys.stderr)
        return 2
    if args.input and args.resume:
        print("error: --resume takes its input from --work-dir; drop the input "
              "argument", file=sys.stderr)
        return 2
    try:
        args.squeeze_factors = sorted(float(v) for v in
                                      args.squeeze_factors.split(",") if v.strip())
    except ValueError:
        print("error: --squeeze-factors must be a comma-separated list of numbers",
              file=sys.stderr)
        return 2
    if args.choose_move and not LIB.joinpath("radial_squeeze.py").is_file():
        print(f"error: --choose-move needs {LIB}/radial_squeeze.py", file=sys.stderr)
        return 2
    try:
        args.vu_ladder = [float(v) for v in args.vu_ladder.split(",") if v.strip()]
    except ValueError:
        print("error: --vu-ladder must be a comma-separated list of numbers",
              file=sys.stderr)
        return 2
    if args.vu_ladder and not LIB.joinpath("refine_xyz.py").is_file():
        print(f"error: --vu-ladder needs {LIB}/refine_xyz.py", file=sys.stderr)
        return 2
    if args.plateau_on.strip().lower() != "auto":
        valid = {"minrad", "struts", "ropelength"}
        args.plateau_on = frozenset(v.strip().lower()
                                    for v in args.plateau_on.split(",") if v.strip())
        if not args.plateau_on or not args.plateau_on <= valid:
            print(f"error: --plateau-on must be 'auto' or a comma-separated "
                  f"subset of {sorted(valid)}", file=sys.stderr)
            return 2
    else:
        args.plateau_on = "auto"
    if args.sym_group and not args.sym_order:
        print("error: --sym-group needs --sym-order", file=sys.stderr)
        return 2
    for binary in ("ropelength", "struts", "residual"):
        if shutil.which(binary) is None:
            print(f"error: {binary} not on PATH (octrope/ridgerunner not installed?)",
                  file=sys.stderr)
            return 2
    if not LIB.joinpath("contract_clusters.py").is_file():
        print(f"error: {LIB}/contract_clusters.py not found", file=sys.stderr)
        return 2
    if args.gui:                      # reachable only via an explicit call
        from tighten_gui import launch
        return launch(parser)
    try:
        return run_cycle(args)
    except CycleError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
