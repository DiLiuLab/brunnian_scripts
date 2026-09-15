#!/usr/bin/env python3
"""Tighten a knot or link stored as an .xyz file using RidgeRunner.

RidgeRunner minimizes the ropelength of a polygonal curve by constrained
gradient descent, which is what turns a hand-drawn or layout-generated link
into an "ideal" (tight) configuration. It speaks Geomview VECT rather than
xyz, so this script does the round trip:

    input.xyz -> input.vect -> ridgerunner -> a saved .vect -> output.xyz

The .xyz convention used here is the one written by the BL drawing scripts:
one "x y z" per line, with a blank line between components. Every component
is treated as a closed loop, which is what VECT encodes with a negative
vertex count.

RidgeRunner resamples the curve while it runs, so the tightened output
usually has a different number of vertices per component than the input.

The last step is not always the best one. A run can reach a nearly
ropelength-critical configuration and then lose it, usually when the
constraint matrix degenerates and the strut set collapses. So by default
this script watches for that failure, stops the run when it happens, and
converts the best saved configuration rather than the final one, saying
which it used (--select final keeps the last one instead).

Install RidgeRunner with:

    brew tap designbynumbers/cantarellalab
    brew install ridgerunner

or build libtsnnls, libplcurve and ridgerunner from source; see
https://www.jasoncantarella.com/wordpress/software/ridgerunner/
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path
from typing import Callable

Point = tuple[float, float, float]
Component = list[Point]

# Directories to search for the ridgerunner binary when it is not on PATH.
RIDGERUNNER_FALLBACK_DIRS = (
    "~/.local/bin",
    "/opt/homebrew/bin",
    "/usr/local/bin",
)

# One RGBA color per component, reused cyclically. These only affect how the
# intermediate VECT files look in Geomview; they are dropped on the way back
# out to xyz.
COMPONENT_COLORS = (
    (0.90, 0.20, 0.20, 1.0),
    (0.20, 0.45, 0.90, 1.0),
    (0.15, 0.70, 0.30, 1.0),
    (0.95, 0.70, 0.10, 1.0),
    (0.60, 0.30, 0.80, 1.0),
    (0.10, 0.75, 0.80, 1.0),
)

# Distance below which the last vertex is treated as a repeat of the first.
CLOSING_VERTEX_TOL = 1e-9


class TightenError(RuntimeError):
    """Raised for malformed input files or a failed ridgerunner run."""


class RunCancelled(TightenError):
    """Raised when ridgerunner is killed by a signal, e.g. the GUI Stop button.

    RidgeRunner has no signal handler of its own, so a terminated run dies
    where it stands. It rewrites <stem>.final.vect periodically rather than
    every step, so a run stopped early may have saved a configuration or may
    have saved nothing at all.
    """

    def __init__(self, message: str, captured: list[str] | None = None) -> None:
        super().__init__(message)
        self.captured = captured or []
        # Set to a short reason when the run was stopped because it
        # degenerated, rather than by the user.
        self.degenerate: str | None = None


def read_xyz(path: Path) -> list[Component]:
    """Read blank-line-separated components of "x y z" rows.

    Also tolerates a leading element/label column ("C 1.0 2.0 3.0"), which is
    what a standard molecular .xyz file looks like, and skips the two-line
    atom-count header such files carry.
    """
    components: list[Component] = []
    current: Component = []

    lines = list(enumerate(path.read_text().splitlines(), start=1))

    # A molecular .xyz file opens with an atom count and a free-form comment
    # line. Drop that pair so both conventions can be read.
    first = next((index for index, (_, raw) in enumerate(lines) if raw.strip()), None)
    if first is not None:
        fields = lines[first][1].split()
        if len(fields) == 1 and fields[0].lstrip("+-").isdigit():
            del lines[first : first + 2]

    for lineno, raw in lines:
        line = raw.split("#", 1)[0].strip()
        if not line:
            # Blank line closes the current component.
            if current:
                components.append(current)
                current = []
            continue

        fields = line.split()
        if len(fields) >= 4:
            # Drop a leading label column, e.g. "C x y z".
            fields = fields[1:4]
        if len(fields) != 3:
            raise TightenError(
                f"{path}:{lineno}: expected 3 coordinates, got {len(fields)}: {raw!r}"
            )

        try:
            current.append((float(fields[0]), float(fields[1]), float(fields[2])))
        except ValueError as exc:
            raise TightenError(f"{path}:{lineno}: {exc}") from exc

    if current:
        components.append(current)

    if not components:
        raise TightenError(f"{path}: no coordinates found")

    return [drop_repeated_closing_vertex(comp) for comp in components]


def drop_repeated_closing_vertex(component: Component) -> Component:
    """Remove a final vertex that just repeats the first one.

    VECT marks a component as closed with a negative vertex count and must
    not also repeat the starting vertex, otherwise ridgerunner sees a
    zero-length edge.
    """
    if len(component) > 1 and math.dist(component[0], component[-1]) < CLOSING_VERTEX_TOL:
        return component[:-1]
    return component


def write_vect(path: Path, components: list[Component]) -> None:
    """Write components as a Geomview VECT file of closed polylines."""
    total_vertices = sum(len(comp) for comp in components)
    counts = " ".join(str(-len(comp)) for comp in components)  # negative == closed
    color_counts = " ".join("1" for _ in components)

    lines = [
        "VECT",
        f"{len(components)} {total_vertices} {len(components)}",
        counts,
        color_counts,
        "",
    ]
    for comp in components:
        lines.extend(f"{x:.17g} {y:.17g} {z:.17g}" for x, y, z in comp)
        lines.append("")
    for index in range(len(components)):
        r, g, b, a = COMPONENT_COLORS[index % len(COMPONENT_COLORS)]
        lines.append(f"{r} {g} {b} {a}")

    path.write_text("\n".join(lines) + "\n")


def read_vect(path: Path) -> list[Component]:
    """Read the polyline vertices out of a Geomview VECT file.

    VECT is whitespace-delimited rather than line-oriented, so this walks a
    flat token stream: the VECT keyword, three counts, per-polyline vertex
    counts, per-polyline color counts, then the vertices themselves.
    """
    text = path.read_text()
    tokens: list[str] = []
    for line in text.splitlines():
        tokens.extend(line.split("#", 1)[0].split())

    if not tokens:
        raise TightenError(f"{path}: empty VECT file")

    position = 0

    def take(count: int) -> list[str]:
        nonlocal position
        chunk = tokens[position : position + count]
        if len(chunk) != count:
            raise TightenError(f"{path}: VECT file ended early")
        position += count
        return chunk

    keyword = take(1)[0]
    if keyword.upper() not in {"VECT", "4VECT", "NVECT"}:
        raise TightenError(f"{path}: expected a VECT header, got {keyword!r}")

    n_polylines, n_vertices, n_colors = (int(value) for value in take(3))
    vertex_counts = [abs(int(value)) for value in take(n_polylines)]
    take(n_polylines)  # per-polyline color counts, not needed here

    if sum(vertex_counts) != n_vertices:
        raise TightenError(
            f"{path}: vertex counts {sum(vertex_counts)} disagree with header {n_vertices}"
        )

    flat = [float(value) for value in take(3 * n_vertices)]

    components: list[Component] = []
    offset = 0
    for count in vertex_counts:
        comp: Component = []
        for index in range(count):
            base = 3 * (offset + index)
            comp.append((flat[base], flat[base + 1], flat[base + 2]))
        components.append(comp)
        offset += count

    del n_colors  # colors trail the vertices; xyz output does not carry them
    return components


def write_xyz(path: Path, components: list[Component], decimals: int) -> None:
    """Write components back out in the blank-line-separated xyz convention."""
    blocks = []
    for comp in components:
        blocks.append(
            "\n".join(
                f"{x:.{decimals}f} {y:.{decimals}f} {z:.{decimals}f}" for x, y, z in comp
            )
        )
    path.write_text("\n\n".join(blocks) + "\n")


def find_ridgerunner(explicit: str | None) -> str:
    """Locate the ridgerunner binary, or explain how to install it."""
    # A binary the user named explicitly is either used or reported. Falling
    # back to a different one would silently run something they did not ask for.
    for value, source in (
        (explicit, "--ridgerunner"),
        (os.environ.get("RIDGERUNNER"), "$RIDGERUNNER"),
    ):
        if not value:
            continue
        resolved = shutil.which(value)
        if resolved:
            return resolved
        raise TightenError(f"{source} is set to {value!r}, which is not an executable.")

    resolved = shutil.which("ridgerunner")
    if resolved:
        return resolved

    for directory in RIDGERUNNER_FALLBACK_DIRS:
        candidate = Path(directory).expanduser() / "ridgerunner"
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)

    raise TightenError(
        "ridgerunner not found. Install it with:\n"
        "    brew tap designbynumbers/cantarellalab\n"
        "    brew install ridgerunner\n"
        "or pass --ridgerunner /path/to/ridgerunner."
    )


def build_ridgerunner_command(binary: str, vect_name: str, args) -> list[str]:
    """Assemble the ridgerunner argument list from the parsed options."""
    command = [binary]

    if args.autoscale:
        # Scale to thickness 0.501 first so the run rebuilds its own set of
        # self-contacts. Recommended even for already-tight input.
        command.append("--Autoscale")
    if args.eq:
        command.append("--EqOn")
    if args.steps is not None:
        command.append(f"--StopSteps={args.steps}")
    if args.stop_res is not None:
        command.append(f"--StopResidual={args.stop_res}")
    if args.stop_time is not None:
        command.append(f"--StopTime={args.stop_time}")
    if args.stop20 is not None:
        command.append(f"--Stop20={args.stop20}")
    if args.tube_radius is not None:
        command.append(f"--TubeRadius={args.tube_radius}")
    if args.stiffness is not None:
        command.append(f"--Lambda={args.stiffness}")
    if args.symmetry:
        command.append(f"--Symmetry={args.symmetry}")
    if args.snapshot_interval is not None:
        command.append(f"--SnapshotInterval={args.snapshot_interval}")
    if not args.keep_snapshots:
        # Skip the vectfiles/ movie frames and the (large) lsqr log.
        command.extend(["--NoOutputFiles", "--NoLsqrLog"])
    if not args.png:
        # The final-snapshot render needs POVRAY, which is usually absent.
        command.append("--NoPNGOutput")

    # Split so one --rr-arg can carry several flags; argparse needs the '='
    # form anyway, since a bare '--Foo' value looks like an option to it.
    for extra in args.rr_arg:
        command.extend(extra.split())
    command.append(vect_name)
    return command


# ridgerunner writes this to its own logfile (via logprintf, never to stdout)
# when the constraint matrix becomes too degenerate to solve. The run does not
# recover: the strut set collapses and the remaining steps are wasted.
DEGENERATE_MARKER = "Linear algebra failure"


# A run that sheds most of its self-contacts has come apart and does not get
# them back. The raw per-step count is useless for detecting this: it spikes to
# near zero constantly even in healthy runs, so a rule on raw values fires
# within the first hundred steps of everything. These thresholds work on the
# median of each COLLAPSE_BIN-step block instead.
#
# Calibrated against seven completed runs (three that collapsed, four that did
# not). Only two settings separated them at all, so the margin is thin and this
# check is opt-in via --stop-on-collapse rather than on by default: a false
# positive truncates a healthy run, whereas a missed collapse only wastes time
# that --select best already recovers from.
COLLAPSE_BIN = 500
COLLAPSE_FRACTION = 0.2
COLLAPSE_MIN_PEAK = 50
COLLAPSE_CONFIRM_BINS = 10


def _tail_new_lines(path: Path, state: dict) -> list[str]:
    """Return complete lines added to path since the last call."""
    try:
        if not path.is_file():
            return []
        with path.open("r", errors="ignore") as handle:
            handle.seek(state.get("offset", 0))
            chunk = handle.read()
            state["offset"] = handle.tell()
    except OSError:
        return []

    text = state.get("pending", "") + chunk
    lines = text.split("\n")
    state["pending"] = lines.pop()  # the last piece may be half-written
    return lines


def watch_for_degeneracy(
    run_dir: Path,
    stem: str,
    process: subprocess.Popen,
    flag: dict,
    stop: "threading.Event",
    watch_collapse: bool = False,
) -> None:
    """Stop the run when it degenerates, by either of two signals.

    ridgerunner logs a linear-algebra failure when its constraint matrix
    becomes unsolvable. A run can also lose its contact set without logging
    anything at all; watching the strut count catches that silent case, but
    only approximately, so it is opt-in.
    """
    log_state: dict = {}
    strut_state: dict = {}
    log_path = run_dir / f"{stem}.log"
    strut_path = run_dir / "logfiles" / "strutcount.dat"

    current_bin = None
    bin_values: list[float] = []
    peak = 0.0
    below = 0

    def collapsed(median: float) -> bool:
        """Update the running peak and say whether the run has come apart."""
        nonlocal peak, below
        peak = max(peak, median)
        if peak < COLLAPSE_MIN_PEAK:
            return False
        if median < COLLAPSE_FRACTION * peak:
            below += 1
            return below >= COLLAPSE_CONFIRM_BINS
        below = 0
        return False

    while not stop.wait(2.0):
        for line in _tail_new_lines(log_path, log_state):
            if DEGENERATE_MARKER in line:
                flag["reason"] = "ridgerunner reported a linear-algebra failure"
                if process.poll() is None:
                    process.terminate()
                return

        if not watch_collapse:
            continue

        for line in _tail_new_lines(strut_path, strut_state):
            fields = line.split()
            if len(fields) < 2:
                continue
            try:
                step, struts = int(fields[0]), float(fields[1])
            except ValueError:
                continue

            index = step // COLLAPSE_BIN
            if current_bin is None:
                current_bin = index
            if index != current_bin and bin_values:
                bin_values.sort()
                median = bin_values[len(bin_values) // 2]
                if collapsed(median):
                    flag["reason"] = (
                        f"the contact set collapsed: {median:.0f} struts against a peak "
                        f"of {peak:.0f}, sustained over "
                        f"{COLLAPSE_CONFIRM_BINS * COLLAPSE_BIN} steps"
                    )
                    if process.poll() is None:
                        process.terminate()
                    return
                bin_values = []
                current_bin = index
            bin_values.append(struts)


def run_ridgerunner(
    command: list[str],
    work_dir: Path,
    quiet: bool,
    on_line: Callable[[str], None] | None = None,
    on_process: Callable[[subprocess.Popen], None] | None = None,
    degenerate_watch: tuple[Path, str] | None = None,
    watch_collapse: bool = False,
) -> list[str]:
    """Run ridgerunner in work_dir, streaming its progress, and return output.

    on_line receives each output line as it arrives, which is how the GUI
    follows a run. on_process hands the caller the live Popen object so a
    long run can be cancelled.
    """
    captured: list[str] = []
    process = subprocess.Popen(
        command,
        cwd=work_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    if on_process is not None:
        on_process(process)

    flag: dict = {"reason": None}
    stop = threading.Event()
    watcher = None
    if degenerate_watch is not None:
        watcher = threading.Thread(
            target=watch_for_degeneracy,
            args=(*degenerate_watch, process, flag, stop, watch_collapse),
            daemon=True,
        )
        watcher.start()

    # Tee ridgerunner's stdout to a file in the run directory. It was
    # previously kept only in `captured`, which is surfaced on FAILURE and
    # parsed for ropelength, but never written anywhere on a successful quiet
    # run -- so every diagnostic it prints was silently lost. That is not
    # cosmetic: the --Symmetry warnings ("turning off EqOn", and the eq FORCE
    # it turns ON in exchange), the "Built <group> ... Error before
    # symmetrizing X, after Y" line, and the half-edge-threshold WARNING all go
    # to stdout, not into ridgerunner's own .rr log. Losing them is why a run
    # could silently use a different stepper and a different eq regime than the
    # command line implied, and nobody noticed for five runs.
    stdout_log = None
    try:
        stdout_log = (Path(work_dir) / "ridgerunner_stdout.log").open("w")
    except Exception:
        stdout_log = None
    assert process.stdout is not None
    for line in process.stdout:
        captured.append(line.rstrip("\n"))
        if stdout_log is not None:
            stdout_log.write(line)
            stdout_log.flush()
        if on_line is not None:
            on_line(line.rstrip("\n"))
        elif not quiet:
            sys.stderr.write(line)
    if stdout_log is not None:
        stdout_log.close()
    returncode = process.wait()
    stop.set()
    if watcher is not None:
        watcher.join(timeout=5)

    if returncode < 0:
        cancelled = RunCancelled(
            f"ridgerunner was stopped by signal {-returncode}.", captured
        )
        cancelled.degenerate = flag["reason"]
        raise cancelled
    if returncode != 0:
        tail = "\n".join(captured[-25:])
        raise TightenError(
            f"ridgerunner exited with status {returncode}. Last output:\n{tail}"
        )
    return captured


def latest_saved_vect(run_dir: Path, stem: str) -> Path | None:
    """Find the newest configuration ridgerunner has written to disk.

    Despite what the documentation says, <stem>.final.vect is only written
    when a run finishes, so a run that was stopped early has to be salvaged
    from the periodic snapshots/<stem>.<step>.vect files instead (written
    every --SnapshotInterval steps, 10000 by default).
    """
    final = run_dir / f"{stem}.final.vect"
    if final.is_file():
        return final

    snapshots: list[tuple[int, Path]] = []
    prefix = f"{stem}."
    for path in (run_dir / "snapshots").glob(f"{stem}.*.vect"):
        # Keep <stem>.<step>.vect; skip the .dlen/.dVdt/.struts companions.
        middle = path.name[len(prefix) : -len(".vect")]
        if middle.isdigit():
            snapshots.append((int(middle), path))

    return max(snapshots)[1] if snapshots else None


# Steps either side of a snapshot averaged when scoring it. Residual and strut
# count swing wildly from step to step, so a single step's value is noise.
SELECTION_WINDOW = 500

# An earlier snapshot has to beat the final residual by more than this to be
# worth rolling back to; otherwise the run keeps its last configuration.
ROLLBACK_MARGIN = 0.8

# Ropelength is the quantity being minimised, so a rollback that gives it up is
# only worth making when the earlier configuration is better in every other way
# too. Outside that case, refuse to trade away more than this fraction of it.
# Measured cost of not having this guard: three runs in one afternoon whose
# written output was 3 to 11 units worse than the final configuration they
# discarded, because the rollback was decided on residual alone.
ROLLBACK_ROPE_TOLERANCE = 0.02


def read_metric_log(run_dir: Path, name: str, column: int = 1) -> dict[int, float]:
    """Read logfiles/<name>.dat into {step: value}."""
    path = run_dir / "logfiles" / f"{name}.dat"
    if not path.is_file():
        return {}

    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    if lines and not text.endswith("\n"):
        lines.pop()  # a live run can leave the final line half-written

    values: dict[int, float] = {}
    for line in lines:
        fields = line.split()
        if len(fields) <= column:
            continue
        try:
            values[int(fields[0])] = float(fields[column])
        except ValueError:
            continue
    return values


def nearest_logged(series: dict[int, float], step: int) -> float | None:
    """Value at the logged step closest to `step`.

    RidgeRunner decimates its logfiles as a run grows -- a long run logs only
    every other step, or coarser -- so a snapshot's step number is usually
    absent from them. Looking it up exactly and falling back to the last
    logged step would describe the wrong configuration entirely.
    """
    if not series:
        return None
    return series[min(series, key=lambda logged: abs(logged - step))]


def window_median(series: dict[int, float], centre: int, window: int) -> float | None:
    """Median of a logged quantity within +/- window steps of centre."""
    nearby = [v for step, v in series.items() if abs(step - centre) <= window]
    if not nearby:
        return None
    nearby.sort()
    middle = len(nearby) // 2
    if len(nearby) % 2:
        return nearby[middle]
    return 0.5 * (nearby[middle - 1] + nearby[middle])


def configuration_step(path: Path, stem: str) -> int | None:
    """The step a saved configuration belongs to; None for <stem>.final.vect."""
    prefix = f"{stem}."
    if not path.name.startswith(prefix) or not path.name.endswith(".vect"):
        return None
    middle = path.name[len(prefix) : -len(".vect")]
    return int(middle) if middle.isdigit() else None


def snapshot_steps(run_dir: Path, stem: str) -> list[tuple[int, Path]]:
    """List saved configurations as (step, path), excluding the step-0 input."""
    found: list[tuple[int, Path]] = []
    prefix = f"{stem}."
    for path in (run_dir / "snapshots").glob(f"{stem}.*.vect"):
        middle = path.name[len(prefix) : -len(".vect")]
        if middle.isdigit() and int(middle) > 0:
            found.append((int(middle), path))
    return sorted(found)


def choose_configuration(run_dir: Path, stem: str) -> tuple[Path, str] | None:
    """Pick the saved configuration that is closest to ropelength-critical.

    A run without --EqOn can reach a very good state and then degenerate as
    edges collapse, so the last configuration is not always the best one.

    Three quantities matter and they do not always agree: ropelength (the thing
    being minimised), windowed median residual (ridgerunner's own convergence
    measure) and strut count (how much of the curve is actually in contact).
    Ranking on residual alone gets both directions wrong, as observed:

      * rollbacks fired onto snapshots 3.1 to 11.3 ropelength units worse than
        the final they discarded, across four runs, purely because those
        snapshots had a lower residual;
      * and a configuration that beat the final on ropelength (217.19 against
        243.71), residual (0.400 against 0.424) and struts (119 against 88) was
        never even weighed, because the argmin over residual landed on a
        different snapshot entirely (step 192000, ropelength 220.94), whose own
        ratio of 0.934 then failed the margin. Being better on every measure
        bought nothing, because nothing looked.

    So dominance decides first: a candidate better on all three is taken, and
    one worse on all three is refused, both without consulting the margin. Only
    genuine trade-offs fall through to the residual test, and even then a
    rollback may not give up more than ROLLBACK_ROPE_TOLERANCE of ropelength.
    """
    residual = read_metric_log(run_dir, "residual")
    ropelength = read_metric_log(run_dir, "ropelength")
    struts = read_metric_log(run_dir, "strutcount")
    if not residual:
        return None

    final_step = max(residual)
    candidates = snapshot_steps(run_dir, stem)
    final_vect = run_dir / f"{stem}.final.vect"
    if final_vect.is_file():
        candidates.append((final_step, final_vect))
    if not candidates:
        return None

    def metrics(step: int) -> tuple[float, float, float]:
        """(ropelength, residual, struts) - first two lower is better, third higher."""
        return (
            nearest_logged(ropelength, step) or float("inf"),
            window_median(residual, step, SELECTION_WINDOW) or 1.0,
            window_median(struts, step, SELECTION_WINDOW) or 0.0,
        )

    def dominates(a: tuple[float, float, float], b: tuple[float, float, float]) -> bool:
        """`a` is no worse than `b` on all three, and strictly better on one."""
        no_worse = a[0] <= b[0] and a[1] <= b[1] and a[2] >= b[2]
        strictly = a[0] < b[0] or a[1] < b[1] or a[2] > b[2]
        return no_worse and strictly

    def score(step: int) -> tuple[float, float]:
        rope, res, _ = metrics(step)
        return (res, rope)

    best_step, best_path = min(candidates, key=lambda item: score(item[0]))
    final_present = final_vect.is_file()

    def describe_step(step: int) -> str:
        rope, res, active = metrics(step)
        parts = [f"step {step}"]
        if rope != float("inf"):
            parts.append(f"ropelength {rope:.4f}")
        parts.append(f"residual {res:.4f}")
        parts.append(f"struts {active:.0f}")
        return ", ".join(parts)

    if not final_present:
        # A killed or crashed run has no final.vect, so there is no incumbent to
        # test dominance against. Rank the snapshots among themselves instead of
        # falling back to the residual argmin this rule exists to replace: keep
        # only those no other snapshot dominates, and take the shortest of them.
        scored = [(step, path, metrics(step)) for step, path in candidates]
        front = [(m[0], step, path) for step, path, m in scored
                 if not any(dominates(other, m) for _, _, other in scored)]
        if front:
            _, step, path = min(front)
            return path, (
                f"Run did not finish. Of the snapshots no other snapshot beats outright, "
                f"using the shortest: {describe_step(step)}."
            )
        return best_path, (
            f"Run did not finish. Using the best snapshot: {describe_step(best_step)}."
        )

    if best_path == final_vect:
        return final_vect, f"Best configuration is the final one ({describe_step(best_step)})."

    final_metrics = metrics(final_step)

    # A snapshot that beats the final on every measure is taken outright. This
    # is the degenerate-run case the rollback exists for: the strut set
    # collapses, residual goes to 1 and ropelength drifts back up, so the
    # rescuing snapshot wins on all three at once and needs no margin test.
    dominating = [
        (metrics(step)[0], step, path)
        for step, path in candidates
        if path != final_vect and dominates(metrics(step), final_metrics)
    ]
    if dominating:
        rope, step, path = min(dominating)
        return path, (
            f"Rolled back to a configuration that is better in every way: "
            f"{describe_step(step)}.\n"
            f"  The run ended worse on ropelength, residual and struts, at "
            f"{describe_step(final_step)}. "
            "Pass --select final to keep the last configuration instead."
        )

    best_metrics = metrics(best_step)
    if dominates(final_metrics, best_metrics):
        return final_vect, (
            f"Keeping the final configuration ({describe_step(final_step)}): it is "
            f"better in every way than the best snapshot ({describe_step(best_step)})."
        )

    # Neither dominates, so this is a real trade-off. Never buy a lower residual
    # with a materially worse ropelength -- ropelength is the objective.
    if best_metrics[0] > final_metrics[0] * (1.0 + ROLLBACK_ROPE_TOLERANCE):
        return final_vect, (
            f"Keeping the final configuration ({describe_step(final_step)}). "
            f"The best snapshot by residual (step {best_step}, "
            f"ropelength {best_metrics[0]:.4f}) would give up "
            f"{100 * (best_metrics[0] / final_metrics[0] - 1):.1f}% of ropelength "
            "to get there."
        )

    # Residual jitters, so only roll back for a clear improvement. Giving up
    # thousands of steps of progress to chase noise would be a bad trade.
    if best_metrics[1] > ROLLBACK_MARGIN * final_metrics[1]:
        return final_vect, (
            f"Keeping the final configuration ({describe_step(final_step)}). "
            f"The best snapshot (step {best_step}, residual {best_metrics[1]:.4f}) is not "
            "meaningfully better."
        )

    note = (
        f"Rolled back to a better earlier configuration: {describe_step(best_step)}.\n"
        f"  The run ended worse, at {describe_step(final_step)}. "
        "Pass --select final to keep the last configuration instead."
    )
    return best_path, note


def prune_snapshot_matrices(run_dir: Path) -> int:
    """Delete the per-snapshot linear-algebra dumps, which carry no geometry.

    Each snapshot writes its own copy of the constraint matrix, tens of MB at
    a time, and that - not the curves - is what makes run directories huge.
    """
    removed = 0
    snapshots = run_dir / "snapshots"
    if not snapshots.is_dir():
        return 0
    for path in snapshots.iterdir():
        if path.suffix in {".dat", ".mat", ".sparse"} and ".A." in path.name:
            try:
                path.unlink()
                removed += 1
            except OSError:
                pass
    return removed


def ropelength_from_output(lines: list[str]) -> tuple[float, float] | None:
    """Fall back to the "Rop:" values ridgerunner prints for each step."""
    values: list[float] = []
    for line in lines:
        fields = line.split()
        if "Rop:" in fields:
            index = fields.index("Rop:")
            if index + 1 < len(fields):
                try:
                    values.append(float(fields[index + 1]))
                except ValueError:
                    continue
    if not values:
        return None
    return values[0], values[-1]


def info(message: str) -> None:
    """Print progress unbuffered so it interleaves with ridgerunner's output."""
    print(message, flush=True)


def describe(components: list[Component]) -> str:
    counts = "/".join(str(len(comp)) for comp in components)
    return f"{len(components)} component(s), {sum(map(len, components))} vertices ({counts})"


def tighten(
    args,
    emit: Callable[[str], None] | None = None,
    on_process: Callable[[subprocess.Popen], None] | None = None,
) -> int:
    say = emit if emit is not None else info
    binary = find_ridgerunner(args.ridgerunner)
    input_path: Path = args.input.expanduser()
    if not input_path.is_file():
        raise TightenError(f"{input_path}: no such file")
    components = read_xyz(input_path)
    say(f"Read {input_path}: {describe(components)}")

    output_path: Path = (
        args.output.expanduser()
        if args.output
        else input_path.with_name(f"{input_path.stem}_tight.xyz")
    )

    # ridgerunner writes its results into <stem>.rr/ beside the input VECT, so
    # give it a dedicated directory instead of littering the data folder.
    work_dir: Path = (
        args.work_dir.expanduser()
        if args.work_dir
        else output_path.parent / f"{input_path.stem}_rr"
    )
    work_dir.mkdir(parents=True, exist_ok=True)

    stem = input_path.stem
    vect_path = work_dir / f"{stem}.vect"
    write_vect(vect_path, components)
    say(f"Wrote {vect_path}")

    command = build_ridgerunner_command(binary, vect_path.name, args)
    say("Running: " + " ".join(command))
    cancelled = False
    try:
        captured = run_ridgerunner(
            command,
            work_dir,
            args.quiet,
            on_line=emit,
            on_process=on_process,
            degenerate_watch=(
                (work_dir / f"{stem}.rr", stem) if args.stop_on_failure else None
            ),
            watch_collapse=args.stop_on_collapse,
        )
    except RunCancelled as exc:
        # A stopped run is still worth salvaging: convert whatever
        # configuration ridgerunner last wrote out.
        cancelled = True
        captured = exc.captured
        if exc.degenerate:
            say(
                f"Stopped early: {exc.degenerate}. A run does not recover from this, so "
                "the remaining steps would be wasted. Falling back to the best saved "
                "configuration."
            )
        else:
            say(str(exc))

    run_dir = work_dir / f"{stem}.rr"

    saved_vect = None
    if args.select in {"best", "both"}:
        chosen = choose_configuration(run_dir, stem)
        if chosen is not None:
            saved_vect, note = chosen
            say(note)

    if saved_vect is None:
        saved_vect = latest_saved_vect(run_dir, stem)
        if saved_vect is not None and saved_vect.name != f"{stem}.final.vect":
            say(f"No final.vect (the run did not finish); using snapshot {saved_vect.name}.")

    if saved_vect is None:
        if cancelled:
            raise TightenError(
                "Stopped before ridgerunner saved any configuration. Snapshots are only "
                "written every --SnapshotInterval steps (10000 by default), so there is "
                f"nothing to convert yet. Run logs are in {run_dir}."
            )
        raise TightenError(
            f"ridgerunner finished but no {stem}.final.vect was written; inspect {run_dir}"
        )

    tightened = read_vect(saved_vect)
    write_xyz(output_path, tightened, args.decimals)

    # Report the ropelength of the configuration actually written out, which
    # after a rollback is not the last step's.
    ropelength = None
    rope_log = read_metric_log(run_dir, "ropelength")
    if rope_log:
        chosen_step = configuration_step(saved_vect, stem)
        if chosen_step is None:
            chosen_step = max(rope_log)
        end = nearest_logged(rope_log, chosen_step)
        ropelength = (rope_log[min(rope_log)], end) if end is not None else None
    if ropelength is None:
        ropelength = ropelength_from_output(captured)

    if ropelength:
        start, end = ropelength
        change = f"{100.0 * (end - start) / start:+.2f}%" if start else "n/a"
        say(f"Ropelength: {start:.6f} -> {end:.6f} ({change})")
    else:
        say("Ropelength: not reported by ridgerunner")

    say(f"Tightened: {describe(tightened)}")
    say(f"Wrote {output_path}")

    if args.select == "both":
        final_vect = run_dir / f"{stem}.final.vect"
        if not final_vect.is_file():
            say("No final.vect to write alongside: the run did not finish.")
        elif final_vect == saved_vect:
            say("Not writing a separate _final file: the best configuration is the final one.")
        else:
            final_path = output_path.with_name(f"{output_path.stem}_final{output_path.suffix}")
            final_components = read_vect(final_vect)
            write_xyz(final_path, final_components, args.decimals)
            last = max(rope_log) if rope_log else None
            detail = f" (ropelength {rope_log[last]:.4f})" if last is not None else ""
            say(f"Wrote {final_path}: the run's last configuration{detail}, "
                f"{describe(final_components)}")
    if cancelled:
        say("The run was stopped early, so this is a saved configuration rather than a converged one.")

    if cancelled or args.keep_run_dir or args.work_dir or args.keep_snapshots:
        if not args.keep_snapshots:
            pruned = prune_snapshot_matrices(run_dir)
            if pruned:
                say(f"Pruned {pruned} snapshot matrix dumps (no geometry lost).")
        say(f"RidgeRunner run directory kept at {run_dir}")
    else:
        shutil.rmtree(work_dir, ignore_errors=True)

    return 0


def gui_defaults() -> argparse.Namespace:
    """Parse an empty CLI so the GUI starts from the documented defaults."""
    return build_parser().parse_args(["__gui_placeholder__"])


def run_gui() -> None:
    """Launch a Tkinter GUI for setting the tightening parameters."""
    try:
        import queue
        import threading
        import tkinter as tk
        from tkinter import filedialog, messagebox
        from tkinter.scrolledtext import ScrolledText
    except Exception as exc:  # pragma: no cover - only used when Tkinter is absent.
        raise RuntimeError(f"Tkinter GUI is not available in this Python environment: {exc}")

    defaults = gui_defaults()

    root = tk.Tk()
    root.title("tighten_link_xyz: tighten a knot or link with RidgeRunner")
    root.geometry("1000x780")

    # --- parameter variables, seeded from the CLI defaults ------------------
    fields = {
        "input": tk.StringVar(),
        "output": tk.StringVar(),
        "work_dir": tk.StringVar(),
        "steps": tk.StringVar(value=str(defaults.steps)),
        "stop_res": tk.StringVar(value=str(defaults.stop_res)),
        "stop_time": tk.StringVar(),
        "stop20": tk.StringVar(),
        "tube_radius": tk.StringVar(),
        "stiffness": tk.StringVar(),
        "symmetry": tk.StringVar(),
        "decimals": tk.StringVar(value=str(defaults.decimals)),
        "snapshot_interval": tk.StringVar(),
        "rr_arg": tk.StringVar(),
        "ridgerunner": tk.StringVar(),
    }
    flags = {
        "autoscale": tk.BooleanVar(value=defaults.autoscale),
        "eq": tk.BooleanVar(value=defaults.eq),
        "keep_run_dir": tk.BooleanVar(value=defaults.keep_run_dir),
        "keep_snapshots": tk.BooleanVar(value=defaults.keep_snapshots),
        "png": tk.BooleanVar(value=defaults.png),
        "select_best": tk.BooleanVar(value=defaults.select == "best"),
        "stop_on_failure": tk.BooleanVar(value=defaults.stop_on_failure),
        "stop_on_collapse": tk.BooleanVar(value=defaults.stop_on_collapse),
        "also_final": tk.BooleanVar(value=False),
    }

    def add_row(parent, row, label, var, hint="", browse=None):
        tk.Label(parent, text=label, anchor="w").grid(
            row=row, column=0, sticky="w", padx=(8, 6), pady=3
        )
        entry = tk.Entry(parent, textvariable=var, width=52 if browse else 18)
        entry.grid(row=row, column=1, sticky="w", padx=(0, 6), pady=3)
        if browse is not None:
            tk.Button(parent, text="Browse...", command=browse).grid(
                row=row, column=2, sticky="w", padx=(0, 6), pady=3
            )
        if hint:
            tk.Label(parent, text=hint, anchor="w", fg="#555555").grid(
                row=row, column=3, sticky="w", padx=(4, 8), pady=3
            )
        return entry

    def pick_input():
        chosen = filedialog.askopenfilename(
            title="Select an .xyz link file",
            filetypes=[("xyz files", "*.xyz"), ("all files", "*.*")],
        )
        if chosen:
            fields["input"].set(chosen)
            # Mirror the CLI's default output name so the target is visible.
            if not fields["output"].get().strip():
                source = Path(chosen)
                fields["output"].set(str(source.with_name(f"{source.stem}_tight.xyz")))

    def pick_output():
        chosen = filedialog.asksaveasfilename(
            title="Tightened .xyz output", defaultextension=".xyz",
            filetypes=[("xyz files", "*.xyz"), ("all files", "*.*")],
        )
        if chosen:
            fields["output"].set(chosen)

    def pick_work_dir():
        chosen = filedialog.askdirectory(title="Directory for the VECT files and .rr tree")
        if chosen:
            fields["work_dir"].set(chosen)

    files_frame = tk.LabelFrame(root, text="Files")
    files_frame.pack(fill="x", padx=10, pady=(10, 4))
    add_row(files_frame, 0, "Input .xyz", fields["input"], browse=pick_input)
    add_row(files_frame, 1, "Output .xyz", fields["output"], browse=pick_output,
            hint="blank = <input>_tight.xyz")
    add_row(files_frame, 2, "Work directory", fields["work_dir"], browse=pick_work_dir,
            hint="blank = <output dir>/<stem>_rr")

    stop_frame = tk.LabelFrame(root, text="Stopping criteria (the run halts at the first one met)")
    stop_frame.pack(fill="x", padx=10, pady=4)
    add_row(stop_frame, 0, "Max steps", fields["steps"], hint="-s / --StopSteps")
    add_row(stop_frame, 1, "Stop residual", fields["stop_res"],
            hint="--StopResidual; 1e-4 is very tight, 1.0 means nothing is constrained")
    add_row(stop_frame, 2, "Stop time (whole minutes)", fields["stop_time"],
            hint="--StopTime; blank = no wall-clock limit")
    add_row(stop_frame, 3, "Stop20 (rope change)", fields["stop20"],
            hint="--Stop20; blank = off")

    geom_frame = tk.LabelFrame(root, text="Tube geometry")
    geom_frame.pack(fill="x", padx=10, pady=4)
    add_row(geom_frame, 0, "Tube radius", fields["tube_radius"], hint="blank = ridgerunner default 0.5")
    add_row(geom_frame, 1, "Stiffness (Lambda)", fields["stiffness"], hint="min radius of curvature")
    add_row(geom_frame, 2, "Symmetry group", fields["symmetry"], hint="e.g. Z/4Z, D2, cplanes; blank = none")

    run_frame = tk.LabelFrame(root, text="Run behaviour")
    run_frame.pack(fill="x", padx=10, pady=4)
    check_row = tk.Frame(run_frame)
    check_row.grid(row=0, column=0, columnspan=4, sticky="w", padx=8, pady=(4, 2))
    tk.Checkbutton(check_row, text="Autoscale to thickness 0.501 (-a)",
                   variable=flags["autoscale"]).pack(side="left", padx=(0, 12))
    tk.Checkbutton(check_row, text="Equilateralize (--EqOn)",
                   variable=flags["eq"]).pack(side="left", padx=(0, 12))
    tk.Checkbutton(check_row, text="Keep .rr run directory",
                   variable=flags["keep_run_dir"]).pack(side="left", padx=(0, 12))
    tk.Checkbutton(check_row, text="Keep snapshots + lsqr log",
                   variable=flags["keep_snapshots"]).pack(side="left", padx=(0, 12))
    tk.Checkbutton(check_row, text="POVRAY snapshot",
                   variable=flags["png"]).pack(side="left")
    select_row = tk.Frame(run_frame)
    select_row.grid(row=4, column=0, columnspan=4, sticky="w", padx=8, pady=(0, 4))
    tk.Checkbutton(
        select_row,
        text="Roll back to the best saved configuration if the run degenerates",
        variable=flags["select_best"],
    ).pack(side="left")
    tk.Checkbutton(
        select_row,
        text="Stop on linear-algebra failure",
        variable=flags["stop_on_failure"],
    ).pack(side="left", padx=(12, 0))
    select_row2 = tk.Frame(run_frame)
    select_row2.grid(row=6, column=0, columnspan=4, sticky="w", padx=8, pady=(0, 4))
    tk.Checkbutton(
        select_row2,
        text="Also write the final configuration (_final.xyz)",
        variable=flags["also_final"],
    ).pack(side="left")
    tk.Checkbutton(
        select_row2,
        text="Stop when the contact set collapses (approximate)",
        variable=flags["stop_on_collapse"],
    ).pack(side="left", padx=(12, 0))
    add_row(run_frame, 1, "Output decimals", fields["decimals"])
    add_row(run_frame, 2, "Extra ridgerunner args", fields["rr_arg"],
            hint="passed through verbatim, e.g. --Timewarp")
    add_row(run_frame, 3, "ridgerunner binary", fields["ridgerunner"],
            hint="blank = $RIDGERUNNER, then PATH")
    add_row(run_frame, 5, "Snapshot interval (steps)", fields["snapshot_interval"],
            hint="blank = ridgerunner default 10000; smaller = finer rollback points")

    button_frame = tk.Frame(root)
    button_frame.pack(fill="x", padx=10, pady=(6, 2))
    run_button = tk.Button(button_frame, text="Tighten")
    run_button.pack(side="left")
    stop_button = tk.Button(button_frame, text="Stop", state="disabled")
    stop_button.pack(side="left", padx=(8, 0))
    tk.Button(button_frame, text="Quit", command=root.destroy).pack(side="right")

    status_var = tk.StringVar(value="")
    tk.Label(root, textvariable=status_var, anchor="w").pack(fill="x", padx=10, pady=(2, 2))

    log_widget = ScrolledText(root, height=18, wrap="none", font=("Courier New", 10))
    log_widget.pack(fill="both", expand=True, padx=10, pady=(0, 10))

    # --- run plumbing -------------------------------------------------------
    log_queue: queue.Queue = queue.Queue()
    state: dict = {"process": None, "cancelled": False, "run_dir": None}

    def append_log(text: str) -> None:
        log_widget.insert("end", text + "\n")
        log_widget.see("end")

    def optional_number(key: str, label: str, cast):
        text = fields[key].get().strip()
        if not text:
            return None
        try:
            return cast(text)
        except ValueError:
            raise TightenError(f"{label}: expected a number, got {text!r}")

    def build_namespace() -> argparse.Namespace:
        source = fields["input"].get().strip()
        if not source:
            raise TightenError("Choose an input .xyz file first.")
        output = fields["output"].get().strip()
        work_dir = fields["work_dir"].get().strip()

        # Start from the parser's own defaults and override only what the GUI
        # exposes, so an option added to the CLI later cannot be missing here.
        settings = gui_defaults()
        settings.input = Path(source)
        settings.output = Path(output) if output else None
        settings.work_dir = Path(work_dir) if work_dir else None
        settings.steps = optional_number("steps", "Max steps", int)
        settings.stop_res = optional_number("stop_res", "Stop residual", float)
        settings.stop_time = optional_number("stop_time", "Stop time", int)
        settings.stop20 = optional_number("stop20", "Stop20", float)
        settings.tube_radius = optional_number("tube_radius", "Tube radius", float)
        settings.stiffness = optional_number("stiffness", "Stiffness", float)
        settings.symmetry = fields["symmetry"].get().strip() or None
        settings.decimals = (
            optional_number("decimals", "Output decimals", int) or defaults.decimals
        )
        settings.rr_arg = fields["rr_arg"].get().split()
        settings.ridgerunner = fields["ridgerunner"].get().strip() or None
        settings.autoscale = flags["autoscale"].get()
        settings.eq = flags["eq"].get()
        settings.keep_run_dir = flags["keep_run_dir"].get()
        settings.keep_snapshots = flags["keep_snapshots"].get()
        settings.png = flags["png"].get()
        settings.stop_on_failure = flags["stop_on_failure"].get()
        settings.stop_on_collapse = flags["stop_on_collapse"].get()
        settings.snapshot_interval = optional_number(
            "snapshot_interval", "Snapshot interval", int
        )
        if not flags["select_best"].get():
            settings.select = "final"
        else:
            settings.select = "both" if flags["also_final"].get() else "best"
        settings.quiet = True  # output is streamed into the log pane instead
        return settings

    def finish(error: str | None) -> None:
        run_button.configure(state="normal")
        stop_button.configure(state="disabled")
        if state["cancelled"]:
            if error:
                # RidgeRunner only rewrites <stem>.final.vect periodically, so a
                # run stopped in its first seconds has nothing to convert.
                status_var.set("Stopped before any configuration was saved; no output written.")
                append_log(f"--- stopped: {error} ---")
            else:
                status_var.set("Stopped early. Wrote the last configuration RidgeRunner had saved.")
                append_log("--- stopped by user ---")
        elif error:
            status_var.set("Failed.")
            append_log(f"--- error ---\n{error}")
            messagebox.showerror("Tightening failed", error)
        else:
            status_var.set("Done.")
        state["process"] = None

    def drain() -> None:
        try:
            while True:
                kind, payload = log_queue.get_nowait()
                if kind == "line":
                    append_log(payload)
                elif kind == "done":
                    finish(None)
                elif kind == "error":
                    finish(payload)
        except queue.Empty:
            pass
        root.after(120, drain)

    def worker(namespace: argparse.Namespace) -> None:
        try:
            tighten(
                namespace,
                emit=lambda line: log_queue.put(("line", line)),
                on_process=lambda proc: state.__setitem__("process", proc),
            )
            log_queue.put(("done", None))
        except TightenError as exc:
            log_queue.put(("error", str(exc)))
        except Exception as exc:  # pragma: no cover - unexpected failure path
            log_queue.put(("error", f"{type(exc).__name__}: {exc}"))

    def start_run() -> None:
        try:
            namespace = build_namespace()
        except TightenError as exc:
            messagebox.showerror("Check the parameters", str(exc))
            return

        log_widget.delete("1.0", "end")
        state["cancelled"] = False
        run_button.configure(state="disabled")
        stop_button.configure(state="normal")
        status_var.set("Running. RidgeRunner runs can take minutes to hours; Stop keeps what it has.")
        threading.Thread(target=worker, args=(namespace,), daemon=True).start()

    def stop_run() -> None:
        process = state["process"]
        if process is not None and process.poll() is None:
            state["cancelled"] = True
            status_var.set("Stopping...")
            process.terminate()

    run_button.configure(command=start_run)
    stop_button.configure(command=stop_run)

    # Report up front whether ridgerunner can be found, rather than at run time.
    try:
        status_var.set(f"Found ridgerunner: {find_ridgerunner(None)}")
    except TightenError as exc:
        status_var.set("ridgerunner not found - set the binary path below.")
        append_log(str(exc))

    root.after(120, drain)
    root.mainloop()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "input",
        type=Path,
        nargs="?",
        help="Input .xyz file; blank lines separate components. Omit it to open the GUI.",
    )
    parser.add_argument("-o", "--output", type=Path, help="Output .xyz path (default: <input>_tight.xyz).")
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Open the Tkinter GUI. This is also the default when no arguments are provided.",
    )

    stopping = parser.add_argument_group("stopping criteria (the run halts at the first one met)")
    stopping.add_argument("-s", "--steps", type=int, default=1000, help="Max gradient-descent steps (default: 1000).")
    stopping.add_argument(
        "--stop-res",
        type=float,
        default=0.01,
        help="Stop when the residual drops below this; 1e-4 is very tight (default: 0.01).",
    )
    stopping.add_argument("--stop-time", type=int, metavar="MINUTES", help="Stop after this many whole minutes of wall clock time.")
    stopping.add_argument("--stop20", type=float, metavar="DELTA", help="Stop when the ropelength change over the last 20 steps is below DELTA.")

    geometry = parser.add_argument_group("tube geometry")
    geometry.add_argument("--tube-radius", type=float, help="Tube radius (ridgerunner default: 0.5).")
    geometry.add_argument("--stiffness", type=float, help="Minimum radius of curvature, ridgerunner's --Lambda.")
    geometry.add_argument(
        "--symmetry",
        metavar="GROUP",
        help="Enforce a symmetry group about the z-axis, e.g. Z/4Z, D2 or cplanes (ridgerunner --Symmetry).",
    )

    run = parser.add_argument_group("run behaviour")
    run.add_argument(
        "--no-autoscale",
        dest="autoscale",
        action="store_false",
        help="Do not rescale to thickness 0.501 before stepping. Autoscaling is on by default because "
        "it clears the initial self-contacts that symmetric layouts tend to have.",
    )
    run.add_argument(
        "--eq",
        action="store_true",
        help="Re-equilateralize during the run (ridgerunner --EqOn). Off by default: it holds edges even "
        "and keeps a configuration stable, but it re-splines the curve, which repeatedly spikes ropelength "
        "and caps how low the residual gets. Without it a run reaches a much better residual but can later "
        "degenerate, which is what --select best is for.",
    )
    run.add_argument(
        "--stop-on-collapse",
        action="store_true",
        help="Also stop the run when its strut count collapses, which is how a run degenerates "
        "without ridgerunner logging anything. Off by default: the threshold was calibrated on "
        "only seven runs and the margin is thin, so a false positive would truncate a healthy "
        "run, while a missed collapse only wastes time that --select best recovers from.",
    )
    run.add_argument(
        "--no-stop-on-failure",
        dest="stop_on_failure",
        action="store_false",
        help="Keep stepping after ridgerunner reports a linear-algebra failure. By default the run "
        "is stopped there: the strut set collapses and does not recover, so the remaining steps are "
        "wasted. The best saved configuration is used either way.",
    )
    run.add_argument(
        "--select",
        choices=("best", "final", "both"),
        default="best",
        help="Which saved configuration to convert. 'best' (default) scores every snapshot on its "
        "windowed median residual and rolls back if the run degenerated after its best state; "
        "'final' always keeps the last one; 'both' writes the best to the output path and the "
        "last one alongside it as <output stem>_final.xyz, for comparison.",
    )
    run.add_argument(
        "--snapshot-interval",
        type=int,
        metavar="N",
        help="Save a configuration every N steps (ridgerunner default 10000). Smaller values give "
        "--select best a finer choice of rollback points, at the cost of disk during the run.",
    )
    run.add_argument("--keep-snapshots", action="store_true", help="Keep intermediate VECT frames, the lsqr log, and the snapshot matrix dumps.")
    run.add_argument("--png", action="store_true", help="Let ridgerunner render a final snapshot; needs POVRAY installed.")
    run.add_argument("--keep-run-dir", action="store_true", help="Keep the .rr run directory and its logfiles.")
    run.add_argument("--work-dir", type=Path, help="Directory for the VECT files and the .rr tree (implies --keep-run-dir).")
    run.add_argument("--rr-arg", action="append", default=[], metavar="ARG", help="Extra flags passed straight to ridgerunner; repeatable. Use the = form, since argparse reads a bare --Foo as an option: --rr-arg=--Timewarp, or --rr-arg='--MangleMode --AnimationStepper'.")
    run.add_argument("--ridgerunner", help="Path to the ridgerunner binary (default: $RIDGERUNNER, then PATH).")
    run.add_argument("--decimals", type=int, default=9, help="Decimal places in the output xyz (default: 9).")
    run.add_argument("-q", "--quiet", action="store_true", help="Do not echo ridgerunner's per-step progress.")
    return parser


def main(argv: list[str] | None = None) -> int:
    """Entry point for CLI and GUI modes."""
    args_in = sys.argv[1:] if argv is None else list(argv)
    parser = build_parser()
    args = parser.parse_args(args_in)

    if not args_in or args.gui:
        try:
            run_gui()
            return 0
        except Exception as exc:
            print(f"GUI could not be started: {exc}", file=sys.stderr)
            print("Use CLI mode, for example: python3 tighten_link_xyz.py link.xyz -s 1000", file=sys.stderr)
            return 1

    if args.input is None:
        parser.error("CLI mode requires an input .xyz file, unless --gui is used.")

    try:
        return tighten(args)
    except TightenError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
