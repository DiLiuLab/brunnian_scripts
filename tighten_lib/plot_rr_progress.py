#!/usr/bin/env python3
"""Plot ropelength, residual and strut count against step for RidgeRunner runs.

Runs that branch from a common configuration share a step origin; a run that
continues another is offset onto a cumulative axis, so a chain reads left to
right.

    python3 plot_rr_progress.py --set chain4

Add a new comparison by adding an entry to RUN_SETS.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUNS = Path.home() / "rr_runs"
OUT_DIR = Path("/Users/diliu/Documents/Claude_tmp")

# Reference categorical palette, slots 1-3 (validated all-pairs, light mode).
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8a8984"
SLOTS = ("#2a78d6", "#eb6834", "#1baf7a")

# Each run: key, run directory under ~/rr_runs, legend label, the key of the
# run it continues (None if it branches from the origin), and the step whose
# configuration was finally kept (None to leave unmarked).
RUN_SETS = {
    "chain3": {
        "title": "RidgeRunner tightening of link_sphere_4BL_wider (4-component Brunnian link)",
        "subtitle": (
            "noeq and eq both branch from the same 1 h configuration (ropelength 126.25); "
            "chain3 continues from eq.\nResidual returns to 1 at each run's start: a fresh "
            "run has no active constraints until it rebuilds its contact set."
        ),
        "runs": [
            ("noeq", "noeq/link_sphere_4BL_wider_tight.rr", "5 h  --no-eq", None, None),
            ("eq", "eq/link_sphere_4BL_wider_tight.rr", "5 h  --eq", None, None),
            ("chain3", "chain3/link_sphere_4BL_wider_5h_eq.rr", "5 h  --eq (chained)", "eq", None),
        ],
    },
    "sixhour": {
        "base": Path("/Users/diliu/Documents/Claude_tmp/6h_runs"),
        "title": "Six-hour --no-eq runs for the three other_diagrams_4BL links",
        "subtitle": (
            "Three separate links, each chained from its own 1 h --eq result, so they share a step "
            "origin rather than a lineage.\nRopelength is scale-invariant, so the three are directly "
            "comparable. Dots mark the configuration each run kept."
        ),
        "runs": [
            ("2_BC", "link_sphere_2_BC/link_sphere_2_BC_tight.rr", "link_sphere_2_BC", None, 48000),
            ("3_OC", "link_sphere_3_OC/link_sphere_3_OC_tight.rr", "link_sphere_3_OC", None, 118000),
            ("4_LC", "link_sphere_4_LC/link_sphere_4_LC_tight.rr", "link_sphere_4_LC", None, None),
        ],
    },
    "chain4": {
        "title": "Continuing from the recovered step-120000 state (ropelength 112.43, residual 0.022)",
        "subtitle": (
            "Both arms start from the same rolled-back configuration. direct runs --no-eq straight "
            "from it; repair spends 30 min on\n--eq to even the edges first, then runs --no-eq. "
            "Dots mark the configuration each run finally kept."
        ),
        "runs": [
            ("direct", "chain4_direct/link_sphere_4BL_wider_best120k.rr",
             "5 h  --no-eq  (direct)", None, 38000),
            ("repair-eq", "chain4_repair/stage1/link_sphere_4BL_wider_best120k.rr",
             "30 min  --eq  (repair 1)", None, None),
            ("repair-noeq", "chain4_repair/stage2/link_sphere_4BL_wider_repaired.rr",
             "5 h  --no-eq  (repair 2)", "repair-eq", 40386),
        ],
    },
}

BIN_STEPS = 400  # one bin width for every series, so spreads are comparable


def read_log(path: Path, column: int = 1) -> tuple[list[int], list[float]]:
    """Read a "step value ..." logfile, skipping torn mid-write rows.

    A live run's logfile can end mid-line, which parses as a valid but wildly
    wrong number (a truncated "113.38" reads as "11"). An unterminated final
    line is the tell, so drop it.
    """
    if not path.is_file():
        return [], []
    text = path.read_text(errors="ignore")
    lines = text.splitlines()
    if lines and not text.endswith("\n"):
        lines.pop()

    steps: list[int] = []
    values: list[float] = []
    for line in lines:
        fields = line.split()
        if len(fields) <= column:
            continue
        try:
            step = int(fields[0])
            value = float(fields[column])
        except ValueError:
            continue
        # Steps are strictly increasing; anything else is a torn row.
        if steps and step <= steps[-1]:
            continue
        steps.append(step)
        values.append(value)
    return steps, values


def envelope(xs, ys, width=BIN_STEPS):
    """Bin into fixed-width step windows, returning min/max/median per bin.

    Residual and strut count swing over most of their range from one step to
    the next, so a raw 100k-point line renders as a solid block. The median
    carries the trend and the min-max band shows the real spread.
    """
    import numpy as np

    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    bins = ((xs - xs.min()) // width).astype(int)
    centre, low, high, mid = [], [], [], []
    for b in np.unique(bins):
        mask = bins == b
        chunk = ys[mask]
        centre.append(xs[mask].mean())
        low.append(chunk.min())
        high.append(chunk.max())
        mid.append(float(np.median(chunk)))
    return centre, low, high, mid


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--set", dest="which", choices=sorted(RUN_SETS), default="chain4")
    parser.add_argument("--out", type=Path, help="Output PNG (default: <set>_progress.png).")
    args = parser.parse_args()

    spec = RUN_SETS[args.which]
    out = args.out or OUT_DIR / f"{args.which}_progress.png"

    base = spec.get("base", RUNS)
    data = {}
    for key, rel, _label, _after, _kept in spec["runs"]:
        logs = base / rel / "logfiles"
        data[key] = {
            "rope": read_log(logs / "ropelength.dat"),
            "residual": read_log(logs / "residual.dat"),
            "struts": read_log(logs / "strutcount.dat", column=1),
        }

    # A run that continues another starts where that one ended.
    offsets: dict[str, int] = {}
    for key, _rel, _label, after, _kept in spec["runs"]:
        if after and data[after]["rope"][0]:
            offsets[key] = offsets.get(after, 0) + max(data[after]["rope"][0])
        else:
            offsets[key] = 0

    colors = {key: SLOTS[i % len(SLOTS)] for i, (key, *_rest) in enumerate(spec["runs"])}

    fig, axes = plt.subplots(3, 1, figsize=(11, 9.5), sharex=True)
    fig.patch.set_facecolor(SURFACE)

    panels = [
        ("rope", "Ropelength", "linear", False),
        ("residual", "Residual  (0 = ropelength-critical)", "log", True),
        ("struts", "Active struts  (self-contacts)", "linear", True),
    ]

    for ax, (metric, title, scale, smooth) in zip(axes, panels):
        ax.set_facecolor(SURFACE)
        for key, _rel, label, _after, kept in spec["runs"]:
            steps, values = data[key][metric]
            if not steps:
                continue
            raw = dict(zip(steps, values))
            xs = [s + offsets[key] for s in steps]
            if scale == "log":
                pairs = [(x, v) for x, v in zip(xs, values) if v > 0]
                xs = [p[0] for p in pairs]
                values = [p[1] for p in pairs]
            if smooth:
                cx, lo, hi, mid = envelope(xs, values)
                ax.fill_between(cx, lo, hi, color=colors[key], alpha=0.16, linewidth=0)
                ax.plot(cx, mid, color=colors[key], linewidth=1.6, label=label)
                xs, values = cx, mid
            else:
                ax.plot(xs, values, color=colors[key], linewidth=1.0, alpha=0.95, label=label)

            # Direct label at the line end: identity is never colour-alone.
            ax.annotate(key, xy=(xs[-1], values[-1]), xytext=(6, 0),
                        textcoords="offset points", color=INK_2, fontsize=8.5,
                        va="center", clip_on=False)

            # Mark the configuration the run actually kept. Place it on the
            # plotted series -- on a smoothed panel a single step's raw value
            # is noise and would sit far off the median line.
            if kept is not None and xs:
                target = kept + offsets[key]
                idx = min(range(len(xs)), key=lambda i: abs(xs[i] - target))
                ax.plot([xs[idx]], [values[idx]], marker="o", markersize=7,
                        color=colors[key], markeredgecolor=SURFACE, markeredgewidth=2, zorder=5)
            del raw

        ax.set_yscale(scale)
        ax.set_title(title, color=INK, fontsize=11.5, loc="left", pad=8)
        ax.grid(True, color="#e6e5e1", linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#d5d4cf")
        ax.tick_params(colors=INK_2, labelsize=9)
        for key, _rel, _label, after, _kept in spec["runs"]:
            if after:
                ax.axvline(offsets[key], color=MUTED, linewidth=0.8, linestyle=(0, (4, 3)))

    axes[1].axhline(1e-4, color=MUTED, linewidth=0.8, linestyle=(0, (1, 2)))
    axes[1].annotate("--stop-res target 1e-4", xy=(0, 1e-4), xytext=(4, 4),
                     textcoords="offset points", color=MUTED, fontsize=8)

    axes[-1].set_xlabel("Cumulative gradient-descent step", color=INK, fontsize=10.5)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=9, labelcolor=INK_2,
               ncols=len(labels), loc="upper left", bbox_to_anchor=(0.05, 0.935))
    fig.suptitle(spec["title"], color=INK, fontsize=12.5, x=0.055, ha="left", y=0.985)
    fig.text(0.055, 0.962, spec["subtitle"], color=INK_2, fontsize=9, ha="left", va="top")
    fig.tight_layout(rect=(0, 0, 1, 0.895))
    fig.savefig(out, dpi=170, facecolor=SURFACE)
    print(f"wrote {out}")

    for key, _rel, _label, _after, _kept in spec["runs"]:
        rope = data[key]["rope"][1]
        res = data[key]["residual"][1]
        struts = data[key]["struts"][1]
        if not rope:
            print(f"{key:>12}: no data")
            continue
        best = min((v for v in res if v > 0), default=float("nan"))
        print(f"{key:>12}: steps {len(rope):>6}  rope {rope[0]:.3f} -> {rope[-1]:.3f}  "
              f"resid final {res[-1]:.4f} best {best:.4f}  "
              f"struts final {struts[-1]:.0f} peak {max(struts):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
