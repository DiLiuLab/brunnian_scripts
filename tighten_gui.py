#!/usr/bin/env python3
"""A tkinter front end for tighten_cycle.py, built from its own argument parser.

Every label, help string, default and type here is READ OUT OF THE PARSER at
run time. Nothing is duplicated: add an argument to tighten_cycle.build_parser()
and it appears in this window with its help text, in its own group, with the
right widget. The one thing this module adds is EXAMPLES, which argparse has
nowhere to put -- and only for the arguments where an example teaches something
a sentence cannot.

Why it does not run the job in-process. A cycle is an hours-long unattended
job. Running it inside the GUI would tie it to the window: close the window, or
let the laptop's Python get killed, and hours of descent go with it. So "Run in
background" spawns a detached process writing to a log file and hands back the
path and the command to watch it. The window is a launcher, not a host.

tkinter is in the standard library, so this adds no dependency to a project
that a volunteer has to install on their own machine.

Entered automatically when tighten_cycle.py is run with no arguments, or
explicitly with --gui. tighten_lib/slp_tighten.py opens the same window for the
SLP polish the same way (no arguments, or --gui), passing its own blurb,
examples, file pickers and log location; App knows nothing about either tool.
"""
from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path

import tkinter as tk
from tkinter import filedialog, ttk

HELP_BLUE = "#bfe3f5"      # the "?" chip; light enough for black text on top
HELP_BLUE_HOT = "#94cfeb"
INK = "#1a1a1a"

REPO = Path(__file__).resolve().parent

# Examples only where one earns its place. argparse has no slot for these, and
# a bad example is worse than none, so most arguments deliberately have none
# and fall back to their help text alone.
EXAMPLES: dict[str, str] = {
    "input": "~/tightening/my_link/raw.xyz\n\n"
             "Blank lines separate components. Leave EMPTY when using --resume: "
             "the input then comes from the work directory.",
    "work_dir": "~/tightening/my_link/cycle\n\n"
                "Keep it OUT of Dropbox or any syncing folder. RidgeRunner "
                "rewrites a multi-megabyte constraint matrix continuously and a "
                "sync client will fight it for the file.",
    "resume": "Tick this and leave the input blank to continue a cycle whose "
              "driver died (crash, reboot, kill).\n\n"
              "It restarts at the first round that produced no output and "
              "rebuilds the topology reference from round 0's input, refusing "
              "to continue if that disagrees with the recorded one. A descent "
              "cut off part-way is NOT resumed: its round is redone.",
    "symmetry": "Z/5Z for a 5-loop link, Z/7Z for a 7-loop one.\n\n"
                "This is RidgeRunner's own --Symmetry, held during the descent. "
                "Leave blank for a link with no usable symmetry.",
    "sym_group": "Cn\n\nUsed by our own symmetry tool, together with "
                 "--sym-order. Pair 'Cn' with the loop count.",
    "sym_order": "5 for a 5-loop link, 7 for a 7-loop one.\n\n"
                 "Must match the n in --symmetry Z/nZ.",
    "sym_axis": "0,0,1\n\nThe rotation axis. Files from our drawing scripts put "
                "it on z, which is this default.",
    "no_permute": "Leave this OFF for our standard links.\n\n"
                  "A Cn rotation on an n-component link CYCLES the components, "
                  "so forbidding the permutation asks for a symmetry the link "
                  "does not have.",
    "squeeze_factors": "0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9\n\n"
                       "f is the SLOPE of the radial map, so SMALLER f is the "
                       "MORE aggressive squeeze. f=1.0 would be no squeeze at "
                       "all. Add gentler values like 0.95 to explore near the "
                       "identity.",
    "plateau_on": "auto   (the default)\n"
                  "minrad,struts\n"
                  "minrad,struts,ropelength\n\n"
                  "'auto' keys the choice on the round's own move. Use the full "
                  "three for any round whose number you will REPORT: stopping "
                  "early propagates, because the next round then starts from a "
                  "less converged configuration.",
    "vu_ladder": "2,4,8   (the author's rungs)\n\n"
                 "Leave EMPTY to do no refinement. Refinement is not free and "
                 "the two modes differ sharply -- see --refine-mode.",
    "refine_mode": "spline reallocates vertices by arclength and evens the edge "
                   "lengths out, which can RAISE minRad.\n\n"
                   "subdivide inserts midpoints, which divides minRad by the "
                   "subdivision factor. Measured on one link: spline took "
                   "minRad 0.5098 -> 0.6671 while subdivide on another took "
                   "0.5115 -> 0.4474. They are not interchangeable.",
    "find_axis": "Leave OFF for any link with a usable symmetry.\n\n"
                 "On a Cn link the squeeze axis MUST be the symmetry axis. A "
                 "searched axis came back 6.4 degrees off on a 7-fold link, "
                 "which broke C7 and aborted the round. Turn this on only when "
                 "the structure has no symmetry to work with.",
    "dry_run": "Strongly recommended before committing hours.\n\n"
               "Measures everything and prints exactly what it WOULD do -- "
               "which move, which factor, and why -- without modifying a file "
               "or spending compute.",
    "rr_arg": "--Timewarp\n\nPassed straight through to ridgerunner. Separate "
              "several with spaces; each becomes its own --rr-arg.",
    "min_gain": "0.5\n\nGive up when a round's predicted gain is smaller than "
                "this, rather than grinding on for fractions.",
}


# --------------------------------------------------------------------------- #
# widgets
# --------------------------------------------------------------------------- #
class HelpChip(tk.Label):
    """A '?' that opens the argument's own help text in a popup.

    A tk.Label rather than a tk.Button on purpose: the macOS aqua theme ignores
    a Button's background colour, so a Button here would render grey and the
    chips would be invisible as a group. A Label honours bg everywhere.
    """

    def __init__(self, master, title, helptext, default, example=None):
        super().__init__(master, text="?", bg=HELP_BLUE, fg=INK, width=2,
                         relief="raised", borderwidth=1, cursor="hand2",
                         font=("Helvetica", 11, "bold"))
        self._args = (title, helptext, default, example)
        self.bind("<Button-1>", self._open)
        self.bind("<Enter>", lambda e: self.config(bg=HELP_BLUE_HOT))
        self.bind("<Leave>", lambda e: self.config(bg=HELP_BLUE))

    def _open(self, _event=None):
        title, helptext, default, example = self._args
        win = tk.Toplevel(self)
        win.title(title)
        win.transient(self.winfo_toplevel())
        win.configure(bg="white")
        frm = tk.Frame(win, bg="white", padx=16, pady=14)
        frm.pack(fill="both", expand=True)
        tk.Label(frm, text=title, bg="white", fg="#1f4e79",
                 font=("Helvetica", 13, "bold")).pack(anchor="w")
        tk.Label(frm, text=helptext or "(no help text)", bg="white", fg=INK,
                 wraplength=460, justify="left",
                 font=("Helvetica", 11)).pack(anchor="w", pady=(8, 0))
        if default not in (None, "", False, []):
            tk.Label(frm, text=f"default:  {default}", bg="white", fg="#555",
                     font=("Helvetica", 10, "italic")).pack(anchor="w", pady=(8, 0))
        if example:
            tk.Label(frm, text="example", bg="white", fg="#1f4e79",
                     font=("Helvetica", 11, "bold")).pack(anchor="w", pady=(12, 2))
            box = tk.Text(frm, height=max(3, example.count("\n") + 2), width=56,
                          wrap="word", bg="#f5f9fc", fg=INK, relief="flat",
                          font=("Menlo", 10), padx=8, pady=6)
            box.insert("1.0", example)
            box.config(state="disabled")
            box.pack(anchor="w", fill="x")
        tk.Button(frm, text="Close", command=win.destroy).pack(anchor="e",
                                                               pady=(14, 0))
        win.bind("<Escape>", lambda e: win.destroy())
        win.update_idletasks()
        # place it near the chip rather than at the screen origin
        x = self.winfo_rootx() + 30
        y = max(0, self.winfo_rooty() - 40)
        win.geometry(f"+{x}+{y}")


class Row:
    """One argparse action rendered as a labelled widget plus a help chip."""

    def __init__(self, parent, action, row, browse=None, examples=None):
        self.action = action
        self.browse = browse if browse is not None else {"input": "open", "work_dir": "dir"}
        self.examples = examples if examples is not None else EXAMPLES
        self.dest = action.dest
        self.is_flag = action.nargs == 0
        name = action.option_strings[0] if action.option_strings else action.dest

        tk.Label(parent, text=name, anchor="w", width=22,
                 font=("Menlo", 11)).grid(row=row, column=0, sticky="w",
                                          padx=(10, 6), pady=3)

        if self.is_flag:
            self.var = tk.BooleanVar(value=bool(action.default))
            self.widget = tk.Checkbutton(parent, variable=self.var)
            self.widget.grid(row=row, column=1, sticky="w")
        elif action.choices:
            self.var = tk.StringVar(value=str(action.default))
            self.widget = ttk.Combobox(parent, textvariable=self.var, width=28,
                                       values=[str(c) for c in action.choices],
                                       state="readonly")
            self.widget.grid(row=row, column=1, sticky="w")
        else:
            default = "" if action.default in (None, []) else str(action.default)
            self.var = tk.StringVar(value=default)
            self.widget = tk.Entry(parent, textvariable=self.var, width=31,
                                   font=("Menlo", 11))
            self.widget.grid(row=row, column=1, sticky="w")

        col = 2
        if self.dest in self.browse:
            tk.Button(parent, text="Browse",
                      command=self._browse).grid(row=row, column=col, padx=4)
            col += 1
        HelpChip(parent, name, action.help, action.default,
                 self.examples.get(self.dest)).grid(row=row, column=col, padx=(6, 10))

    def _browse(self):
        kind = self.browse.get(self.dest, "open")
        title = {"gate_ref": "Reference link for the HOMFLY gate"}.get(self.dest, "Input link")
        if kind == "dir":
            got = filedialog.askdirectory(title="Work directory (keep out of Dropbox)")
        elif kind == "save":
            got = filedialog.asksaveasfilename(
                title="Output file", defaultextension=".xyz",
                filetypes=[("link files", "*.xyz"), ("all files", "*")])
        else:
            got = filedialog.askopenfilename(
                title=title, filetypes=[("link files", "*.xyz *.vect"),
                                               ("all files", "*")])
        if got:
            self.var.set(got)

    def argv(self) -> list[str]:
        """This row's contribution to the command line, empty when unchanged."""
        a = self.action
        if self.is_flag:
            return [a.option_strings[0]] if self.var.get() else []
        val = self.var.get().strip()
        if not val:
            return []
        if not a.option_strings:                    # the positional
            return [val]
        if a.default is not None and val == str(a.default):
            return []                               # leave defaults implicit
        if getattr(a, "nargs", None) is None and a.dest == "rr_arg":
            # opt=value, so a piece that starts with '-' (--Timewarp) is not read as an option
            return [f"{a.option_strings[0]}={piece}" for piece in val.split()]
        if val.startswith("-"):                     # e.g. --ref -1,0,0 or --slp-arg --no-soc
            return [f"{a.option_strings[0]}={val}"]
        return [a.option_strings[0], val]


# --------------------------------------------------------------------------- #
# window
# --------------------------------------------------------------------------- #
class App:
    """Renders any argparse parser. It imports nothing from this project: the
    parser and the script to launch are both handed in, so the same window
    serves tighten_cycle.py, sono_link_xyz.py or anything else with a parser."""

    def __init__(self, parser, script: Path = REPO / "tighten_cycle.py", *, blurb=None,
                 examples=None, browse=None, log_for=None, quick_run=False):
        """blurb: the header text; examples: dest -> example text; browse: dest ->
        'open' | 'save' | 'dir' (which rows get a Browse button); log_for(app) ->
        (log Path, None) or (None, reason) for Run in background; quick_run, if
        given, adds a button that runs the command here and shows its output; it is
        True or a predicate quick_run(app) -> (ok, reason) that refuses anything but
        the modes that finish in seconds (such as slp_tighten.py --group detect)."""
        self.parser = parser
        self.script = Path(script)
        self.examples = examples if examples is not None else EXAMPLES
        self.browse = browse if browse is not None else {"input": "open", "work_dir": "dir"}
        self.log_for = log_for if log_for is not None else App._cycle_log
        self.quick_ok = quick_run if callable(quick_run) else (lambda app: (True, None))
        self.root = tk.Tk()
        self.root.title(f"{self.script.stem} — parameters")
        self.rows: list[Row] = []

        head = tk.Frame(self.root, bg="white", padx=14, pady=10)
        head.pack(fill="x")
        tk.Label(head, text=self.script.stem, bg="white", fg="#1f4e79",
                 font=("Helvetica", 16, "bold")).pack(anchor="w")
        tk.Label(head, bg="white", fg="#555", justify="left",
                 font=("Helvetica", 11), wraplength=760,
                 text=blurb or ("Fill in what you need and leave the rest at its default. "
                       "Every field has a  ?  that explains it. A cycle runs for "
                       "hours, so Run starts it in the background, detached from "
                       "this window — closing the window will not stop it.")
                 ).pack(anchor="w", pady=(4, 0))

        nb = ttk.Notebook(self.root)
        nb.pack(fill="both", expand=True, padx=10, pady=(6, 0))
        for group in parser._action_groups:
            acts = [a for a in group._group_actions
                    if a.dest not in ("help", "gui")]
            if not acts:
                continue
            tab = tk.Frame(nb)
            nb.add(tab, text=group.title.replace("arguments", "args"))
            for i, action in enumerate(acts):
                self.rows.append(Row(tab, action, i, self.browse, self.examples))

        self.cmd = tk.Text(self.root, height=5, wrap="word", font=("Menlo", 10),
                           bg="#f5f9fc", relief="flat", padx=10, pady=8)
        self.cmd.pack(fill="x", padx=10, pady=(10, 0))
        self._show("Press “Build command” to assemble the command line.")

        bar = tk.Frame(self.root, padx=10, pady=10)
        bar.pack(fill="x")
        tk.Button(bar, text="Build command", command=self.build).pack(side="left")
        tk.Button(bar, text="Copy", command=self.copy).pack(side="left", padx=6)
        tk.Button(bar, text="Run in background",
                  command=self.run).pack(side="left", padx=6)
        if quick_run:
            tk.Button(bar, text="Run now (short jobs)",
                      command=self.run_now).pack(side="left", padx=6)
        tk.Button(bar, text="Quit", command=self.root.destroy).pack(side="right")

    # -- helpers ---------------------------------------------------------- #
    def _show(self, text):
        self.cmd.config(state="normal")
        self.cmd.delete("1.0", "end")
        self.cmd.insert("1.0", text)
        self.cmd.config(state="disabled")

    def _argv(self) -> list[str] | None:
        """Assemble argv and let argparse itself validate it."""
        argv: list[str] = []
        for r in self.rows:
            if not r.action.option_strings:      # positional goes first
                argv = r.argv() + argv
            else:
                argv += r.argv()
        import io
        import contextlib
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                self.parser.parse_args(argv)
        except SystemExit:
            self._popup("That does not parse", err.getvalue().strip()
                        or "argparse rejected these values.")
            return None
        return argv

    def _popup(self, title, body):
        win = tk.Toplevel(self.root)
        win.title(title)
        tk.Label(win, text=body, wraplength=520, justify="left",
                 padx=16, pady=14, font=("Helvetica", 11)).pack()
        tk.Button(win, text="Close", command=win.destroy).pack(pady=(0, 12))

    def command(self) -> list[str] | None:
        argv = self._argv()
        if argv is None:
            return None
        return [sys.executable, str(self.script)] + argv

    # -- buttons ----------------------------------------------------------- #
    def build(self):
        cmd = self.command()
        if cmd:
            self._show(" ".join(shlex.quote(c) for c in cmd))

    def copy(self):
        cmd = self.command()
        if not cmd:
            return
        line = " ".join(shlex.quote(c) for c in cmd)
        self.root.clipboard_clear()
        self.root.clipboard_append(line)
        self._show(line + "\n\n(copied to the clipboard)")

    def value(self, dest: str) -> str:
        return next((r.var.get().strip() for r in self.rows if r.dest == dest), "")

    @staticmethod
    def _cycle_log(app):
        wd = app.value("work_dir")
        if not wd:
            return None, ("Set --work-dir before running. It holds the run tree, "
                          "and it should be outside Dropbox.")
        work = Path(wd).expanduser()
        work.mkdir(parents=True, exist_ok=True)
        return work / "cycle.log", None

    def run_now(self, timeout: float = 180):
        """Run the command here and show its output: for jobs that finish in seconds.
        The window waits while it runs, so the tool's predicate refuses anything longer."""
        cmd = self.command()
        if not cmd:
            return
        ok, why = self.quick_ok(self)
        if not ok:
            self._popup("Use Run in background", why)
            return
        line = " ".join(shlex.quote(c) for c in cmd)
        self._show(line + "\n\nrunning ...")
        self.root.update_idletasks()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
            out = (p.stdout + p.stderr).strip()
            self._show(f"{line}\n\nexit {p.returncode}\n\n{out[-6000:]}")
        except subprocess.TimeoutExpired as exc:
            got = ((exc.stdout or b"") + (exc.stderr or b""))
            got = got.decode(errors="replace") if isinstance(got, bytes) else got
            self._show(line + f"\n\nSTOPPED after {timeout:.0f} s (the job was killed) -- use "
                              f"Run in background for anything longer.\n\n{got[-4000:]}")

    def run(self):
        cmd = self.command()
        if not cmd:
            return
        log, why = self.log_for(self)
        if log is None:
            self._popup("Cannot start in the background", why)
            return
        with log.open("a") as sink:
            subprocess.Popen(cmd, stdout=sink, stderr=subprocess.STDOUT,
                             start_new_session=True)
        self._show(" ".join(shlex.quote(c) for c in cmd)
                   + f"\n\nstarted in the background; it keeps running if you "
                     f"close this window.\nlog:   {log}\nwatch: tail -f {shlex.quote(str(log))}")

    def go(self):
        self.root.mainloop()


def launch(parser, script: Path = REPO / "tighten_cycle.py", **options) -> int:
    """Open the window for any parser; options go to App (blurb, examples,
    browse, log_for, quick_run)."""
    try:
        App(parser, script, **options).go()
    except tk.TclError as exc:
        print(f"error: cannot open a window ({exc}).\n"
              f"       Run {Path(script).name} with --help for the command line "
              f"interface instead.", file=sys.stderr)
        return 2
    return 0
