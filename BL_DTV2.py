#!/usr/bin/env python3
# BL_DTV2.py
# Version: V2_1
# Generate signed DT codes and generalized formulas for four Brunnian-link
# construction families. Inputs: a pattern name and component number n.
# Outputs: the selected generalized formula and the corresponding DT code.
# Example: python BL_DTV2.py --pattern cyclic_larks --n 5

import argparse
import os
import sys
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Component = Tuple[int, ...]
DTCode = List[Component]
VERSION = "V2_1"
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(BASE_DIR, "Assets")
APP_ICON_FILENAME = "cyclic larks.png"

PATTERN_ALIASES = {
    "cyclic_squares": "cyclic_squares",
    "cyclic_square": "cyclic_squares",
    "squares": "cyclic_squares",
    "square": "cyclic_squares",
    "cs": "cyclic_squares",
    "1": "cyclic_squares",
    "cyclic_larks": "cyclic_larks",
    "cyclic_lark": "cyclic_larks",
    "larks": "cyclic_larks",
    "lark": "cyclic_larks",
    "cl": "cyclic_larks",
    "2": "cyclic_larks",
    "cyclic_rubberband": "cyclic_rubberband",
    "cyclic_rubberbands": "cyclic_rubberband",
    "rubberband": "cyclic_rubberband",
    "rubberbands": "cyclic_rubberband",
    "cr": "cyclic_rubberband",
    "3": "cyclic_rubberband",
    "linear_rubberband": "linear_rubberband",
    "linear_rubberbands": "linear_rubberband",
    "linear": "linear_rubberband",
    "lr": "linear_rubberband",
    "4": "linear_rubberband",
}

PATTERN_TITLES = {
    "cyclic_squares": "Cyclic squares",
    "cyclic_larks": "Cyclic larks",
    "cyclic_rubberband": "Cyclic rubberband",
    "linear_rubberband": "Linear rubberband",
}

ORDERED_PATTERNS = [
    "cyclic_squares",
    "cyclic_larks",
    "cyclic_rubberband",
    "linear_rubberband",
]

PATTERN_IMAGE_FILENAMES = {
    "cyclic_squares": "cyclic squares.png",
    "cyclic_larks": "cyclic larks.png",
    "cyclic_rubberband": "cyclic rubberband.png",
    "linear_rubberband": "linear rubberband.png",
}


def normalize_pattern(pattern: str) -> str:
    """Return the canonical pattern key from a user-facing alias."""
    key = pattern.strip().lower().replace("-", "_").replace(" ", "_")
    if key not in PATTERN_ALIASES:
        valid = ", ".join(ORDERED_PATTERNS)
        raise ValueError("Unknown pattern {!r}. Valid patterns: {}".format(pattern, valid))
    return PATTERN_ALIASES[key]


def validate_n(n: int) -> None:
    """Require n >= 3 for the Brunnian-link families used here."""
    if n < 3:
        raise ValueError("n must be at least 3 for these Brunnian-link series.")


def r(value: int, modulus: int) -> int:
    """Reduce value modulo modulus, with residue 0 represented by modulus."""
    reduced = value % modulus
    return modulus if reduced == 0 else reduced


def cyclic_squares_dt(n: int) -> DTCode:
    """Generate DT code for the cyclic-squares Brunnian-link family."""
    validate_n(n)
    modulus = 12 * n
    components: DTCode = []
    for i in range(1, n + 1):
        components.append(
            (
                r(12 * i - 20, modulus),
                r(12 * i - 18, modulus),
                r(12 * i + 2, modulus),
                -r(12 * i + 10, modulus),
                -r(12 * i + 12, modulus),
                -r(12 * i - 16, modulus),
            )
        )
    return components


def cyclic_larks_dt(n: int) -> DTCode:
    """Generate DT code for the cyclic-larks Brunnian-link family."""
    validate_n(n)
    modulus = 12 * n
    components: DTCode = []
    for i in range(1, n + 1):
        components.append(
            (
                r(12 * i - 20, modulus),
                r(12 * i + 12, modulus),
                r(12 * i - 10, modulus),
                -r(12 * i + 10, modulus),
                -r(12 * i - 6, modulus),
                -r(12 * i - 16, modulus),
            )
        )
    return components


def cyclic_rubberband_dt(n: int) -> DTCode:
    """Generate DT code for the cyclic-rubberband Brunnian-link family."""
    validate_n(n)
    modulus = 16 * n
    components: DTCode = []
    for i in range(1, n + 1):
        components.append(
            (
                r(16 * i - 20, modulus),
                -r(16 * i - 28, modulus),
                r(16 * i + 8, modulus),
                -r(16 * i + 2, modulus),
                r(16 * i - 26, modulus),
                -r(16 * i - 18, modulus),
                -r(16 * i + 10, modulus),
                r(16 * i + 16, modulus),
            )
        )
    return components


def linear_rubberband_dt(n: int) -> DTCode:
    """Generate DT code for the linear-rubberband Brunnian-link family."""
    validate_n(n)
    if n == 3:
        return [(12, -8), (-16, -2, 14, 4), (-6, 10)]
    if n == 4:
        return [(16, -10), (28, -24, -2, -18, 22, 4), (8, -14, 30, 12, -6, -32), (-26, 20)]

    components: DTCode = []
    components.append((16, -10))
    components.append((32, -26, -2, -18, 24, 4))

    for i in range(3, n - 1):
        if i == 3:
            left = (8, -14, 12, -6)
        else:
            left = (16 * i - 36, -(16 * i - 44), 16 * i - 42, -(16 * i - 34))

        if i == n - 2:
            right = (16 * n - 42, -(16 * n - 46), -(16 * n - 40), 16 * n - 36)
        else:
            right = (16 * i - 8, -(16 * i - 14), -(16 * i - 6), 16 * i)

        l1, l2, l3, l4 = left
        r1, r2, r3, r4 = right
        components.append((l1, l2, r1, r2, l3, l4, r3, r4))

    components.append(
        (
            16 * n - 52,
            -(16 * n - 60),
            16 * n - 34,
            16 * n - 58,
            -(16 * n - 50),
            -(16 * n - 32),
        )
    )
    components.append((-(16 * n - 38), 16 * n - 44))
    return components


FORMULA_HEADER = """Notation used below
DT_n = [C_1, C_2, ..., C_n]
r_m(x) = x mod m, with residue 0 written as m
"""


def cyclic_squares_formula() -> str:
    return FORMULA_HEADER + """
Cyclic squares, n >= 3, m = 12n

Unified cyclic form, 1 <= i <= n:
C_i = ( r_m(12i - 20),  r_m(12i - 18),  r_m(12i + 2),
       -r_m(12i + 10), -r_m(12i + 12), -r_m(12i - 16) )

Expanded boundary form:
C_1 = (12n - 8, 12n - 6, 14, -22, -24, -(12n - 4))

C_i = (12i - 20, 12i - 18, 12i + 2,
       -(12i + 10), -(12i + 12), -(12i - 16)),  2 <= i <= n - 1

C_n = (12n - 20, 12n - 18, 2, -10, -12, -(12n - 16))
"""


def cyclic_larks_formula() -> str:
    return FORMULA_HEADER + """
Cyclic larks, n >= 3, m = 12n

Unified cyclic form, 1 <= i <= n:
C_i = ( r_m(12i - 20),  r_m(12i + 12),  r_m(12i - 10),
       -r_m(12i + 10), -r_m(12i - 6),  -r_m(12i - 16) )

Expanded boundary form:
C_1 = (12n - 8, 24, 2, -22, -6, -(12n - 4))

C_i = (12i - 20, 12i + 12, 12i - 10,
       -(12i + 10), -(12i - 6), -(12i - 16)),  2 <= i <= n - 1

C_n = (12n - 20, 12, 12n - 10, -10, -(12n - 6), -(12n - 16))
"""


def cyclic_rubberband_formula() -> str:
    return FORMULA_HEADER + """
Cyclic rubberband, n >= 3, m = 16n

Unified cyclic form, 1 <= i <= n:
C_i = ( r_m(16i - 20), -r_m(16i - 28),  r_m(16i + 8),  -r_m(16i + 2),
        r_m(16i - 26), -r_m(16i - 18), -r_m(16i + 10),  r_m(16i + 16) )

Expanded boundary form:
C_1 = (16n - 4, -(16n - 12), 24, -18,
       16n - 10, -(16n - 2), -26, 32)

C_i = (16i - 20, -(16i - 28), 16i + 8, -(16i + 2),
       16i - 26, -(16i - 18), -(16i + 10), 16i + 16),  2 <= i <= n - 1

C_n = (16n - 20, -(16n - 28), 8, -2,
       16n - 26, -(16n - 18), -10, 16)
"""


def linear_rubberband_formula() -> str:
    return FORMULA_HEADER + """
Linear rubberband, n >= 5 for the general formula
Special n = 3 and n = 4 cases are generated separately.

C_1 = (16, -10)

C_2 = (32, -26, -2, -18, 24, 4)

For 3 <= i <= n - 2, define L_i and R_i:

L_i = (8, -14, 12, -6),                                           if i = 3
L_i = (16i - 36, -(16i - 44), 16i - 42, -(16i - 34)),             if i > 3

R_i = (16i - 8, -(16i - 14), -(16i - 6), 16i),                   if i < n - 2
R_i = (16n - 42, -(16n - 46), -(16n - 40), 16n - 36),            if i = n - 2

If L_i = (l_1, l_2, l_3, l_4) and R_i = (r_1, r_2, r_3, r_4), then:
C_i = (l_1, l_2, r_1, r_2, l_3, l_4, r_3, r_4),  3 <= i <= n - 2

C_(n-1) = (16n - 52, -(16n - 60), 16n - 34,
           16n - 58, -(16n - 50), -(16n - 32))

C_n = (-(16n - 38), 16n - 44)
"""


DT_GENERATORS: Dict[str, Callable[[int], DTCode]] = {
    "cyclic_squares": cyclic_squares_dt,
    "cyclic_larks": cyclic_larks_dt,
    "cyclic_rubberband": cyclic_rubberband_dt,
    "linear_rubberband": linear_rubberband_dt,
}

FORMULA_GENERATORS: Dict[str, Callable[[], str]] = {
    "cyclic_squares": cyclic_squares_formula,
    "cyclic_larks": cyclic_larks_formula,
    "cyclic_rubberband": cyclic_rubberband_formula,
    "linear_rubberband": linear_rubberband_formula,
}


def format_component(component: Sequence[int]) -> str:
    """Format one DT component tuple without spaces inside the tuple."""
    return "(" + ",".join(str(value) for value in component) + ")"


def format_dt_code(components: DTCode) -> str:
    """Format a full DT code in the compact style used in the notes."""
    return "DT: [{}]".format(", ".join(format_component(component) for component in components))


def generate_output(pattern: str, n: int) -> str:
    """Generate formula and DT code text for a selected pattern and n."""
    canonical = normalize_pattern(pattern)
    validate_n(n)
    formula = FORMULA_GENERATORS[canonical]().rstrip()
    dt_code = format_dt_code(DT_GENERATORS[canonical](n))
    return "{}\n\n{}\n".format(formula, dt_code)


def write_text(path: str, text: str) -> None:
    """Write text output to a user-selected path."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def default_output_path(pattern: str, n: int) -> str:
    """Create a default text filename for saving output."""
    canonical = normalize_pattern(pattern)
    filename = "BL_DT_{}_{}_n{}.txt".format(VERSION, canonical, n)
    return os.path.abspath(filename)


def asset_path(filename: str) -> str:
    """Return an absolute path for a bundled asset filename."""
    return os.path.join(ASSETS_DIR, filename)


def pattern_image_path(pattern: str) -> str:
    """Return the snapshot path for a canonical pattern."""
    return asset_path(PATTERN_IMAGE_FILENAMES[pattern])


def run_gui() -> None:
    """Launch a Tkinter GUI for selecting pattern and n."""
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
        from tkinter.scrolledtext import ScrolledText
    except Exception as exc:  # pragma: no cover - only used when Tkinter is absent.
        raise RuntimeError("Tkinter GUI is not available in this Python environment: {}".format(exc))

    root = tk.Tk()
    root.title("BL_DTV2 {}: Brunnian-link DT code generator".format(VERSION))
    root.geometry("1180x760")

    try:
        icon_photo = tk.PhotoImage(file=asset_path(APP_ICON_FILENAME))
        root.iconphoto(True, icon_photo)
        root._bl_dtv2_icon = icon_photo
    except Exception:
        root._bl_dtv2_icon = None

    pattern_label = tk.Label(root, text="Pattern")
    pattern_label.grid(row=0, column=0, sticky="w", padx=10, pady=(10, 4))

    pattern_var = tk.StringVar(value="cyclic_squares")
    pattern_combo = ttk.Combobox(root, textvariable=pattern_var, values=ORDERED_PATTERNS, state="readonly", width=26)
    pattern_combo.grid(row=0, column=1, sticky="w", padx=10, pady=(10, 4))

    n_label = tk.Label(root, text="Number of components, n")
    n_label.grid(row=0, column=2, sticky="w", padx=10, pady=(10, 4))

    n_var = tk.StringVar(value="5")
    n_entry = tk.Entry(root, textvariable=n_var, width=10)
    n_entry.grid(row=0, column=3, sticky="w", padx=10, pady=(10, 4))

    formula_label = tk.Label(root, text="Generalized formula (selectable)")
    formula_label.grid(row=1, column=0, columnspan=4, sticky="w", padx=10, pady=(10, 4))

    formula_text = ScrolledText(root, height=20, wrap="word", font=("Courier New", 11))
    formula_text.grid(row=2, column=0, columnspan=4, sticky="nsew", padx=10, pady=(0, 8))

    preview_frame = tk.Frame(root)
    preview_frame.grid(row=1, column=4, rowspan=4, sticky="nsew", padx=(0, 10), pady=(10, 6))

    preview_title_var = tk.StringVar(value="")
    preview_title = tk.Label(preview_frame, textvariable=preview_title_var, anchor="center")
    preview_title.pack(fill="x", pady=(0, 6))

    preview_image_label = tk.Label(preview_frame, bd=1, relief="solid", width=330, height=330)
    preview_image_label.pack(fill="both", expand=True)

    dt_label = tk.Label(root, text="Generated DT code (auto-updates and selectable)")
    dt_label.grid(row=3, column=0, columnspan=4, sticky="w", padx=10, pady=(10, 4))

    dt_text = ScrolledText(root, height=10, wrap="word", font=("Courier New", 11))
    dt_text.grid(row=4, column=0, columnspan=4, sticky="nsew", padx=10, pady=(0, 6))

    status_var = tk.StringVar(value="")
    status_label = tk.Label(root, textvariable=status_var, anchor="w")
    status_label.grid(row=5, column=0, columnspan=4, sticky="ew", padx=10, pady=(0, 6))

    button_frame = tk.Frame(root)
    button_frame.grid(row=6, column=0, columnspan=4, sticky="w", padx=10, pady=(0, 10))

    root.grid_columnconfigure(0, weight=0)
    root.grid_columnconfigure(1, weight=1)
    root.grid_columnconfigure(2, weight=0)
    root.grid_columnconfigure(3, weight=1)
    root.grid_columnconfigure(4, weight=0, minsize=350)
    root.grid_rowconfigure(2, weight=3)
    root.grid_rowconfigure(4, weight=2)

    current_output: Dict[str, str] = {"text": ""}
    pending_update: Dict[str, Optional[str]] = {"job": None}
    preview_state: Dict[str, Optional[tk.PhotoImage]] = {"photo": None}

    def clear_and_insert(widget: ScrolledText, text: str) -> None:
        widget.configure(state="normal")
        widget.delete("1.0", "end")
        widget.insert("1.0", text)
        widget.configure(state="normal")  # Keep selectable/editable for copying.

    def current_n() -> int:
        n_text = n_var.get().strip()
        if not n_text:
            raise ValueError("Enter an integer n >= 3.")
        return int(n_text)

    def load_pattern_photo(pattern: str, max_width: int = 330, max_height: int = 330) -> tk.PhotoImage:
        photo = tk.PhotoImage(file=pattern_image_path(pattern))
        factor = max(
            1,
            (photo.width() + max_width - 1) // max_width,
            (photo.height() + max_height - 1) // max_height,
        )
        if factor > 1:
            photo = photo.subsample(factor, factor)
        return photo

    def update_preview(pattern: str) -> None:
        preview_title_var.set(PATTERN_TITLES[pattern])
        try:
            photo = load_pattern_photo(pattern)
            preview_state["photo"] = photo
            preview_image_label.configure(image=photo, text="")
        except Exception:
            preview_state["photo"] = None
            preview_image_label.configure(image="", text="Snapshot unavailable")

    def update_output(show_popup: bool = False) -> None:
        try:
            pattern = normalize_pattern(pattern_var.get())
            update_preview(pattern)
            formula = FORMULA_GENERATORS[pattern]().rstrip() + "\n"
            clear_and_insert(formula_text, formula)

            n = current_n()
            validate_n(n)
            dt_code = format_dt_code(DT_GENERATORS[pattern](n)) + "\n"
            clear_and_insert(dt_text, dt_code)
            current_output["text"] = "{}\n{}".format(formula, dt_code)
            status_var.set("Updated: {}, n = {}".format(PATTERN_TITLES[pattern], n))
        except Exception as exc:
            current_output["text"] = ""
            clear_and_insert(dt_text, "Error: {}\n".format(exc))
            status_var.set("Fix the input to update the DT code.")
            if show_popup:
                messagebox.showerror("BL_DTV2 error", str(exc))

    def schedule_update(*_args) -> None:
        job = pending_update.get("job")
        if job is not None:
            root.after_cancel(job)
        pending_update["job"] = root.after(180, update_output)

    def save_output() -> None:
        if not current_output["text"]:
            update_output(show_popup=True)
        if not current_output["text"]:
            return
        try:
            n = current_n()
            default_name = os.path.basename(default_output_path(pattern_var.get(), n))
        except Exception:
            default_name = "BL_DTV2_output.txt"
        path = filedialog.asksaveasfilename(
            title="Save BL_DTV2 output",
            initialfile=default_name,
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if path:
            try:
                write_text(path, current_output["text"])
                messagebox.showinfo("BL_DTV2", "Saved to {}".format(path))
            except Exception as exc:
                messagebox.showerror("BL_DTV2 error", str(exc))

    refresh_button = tk.Button(button_frame, text="Refresh", command=lambda: update_output(show_popup=True))
    refresh_button.pack(side="left", padx=(0, 8))

    save_button = tk.Button(button_frame, text="Save output...", command=save_output)
    save_button.pack(side="left")

    pattern_combo.bind("<<ComboboxSelected>>", lambda _event: schedule_update())
    n_entry.bind("<Return>", lambda _event: update_output(show_popup=True))
    n_var.trace_add("write", schedule_update)

    update_output()
    root.mainloop()


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        description="Generate DT codes for four Brunnian-link construction families. Version: {}".format(VERSION)
    )
    parser.add_argument("--version", action="version", version="BL_DTV2 {}".format(VERSION))
    parser.add_argument(
        "--pattern",
        "-p",
        help="Pattern name or alias: cyclic_squares, cyclic_larks, cyclic_rubberband, or linear_rubberband.",
    )
    parser.add_argument("--n", "-n", type=int, help="Number of components. Must be >= 3.")
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Open the Tkinter GUI. This is also the default when no arguments are provided.",
    )
    parser.add_argument("--output", "-o", help="Optional text file to save the formula and DT code.")
    parser.add_argument("--list-patterns", action="store_true", help="List canonical pattern names and exit.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point for CLI and GUI modes."""
    args_in = sys.argv[1:] if argv is None else list(argv)
    parser = build_parser()
    args = parser.parse_args(args_in)

    if args.list_patterns:
        for pattern in ORDERED_PATTERNS:
            print(pattern)
        return 0

    if len(args_in) == 0 or args.gui:
        try:
            run_gui()
            return 0
        except Exception as exc:
            print("GUI could not be started: {}".format(exc), file=sys.stderr)
            print("Use CLI mode, for example: python BL_DTV2.py --pattern cyclic_larks --n 5", file=sys.stderr)
            return 1

    if args.pattern is None or args.n is None:
        parser.error("CLI mode requires both --pattern and --n, unless --gui is used.")

    try:
        output = generate_output(args.pattern, args.n)
        print(output, end="")
        if args.output:
            write_text(args.output, output)
    except Exception as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
