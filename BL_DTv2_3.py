#!/usr/bin/env python3
# BL_DTv2_3.py
# Version: V2_3
# Generate signed DT codes and generalized formulas for eight link-construction
# families: four Brunnian-link series (cyclic squares, cyclic larks, cyclic
# rubberband, linear rubberband), Edwards' Venn (AM_n, Maes-Cerf), Brunn's
# classic (BR_n, Brunn 1892 / Rolfsen), the Fishtail bracelet, and the
# Mirror fishtail (neither of the last two is a Brunnian link, but both are
# closely related and interesting).
#
# V2_3 CORRECTION (July 2026).  The 'fishtail' code emitted by V2_2 was the
# MIRRORED stitch, not the true fishtail.  The decisive test: the two-peg
# stitch at grip span j = 1 must reduce to the classical single chain
# (= cyclic_rubberband).  The same-fold wrap rule passes that test (isometric
# exteriors at n = 5, 6); the mirrored rule fails it (58.0854 vs 57.8144 at
# n = 5).  V2_3 therefore emits the SAME-FOLD code under 'fishtail', and keeps
# the V2_2 code available as the new pattern 'mirror_fishtail' (aliases:
# v2_2_fishtail, twisted_fishtail).  The two families share every unsigned
# entry and differ in exactly 8 of the 16 signs; they are distinct links
# (n = 5 volumes 91.7463 vs 92.7683).
# Inputs: a pattern name and component number n.
# Outputs: the selected generalized formula and the corresponding DT code.
# Example: python BL_DTv2_3.py --pattern fishtail --n 5

import argparse
import os
import sys
from typing import Callable, Dict, List, Optional, Sequence, Tuple

Component = Tuple[int, ...]
DTCode = List[Component]
VERSION = "V2_3"
APP_NAME = "BL_DTv2_3"
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
    "edwards_venn": "edwards_venn",
    "edwards": "edwards_venn",
    "venn": "edwards_venn",
    "am": "edwards_venn",
    "amn": "edwards_venn",
    "ev": "edwards_venn",
    "5": "edwards_venn",
    "brunn_classic": "brunn_classic",
    "brunns_classic": "brunn_classic",
    "brunn": "brunn_classic",
    "brunns": "brunn_classic",
    "br": "brunn_classic",
    "brn": "brunn_classic",
    "classic": "brunn_classic",
    "6": "brunn_classic",
    "fishtail": "fishtail",
    "fishtail_bracelet": "fishtail",
    "fish": "fishtail",
    "ft": "fishtail",
    "7": "fishtail",
    "mirror_fishtail": "mirror_fishtail",
    "mirror": "mirror_fishtail",
    "mirrored_fishtail": "mirror_fishtail",
    "twisted_fishtail": "mirror_fishtail",
    "v2_2_fishtail": "mirror_fishtail",
    "mf": "mirror_fishtail",
    "8": "mirror_fishtail",
}

PATTERN_TITLES = {
    "cyclic_squares": "Cyclic squares",
    "cyclic_larks": "Cyclic larks",
    "cyclic_rubberband": "Cyclic rubberband",
    "linear_rubberband": "Linear rubberband",
    "edwards_venn": "Edwards' Venn (AM_n)",
    "brunn_classic": "Brunn's classic (BR_n)",
    "fishtail": "Fishtail bracelet",
    "mirror_fishtail": "Mirror fishtail (V2_2 'fishtail')",
}

ORDERED_PATTERNS = [
    "cyclic_squares",
    "cyclic_larks",
    "cyclic_rubberband",
    "linear_rubberband",
    "edwards_venn",
    "brunn_classic",
    "fishtail",
    "mirror_fishtail",
]

PATTERN_IMAGE_FILENAMES = {
    "cyclic_squares": "cyclic squares.png",
    "cyclic_larks": "cyclic larks.png",
    "cyclic_rubberband": "cyclic rubberband.png",
    "linear_rubberband": "linear rubberband.png",
    "edwards_venn": "Edward.png",
    "brunn_classic": "Brunn.png",
    "fishtail": "fishtail.png",
    "mirror_fishtail": "mirror fishtail.png",
}

# Default n per pattern: 7 for the four original Brunnian series, 5 for the
# three families added in V2_2.
PATTERN_DEFAULT_N = {
    "cyclic_squares": 7,
    "cyclic_larks": 7,
    "cyclic_rubberband": 7,
    "linear_rubberband": 7,
    "edwards_venn": 5,
    "brunn_classic": 5,
    "fishtail": 5,
    "mirror_fishtail": 5,
}

# Snapshot images in Assets/ show examples with this many components.
PATTERN_PREVIEW_COMPONENTS = dict(PATTERN_DEFAULT_N)

# AM_n and BR_n crossing numbers grow exponentially (2^n - 2 and 3*2^(n-1) - 4);
# cap n so the GUI and CLI stay responsive.
PATTERN_MAX_N = {
    "edwards_venn": 14,
    "brunn_classic": 14,
}

PATTERN_NOTES = {
    "edwards_venn": (
        "Note: AM_n is a Brunnian link whose projection is Edwards' construction\n"
        "of rotationally symmetric Venn diagrams (Maes & Cerf, JKTR 10(1), 2001).\n"
        "The diagram has 2^n - 2 crossings."
    ),
    "brunn_classic": (
        "Note: BR_n is Brunn's classical Brunnian link (Brunn 1892; Rolfsen,\n"
        "Knots and Links, p.67 Exercise 8): n-1 round circles plus one component\n"
        "woven through them, with 3*2^(n-1) - 4 crossings."
    ),
    "fishtail": (
        "Note: the Fishtail bracelet is NOT a Brunnian link -- removing any one\n"
        "band never unlinks it -- but it is closely related and interesting: it\n"
        "obeys a position-graded unlinking law (a removal set unlinks the\n"
        "bracelet exactly when it contains two CYCLICALLY ADJACENT bands; any\n"
        "two bands further apart leave a nontrivial link).\n"
        "V2_3: this is the corrected, same-fold code. Its grip-span-1 reduction\n"
        "is the classical chain (= cyclic_rubberband), which is what identifies\n"
        "it as the physical fishtail. The V2_2 code is now 'mirror_fishtail'."
    ),
    "mirror_fishtail": (
        "Note: the Mirror fishtail is the mirrored-rail stitch that V2_2 emitted\n"
        "under the name 'fishtail'. It is a closely related but DIFFERENT link\n"
        "(same unsigned DT entries, 8 of 16 signs flipped; n=5 volumes 92.7683 vs\n"
        "91.7463 for the true fishtail), and it does NOT reduce to the classical\n"
        "chain at grip span 1. Kept for reproducibility of V2_2 output and\n"
        "because it is an interesting family in its own right."
    ),
}


def normalize_pattern(pattern: str) -> str:
    """Return the canonical pattern key from a user-facing alias."""
    key = pattern.strip().lower().replace("-", "_").replace(" ", "_")
    if key not in PATTERN_ALIASES:
        valid = ", ".join(ORDERED_PATTERNS)
        raise ValueError("Unknown pattern {!r}. Valid patterns: {}".format(pattern, valid))
    return PATTERN_ALIASES[key]


def default_n(pattern: str) -> int:
    """Default component number for a canonical pattern."""
    return PATTERN_DEFAULT_N[normalize_pattern(pattern)]


def validate_n(n: int, pattern: Optional[str] = None) -> None:
    """Require n >= 3, plus per-pattern caps for the exponential families."""
    if n < 3:
        raise ValueError("n must be at least 3 for these link series.")
    if pattern is not None:
        canonical = normalize_pattern(pattern)
        max_n = PATTERN_MAX_N.get(canonical)
        if max_n is not None and n > max_n:
            raise ValueError(
                "n must be at most {} for {} (its crossing number grows "
                "exponentially with n).".format(max_n, PATTERN_TITLES[canonical])
            )


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


# ---------------------------------------------------------------------------
# Braid-based families added in V2_2 (Edwards' Venn AM_n and Brunn's classic
# BR_n, braids from Maes & Cerf, JKTR 10(1), 2001).
# ---------------------------------------------------------------------------


def _inv(word: List[int]) -> List[int]:
    """Inverse of a braid word."""
    return [-x for x in reversed(word)]


def _shift(word: List[int]) -> List[int]:
    """Index-shift map S: s_i -> s_(i+1)."""
    return [x + 1 if x > 0 else x - 1 for x in word]


def edwards_venn_braid(n: int) -> List[int]:
    """Braid word for AM_n (Maes & Cerf, Lemma 2.3): alpha_n = u_n v_n,
    u_3 = (s1^-1 s2)^2, v_3 = s1^-1 s2,
    u_(k+1) = u_k S(v_k^-1) s1^-1 S(u_k^-1),  v_(k+1) = v_k s_k."""
    u, v = [-1, 2, -1, 2], [-1, 2]
    for k in range(3, n):
        u, v = u + _shift(_inv(v)) + [-1] + _shift(_inv(u)), v + [k]
    return u + v


def brunn_classic_braid(n: int) -> List[int]:
    """Braid word for BR_n (Maes & Cerf, Fig. 7): beta_3 = x3 y3,
    x3 = s1^-1 s2 s2 s1 s1 s2^-1, y3 = s2^-1 s1^-1;
    x_(k+1) = x_k s_k^2 x_k^-1 y_k^-1 s_k^-1,  y_(k+1) = s_k^-1 y_k."""
    x, y = [-1, 2, 2, 1, 1, -2], [-2, -1]
    for k in range(3, n):
        x, y = x + [k, k] + _inv(x) + _inv(y) + [-k], [-k] + y
    return x + y


def _braid_closure_dt(word: List[int], n: int) -> DTCode:
    """DT code of the closed-braid diagram of an n-strand braid word.
    Convention: component C_i = braid strand i, all oriented in the braid
    direction; traversal order C_1, ..., C_n; basepoint of C_i at the left
    edge of the braid for odd i, one crossing-pass later for even i; an even
    label is positive iff its pass is an over-pass."""
    pos = list(range(1, n + 1))
    passes: Dict[int, List[Tuple[int, bool]]] = {s: [] for s in range(1, n + 1)}
    for cid, letter in enumerate(word):  # sigma_i: lower strand over
        i = abs(letter)
        a, b = pos[i - 1], pos[i]
        passes[a].append((cid, letter > 0))
        passes[b].append((cid, letter < 0))
        pos[i - 1], pos[i] = b, a
    label = 1
    cross: Dict[int, List[Tuple[int, bool]]] = {}
    ranges: List[Tuple[int, int]] = []
    for s in range(1, n + 1):
        seq = passes[s] if s % 2 else passes[s][1:] + passes[s][:1]
        ranges.append((label, label + len(seq) - 1))
        for cid, over in seq:
            cross.setdefault(cid, []).append((label, over))
            label += 1
    pair: Dict[int, int] = {}
    for (l1, o1), (l2, o2) in cross.values():
        if (l1 + l2) % 2 != 1:
            raise RuntimeError("DT labeling parity failed; unexpected braid word.")
        odd = l1 if l1 % 2 else l2
        even, even_over = (l2, o2) if l1 % 2 else (l1, o1)
        pair[odd] = even if even_over else -even
    return [
        tuple(pair[o] for o in range(a + (a + 1) % 2, b + 1, 2))
        for a, b in ranges
    ]


def edwards_venn_dt(n: int) -> DTCode:
    """DT code of the closed-braid (Venn-diagram) projection of AM_n."""
    validate_n(n, "edwards_venn")
    return _braid_closure_dt(edwards_venn_braid(n), n)


def brunn_classic_dt(n: int) -> DTCode:
    """DT code of the closed-braid (Rolfsen-picture) diagram of BR_n."""
    validate_n(n, "brunn_classic")
    return _braid_closure_dt(brunn_classic_braid(n), n)


# ---------------------------------------------------------------------------
# Fishtail bracelet (added in V2_2). Not a Brunnian link, but closely related.
# ---------------------------------------------------------------------------

# Both fishtail families share these 16 offsets and differ only in signs.
FISHTAIL_OFFSETS = [0, -34, -50, -16, 40, 70, 84, 50, -22, -52, -36, -6, 34, 68, 86, 56]

# V2_3: the TRUE (same-fold) fishtail -- the stitch whose grip-span-1
# reduction is the classical chain.  Verified against the diagram model for
# n = 3..7.
FISHTAIL_SIGNS = [1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, -1]

# The mirrored-rail stitch, emitted as 'fishtail' by V2_2 and earlier.
MIRROR_FISHTAIL_SIGNS = [-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1]


def _fishtail_block_dt(n: int, signs: List[int], pattern: str) -> DTCode:
    """Shared closed form: labels 1..32n; with w = 32i and
    <x> = ((x-1) mod 32n) + 1, band i's block is the fixed 16-entry
    signed-offset pattern, translated by 32 per band."""
    validate_n(n, pattern)
    modulus = 32 * n
    return [
        tuple(
            sign * (((32 * i + offset - 1) % modulus) + 1)
            for sign, offset in zip(signs, FISHTAIL_OFFSETS)
        )
        for i in range(n)
    ]


def fishtail_dt(n: int) -> DTCode:
    """DT code of the canonical 16n-crossing diagram of the perfectly closed
    Rainbow Loom fishtail bracelet with n bands (V2_3 corrected, same-fold)."""
    return _fishtail_block_dt(n, FISHTAIL_SIGNS, "fishtail")


def mirror_fishtail_dt(n: int) -> DTCode:
    """DT code of the mirrored-rail fishtail variant -- the code emitted as
    'fishtail' by BL_DT V2_2 and earlier. Kept for reproducibility."""
    return _fishtail_block_dt(n, MIRROR_FISHTAIL_SIGNS, "mirror_fishtail")


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


def edwards_venn_formula() -> str:
    return FORMULA_HEADER + """
Edwards' Venn AM_n (Maes & Cerf, JKTR 10(1), 2001), 3 <= n <= 14 here
Brunnian link whose projection is Edwards' rotationally symmetric Venn
diagram; the closed-braid diagram has 2^n - 2 crossings.

The DT code is computed from the braid word alpha_n = u_n v_n on n strands,
defined recursively (S shifts every s_i to s_(i+1)):

u_3 = (s1^-1 s2)^2,  v_3 = s1^-1 s2
u_(k+1) = u_k S(v_k^-1) s1^-1 S(u_k^-1)
v_(k+1) = v_k s_k

Convention for the closed braid: component C_i = braid strand i (C_1 is the
circular component), all oriented in the braid direction; traversal order
C_1, ..., C_n; basepoint of C_i at the left edge of the braid for odd i, one
crossing-pass later for even i; an even label is positive iff its pass is an
over-pass.

Block sizes (machine-verified n = 3..8): |C_1| = |C_2| = 2^(n-2),
|C_3| = 3*2^(n-4), |C_4| = 2^(n-3), and |C_(n-1)| = |C_n| = n - 1, summing
to 2^n - 2 crossings. No closed-form per-entry formula is known; the
recursion above is the generalized description.
"""


def brunn_classic_formula() -> str:
    return FORMULA_HEADER + """
Brunn's classic BR_n (Brunn 1892; Rolfsen p.67 Ex. 8), 3 <= n <= 14 here
(n-1) round circles plus one component W woven through them;
the closed-braid diagram has 3*2^(n-1) - 4 crossings.

The DT code is computed from the braid word beta_n = x_n y_n on n strands:

x_3 = s1^-1 s2 s2 s1 s1 s2^-1,  y_3 = s2^-1 s1^-1
x_(k+1) = x_k s_k^2 x_k^-1 y_k^-1 s_k^-1
y_(k+1) = s_k^-1 y_k

Convention: components W = strand 1 (the woven component), then circles
C_1..C_(n-1) = strands 2..n; traversal order W, C_1, ..., C_(n-1); basepoint
at the left edge of the braid for strands 1, 3, 5, ..., one crossing-pass
later for strands 2, 4, ...; an even label is positive iff its pass is an
over-pass.

Closed forms exist for the first circle blocks (nu2 = 2-adic valuation),
machine-verified for n = 3..8:

C_1: signs -,+,-,+,...; |a_1| = 4; steps alternate
     (4 + 2*min(nu2(t) + 1, n - 3), then 4)
C_2: signs +,-,-,+ repeating; |a_1| = 2; steps alternate
     (4, then 4 + 2*min(nu2(t) + 1, n - 3))
C_3 (n >= 5): signs -,-,+,+ repeating; |a_1| = 8; steps alternate
     (10 + 2*min(nu2(t) + 1, n - 4), then 10)

Block sizes (machine-verified n = 3..8): |W| = 3*2^(n-2) - 2,
|C_1| = 2^(n-2), and |C_j| = 2^(n-j) for j = 2, ..., n-1.
The woven block W follows from the braid recursion above.
"""


def fishtail_formula() -> str:
    return FORMULA_HEADER + """
Fishtail bracelet, n >= 3, m = 32n (16n crossings, 16 per band)
V2_3 CORRECTED CODE (same-fold wrap; V2_2 emitted the mirrored variant,
now available as pattern 'mirror_fishtail').

NOT a Brunnian link -- removing any one band never unlinks it -- but
closely related and interesting: a removal set unlinks the bracelet exactly
when it contains two CYCLICALLY ADJACENT bands.

Each band is an unknot; worked off the pegs, its two arcs wrap the spanning
strands of the next TWO bands with the SAME handedness on both rails. The
identifying property: at grip span 1 the same construction reduces to the
classical single chain (= cyclic_rubberband).

Closed form. Write w = 32i and <x> = ((x - 1) mod 32n) + 1. Band i's block
(i = 0, ..., n-1) is:

C_i = (  <w>,     <w-34>, -<w-50>, -<w-16>,
        -<w+40>, -<w+70>,  <w+84>,  <w+50>,
         <w-22>,  <w-52>, -<w-36>, -<w-6>,
         <w+34>,  <w+68>, -<w+86>, -<w+56> )

That is, the entire code is one fixed 16-entry pattern of signed offsets,
translated by 32 from each band to the next, with wrap-around modulo 32n.
"""


def mirror_fishtail_formula() -> str:
    return FORMULA_HEADER + """
Mirror fishtail, n >= 3, m = 32n (16n crossings, 16 per band)
This is the code BL_DT V2_2 (and earlier) emitted under the name
'fishtail'. It is the MIRRORED-rail stitch: the two rails wrap with
OPPOSITE handedness. It is a genuinely different link from the true
fishtail (n = 5 volumes 92.7683 vs 91.7463) and does not reduce to the
classical chain at grip span 1.

Same 16 offsets as the fishtail; 8 of the 16 signs differ.

Closed form. Write w = 32i and <x> = ((x - 1) mod 32n) + 1. Band i's block
(i = 0, ..., n-1) is:

C_i = ( -<w>,    -<w-34>, -<w-50>, -<w-16>,
        -<w+40>, -<w+70>,  <w+84>,  <w+50>,
         <w-22>,  <w-52>,  <w-36>,  <w-6>,
        -<w+34>, -<w+68>,  <w+86>,  <w+56> )
"""


DT_GENERATORS: Dict[str, Callable[[int], DTCode]] = {
    "cyclic_squares": cyclic_squares_dt,
    "cyclic_larks": cyclic_larks_dt,
    "cyclic_rubberband": cyclic_rubberband_dt,
    "linear_rubberband": linear_rubberband_dt,
    "edwards_venn": edwards_venn_dt,
    "brunn_classic": brunn_classic_dt,
    "fishtail": fishtail_dt,
    "mirror_fishtail": mirror_fishtail_dt,
}

FORMULA_GENERATORS: Dict[str, Callable[[], str]] = {
    "cyclic_squares": cyclic_squares_formula,
    "cyclic_larks": cyclic_larks_formula,
    "cyclic_rubberband": cyclic_rubberband_formula,
    "linear_rubberband": linear_rubberband_formula,
    "edwards_venn": edwards_venn_formula,
    "brunn_classic": brunn_classic_formula,
    "fishtail": fishtail_formula,
    "mirror_fishtail": mirror_fishtail_formula,
}

CROSSING_EXPRESSIONS = {
    "cyclic_squares": "6n",
    "cyclic_larks": "6n",
    "cyclic_rubberband": "8n",
    "linear_rubberband": "8n - 16",
    "edwards_venn": "2^n - 2",
    "brunn_classic": "3*2^(n-1) - 4",
    "fishtail": "16n",
    "mirror_fishtail": "16n",
}

CROSSING_COUNTS: Dict[str, Callable[[int], int]] = {
    "cyclic_squares": lambda n: 6 * n,
    "cyclic_larks": lambda n: 6 * n,
    "cyclic_rubberband": lambda n: 8 * n,
    "linear_rubberband": lambda n: 8 * n - 16,
    "edwards_venn": lambda n: 2**n - 2,
    "brunn_classic": lambda n: 3 * 2 ** (n - 1) - 4,
    "fishtail": lambda n: 16 * n,
    "mirror_fishtail": lambda n: 16 * n,
}


def format_component(component: Sequence[int]) -> str:
    """Format one DT component tuple without spaces inside the tuple."""
    return "(" + ",".join(str(value) for value in component) + ")"


def format_dt_code(components: DTCode) -> str:
    """Format a full DT code in the compact style used in the notes."""
    return "DT: [{}]".format(", ".join(format_component(component) for component in components))


def pattern_note(pattern: str) -> Optional[str]:
    """Return the informational note for a canonical pattern, if any."""
    return PATTERN_NOTES.get(normalize_pattern(pattern))


def format_formula(pattern: str, n: Optional[int] = None) -> str:
    """Combine a pattern formula, crossing count, and informational note."""
    canonical = normalize_pattern(pattern)
    formula = FORMULA_GENERATORS[canonical]().rstrip()
    crossing_text = "Number of crossings: c(n) = {}".format(CROSSING_EXPRESSIONS[canonical])
    if n is not None:
        validate_n(n, canonical)
        crossing_text += "\nFor n = {}: c({}) = {}".format(n, n, CROSSING_COUNTS[canonical](n))
    sections = [formula, crossing_text]
    note = PATTERN_NOTES.get(canonical)
    if note:
        sections.append(note)
    return "\n\n".join(sections)


def generate_output(pattern: str, n: int) -> str:
    """Generate formula information followed by DT code for a pattern and n."""
    canonical = normalize_pattern(pattern)
    validate_n(n, canonical)
    formula = format_formula(canonical, n)
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
    root.title("{} {}: link-family DT code generator".format(APP_NAME, VERSION))
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

    n_var = tk.StringVar(value=str(PATTERN_DEFAULT_N["cyclic_squares"]))
    n_entry = tk.Entry(root, textvariable=n_var, width=10)
    n_entry.grid(row=0, column=3, sticky="w", padx=10, pady=(10, 4))

    formula_label = tk.Label(root, text="Generalized formula (selectable)")
    formula_label.grid(row=1, column=0, columnspan=4, sticky="w", padx=10, pady=(10, 4))

    formula_text = ScrolledText(root, height=28, wrap="word", font=("Courier New", 11))
    formula_text.grid(row=2, column=0, columnspan=4, sticky="nsew", padx=10, pady=(0, 8))

    preview_frame = tk.Frame(root)
    preview_frame.grid(row=1, column=4, rowspan=4, sticky="nsew", padx=(0, 10), pady=(10, 6))

    preview_title_var = tk.StringVar(value="")
    preview_title = tk.Label(
        preview_frame,
        textvariable=preview_title_var,
        anchor="center",
        justify="center",
        wraplength=330,
    )
    preview_title.pack(fill="x", pady=(0, 6))

    # Label width/height are measured in text units, so use a frame to reserve
    # a stable pixel-sized preview area when an image is missing.
    preview_image_frame = tk.Frame(preview_frame, bd=1, relief="solid", width=330, height=330)
    preview_image_frame.pack()
    preview_image_frame.pack_propagate(False)

    preview_image_label = tk.Label(preview_image_frame, anchor="center", justify="center", wraplength=290)
    preview_image_label.pack(fill="both", expand=True)

    dt_label = tk.Label(root, text="Generated DT code (auto-updates and selectable)")
    dt_label.grid(row=3, column=0, columnspan=4, sticky="w", padx=10, pady=(10, 4))

    dt_text = ScrolledText(root, height=5, wrap="word", font=("Courier New", 11))
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
    root.grid_rowconfigure(2, weight=5)
    root.grid_rowconfigure(4, weight=1)

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
        preview_n = PATTERN_PREVIEW_COMPONENTS[pattern]
        preview_title_var.set("{} ({}-component example)".format(PATTERN_TITLES[pattern], preview_n))
        try:
            photo = load_pattern_photo(pattern)
            preview_state["photo"] = photo
            preview_image_label.configure(image=photo, text="")
        except Exception:
            preview_state["photo"] = None
            preview_image_label.configure(
                image="",
                text="{}-component snapshot unavailable".format(preview_n),
            )

    def update_output(show_popup: bool = False) -> None:
        try:
            pattern = normalize_pattern(pattern_var.get())
            update_preview(pattern)
            n = current_n()
            validate_n(n, pattern)
            formula = format_formula(pattern, n) + "\n"
            clear_and_insert(formula_text, formula)

            dt_code = format_dt_code(DT_GENERATORS[pattern](n)) + "\n"
            clear_and_insert(dt_text, dt_code)
            current_output["text"] = "{}\n{}".format(formula, dt_code)
            status_var.set("Updated: {}, n = {}".format(PATTERN_TITLES[pattern], n))
        except Exception as exc:
            current_output["text"] = ""
            clear_and_insert(dt_text, "Error: {}\n".format(exc))
            status_var.set("Fix the input to update the DT code.")
            if show_popup:
                messagebox.showerror("{} error".format(APP_NAME), str(exc))

    def schedule_update(*_args) -> None:
        job = pending_update.get("job")
        if job is not None:
            root.after_cancel(job)
        pending_update["job"] = root.after(180, update_output)

    def on_pattern_selected(_event=None) -> None:
        """Reset n to the selected pattern's default, then refresh."""
        try:
            pattern = normalize_pattern(pattern_var.get())
            n_var.set(str(PATTERN_DEFAULT_N[pattern]))
        except Exception:
            pass
        schedule_update()

    def save_output() -> None:
        if not current_output["text"]:
            update_output(show_popup=True)
        if not current_output["text"]:
            return
        try:
            n = current_n()
            default_name = os.path.basename(default_output_path(pattern_var.get(), n))
        except Exception:
            default_name = "{}_output.txt".format(APP_NAME)
        path = filedialog.asksaveasfilename(
            title="Save {} output".format(APP_NAME),
            initialfile=default_name,
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
        )
        if path:
            try:
                write_text(path, current_output["text"])
                messagebox.showinfo(APP_NAME, "Saved to {}".format(path))
            except Exception as exc:
                messagebox.showerror("{} error".format(APP_NAME), str(exc))

    refresh_button = tk.Button(button_frame, text="Refresh", command=lambda: update_output(show_popup=True))
    refresh_button.pack(side="left", padx=(0, 8))

    save_button = tk.Button(button_frame, text="Save output...", command=save_output)
    save_button.pack(side="left")

    pattern_combo.bind("<<ComboboxSelected>>", on_pattern_selected)
    n_entry.bind("<Return>", lambda _event: update_output(show_popup=True))
    n_var.trace_add("write", schedule_update)

    update_output()
    root.mainloop()


def build_parser() -> argparse.ArgumentParser:
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Generate DT codes for eight link-construction families: four "
            "Brunnian-link series plus Edwards' Venn (AM_n), Brunn's classic "
            "(BR_n), the Fishtail bracelet, and the Mirror fishtail. "
            "Version: {}".format(VERSION)
        )
    )
    parser.add_argument("--version", action="version", version="{} {}".format(APP_NAME, VERSION))
    parser.add_argument(
        "--pattern",
        "-p",
        help=(
            "Pattern name or alias: cyclic_squares, cyclic_larks, "
            "cyclic_rubberband, linear_rubberband, edwards_venn, "
            "brunn_classic, fishtail, or mirror_fishtail."
        ),
    )
    parser.add_argument(
        "--n",
        "-n",
        type=int,
        help=(
            "Number of components (>= 3). Default: 7 for the four original "
            "families, 5 for edwards_venn, brunn_classic, fishtail, and "
            "mirror_fishtail."
        ),
    )
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
            print("Use CLI mode, for example: python BL_DTv2_3.py --pattern fishtail --n 5", file=sys.stderr)
            return 1

    if args.pattern is None:
        parser.error("CLI mode requires --pattern, unless --gui is used.")

    try:
        canonical = normalize_pattern(args.pattern)
        n = args.n if args.n is not None else PATTERN_DEFAULT_N[canonical]
        output = generate_output(canonical, n)
        print(output, end="")
        if args.output:
            write_text(args.output, output)
    except Exception as exc:
        print("Error: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
