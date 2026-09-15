# Extra symmetry groups for RidgeRunner

A patch against **ridgerunner 2.3.1** adding inversion, Cnv and the real
dihedral groups to `--Symmetry`, plus an arbitrary symmetry axis.

`libplCurve` is **not** modified. Every group is built in
`src/ridgerunner_main.c` from public API already exported by the installed
library, so a patched RidgeRunner links against an unmodified plCurve.

## Build

```bash
./build_ridgerunner_sym.sh
```

Unpacks the cached source tarball, applies the patch, configures, builds, and
installs the result as `~/.local/bin/ridgerunner-sym`. The stock
`~/.local/bin/ridgerunner` is left untouched.

Point tooling at the patched build explicitly:

```bash
python3 tighten_link_xyz.py link.xyz --ridgerunner ~/.local/bin/ridgerunner-sym
export RIDGERUNNER=~/.local/bin/ridgerunner-sym
```

Explicitly, and not by relying on `PATH`: `/opt/homebrew/bin` sits ahead of
`~/.local/bin`, so a later `brew install ridgerunner` would shadow a patched
binary of the same name without deleting it, silently giving stock behaviour.

Keep the patch, not just the binary. The binary links Homebrew dylibs by
versioned path (`libgsl.28` in particular), so a `brew upgrade` can stop it
loading; re-run the build script when that happens.

## New groups

| `--Symmetry=` | Group | Order | Needs |
|---|---|---|---|
| `Ci` | inversion, `{E, -I}` | 2 | — |
| `Cs` | one mirror | 2 | axis = mirror normal |
| `Cpv` (`C2v`, `C3v`, …) | p-fold rotation + p mirrors **containing** the axis | 2p | axis, ref |
| `RDp` (`RD2`, `RD3`, …) | p-fold axis + p 2-fold axes **perpendicular** to it | 2p | axis, ref |
| `Z/pZ` or `Cp` | p-fold rotation — the two spellings are the same group | p | axis |
| `cplanes` | reflections in all three coordinate planes (pre-existing) | 4 | — |
| `D2` | deprecated alias for `Cs` | 2 | axis |

```
--SymmetryAxis=x,y,z    default 0,0,1
--SymmetryRef=x,y,z     default 1,0,0
```

`Cp` is an alias for `Z/pZ`, verified numerically identical (same error, same
ropelength at step 60 on the Borromean rings). It exists because `C3` is the
natural Schoenflies spelling and must not be confusable with `C3v`, which differs
by three mirror planes.

`--SymmetryAxis` applies to `Z/pZ` and `D2` as well, which previously always
used z and so required the input to be pre-aligned.

**One backward-compatibility regression.** argtable2 accepts unambiguous option
prefixes, and every prefix of `--Symmetry` is also a prefix of
`--SymmetryAxis`/`--SymmetryRef`, so abbreviated spellings that worked on stock
(`--Sym`, `--Symm`, `--Symmetr`) now fail. The full `--Symmetry=` spelling — what
every script and `tighten_link_xyz.py` emit — is unaffected. `--SymmetryRef` fixes where
the first mirror plane (`Cpv`) or first 2-fold axis (`RDp`) sits, and is
projected perpendicular to the axis; nothing else uses it.

`symmetrize_link_xyz.py --group Cnv` aligns its output to exactly these
defaults — rotation axis on z, first mirror normal on x — so the two tools
compose without any flags.

### On the name RDn

RidgeRunner's own `D2` calls `plc_reflection_group` and yields a **single
mirror**: Cs, order 2. That is not the D2 of molecular symmetry, which is three
mutually perpendicular 2-fold axes, order 4, no mirrors. `RD2` is the genuine
one, hence the `R`. `D2` still works, and now prints a note saying what it
actually builds.

## Other fixes in the same patch

**A `printf` format bug.** The old `D2` and `cplanes` branches passed an
uninitialised `int p` where the format string expected one `%g`, so the double
was consumed as an extra argument and the int was reinterpreted as a double.
That is why stock RidgeRunner reports `initial error 1.16911e-315` for those
groups — a denormal from misread bits, not a measurement. Both now report the
real error.

**A segfault on degenerate orders.** `plc_symmetry_group_new` returns NULL for
`n <= 0`, and `plc_rotation_group` then writes `build->sym[0]` without checking,
so **stock RidgeRunner segfaults on `--Symmetry=Z/0Z`**. The real fix belongs in
plCurve; the patch validates the order first, turning the crash into a message.

**`symmetryerror.dat` was never written.** The logfile was created from the
name table but nothing ever wrote a line to it, so a symmetry run left no record
of whether the symmetry held. It is now logged every logged step, which costs
nothing when no group is set (`plc_symmetry_group_check` returns 0 immediately).

**Malformed group specs are rejected.** `--Symmetry=Z/2Z,D2` used to parse as
`Z/2Z` and silently drop the rest, reading as though both had been applied.
Specs are now matched with `%n` and required to consume the whole string, so
`Z/2`, `Z/2Z,D2`, `C3` (which would otherwise build C3v), `C3h`, `RD2x` and
`C2vx` are all rejected. This is the one intentional behaviour change for
previously-accepted input.

A trailing `%c` is *not* sufficient here, and was the first attempt: on `C3` the
`%d` is assigned and only then does the literal `v` fail to match, so `sscanf`
still returns 1 and a naive `== 1` test accepts the prefix. `%n` is not counted
in the return value and is reached only if every preceding literal matched, so
requiring it to land on the terminating NUL rejects a short prefix as well as a
long tail.

**A guard against an upstream validation bug.**
`plc_coordplanes_reflection_group` checks only the first two of its four
elements for a failed build (its loop runs `i<2` over a group of order 4), so it
can return a group carrying NULL entries that crash whatever walks it next.
Every group is now checked properly before use.

**An unbounded order parameter.** `--Symmetry=C1000000000v` computed `2*p`
and walked into plCurve's unchecked `calloc`s, hanging the process until killed.
The order is now required to be no larger than the vertex count, since a group
cannot permute more elements than the curve has vertices.

**Non-finite axis components.** `--SymmetryAxis=inf,0,1` passed the magnitude
check (`plc_norm(inf) > 1e-12`) and produced a NaN transform that failed later,
blaming the wrong thing. Now rejected at parse time.

## Two things worth knowing before trusting a run

**`plc_build_symmetry` has no distance tolerance.** It maps each vertex with
the bare matrix, takes the nearest vertex to the image, and fails *only* when
two vertices claim the same target. A curve far from symmetric can therefore
produce a perfectly injective but entirely wrong vertex map, which
`plc_symmetrize` then bakes in. A successful build is not evidence the symmetry
is real. The patch measures the error **before** symmetrizing and warns when it exceeds
half the mean edge length, which is the validity condition plCurve's own comment
states. Measuring after is useless, and was the first attempt: `plc_symmetrize`
forces every orbit into agreement, so the post-symmetrization error is ~1e-16 by
construction and the warning could never fire. Both numbers are now printed —
`Error before symmetrizing X, after Y` — and only X carries information.

The danger band is real and narrow. On a 200-vertex test curve whose C2 symmetry
was broken smoothly, a break of 0.02–0.06 still *built* a group (no vertex
collision) while being off symmetry by 0.04–0.12, well past half an edge; only at
0.12 did the build finally fail. Inside that band the warning is the only thing
standing between you and a silently wrong result. Symmetrize the input first with
`symmetrize_link_xyz.py`.

**`--Symmetry` forces `--AnimationStepper`**, which the man page says "is not as
effective as the default stepper in producing configurations with very low
ropelength and residual". That coupling predates this patch and is untouched by
it. Treat symmetry-constrained runs as a way to get a well-conditioned
configuration, and re-run unconstrained with the default stepper for a
publishable ropelength.

Also: `Ci` inverts through the **origin**, not the centroid, because
`plc_build_symmetry` applies the bare matrix. Centre the input first.
Autoscaling is safe — `plc_scale` multiplies coordinates, so it scales about
the origin and a centred curve stays centred.

## Verified

Built from a clean tarball; patch applies with no fuzz. Behaviour on the test
set:

| Case | Result |
|---|---|
| `Ci` on a Ci-symmetrized 4-component link | order 2, error 0 |
| `RD2` on a curve with exact RD2 | order 4, error 6.7e-16 |
| `C2v` on a curve with exact C2v | order 4, error 6.3e-16 |
| `Cs` on a Cs-symmetrized link | order 2, error 0 |
| `Z/3Z` on Borromean rings, axis `1,1,1` | order 3, error 2.1e-15 |
| `Z/3Z` on the same, default z axis | correctly refuses |
| `RD2` with the ref 45° off | correctly refuses |
| `C3v` on a curve with only 2 mirrors | correctly refuses |
| `Ci` on a structure with no inversion | correctly refuses |
| `Z/0Z`, `Z/-2Z` | clean error (segfault in stock) |
| `C0v`, `C1v`, `RD0`, `RD1` | clean error |
| unknown name; `Z/2`, `C3`, `C3h`, `RD2x`, `C2vx`, `Z/2Z,D2` | clean error listing what is accepted |
| smoothly-broken C2, break 0.02 and 0.06 | builds, and the half-edge warning fires |
| `C3` vs `Z/3Z` on the Borromean rings | identical error and ropelength |
| `C1000000000v` (previously hung) | clean error naming the vertex count |
| `--SymmetryAxis=inf,0,1` / `nan,0,1` | rejected at parse time |

600-step runs under `Ci` (order 2) and `C2v` (order 4) held their symmetry for
the whole run — `symmetryerror` stayed at 0 and ~3e-16 respectively — and the
tightened output was independently re-checked at 7.6e-17 with the component
swap map preserved.

---

# Handoff prompt

Paste this into a fresh Claude Code session working in this repo.

```text
I have a patched RidgeRunner with extra symmetry groups. Context you need:

WHAT EXISTS
- ridgerunner_patches/0001-extra-symmetry-groups.patch — the diff against
  ridgerunner 2.3.1 (src/ridgerunner_main.c, src/stepper.c, src/ridgerunner.h)
- ridgerunner_patches/build_ridgerunner_sym.sh — applies it and installs the
  result as ~/.local/bin/ridgerunner-sym
- ridgerunner_patches/README.md — full notes, caveats and test results
- ~/.local/bin/ridgerunner-sym — the patched binary (built and tested)
- ~/.local/bin/ridgerunner — the STOCK binary, deliberately untouched
- libplCurve is NOT patched. Every group is built in ridgerunner_main.c from
  public plCurve API, so the patched binary links an unmodified library.

WHAT THE PATCH ADDS
  --Symmetry=Cp          p-fold rotation — IDENTICAL to Z/pZ (C3 == Z/3Z)
  --Symmetry=Ci          inversion {E,-I}, order 2
  --Symmetry=Cpv         p-fold axis + p mirrors CONTAINING it, order 2p (C2v, C3v, ...)
  --Symmetry=RDp         p-fold axis + p PERPENDICULAR 2-fold axes, order 2p (RD2, RD3, ...)
  --Symmetry=Cs          one mirror, order 2
  --SymmetryAxis=x,y,z   default 0,0,1; also frees the pre-existing Z/pZ and D2
  --SymmetryRef=x,y,z    default 1,0,0; fixes the azimuth of the first mirror
                         (Cpv) or first 2-fold axis (RDp); nothing else uses it
Pre-existing Z/pZ, D2 and cplanes are unchanged in behaviour — stock and patched
were verified to give identical ropelength on the same input. ONE regression:
argtable2 prefix abbreviations of --Symmetry (--Sym, --Symm, --Symmetr) now fail,
because every prefix of --Symmetry is also a prefix of --SymmetryAxis. The full
--Symmetry= spelling is unaffected.

Cp and Z/pZ are the same group. C3 used to silently build C3v (order 6, three
mirror planes) because sscanf returns 1 on a prefix match; specs are now matched
with %n and must consume the whole string.

WHY "RDp" AND NOT "Dp"
RidgeRunner's own --Symmetry=D2 calls plc_reflection_group and yields a SINGLE
MIRROR: Cs, order 2. That is not molecular D2 (three perpendicular 2-fold axes,
order 4, no mirrors). RD2 is the real one. D2 still works and now prints a note
saying what it actually builds.

CRITICAL CAVEATS — read before trusting any symmetry run
1. plc_build_symmetry has NO distance tolerance. It maps each vertex with the
   bare matrix, takes the NEAREST vertex to the image, and fails ONLY when two
   vertices claim the same target. A curve far from symmetric can produce a
   perfectly injective but entirely WRONG vertex map, which plc_symmetrize then
   bakes in. A successful build is NOT evidence the symmetry is real. The patch
   prints "Error before symmetrizing X, after Y" and warns when X exceeds half
   the mean edge length. Only X carries information; Y is ~1e-16 by construction.
   There is a real band (tested: symmetry break 0.02-0.06 on a 200-vertex curve)
   where it builds a wrong map and only that warning catches it.
2. --Symmetry forces --AnimationStepper and silently disables --EqOn. Per the
   man page the animation stepper "is not as effective as the default stepper in
   producing configurations with very low ropelength and residual". This coupling
   predates the patch. Use symmetry runs to get a well-conditioned configuration;
   re-run unconstrained with the default stepper for a publishable ropelength.
3. Ci inverts through the ORIGIN, not the centroid, because plc_build_symmetry
   applies the bare matrix. Centre the input first. Autoscale is safe: plc_scale
   multiplies coordinates, so it scales about the origin.

COMPANION TOOL
symmetrize_link_xyz.py projects an .xyz link onto Cs / Ci / Cn / Cnv and writes a
canonical frame: rotation axis on z, and for Cnv the first mirror normal on x.
That matches the --SymmetryAxis/--SymmetryRef defaults exactly, so its output
feeds ridgerunner-sym with no flags. Use --no-permute when the symmetry maps each
component to itself, and NOT when it exchanges components (inversion often does;
two components with near-equal arclength are the tell).

HOW TO RUN
  python3 tighten_link_xyz.py link.xyz --symmetry Ci \
      --ridgerunner ~/.local/bin/ridgerunner-sym
tighten_link_xyz.py also has --symmetry-axis and --symmetry-ref, which it only
passes through when set (stock ridgerunner would reject them).

HOW TO REBUILD (needed after a brew upgrade breaks the dylib links)
  ./ridgerunner_patches/build_ridgerunner_sym.sh
OpenBLAS is keg-only so PKG_CONFIG_PATH needs its keg; plCurve and tsnnls live in
~/.local and ship no pkg-config files. The script handles both. Never `make
install` into ~/.local — that would overwrite the stock binary. Address the
patched build explicitly rather than via PATH: /opt/homebrew/bin sits ahead of
~/.local/bin, so a later `brew install ridgerunner` would shadow a same-named
patched binary and silently give stock behaviour.

UPSTREAM BUGS THE PATCH WORKS AROUND (all confirmed against the real source)
- The old D2/cplanes branches passed an uninitialised int where the format
  wanted %g, so their reported "initial error" was a misread bit pattern
  (always 1.16911e-315). Fixed.
- plc_symmetry_group_new returns NULL for n<=0 and plc_rotation_group writes
  build->sym[0] unchecked, so STOCK ridgerunner SEGFAULTS on --Symmetry=Z/0Z.
  The patch validates the order first. A proper fix belongs in plCurve.
- plc_coordplanes_reflection_group validates only 2 of its 4 elements (loop runs
  i<2 over order 4), so it can return a group with NULL entries. The patch
  checks every group properly before use. That element set is also NOT closed —
  the product of two of its mirrors is a half turn outside the set — so
  averaging over it is not a projection. Prefer RD2 or C2v.
- symmetryerror.dat was created from the logfile name table but never written.
  Now logged every logged step; free when no group is set. Caveat: it is sampled
  just after plc_symmetrize, so it reads ~0 throughout and confirms the group is
  being applied rather than measuring drift within a step.
- Unbounded p hung the process (C1000000000v); now bounded by the vertex count.
- Non-finite --SymmetryAxis components produced a NaN transform; now rejected.

GOTCHA WHEN REINSTALLING
rm the target before cp-ing a rebuilt binary over it. Copying over a recently
executed binary can leave macOS holding the old mapping, and the "updated"
binary then runs as a no-op: no output, exit 0, nothing to tell you. The build
script does this.

STATUS
Nothing is committed. ridgerunner_patches/ and symmetrize_link_xyz.py are
untracked; tighten_link_xyz.py is modified.
```
