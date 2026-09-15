#!/bin/zsh
#
# Build a patched ridgerunner with the extra symmetry groups and install it as
# ~/.local/bin/ridgerunner-sym, leaving the stock ~/.local/bin/ridgerunner in
# place.
#
# Keeping the patch as a diff rather than only as a binary is the point: the
# binary links against Homebrew dylibs by versioned path (libgsl.28 in
# particular), so a brew upgrade can stop it loading. When that happens, or
# after a macOS or Xcode change, re-run this script.
#
# Why a separate name: /opt/homebrew/bin sits ahead of ~/.local/bin on PATH,
# so a later `brew install ridgerunner` would shadow a patched binary of the
# same name without deleting it -- silently giving stock behaviour. Address the
# patched build explicitly instead:
#
#     tighten_link_xyz.py link.xyz --ridgerunner ~/.local/bin/ridgerunner-sym
#     export RIDGERUNNER=~/.local/bin/ridgerunner-sym
#
set -e

VERSION=2.3.1
PATCH="${0:a:h}/0001-extra-symmetry-groups.patch"
WORK="${TMPDIR:-/tmp}/rr-sym-build"
TARBALL="$HOME/Library/Caches/Homebrew/ridgerunner--${VERSION}.tar.gz"

if [[ ! -f "$TARBALL" ]]; then
  echo "Source tarball not found at:"
  echo "  $TARBALL"
  echo "Fetch ridgerunner ${VERSION} from https://www.jasoncantarella.com/ and put it there,"
  echo "or run: brew fetch designbynumbers/cantarellalab/ridgerunner"
  exit 1
fi

rm -rf "$WORK"
mkdir -p "$WORK"
tar xzf "$TARBALL" -C "$WORK"
cd "$WORK/ridgerunner-${VERSION}"

patch -p1 < "$PATCH"

# OpenBLAS is keg-only so it is off the default pkg-config path; plCurve and
# tsnnls live in ~/.local and ship no pkg-config files. Both points are
# documented in the tarball's own configsonoma.sh.
export PKG_CONFIG_PATH="$(brew --prefix openblas)/lib/pkgconfig:$PKG_CONFIG_PATH"
./configure --prefix="$WORK/install" \
  CPPFLAGS="-I$HOME/.local/include -I$(brew --prefix)/include" \
  LDFLAGS="-L$HOME/.local/lib -L$(brew --prefix)/lib"
make -j4

# Plain copy under a new name. Never `make install` into ~/.local here: that
# would overwrite the stock ridgerunner.
mkdir -p "$HOME/.local/bin"
# rm first. Copying over a binary that has recently been executed can leave
# macOS holding the old mapping, and the "updated" binary then runs as a
# no-op: no output, exit 0, and nothing to tell you the copy did not take.
rm -f "$HOME/.local/bin/ridgerunner-sym"
cp ./ridgerunner "$HOME/.local/bin/ridgerunner-sym"
chmod +x "$HOME/.local/bin/ridgerunner-sym"

echo
echo "Installed $HOME/.local/bin/ridgerunner-sym"
"$HOME/.local/bin/ridgerunner-sym" --help 2>&1 | grep -E "Symmetry=|SymmetryAxis|SymmetryRef" || true
