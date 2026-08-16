#!/bin/sh
# Compile the exact-basis / Peschel helpers with Nuitka.
# Precision is unchanged (mpmath). The .py files stay the record.
# Output is local-only: scripts/_nuitka/ (gitignored *.so).
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
OUT="$HERE/_nuitka"
mkdir -p "$OUT"
PY=${PYTHON:-python3}
echo "nuitka $($PY -m nuitka --version | head -1)"
echo "out $OUT"
# Compile mpmath into the helper modules. Do not follow numpy
# (already compiled).
for mod in m9_60_exact.py m9_62_pair_ent.py; do
    echo "compile $mod"
    $PY -m nuitka \
        --module \
        --output-dir="$OUT" \
        --include-package=mpmath \
        --nofollow-import-to=numpy \
        --nofollow-import-to=m9_60_exact \
        --nofollow-import-to=m9_62_pair_ent \
        --assume-yes-for-downloads \
        --remove-output \
        "$HERE/$mod"
done
echo "done"
ls -la "$OUT"
