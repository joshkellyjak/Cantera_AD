#!/usr/bin/env bash
# Build and run the AD gradient validation against forward- and reverse-mode
# installs of Cantera.
#
# Prerequisite (from the repository root), for MODE in forward and reverse:
#   scons build ad_mode=$MODE prefix=$PWD/build/install-ad-$MODE
#   scons install ad_mode=$MODE prefix=$PWD/build/install-ad-$MODE
#
# Usage: test/ad_validation/run_validation.sh [repo-root]
set -euo pipefail

ROOT="$(cd "${1:-$(dirname "$0")/../..}" && pwd)"
HERE="$ROOT/test/ad_validation"
OUT="$ROOT/build/ad_validation"
CXX="${CXX:-clang++}"
# CoDiPack headers are not installed with Cantera; use the same copy it was built with
CODI="$ROOT/ext/codipack/include"
mkdir -p "$OUT"

status=0
for MODE in forward reverse; do
    PREFIX="$ROOT/build/install-ad-$MODE"
    if [ ! -f "$PREFIX/include/cantera/base/config.h" ]; then
        echo "Missing $PREFIX; build and install Cantera with ad_mode=$MODE first" >&2
        exit 2
    fi
    # No -DAD_* flags on purpose: the AD mode must come from the installed config.h
    "$CXX" -std=c++20 -O2 -isystem "$PREFIX/include" -isystem "$CODI" \
        -isystem /opt/homebrew/include \
        "$HERE/gradient_validation.cpp" -o "$OUT/gradient_validation_$MODE" \
        -L"$PREFIX/lib" -lcantera_shared -Wl,-rpath,"$PREFIX/lib"

    export CANTERA_DATA="$HERE:$PREFIX/share/cantera/data"
    echo "################ $MODE mode ################"
    if [ "$MODE" = forward ]; then
        "$OUT/gradient_validation_$MODE" --write "$OUT/jacobians_forward.txt" \
            | tee "$OUT/report_$MODE.txt" || status=1
    else
        "$OUT/gradient_validation_$MODE" --compare "$OUT/jacobians_forward.txt" \
            | tee "$OUT/report_$MODE.txt" || status=1
    fi
    echo
done
exit $status
