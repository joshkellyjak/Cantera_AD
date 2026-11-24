#!/usr/bin/env bash
# Mechanical part of Cantera's AD conversion: replace the C++ type 'double' with
# 'CanteraDouble' in Cantera's own (tracked) C++ sources and headers.
#
# - Only whole words are replaced; 'doublereal', 'CanteraDouble', etc. are untouched,
#   so running the script twice changes nothing.
# - String and character literals are left alone (e.g. the signature strings in
#   Delegator.h such as "void(double)" are matched at run time).
# - Lines that are only a comment, and trailing // comments, are left alone.
# - The two headers that define CanteraDouble in terms of double, and the
#   third-party headers that SCons copies into include/cantera/ext, are excluded.
#
# Run from anywhere inside the repository, on a clean working tree.
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

files=$(git ls-files -- 'src/*.cpp' 'src/*.h' 'include/cantera/*.h' 'include/cantera/*.hpp' \
    ':!include/cantera/ext/*' \
    ':!include/cantera/base/ct_typedefs.h' \
    ':!include/cantera/base/ct_ad_functionality.h')

before=$(git diff --shortstat | tr -d '\n')
printf '%s\n' $files | xargs perl -pi -e '
    next if m{^\s*(//|/\*|\*(?:\s|/|\*|$))};  # comment lines (not e.g. "*this = ...")
    s{(\x27(?:[^\x27\\]|\\.)\x27|"(?:[^"\\]|\\.)*"|//.*$)|\bdouble\b}{defined $1 ? $1 : "CanteraDouble"}ge;
'
echo "rename_double.sh: $(printf '%s\n' $files | wc -l | tr -d ' ') files scanned"
git diff --shortstat
[ -z "$before" ] || echo "(note: the tree was not clean before running: $before)"
