#!/usr/bin/env python3
"""Make Cantera's source generator (CLib / Julia / C# bindings) work after
rename_double.sh has renamed the C++ type 'double' to 'CanteraDouble'.

1. In the C++ -> C type crosswalks (clib/config.yaml, headers/config.yaml), add a
   'CanteraDouble' twin for every C++ type containing 'double'. The C-side types
   stay 'double'.
2. In the header specifications (headers/*.yaml), rename 'double' in the fields
   that quote C++ signatures ('wraps:' and the items under 'uses:'). The
   'declaration:' fields are C declarations and must keep 'double'.
3. Two checks in the generators that test a C++ type name for 'double' (the
   CLib error sentinel DERR, and 'const' handling of primitive return types)
   also have to accept 'CanteraDouble'.

Idempotent. Run from the repository root:
    python3 sourcegen_ad.py
"""
import re
from pathlib import Path

ROOT = Path("interfaces/sourcegen/src/sourcegen")
WORD = re.compile(r"\bdouble\b")


def add_crosswalk_twins(path):
    out, section, added = [], None, 0
    lines = path.read_text().splitlines(keepends=True)
    existing = {l.split(":")[0].strip() for l in lines if ":" in l}
    for line in lines:
        m = re.match(r"^(\w+):", line)
        if m:
            section = m.group(1)
        out.append(line)
        k = re.match(r"^(\s+)([^#:\s][^:]*):(.*\n?)$", line)
        if section in ("ret_type_crosswalk", "par_type_crosswalk") and k \
                and WORD.search(k.group(2)):
            key = WORD.sub("CanteraDouble", k.group(2))
            if key.strip() not in existing:
                out.append(k.group(1) + key + ":" + k.group(3))
                added += 1
    path.write_text("".join(out))
    return added


def rename_cxx_signatures(path):
    out, in_uses, changed = [], False, 0
    for line in path.read_text().splitlines(keepends=True):
        stripped = line.lstrip()
        if re.match(r"^\s*uses:\s*$", line):
            in_uses = True
        elif re.match(r"^\s*-?\s*\w[\w-]*:", line) and not stripped.startswith("- "):
            in_uses = False
        is_wraps = re.match(r"^\s*-?\s*wraps:", line)
        is_use_item = in_uses and stripped.startswith("- ")
        if (is_wraps or is_use_item) and WORD.search(line):
            line = WORD.sub("CanteraDouble", line)
            changed += 1
        out.append(line)
    path.write_text("".join(out))
    return changed


GENERATOR_PATCHES = [
    (ROOT / "clib/generator.py",
     'cxx_rtype.endswith("double")',
     'cxx_rtype.endswith(("double", "CanteraDouble"))'),
    (ROOT / "headers/generator.py",
     '{"bool", "size_t", "int", "double"}',
     '{"bool", "size_t", "int", "double", "CanteraDouble"}'),
]


def patch_generators():
    done = 0
    for path, old, new in GENERATOR_PATCHES:
        text = path.read_text()
        if new in text:
            continue
        if text.count(old) != 1:
            raise SystemExit(f"{path}: expected exactly one '{old}'; update this script")
        path.write_text(text.replace(old, new))
        done += 1
    return done


if __name__ == "__main__":
    for cfg in ("clib/config.yaml", "headers/config.yaml"):
        print(f"{cfg}: {add_crosswalk_twins(ROOT / cfg)} crosswalk entries added")
    total = sum(rename_cxx_signatures(p) for p in sorted((ROOT / "headers").glob("ct*.yaml")))
    print(f"headers/ct*.yaml: {total} C++ signature lines renamed")
    print(f"generators: {patch_generators()} type checks updated")
