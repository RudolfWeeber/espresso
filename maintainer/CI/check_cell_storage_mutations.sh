#!/usr/bin/env bash
# Guard for the ParticleStore migration (see docs/superpowers/specs/
# 2026-07-03-array-based-particle-storage-design.md, section 3, phase 1):
# cell particle storage may only be mutated through the primitives in
# src/core/cell_system/ParticleListOperations.hpp. This script fails if a
# direct mutation pattern reappears in src/core (unit tests excluded:
# they construct fixtures directly).
#
# The scan is delegated to an embedded Python program because the mutation
# call can be split across several physical lines by clang-format, e.g.
#     object.get_local_cells()[index]
#         ->particles()
#         .insert(std::move(p));
# A line-based grep misses this; the Python scanner reads each file as one
# string and matches the pattern across newlines.
#
# ALIAS LIMITATION: this check cannot see mutations performed through a
# local reference to a cell's storage, e.g.
#     auto &parts = cell->particles();
#     parts.erase(it);
# because the mutating call no longer names particles(). A green run is a
# tripwire that catches the common, direct form; it is NOT a proof that no
# direct cell-storage mutation exists. Reviewers must still watch for
# aliased mutations by hand.
#
# EXCEPTION (decomposition swap/teardown): CellStructure::set_particle_
# decomposition re-inserts particles through the routed add_particle and
# then destroys the old ParticleDecomposition wholesale; the old cells'
# particle lists are torn down by their destructors, not by these
# primitives. This bulk teardown is deliberately outside the hook and is
# handled by migration phase 2 via a full-store rebuild (mark-dirty), not
# per-row hooks. See ParticleListOperations.hpp for details.
set -u
cd "$(dirname "$0")/../.."
exec python3 - "$@" << 'PYEOF'
import os
import re
import sys

ROOT = "src/core"
EXCLUDED_DIRECTORIES = ("src/core/unit_tests/",)
EXCLUDED_FILES = ("src/core/cell_system/ParticleListOperations.hpp",)

# Match particles() followed by a mutating member call, allowing arbitrary
# whitespace (including newlines) around the '.' or '->' operator so that
# clang-format-wrapped calls are caught.
PATTERN = re.compile(
    r"particles\(\)\s*(?:->|\.)\s*(?:insert|erase|clear|resize)\s*\(")


def normalized(path):
    return path.replace(os.sep, "/")


def offending_files():
    findings = []
    for dirpath, _, filenames in os.walk(ROOT):
        for filename in filenames:
            if not filename.endswith((".cpp", ".hpp")):
                continue
            path = normalized(os.path.join(dirpath, filename))
            if path.startswith(EXCLUDED_DIRECTORIES):
                continue
            if path in EXCLUDED_FILES:
                continue
            with open(path, encoding="utf-8", errors="replace") as handle:
                text = handle.read()
            for match in PATTERN.finditer(text):
                line_number = text.count("\n", 0, match.start()) + 1
                snippet = " ".join(match.group(0).split())
                findings.append((path, line_number, snippet))
    return findings


def main():
    findings = offending_files()
    if findings:
        print("ERROR: direct cell-storage mutation outside "
              "CellParticleStorage:")
        for path, line_number, snippet in findings:
            print(f"  {path}:{line_number}: {snippet}")
        return 1
    print("OK: no direct cell-storage mutations found.")
    return 0


sys.exit(main())
PYEOF
