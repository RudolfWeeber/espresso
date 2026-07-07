#!/usr/bin/env bash
# Guard for the ParticleStore migration (see docs/superpowers/specs/
# 2026-07-03-array-based-particle-storage-design.md, section 3):
# cell particle storage may only be mutated through the primitives in
# src/core/cell_system/ParticleListOperations.hpp. This script fails if a
# direct mutation pattern reappears in src/core (unit tests excluded:
# they construct fixtures directly).
#
# Phase 7a: cells no longer own a Bag<Particle>; a Cell holds ParticleStore ROW
# indices (Cell::rows()) plus a staging buffer of not-yet-committed detached
# particles (Cell::staged()). The choke points now stage inserts, snapshot +
# remove rows on extract, clear both, and stage ghost slots on resize. Direct
# mutation of rows()/staged() (and, defensively, the legacy particles()) outside
# the choke points and the sanctioned store rebuild is the out-of-band-mutation
# hazard this guard catches.
#
# The scan is delegated to an embedded Python program because the mutation
# call can be split across several physical lines by clang-format, e.g.
#     object.get_local_cells()[index]
#         ->rows()
#         .insert(row);
# A line-based grep misses this; the Python scanner reads each file as one
# string and matches the pattern across newlines.
#
# ALIAS LIMITATION: this check cannot see mutations performed through a
# local reference to a cell's storage, e.g.
#     auto &rows = cell->rows();
#     rows.erase(it);
# because the mutating call no longer names rows(). A green run is a tripwire
# that catches the common, direct form; it is NOT a proof that no direct
# cell-storage mutation exists. Reviewers must still watch for aliased mutations
# by hand.
#
# EXCEPTIONS: (1) ParticleListOperations.hpp defines the choke points.
# (2) CellStructure.cpp's ensure_particle_store_synchronized owns the store
# rebuild and legitimately clears + refills every cell's row bag as it renumbers
# rows. Both files are excluded below.
set -u
cd "$(dirname "$0")/../.."
exec python3 - "$@" << 'PYEOF'
import os
import re
import sys

ROOT = "src/core"
EXCLUDED_DIRECTORIES = ("src/core/unit_tests/",)
# Choke points (ParticleListOperations.hpp) legitimately mutate the row bag /
# staging area; the store rebuild in CellStructure.cpp legitimately clears and
# refills every cell's row bag as it renumbers rows. Both are the sanctioned
# owners of cell storage and are exempt.
EXCLUDED_FILES = (
    "src/core/cell_system/ParticleListOperations.hpp",
    "src/core/cell_system/CellStructure.cpp",
)

# Phase 7a: cells hold ParticleStore ROW indices (Cell::rows()) plus a staging
# buffer of not-yet-committed particles (Cell::staged()); they no longer own a
# Bag<Particle>. Direct mutation of either -- outside the CellParticleStorage
# choke points and the store rebuild -- is the new out-of-band-mutation hazard.
# Match rows()/staged() (and, defensively, the legacy particles()) followed by a
# mutating member call, allowing arbitrary whitespace (including newlines) around
# the '.'/'->' operator so clang-format-wrapped calls are caught.
PATTERN = re.compile(
    r"(?:particles|rows|staged)\(\)\s*(?:->|\.)\s*"
    r"(?:insert|erase|clear|resize|push_back|emplace_back|pop_back|emplace)"
    r"\s*\(")


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
