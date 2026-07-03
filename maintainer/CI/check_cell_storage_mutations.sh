#!/usr/bin/env bash
# Guard for the ParticleStore migration (see docs/superpowers/specs/
# 2026-07-03-array-based-particle-storage-design.md, section 3, phase 1):
# cell particle storage may only be mutated through the primitives in
# src/core/cell_system/ParticleListOperations.hpp. This script fails if a
# direct mutation pattern reappears in src/core (unit tests excluded:
# they construct fixtures directly).
set -u
cd "$(dirname "$0")/../.."
matches=$(grep -rn --include='*.cpp' --include='*.hpp' \
    -e 'particles()\.insert(' \
    -e 'particles()\.erase(' \
    -e 'particles()\.clear(' \
    -e 'particles()\.resize(' \
    src/core \
    | grep -v 'src/core/unit_tests/' \
    | grep -v 'src/core/cell_system/ParticleListOperations.hpp')
if [ -n "${matches}" ]; then
    echo "ERROR: direct cell-storage mutation outside CellParticleStorage:"
    echo "${matches}"
    exit 1
fi
echo "OK: no direct cell-storage mutations found."
exit 0
