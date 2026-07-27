#!/usr/bin/env bash
# Usage: ghost_ab.sh <label>   (run from the build dir)
# Runs lj.py and p3m.py on 1/2/4/8 ranks after checking the machine is idle.
set -euo pipefail
label="${1:?need a label, e.g. baseline or async}"
build="$(pwd)"

# --- machine-idle gate (per project constraint) ---
read -r load1 _ < <(awk '{print $1}' /proc/loadavg | tr '\n' ' ')
ncpu="$(nproc)"
others="$(who | awk '{print $1}' | sort -u | grep -v "^${USER}$" | wc -l)"
if (( $(echo "${load1} > 2.0" | bc -l) )) || (( others > 0 )); then
  echo "REFUSING to benchmark: load1=${load1}, other users=${others}, ncpu=${ncpu}." >&2
  echo "Re-run when the machine is idle." >&2
  exit 3
fi

out="${build}/ghost_ab_${label}.csv"
: > "${out}"
for ranks in 1 2 4 8; do
  mpiexec --bind-to core -n "${ranks}" ./pypresso ../maintainer/benchmarks/lj.py \
    --particles_per_core=1000 --volume_fraction=0.50 --output="${out}"
  mpiexec --bind-to core -n "${ranks}" ./pypresso ../maintainer/benchmarks/p3m.py \
    --particles_per_core=1000 --volume_fraction=0.25 --output="${out}"
done
echo "wrote ${out}"
