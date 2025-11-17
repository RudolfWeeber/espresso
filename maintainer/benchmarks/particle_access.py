#
# Copyright (C) 2025 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import espressomd
import numpy as np
import argparse
import time

parser = argparse.ArgumentParser(
    description="Benchmark particle creation and property access.")
parser.add_argument("--particles_per_core", metavar="N", action="store",
                    type=int, default=1000, required=False,
                    help="Number of particles per core (default: 1000)")
parser.add_argument("--additional_properties", metavar="PROPS", action="store",
                    type=str, default="", required=False,
                    help="Comma-separated list of additional properties to test (e.g., 'mass,type,v')")
parser.add_argument("--n_iter", metavar="N", action="store",
                    type=int, default=100, required=False,
                    help="Number of iterations for timing (default: 10)")

args = parser.parse_args()

# Parse additional properties
additional_props = []
if args.additional_properties:
    additional_props = [p.strip()
                        for p in args.additional_properties.split(',')]

# System setup
system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.time_step = 0.01
system.cell_system.skin = 0.5

n_proc = system.cell_system.get_state()["n_nodes"]
n_part = n_proc * args.particles_per_core

print(f"Benchmark Configuration:")
print(f"  Particles per core: {args.particles_per_core}")
print(f"  Total particles: {n_part}")
print(f"  MPI ranks: {n_proc}")
print(f"  Iterations: {args.n_iter}")
print(f"  Additional properties: {
      additional_props if additional_props else 'None'}")
print()

# Test 1: Bulk particle creation without providing IDs
print("=" * 60)
print("Test 1: Bulk particle creation without IDs")
print("=" * 60)
timings = []
for i in range(args.n_iter):
    system.part.clear()
    positions = np.random.random((n_part, 3)) * system.box_l
    tick = time.time()
    system.part.add(pos=positions)
    tock = time.time()
    t = tock - tick
    timings.append(t)

avg = np.mean(timings)
std = np.std(timings)
print(f"Average: {avg:.6f} s ± {std:.6f} s ({n_part / avg:.0f} particles/s)")
print()

# Test 2: Bulk particle creation with linear IDs
print("=" * 60)
print("Test 2: Bulk particle creation with linear IDs")
print("=" * 60)
timings = []
for i in range(args.n_iter):
    system.part.clear()
    positions = np.random.random((n_part, 3)) * system.box_l
    ids = np.arange(n_part, dtype=int)
    tick = time.time()
    system.part.add(pos=positions, id=ids)
    tock = time.time()
    t = tock - tick
    timings.append(t)

avg = np.mean(timings)
std = np.std(timings)
print(f"Average: {avg:.6f} s ± {std:.6f} s ({n_part / avg:.0f} particles/s)")
print()

# Setup particles for property access tests
system.part.clear()
positions = np.random.random((n_part, 3)) * system.box_l
system.part.add(pos=positions)

# Test 3: Single particle property access
print("=" * 60)
print("Test 3: Single particle property access")
print("=" * 60)

# Base properties to test
base_properties_read = ['pos', 'pos_folded', 'q', 'f']
base_properties_write = ['pos', 'q']

# Add additional properties
all_props_read = base_properties_read + additional_props
all_props_write = base_properties_write + \
    [p for p in additional_props if p not in ['f', 'pos_folded']]

# Test single particle read
print("\nSingle particle property READ:")
for prop in all_props_read:
    timings = []
    for i in range(args.n_iter):
        tick = time.time()
        _ = getattr(system.part.by_id(0), prop)
        tock = time.time()
        timings.append(tock - tick)
    avg = np.mean(timings)
    std = np.std(timings)
    print(f"  {prop}: {avg * 1e6:.3f} µs ± {std *
          1e6:.3f} µs, {1 / avg} particles/s")

# Test single particle write
print("\nSingle particle property WRITE:")
for prop in all_props_write:
    timings = []
    # Read valid value from particle 0
    test_value = getattr(system.part.by_id(0), prop)

    for i in range(args.n_iter):
        tick = time.time()
        setattr(system.part.by_id(0), prop, test_value)
        tock = time.time()
        timings.append(tock - tick)
    avg = np.mean(timings)
    std = np.std(timings)
    print(f"  {prop}: {avg * 1e6:.3f} µs ± {std *
          1e6:.3f} µs, {1 / avg} particles/s")

print()

# Test 4: Slice property access (re-instantiating slice)
print("=" * 60)
print("Test 4: Slice property access (re-instantiating slice)")
print("=" * 60)

# Test slice read - re-instantiating slice on every access
print("\nSlice property READ (re-instantiating slice):")
for prop in all_props_read:
    timings = []
    for i in range(args.n_iter):
        tick = time.time()
        _ = getattr(system.part.all(), prop)
        tock = time.time()
        timings.append(tock - tick)
    avg = np.mean(timings)
    std = np.std(timings)
    print(f"  {prop}: {avg * 1e3:.3f} ms ± {std *
          1e3:.3f} ms ({n_part / avg:.0f} particles/s)")

# Test slice write - re-instantiating slice on every access
print("\nSlice property WRITE (re-instantiating slice):")
for prop in all_props_write:
    timings = []
    # Read valid value from particle 0 to get correct type/shape
    test_value = getattr(system.part.by_id(0), prop)
    # Prepare test array
    if prop in ['pos', 'v']:
        test_array = np.tile(test_value, (n_part, 1))
    elif prop in ['q', 'mass', 'charge', 'type']:
        test_array = np.full(n_part, test_value)
    else:
        test_array = np.full(n_part, test_value)

    for i in range(args.n_iter):
        tick = time.time()
        setattr(system.part.all(), prop, test_array)
        tock = time.time()
        timings.append(tock - tick)
    avg = np.mean(timings)
    std = np.std(timings)
    print(f"  {prop}: {avg * 1e3:.3f} ms ± {std *
          1e3:.3f} ms ({n_part / avg:.0f} particles/s)")

print()

# Test 5: Slice property access (cached slice)
print("=" * 60)
print("Test 5: Slice property access (cached slice)")
print("=" * 60)

# Test slice read - cached slice
print("\nSlice property READ (cached slice):")
for prop in all_props_read:
    timings = []
    for i in range(args.n_iter):
        particle_slice = system.part.all()
        # Warm-up access
        _ = getattr(particle_slice, prop)
        # Timed access
        tick = time.time()
        _ = getattr(particle_slice, prop)
        tock = time.time()
        timings.append(tock - tick)
    avg = np.mean(timings)
    std = np.std(timings)
    print(f"  {prop}: {avg * 1e3:.3f} ms ± {std *
          1e3:.3f} ms ({n_part / avg:.0f} particles/s)")

# Test slice write - cached slice
print("\nSlice property WRITE (cached slice):")
for prop in all_props_write:
    timings = []
    # Read valid value from particle 0 to get correct type/shape
    test_value = getattr(system.part.by_id(0), prop)
    # Prepare test array
    if prop in ['pos', 'v']:
        test_array = np.tile(test_value, (n_part, 1))
    elif prop in ['q', 'mass', 'charge', 'type']:
        test_array = np.full(n_part, test_value)
    else:
        test_array = np.full(n_part, test_value)

    for i in range(args.n_iter):
        particle_slice = system.part.all()
        # Warm-up write
        setattr(particle_slice, prop, test_array)
        # Timed write
        tick = time.time()
        setattr(particle_slice, prop, test_array)
        tock = time.time()
        timings.append(tock - tick)
    avg = np.mean(timings)
    std = np.std(timings)
    print(f"  {prop}: {avg * 1e3:.3f} ms ± {std *
          1e3:.3f} ms ({n_part / avg:.0f} particles/s)")

print()
print("=" * 60)
print("Benchmark completed")
print("=" * 60)
