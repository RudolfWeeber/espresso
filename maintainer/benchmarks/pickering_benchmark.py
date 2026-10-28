#
# Copyright (C) 2026 The ESPResSo project
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

"""
Performance benchmark: magnetic Pickering emulsion in a two-component color gradient LB fluid.

A droplet (radius 12.5 lattice units) sits in the centre of a cubic box. 312 magnetic
particles sit on its interface in a homogeneous magnetic field. Each particle is a rigid
tripod (anchor + two solvation virtual sites straddling the interface) carrying a
free-spinning dipole virtual site. The droplet is the same for every box size, so a larger
box only adds fluid.

After 1000 warm-up steps without field, the field is switched on and 10 blocks of 100 steps
are timed. The mean time per step with its 95% confidence interval and the LB throughput
(MLUPS, million lattice node updates per second) are printed and appended to
benchmarks.csv in the current directory (same columns as the ESPResSo benchmarks in
maintainer/benchmarks; the label column says cpu or gpu).

Requirements
------------
ESPResSo built from the color gradient branch
    https://github.com/RebeccaStephan/espresso  (branch color_gradient, tested at 82e6d3d92)
with the features WALBERLA, DIPOLES, ROTATION, VIRTUAL_SITES_RELATIVE and WCA,
plus CUDA for --gpu. Only numpy is needed besides ESPResSo.

How to run
----------
LB on the CPU, with N MPI ranks (OpenMP fixed to one thread):
    OMP_NUM_THREADS=1 mpiexec -n N /path/to/espresso/build/pypresso pickering_benchmark.py --box_l 50

LB on the GPU (needs a CUDA build; one rank is fastest, the particles stay on the CPU):
    OMP_NUM_THREADS=1 /path/to/espresso/build_cuda/pypresso pickering_benchmark.py --box_l 50 --gpu

--box_l sets the edge length of the cubic box in lattice units (default 50, minimum 35.7).
--single_precision runs the LB fluid in single precision.
Each run takes about 1 minute on the GPU at box size 50. The timings are sensitive to
other load on the machine, also for GPU runs.
"""

import argparse
import os
import sys
import time

import numpy as np
import espressomd
import espressomd.constraints
import espressomd.lb
import espressomd.magnetostatics
import espressomd.propagation

parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
parser.add_argument("--box_l", type=int, default=50,
                    help="edge length of the cubic box in lattice units (default: 50); "
                         "the droplet size stays fixed")
parser.add_argument("--gpu", action="store_true",
                    help="run the LB fluid on the GPU (needs a CUDA build of ESPResSo)")
parser.add_argument("--single_precision", action="store_true",
                    help="run the LB fluid in single precision")
args = parser.parse_args()
DEVICE = "gpu" if args.gpu else "cpu"
espressomd.assert_features(["WALBERLA", "DIPOLES", "ROTATION", "VIRTUAL_SITES_RELATIVE", "WCA"]
                           + (["CUDA"] if args.gpu else []))

#########################################################################
# Benchmark parameters
#########################################################################

WARMUP_STEPS = 1000
MEASUREMENT_STEPS = 100
N_ITERATIONS = 10

#########################################################################
# Simulation parameters
#########################################################################

DOMAIN_SIZE = args.box_l
AGRID = 1.0
TAU = 0.01
RHO_0 = 1.0
RADIUS = 12.5  # droplet radius, independent of the box size
# initial tanh interface width (equilibrium width at BETA=0.9)
SMOOTHING_WIDTH = 3.76
SIGMA = 0.05
GAMMA = 0.5
VISCOSITY = 1.0 / 6.0 * AGRID**2 / TAU
KT = 1.0
EPSILON = 1e-6
BETA = 0.9
# rotational friction/noise of the free-spinning dipoles (the LB thermostat only
# thermalizes translation)
GAMMA_ROTATION = 0.5

# Tripods: an anchor particle on the interface (WCA, LB coupled) and two LB coupled virtual
# sites at +-VS_OFFSET along the interface normal with solvation energy +-DELTA_MU, which pin
# the anchor to the interface and keep the tripod upright. A fourth virtual site at the
# anchor position carries the dipole and rotates freely.
particle_radius = 1
DELTA_MU = 100
VS_OFFSET = 1.6 * particle_radius
COVERAGE = 0.5  # surface coverage N*pi*r^2/(4*pi*R^2) of the anchors
N_PART = round(COVERAGE * 4 * (RADIUS / particle_radius)**2)
PREFACTOR = 1.0  # prefactor of the dipolar direct sum
# dipole moment, dipolar coupling lambda = PREFACTOR*DIPM^2/(KT*(2*particle_radius)^3)
DIPM = 1.0
ALPHA = 10.0  # Langevin parameter alpha = mu_0*DIPM*H/KT
H = ALPHA * KT / DIPM  # applied field magnitude (along x), mu_0 = 1

# the droplet with its tripods and interface has to fit into the box
min_box_l = 2 * (RADIUS + VS_OFFSET + SMOOTHING_WIDTH)
if DOMAIN_SIZE < min_box_l:
    parser.error(
        f"--box_l must be at least {min_box_l:.1f} to fit the droplet")

#########################################################################
# Initial configuration helpers
#########################################################################


def droplet_densities(grid_size, radius, smoothing_width, rho_0, epsilon):
    """Densities of component a (outside) and b (droplet) for a droplet in the box centre,
    with a tanh profile across the interface; the minority phase has density epsilon*rho_0."""
    x = np.arange(grid_size) + 0.5
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    center = grid_size / 2.0
    signed = np.sqrt((xx - center)**2 + (yy - center)
                     ** 2 + (zz - center)**2) - radius
    profile = (np.tanh(signed / smoothing_width * 2.64665) + 1.0) / \
        2.0  # 0 inside, 1 outside
    rho_a = epsilon * rho_0 + (rho_0 - epsilon * rho_0) * profile
    rho_b = rho_0 + (epsilon * rho_0 - rho_0) * profile
    return rho_a, rho_b


def fibonacci_sphere(n_points):
    """Quasi-uniformly distributed unit vectors on a sphere (Fibonacci lattice)."""
    indices = np.arange(n_points) + 0.5
    phi = np.arccos(1.0 - 2.0 * indices / n_points)
    theta = np.pi * (3.0 - np.sqrt(5.0)) * indices
    return np.stack([np.sin(phi) * np.cos(theta), np.sin(phi) * np.sin(theta), np.cos(phi)],
                    axis=-1)


def random_dipoles(n_particles):
    """Random unit dipole orientations (fixed seed, so every run starts identically)."""
    np.random.seed(seed=1)
    dip_phi = 2. * np.pi * np.random.random((n_particles, 1))
    dip_cos_theta = 2 * np.random.random((n_particles, 1)) - 1
    dip_sin_theta = np.sin(np.arccos(dip_cos_theta))
    return np.hstack((dip_sin_theta * np.sin(dip_phi),
                      dip_sin_theta * np.cos(dip_phi),
                      dip_cos_theta))


#########################################################################
# System setup
#########################################################################

system = espressomd.System(box_l=[DOMAIN_SIZE] * 3)
system.time_step = TAU
system.cell_system.skin = 0.4
system.periodicity = [True, True, True]
# with MPI the solvation virtual sites must lie within the minimum global cutoff of their
# anchor, otherwise vs_auto_relate_to raises
system.min_global_cut = VS_OFFSET + 0.1

lbf = espressomd.lb.LBFluid(agrid=AGRID,
                            density=RHO_0,
                            tau=system.time_step,
                            kT=KT,
                            kinematic_viscosity=[VISCOSITY, VISCOSITY],
                            sigma=SIGMA,
                            beta=BETA,
                            single_precision=args.single_precision,
                            gpu=args.gpu)
system.lb = lbf

rho_a, rho_b = droplet_densities(
    DOMAIN_SIZE/AGRID, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)
lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
lbf.init_two_component()
system.integrator.run(100)

#########################################################################
# Tripods on the interface and thermostats
#########################################################################

system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42, kT=KT)
# translational gamma=0: LB already provides translational friction/noise; this Langevin
# thermostat only adds gamma_rotation for the free-spinning dipoles (ROT_LANGEVIN)
system.thermostat.set_langevin(
    kT=KT, gamma=0.0, gamma_rotation=GAMMA_ROTATION, seed=43)

# anchor: LB coupled translation, rotates only by the torque of its virtual sites
ANCHOR_PROPAGATION = (espressomd.propagation.Propagation.TRANS_LB_MOMENTUM_EXCHANGE |
                      espressomd.propagation.Propagation.ROT_EULER)
# dipole virtual site: follows the anchor position, rotates freely (no ROT_VS_* flag, so
# the dip set at creation survives the zero-distance vs_auto_relate_to)
CORE_PROPAGATION = (espressomd.propagation.Propagation.TRANS_VS_RELATIVE |
                    espressomd.propagation.Propagation.ROT_LANGEVIN)

directions = fibonacci_sphere(N_PART)
center = np.array(system.box_l) / 2.0
dipoles = DIPM * random_dipoles(N_PART)

anchors = system.part.add(pos=center + RADIUS * directions,
                          type=[0] * N_PART,
                          rotation=N_PART * [(True, True, True)])
cores = system.part.add(pos=center + RADIUS * directions,
                        type=[2] * N_PART,
                        dip=dipoles,
                        rotation=N_PART * [(True, True, True)])
# positive delta_mu is attracted to phase b (the droplet): +DELTA_MU on the inner site,
# -DELTA_MU on the outer site (towards the continuous phase)
vs_inner = system.part.add(pos=center + (RADIUS - VS_OFFSET) * directions,
                           type=[1] * N_PART,
                           solvation_delta_mu=[DELTA_MU] * N_PART)
vs_outer = system.part.add(pos=center + (RADIUS + VS_OFFSET) * directions,
                           type=[1] * N_PART,
                           solvation_delta_mu=[-DELTA_MU] * N_PART)
for anchor, core, inner, outer in zip(anchors, cores, vs_inner, vs_outer):
    anchor.propagation = ANCHOR_PROPAGATION
    core.vs_auto_relate_to(anchor)
    core.propagation = CORE_PROPAGATION
    inner.vs_auto_relate_to(anchor, couple_to_lb=True)
    outer.vs_auto_relate_to(anchor, couple_to_lb=True)

system.non_bonded_inter[0, 0].wca.set_params(
    epsilon=1.0, sigma=2*particle_radius)
system.magnetostatics.solver = espressomd.magnetostatics.DipolarDirectSum(
    prefactor=PREFACTOR)

state = system.cell_system.get_state()
n_ranks, n_threads = state["n_nodes"], state["omp_num_threads"]
n_lb_nodes = round(DOMAIN_SIZE / AGRID)**3
precision = "single" if args.single_precision else "double"
node_grid = [int(n) for n in state["node_grid"]]
print(f"LB on the {DEVICE.upper()} in {precision} precision, {n_ranks} MPI rank(s) "
      f"(node grid {node_grid}), {n_threads} OpenMP thread(s)")
print(f"{n_lb_nodes} LB nodes, {N_PART} tripods ({4 * N_PART} particles)")

#########################################################################
# Warm-up, then switch on the field
#########################################################################

print(f"Warm-up: {WARMUP_STEPS} steps")
system.integrator.run(WARMUP_STEPS)
system.constraints.add(
    espressomd.constraints.HomogeneousMagneticField(H=[H, 0, 0]))
system.integrator.run(1)

#########################################################################
# Timing
#########################################################################

print(f"Timing {N_ITERATIONS} x {MEASUREMENT_STEPS} steps")
timings = []
for i in range(N_ITERATIONS):
    tick = time.perf_counter()
    system.integrator.run(MEASUREMENT_STEPS)
    t = (time.perf_counter() - tick) / MEASUREMENT_STEPS
    timings.append(t)
    print(f"iteration {i}: {1000 * t:.2f} ms/step")
timings = np.array(timings)

avg = timings.mean()
ci = 1.96 * timings.std() / np.sqrt(len(timings) - 1)
print(f"average: {1000 * avg:.2f} +/- {1000 * ci:.2f} ms/step (95% C.I.), "
      f"{n_lb_nodes / avg / 1e6:.2f} MLUPS")

report = "benchmarks.csv"
header = '"script","arguments","ranks","threads","mean","ci","nsteps","duration","label"\n'
with open(report, "a") as f:
    if f.tell() == 0:
        f.write(header)
    f.write(f'"{os.path.basename(__file__)}","{" ".join(sys.argv[1:])}",{n_ranks},{n_threads},'
            f'{avg:.3e},{ci:.3e},{MEASUREMENT_STEPS},{timings.sum():.1f},"{DEVICE}"\n')
print(f"Appended result to {os.path.abspath(report)}")
