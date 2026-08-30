#
# Copyright (C) 2013-2026 The ESPResSo project
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

import sys
import espressomd
import espressomd.lees_edwards
import benchmarks
import numpy as np
import argparse

DPD_GAMMA = 4.5
DPD_R_CUT = 1.5
KT = 1.0
SEED = 42
DEFAULT_PARTICLES_PER_CORE = 1000
INITIAL_SKIN = 0.4
SHEAR_DIRECTION = "x"
SHEAR_PLANE_NORMAL = "y"

# Overridable by the importlib smoke test (see testsuite/benchmarks).
measurement_steps = None
n_iterations = 30
min_skin = 0.2
max_skin = 1.0

parser = argparse.ArgumentParser(
    description="Benchmark a sheared DPD fluid (Lees-Edwards). "
                "Save the results to a CSV file.")
benchmarks.add_common_args(parser, DEFAULT_PARTICLES_PER_CORE)
parser.add_argument("--density", metavar="RHO", action="store", type=float,
                    default=0.5, required=False,
                    help="Number density of DPD particles (default: 0.5)")
parser.add_argument("--shear_velocity", metavar="V", action="store",
                    type=float, default=1.0, required=False,
                    help="Lees-Edwards shear velocity (default: 1.0)")
parser.add_argument("--fully_connected", action="store_true",
                    help="Use the legacy fully_connected_boundary decomposition "
                         "(reference). Default: dynamic sheared halo.")
parser.add_argument("--retune_skin_after", metavar="N", action="store",
                    type=int, default=None, required=False,
                    help="Retune skin every N timing iterations "
                    "(0 disables, default: 5).")
parser.add_argument("--output", metavar="FILEPATH", action="store",
                    type=str, required=False, default="benchmarks.csv",
                    help="Output file (default: benchmarks.csv)")

args = parser.parse_args()
benchmarks.validate_mode(args)

espressomd.assert_features(["DPD"])
np.random.seed(SEED)

system = espressomd.System(box_l=[1, 1, 1])
system.time_step = 0.01


def configure_dpd(system):
    system.non_bonded_inter[0, 0].dpd.set_params(
        weight_function=1, gamma=DPD_GAMMA, r_cut=DPD_R_CUT,
        trans_weight_function=1, trans_gamma=DPD_GAMMA, trans_r_cut=DPD_R_CUT)


def configure_decomposition(system, args):
    """Select the cell decomposition for LE shear."""
    if args.fully_connected:
        normal_axis = {"x": 0, "y": 1, "z": 2}[SHEAR_PLANE_NORMAL]
        n_nodes = system.cell_system.get_state()["n_nodes"]
        system.cell_system.node_grid = [
            n_nodes if i == normal_axis else 1 for i in range(3)]
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True,
            fully_connected_boundary={"boundary": SHEAR_PLANE_NORMAL,
                                      "direction": SHEAR_DIRECTION})
    else:
        system.cell_system.set_regular_decomposition(
            use_verlet_lists=True, fully_connected_boundary=None)


def configure_lees_edwards(system, args):
    protocol = espressomd.lees_edwards.LinearShear(
        initial_pos_offset=0.0, shear_velocity=args.shear_velocity)
    system.lees_edwards.set_boundary_conditions(
        shear_direction=SHEAR_DIRECTION,
        shear_plane_normal=SHEAR_PLANE_NORMAL,
        protocol=protocol)


def resolve_retune(args, meta=None):
    if args.retune_skin_after is not None:
        value = args.retune_skin_after
    elif args.skin is not None:
        value = 0
    elif meta is not None:
        value = int(meta["retune_skin_after"])
    else:
        value = 5
    return None if value <= 0 else value


def build_and_tune(system, args):
    n_part = benchmarks.resolve_n_part(system, args)
    steps = measurement_steps
    if steps is None:
        steps = max(20, int(np.round(
            5e6 / n_part * system.cell_system.get_state()["n_nodes"], -2)))
    box_l = (n_part / args.density)**(1. / 3.)
    system.box_l = 3 * (box_l,)
    system.cell_system.skin = INITIAL_SKIN if args.skin is None else args.skin
    configure_dpd(system)
    system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)

    benchmarks.minimize(system, 1 / system.time_step)
    system.integrator.set_vv()
    system.thermostat.set_dpd(kT=KT, seed=SEED)

    print("Equilibration")
    system.time_step /= 10.
    system.integrator.run(100)
    system.time_step *= 10.
    system.integrator.run(min(5 * steps, 60000))
    print("Tune skin: {:.3f}".format(benchmarks.tune_skin_unless_fixed(
        system, args, min_skin, max_skin, tol=0.025, int_steps=200)))

    configure_decomposition(system, args)
    configure_lees_edwards(system, args)
    print("Equilibration")
    system.integrator.run(min(10 * steps, 60000))

    return {"n_part": n_part, "measurement_steps": steps}


def save_dpd_state(system, args, ctx):
    resolved_retune = resolve_retune(args)
    p = system.part.all()
    meta = {
        "skin": float(system.cell_system.skin),
        "box_l": [float(x) for x in system.box_l],
        "n_part": int(ctx["n_part"]),
        "time_step": float(system.time_step),
        "measurement_steps": int(ctx["measurement_steps"]),
        "n_iterations": n_iterations,
        "density": float(args.density),
        "shear_velocity": float(args.shear_velocity),
        "fully_connected": bool(args.fully_connected),
        "retune_skin_after": 0 if resolved_retune is None else resolved_retune,
        "kT": KT, "gamma": DPD_GAMMA, "r_cut": DPD_R_CUT, "seed": SEED,
    }
    meta.update(benchmarks.topology_meta(system))
    benchmarks.save_state(
        args.state_file, meta,
        pos=np.copy(p.pos), vel=np.copy(p.v))
    print(f"Saved state to {args.state_file}")


def run_from_state(system, args):
    meta, handle = benchmarks.load_state(args.state_file)
    benchmarks.verify_topology(system, meta)
    system.box_l = meta["box_l"]
    system.time_step = meta["time_step"]
    system.cell_system.skin = args.skin if args.skin is not None else meta["skin"]
    configure_dpd(system)
    system.part.add(pos=handle["pos"], v=handle["vel"])
    system.integrator.set_vv()
    system.thermostat.set_dpd(kT=meta["kT"], seed=int(meta["seed"]))
    args.fully_connected = bool(meta["fully_connected"])
    args.shear_velocity = float(meta["shear_velocity"])
    configure_decomposition(system, args)
    configure_lees_edwards(system, args)
    retune = resolve_retune(args, meta)
    return meta["measurement_steps"], meta["n_iterations"], retune


def time_and_report(system, args, measurement_steps, n_iterations, retune):
    timings = benchmarks.get_timings(
        system, measurement_steps, n_iterations,
        retune_skin_after_steps=retune)
    avg, ci = benchmarks.get_average_time(timings)
    print(f"average: {avg:.3e} +/- {ci:.3e} (95% C.I.)")
    n_proc = system.cell_system.get_state()["n_nodes"]
    benchmarks.write_report(args.output, n_proc, timings, measurement_steps,
                            n_threads=benchmarks.get_omp_num_threads(system))


if args.mode == "run":
    measurement_steps, n_iterations, retune = run_from_state(system, args)
    time_and_report(system, args, measurement_steps, n_iterations, retune)
    sys.exit(0)

ctx = build_and_tune(system, args)

if args.state_file:
    save_dpd_state(system, args, ctx)

if args.mode == "tune":
    sys.exit(0)

time_and_report(system, args, ctx["measurement_steps"], n_iterations,
                resolve_retune(args))
