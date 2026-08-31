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

import argparse
import itertools
import numpy as np
import espressomd
import espressomd.lees_edwards

SHEAR_DIRECTION = "x"
SHEAR_PLANE_NORMAL = "y"
AXIS = {"x": 0, "y": 1, "z": 2}
DPD_GAMMA = 4.5
DPD_R_CUT = 1.5
KT = 1.0
SEED = 42


def build_system(node_grid, density, box_l, shear_velocity,
                 initial_pos_offset):
    system = espressomd.System(box_l=3 * [box_l])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    system.non_bonded_inter[0, 0].dpd.set_params(
        weight_function=1, gamma=DPD_GAMMA, r_cut=DPD_R_CUT,
        trans_weight_function=1, trans_gamma=DPD_GAMMA, trans_r_cut=DPD_R_CUT)
    # Deterministic identical initial condition on every decomposition.
    rng = np.random.default_rng(SEED)
    n_part = int(round(density * box_l**3))
    pos = rng.random((n_part, 3)) * box_l
    vel = (rng.random((n_part, 3)) - 0.5)
    system.part.add(pos=pos, v=vel)
    system.integrator.set_vv()
    system.thermostat.set_dpd(kT=KT, seed=SEED)

    system.cell_system.node_grid = list(node_grid)
    system.cell_system.set_regular_decomposition(use_verlet_lists=True)
    protocol = espressomd.lees_edwards.LinearShear(
        initial_pos_offset=initial_pos_offset, shear_velocity=shear_velocity)
    system.lees_edwards.set_boundary_conditions(
        shear_direction=SHEAR_DIRECTION,
        shear_plane_normal=SHEAR_PLANE_NORMAL, protocol=protocol)
    return system


def trace_selfcheck(system):
    """Within-cutoff pairs the core reports vs Python brute force on the SAME
    folded positions. Returns (ok, within_cut_pairs, n_missing, n_extra)."""
    system.integrator.run(0)
    folded = {p.id: p.pos_folded for p in system.part.all()}
    core = set()
    for id1, id2, _p1, _p2, vec21, _node in \
            system.cell_system.non_bonded_loop_trace():
        if np.linalg.norm(vec21) <= DPD_R_CUT:
            core.add((id1, id2) if id1 < id2 else (id2, id1))
    py = set()
    for i, j in itertools.combinations(folded.keys(), 2):
        if np.linalg.norm(system.distance_vec(folded[i], folded[j])) <= DPD_R_CUT:
            py.add((i, j) if i < j else (j, i))
    missing = py - core
    extra = core - py
    return (len(missing) == 0 and len(extra) == 0), core, len(missing), len(extra)


def shear_stress(system):
    t = system.analysis.pressure_tensor()["total"]
    return float(t[AXIS[SHEAR_DIRECTION]][AXIS[SHEAR_PLANE_NORMAL]])


def dump(args):
    system = build_system(args.node_grid, args.density,
                          args.box_l, args.shear_velocity,
                          args.initial_pos_offset)
    # Step-0 within-cutoff self-check on the identical deterministic config:
    # this pair set is the cross-decomposition criterion-2 proof.
    ok0, pairs0, miss0, extra0 = trace_selfcheck(system)
    assert ok0, f"step-0 trace mismatch: missing={miss0} extra={extra0}"
    stresses = []
    for _ in range(args.steps):
        system.integrator.run(1)
        stresses.append(shear_stress(system))
    # Per-run correctness on the (now diverged) final config.
    okN, _pairsN, missN, extraN = trace_selfcheck(system)
    stresses = np.array(stresses)
    p = system.part.all()
    order = np.argsort(p.id)
    np.savez(f"{args.out}/{args.tag}.npz",
             pos=np.copy(p.pos_folded)[order],
             pairs0=np.array(sorted(pairs0), dtype=int),
             stress_mean=float(np.mean(stresses)),
             stress_std=float(np.std(stresses)),
             stress_n=int(len(stresses)),
             trace_ok=bool(ok0 and okN))
    print(f"dumped {args.tag}: step0_within_cut={len(pairs0)} "
          f"trace_ok={ok0 and okN} finalcheck(missing={missN},extra={extraN})")


def compare(args):
    a = np.load(args.a, allow_pickle=True)
    b = np.load(args.b, allow_pickle=True)
    assert bool(a["trace_ok"]) and bool(b["trace_ok"]), \
        "a per-run within-cutoff trace self-check failed"
    print("trace: both runs find exactly the within-cutoff pairs (vs Python)")
    sa = {tuple(r) for r in a["pairs0"]}
    sb = {tuple(r) for r in b["pairs0"]}
    assert sa == sb, \
        f"step-0 within-cutoff pair sets differ: {len(sa ^ sb)} symmetric-diff"
    print(f"trace: identical step-0 within-cutoff pair sets ({len(sa)})")
    if args.bitwise:
        assert np.array_equal(
            a["pos"], b["pos"]), "positions not bitwise equal"
        print("pos: bitwise identical")
    ma, mb = float(a["stress_mean"]), float(b["stress_mean"])
    sea = float(a["stress_std"]) / max(1, int(a["stress_n"]))**0.5
    seb = float(b["stress_std"]) / max(1, int(b["stress_n"]))**0.5
    denom = (sea**2 + seb**2)**0.5
    z = abs(ma - mb) / denom if denom > 0 else 0.0
    print(f"shear stress: {ma:.6e} vs {mb:.6e}  z={z:.2f} "
          f"(se_a={sea:.2e}, se_b={seb:.2e})")
    assert z < args.stress_z, f"shear stress differs beyond {
        args.stress_z} sigma"


p = argparse.ArgumentParser()
sub = p.add_subparsers(dest="cmd", required=True)
d = sub.add_parser("dump")
d.add_argument("--tag", required=True)
d.add_argument("--out", default="/tmp")
d.add_argument("--node_grid", type=int, nargs=3, default=[1, 1, 1])
d.add_argument("--steps", type=int, default=200)
d.add_argument("--density", type=float, default=2.0)
d.add_argument("--box_l", type=float, default=10.0)
d.add_argument("--shear_velocity", type=float, default=1.0)
d.add_argument("--initial_pos_offset", type=float, default=1.3)
d.set_defaults(func=dump)
c = sub.add_parser("compare")
c.add_argument("--a", required=True)
c.add_argument("--b", required=True)
c.add_argument("--bitwise", action="store_true")
c.add_argument("--stress_z", type=float, default=5.0)
c.set_defaults(func=compare)
args = p.parse_args()
args.func(args)
