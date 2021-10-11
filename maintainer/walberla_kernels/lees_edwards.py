from pystencils.astnodes import LoopOverCoordinate
#from pystencils.field import fields
from pystencils.data_types import type_all_numbers
from pystencils.data_types import TypedSymbol
from pystencils import Assignment

#from lbmpy.creationfunctions import create_lb_update_rule, create_mrt_orthogonal, force_model_from_string
#from lbmpy.stencils import get_stencil
#from lbmpy.updatekernels import create_stream_pull_with_output_kernel
#from lbmpy.macroscopic_value_kernels import macroscopic_values_setter
#from pystencils_walberla import CodeGeneration, generate_sweep, generate_info_header
#from lbmpy_walberla import RefinementScaling, generate_boundary, generate_lb_pack_info
#import relaxation_rates

import sympy as sp


def velocity_shift():
    dim = 3
    counters = [LoopOverCoordinate.get_loop_counter_symbol(
        i) for i in range(dim)]
    points_up = TypedSymbol('points_up', bool)
    points_down = TypedSymbol('points_down', bool)
    shear_velocity = TypedSymbol('shear_velocity', float)
    grid_size = TypedSymbol('grid_size', int)

    u_p = sp.Piecewise((1, sp.And(type_all_numbers(counters[1] <= 0, 'int'), points_down)),
                       (-1,
                        sp.And(type_all_numbers(counters[1] >= grid_size - 1,
                                                'int'),
                               points_up)),
                       (0, True)) * shear_velocity
    return u_p


def modify_method(collision, shear_dir, shear_dir_normal):
    points_up = TypedSymbol('points_up', bool)
    points_down = TypedSymbol('points_down', bool)

    to_insert = [s.lhs for s in collision.subexpressions
                 if collision.method.first_order_equilibrium_moment_symbols[shear_dir]
                 in s.free_symbols]
    for s in to_insert:
        collision = collision.new_with_inserted_subexpression(s)
    ma = []
    for a, c in zip(collision.main_assignments, collision.method.stencil):
        if c[shear_dir_normal] == -1:
            b = (True, False)
        elif c[shear_dir_normal] == 1:
            b = (False, True)
        else:
            b = (False, False)
        a = Assignment(a.lhs, a.rhs.replace(points_down, b[0]))
        a = Assignment(a.lhs, a.rhs.replace(points_up, b[1]))
        ma.append(a)
    collision.main_assignments = ma
    return collision
