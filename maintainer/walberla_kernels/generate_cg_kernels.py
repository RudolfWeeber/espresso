#
# Copyright (C) 2020-2024 The ESPResSo project
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
import packaging.specifiers

import pystencils as ps

import lbmpy

kernel_codes = "cg_init cg_stream cg_color_gradient cg_collide".split()
parser = argparse.ArgumentParser(description="Generate the waLBerla kernels.")
parser.add_argument("--single-precision", action="store_true", required=False,
                    help="Use single-precision")
parser.add_argument("--gpu", action="store_true")
parser.add_argument("--kernels", nargs="+", type=str, default="all",
                    choices=["all"] + kernel_codes,
                    help="Which kernels to generate")
args = parser.parse_args()

# Make sure we have the correct versions of the required dependencies
for module, requirement in [(ps, "==1.4.0"), (lbmpy, "==1.4.0")]:
    assert packaging.specifiers.SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} " \
        f"doesn't match requirement {requirement}"

import pystencils_walberla
import pystencils_espresso
import lbmpy.stencils
import lbmpy.enums

import color_gradient
import code_generation_context
from kernel_generation_utils import paramlist, get_ext_source, patch_openmp_kernels

if args.gpu:
    target = ps.Target.GPU
else:
    target = ps.Target.CPU
if args.kernels == "all" or args.kernels == ["all"]:
    args.kernels = kernel_codes

# vectorization parameters
parameters = {}
if target == ps.Target.GPU:
    default_key = "GPU"
    parameters["GPU"] = ({"target": target}, "CUDA")
else:
    default_key = "CPU"
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": False}
    parameters["CPU"] = ({"target": target,
                          "cpu_openmp": True}, "")
    parameters["CPU_linear"] = ({"target": target}, "")
    parameters["AVX"] = ({"target": target,
                          "cpu_openmp": True,
                          "cpu_vectorize_info": cpu_vectorize_info}, "AVX")

# global parameters
stencil = lbmpy.stencils.LBStencil(lbmpy.enums.Stencil.D3Q19)


def generate_cg_init_kernels(ctx, cg_fields, cg_methods):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    init = color_gradient.create_init_operator(cg_fields, cg_methods)
    for params, target_suffix in paramlist(parameters, (default_key,)):
        stem = f"ColorGradientInitialPDFsSetter{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(ctx, stem, init, **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


def generate_cg_color_gradient_kernels(ctx, cg_fields):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    cg_op = color_gradient.create_color_gradient_operator(cg_fields)
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"ColorGradientSweep{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(ctx, stem, cg_op, **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


def generate_cg_collide_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    collide = color_gradient.create_collide_perturb_recolor_operator(cg_fields, cg_methods, cg_configs, cg_opts)
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"ColorGradientCollideSweep{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(ctx, stem, collide, **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


def generate_cg_stream_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    stream = color_gradient.create_stream_operator(cg_fields, cg_methods, cg_configs, cg_opts)
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"ColorGradientStreamSweep{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(
            ctx, stem, stream,
            field_swaps=[
                (cg_fields["pdfs_a"], cg_fields["pdfs_a_tmp"]),
                (cg_fields["pdfs_b"], cg_fields["pdfs_b_tmp"]),
            ],
            **params)
        ctx.patch_file(stem, get_ext_source(target_suffix), patch_openmp_kernels)


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    if target == ps.Target.CPU:
        ctx.openmp = True
    if target == ps.Target.GPU:
        ctx.gpu = True
        ctx.cuda = True

    # codegen configuration
    config = pystencils_espresso.generate_config(
        ctx, parameters[default_key][0])
    data_type = "float64" if ctx.double_accuracy else "float32"
    
    # CG fields    
    cg_fields = color_gradient.generate_fields(stencil, data_type)

    # CG Method definition
    cg_methods = color_gradient.create_methods(stencil, cg_fields)
    
    cg_configs = color_gradient.create_configs(cg_fields, cg_methods)
    cg_opts = color_gradient.create_opts(cg_fields)

    if "cg_init" in args.kernels:
        generate_cg_init_kernels(ctx, cg_fields, cg_methods)
    if "cg_stream" in args.kernels:
        generate_cg_stream_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts)
    if "cg_color_gradient" in args.kernels:
        generate_cg_color_gradient_kernels(ctx, cg_fields)
    if "cg_collide" in args.kernels:
        generate_cg_collide_kernels(ctx, cg_fields, cg_methods, cg_configs, cg_opts)