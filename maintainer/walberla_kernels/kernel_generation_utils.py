import re

def paramlist(parameters, keys):
    for key in keys:
        if key in parameters:
            yield parameters[key]


def get_ext_header(target_suffix):
    return {"CUDA": "h"}.get(target_suffix, "h")


def get_ext_source(target_suffix):
    return {"CUDA": "cu"}.get(target_suffix, "cpp")


def patch_openmp_kernels(content):
    # surrounds omp pragmas with ifdefs
    content = re.sub("^( *#pragma omp .*)$",
                     r"#ifdef _OPENMP\n\1\n#endif", content, flags=re.MULTILINE)
    return content