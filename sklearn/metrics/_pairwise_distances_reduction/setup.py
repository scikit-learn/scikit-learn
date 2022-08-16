import os

import numpy as np
from numpy.distutils.misc_util import Configuration

from sklearn._build_utils import gen_from_templates


def configuration(parent_package="", top_path=None):
    config = Configuration("_pairwise_distances_reduction", parent_package, top_path)
    libraries = []
    if os.name == "posix":
        libraries.append("m")

    templates = [
        "sklearn/metrics/_pairwise_distances_reduction/_datasets_pair.pyx.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_datasets_pair.pxd.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_gemm_term_computer.pyx.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_gemm_term_computer.pxd.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_base.pyx.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_base.pxd.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_argkmin.pyx.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_argkmin.pxd.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_radius_neighborhood.pyx.tp",
        "sklearn/metrics/_pairwise_distances_reduction/_radius_neighborhood.pxd.tp",
    ]

    gen_from_templates(templates)

    cython_sources = [
        "_datasets_pair.pyx",
        "_gemm_term_computer.pyx",
        "_base.pyx",
        "_argkmin.pyx",
        "_radius_neighborhood.pyx",
    ]

    for source_file in cython_sources:
        private_extension_name = source_file.replace(".pyx", "")
        config.add_extension(
            name=private_extension_name,
            sources=[source_file],
            include_dirs=[np.get_include()],
            language="c++",
            libraries=libraries,
            extra_compile_args=["-std=c++11"],
        )

    return config


if __name__ == "__main__":
    from numpy.distutils.core import setup

    setup(**configuration().todict())
