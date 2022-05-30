import os


def configuration(parent_package="", top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration("feature_extraction", parent_package, top_path)
    libraries = []
    if os.name == "posix":
        libraries.append("m")

    config.add_extension(
        "_hashing_fast",
        sources=["_hashing_fast.pyx"],
        include_dirs=[numpy.get_include()],
        language="c++",
        libraries=libraries,
    )
    config.add_subpackage("tests")

    return config
