def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration

    config = Configuration('neighbors', parent_package, top_path)

    config.add_extension('ball_tree',
         sources=['ball_tree.c'],
         extra_compile_args=["-Wno-unused-function", "-Wno-unused-but-set-variable"],
         include_dirs=[numpy.get_include()])

    return config
