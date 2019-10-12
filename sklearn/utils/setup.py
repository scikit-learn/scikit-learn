import os
from os.path import join


def configuration(parent_package='', top_path=None):
    import numpy
    from numpy.distutils.misc_util import Configuration
    from Cython import Tempita

    config = Configuration('utils', parent_package, top_path)

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('sparsefuncs_fast',
                         sources=['sparsefuncs_fast.pyx'],
                         libraries=libraries)

    config.add_extension('_cython_blas',
                         sources=['_cython_blas.pyx'],
                         libraries=libraries)

    config.add_extension('arrayfuncs',
                         sources=['arrayfuncs.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension('murmurhash',
                         sources=['murmurhash.pyx', join(
                             'src', 'MurmurHash3.cpp')],
                         include_dirs=['src'])

    config.add_extension('graph_shortest_path',
                         sources=['graph_shortest_path.pyx'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('fast_dict',
                         sources=['fast_dict.pyx'],
                         language="c++",
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    # generate files from a template
    pyx_templates = ['sklearn/utils/seq_dataset.pyx.tp',
                     'sklearn/utils/seq_dataset.pxd.tp']

    for pyxfiles in pyx_templates:
        outfile = pyxfiles.replace('.tp', '')
        # if .pyx.tp is not updated, no need to output .pyx
        if (os.path.exists(outfile) and
                os.stat(pyxfiles).st_mtime < os.stat(outfile).st_mtime):
            continue

        with open(pyxfiles, "r") as f:
            tmpl = f.read()
        pyxcontent = Tempita.sub(tmpl)

        with open(outfile, "w") as f:
            f.write(pyxcontent)

    config.add_extension('seq_dataset',
                         sources=['seq_dataset.pyx'],
                         include_dirs=[numpy.get_include()])

    config.add_extension('_weight_vector',
                         sources=['_weight_vector.pyx'],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension("_random",
                         sources=["_random.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_extension("_logistic_sigmoid",
                         sources=["_logistic_sigmoid.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=libraries)

    config.add_subpackage('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
