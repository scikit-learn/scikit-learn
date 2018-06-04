import os


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('skimage', parent_package, top_path)

    config.add_subpackage('_shared')
    config.add_subpackage('color')
    config.add_subpackage('data')
    config.add_subpackage('draw')
    config.add_subpackage('exposure')
    config.add_subpackage('feature')
    config.add_subpackage('restoration')
    config.add_subpackage('filters')
    config.add_subpackage('future')
    config.add_subpackage('graph')
    config.add_subpackage('io')
    config.add_subpackage('measure')
    config.add_subpackage('morphology')
    config.add_subpackage('transform')
    config.add_subpackage('util')
    config.add_subpackage('segmentation')
    config.add_subpackage('external')

    def add_test_directories(arg, dirname, fnames):
        if dirname.split(os.path.sep)[-1] == 'tests':
            config.add_data_dir(dirname)

    # Add test directories
    from os.path import isdir, dirname, join
    rel_isdir = lambda d: isdir(join(curpath, d))

    curpath = join(dirname(__file__), './')
    subdirs = [join(d, 'tests') for d in os.listdir(curpath) if rel_isdir(d)]
    subdirs = [d for d in subdirs if rel_isdir(d)]
    for test_dir in subdirs:
        config.add_data_dir(test_dir)
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup

    config = configuration(top_path='').todict()
    setup(**config)
