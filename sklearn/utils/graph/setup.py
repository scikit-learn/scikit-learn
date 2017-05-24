import numpy


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('graph', parent_package, top_path)

    config.add_extension('_shortest_path',
                         sources=['_shortest_path.pyx'],
                         include_dirs=[numpy.get_include()])

    config.add_subpackage('tests')

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
