import os


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('feature_extraction', parent_package, top_path)

    libraries = []
    if os.name == 'posix':
        libraries.append('m')

    config.add_extension('_term_counter',
         sources=['_term_counter.c'])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
