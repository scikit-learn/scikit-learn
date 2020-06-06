from __future__ import division, print_function


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('fft', parent_package, top_path)

    config.add_data_dir('tests')

    # Configure pocketfft_internal
    config.add_extension('_pocketfft_internal',
                         sources=['_pocketfft.c']
                         )

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(configuration=configuration)
