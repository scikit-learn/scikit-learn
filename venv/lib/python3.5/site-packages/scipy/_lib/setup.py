from __future__ import division, print_function, absolute_import

import os


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('_lib', parent_package, top_path)
    config.add_data_files('tests/*.py')

    include_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'src'))
    depends = [os.path.join(include_dir, 'ccallback.h')]

    config.add_extension("_ccallback_c",
                         sources=["_ccallback_c.c"],
                         depends=depends,
                         include_dirs=[include_dir])

    config.add_extension("_test_ccallback",
                         sources=["src/_test_ccallback.c"],
                         depends=depends,
                         include_dirs=[include_dir])

    config.add_extension("_fpumode",
                         sources=["_fpumode.c"])

    def get_messagestream_config(ext, build_dir):
        # Generate a header file containing defines
        config_cmd = config.get_config_cmd()
        defines = []
        if config_cmd.check_func('open_memstream', decl=True, call=True):
            defines.append(('HAVE_OPEN_MEMSTREAM', '1'))
        target = os.path.join(os.path.dirname(__file__), 'src',
                              'messagestream_config.h')
        with open(target, 'w') as f:
            for name, value in defines:
                f.write('#define {0} {1}\n'.format(name, value))

    depends = [os.path.join(include_dir, 'messagestream.h')]
    config.add_extension("messagestream",
                         sources=["messagestream.c"] + [get_messagestream_config],
                         depends=depends,
                         include_dirs=[include_dir])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup

    setup(**configuration(top_path='').todict())
