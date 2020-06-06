
def pre_build_hook(build_ext, ext):
    from scipy._build_utils.compiler_helper import (
        get_cxx_std_flag, has_flag, try_add_flag)
    cc = build_ext._cxx_compiler
    args = ext.extra_compile_args

    std_flag = get_cxx_std_flag(cc)
    if std_flag is not None:
        args.append(std_flag)

    if cc.compiler_type == 'msvc':
        args.append('/EHsc')
    else:
        try_add_flag(args, cc, '-fvisibility=hidden')

        min_macos_flag = '-mmacosx-version-min=10.9'
        import sys
        if sys.platform == 'darwin' and has_flag(cc, min_macos_flag):
            args.append(min_macos_flag)
            ext.extra_link_args.append(min_macos_flag)


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('_uarray', parent_package, top_path)
    config.add_data_files('LICENSE')
    ext = config.add_extension('_uarray',
                               sources=['_uarray_dispatch.cxx'],
                               language='c++')
    ext._pre_build_hook = pre_build_hook
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
