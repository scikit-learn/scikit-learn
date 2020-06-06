
def pre_build_hook(build_ext, ext):
    from scipy._build_utils.compiler_helper import (
        get_cxx_std_flag, try_add_flag, try_compile, has_flag)
    cc = build_ext._cxx_compiler
    args = ext.extra_compile_args

    std_flag = get_cxx_std_flag(build_ext._cxx_compiler)
    if std_flag is not None:
        args.append(std_flag)

    if cc.compiler_type == 'msvc':
        args.append('/EHsc')
    else:
        try_add_flag(args, cc, '-fvisibility=hidden')

        has_pthreads = try_compile(cc, code='#include <pthread.h>\n'
                                   'int main(int argc, char **argv) {}')
        if has_pthreads:
            ext.define_macros.append(('POCKETFFT_PTHREADS', None))

        min_macos_flag = '-mmacosx-version-min=10.9'
        import sys
        if sys.platform == 'darwin' and has_flag(cc, min_macos_flag):
            args.append(min_macos_flag)
            ext.extra_link_args.append(min_macos_flag)


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    import pybind11
    include_dirs = [pybind11.get_include(True), pybind11.get_include(False)]

    config = Configuration('_pocketfft', parent_package, top_path)
    ext = config.add_extension('pypocketfft',
                               sources=['pypocketfft.cxx'],
                               depends=['pocketfft_hdronly.h'],
                               include_dirs=include_dirs,
                               language='c++')
    ext._pre_build_hook = pre_build_hook

    config.add_data_files('LICENSE.md')
    config.add_data_dir('tests')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
