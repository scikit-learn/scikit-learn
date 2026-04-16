try:
    # python3 vs. python2
    from configparser import ConfigParser
except ImportError:
    from ConfigParser import SafeConfigParser as ConfigParser
import io
import logging
import os
import re
from shlex import split as shsplit
import sys

import numpy

logger = logging.getLogger('pythran')


def get_include():
    # using / as separator as advised in the distutils doc
    return (os.path.dirname(os.path.dirname(__file__)) or '.') + '/pythran'


class silent(object):
    '''
    Silent sys.stderr at the system level
    '''

    def __enter__(self):
        try:
            self.prevfd = os.dup(sys.stderr.fileno())
            os.close(sys.stderr.fileno())
        except io.UnsupportedOperation:
            self.prevfd = None

        self.prevstream = sys.stderr
        sys.stderr = open(os.devnull, 'r')

    def __exit__(self, exc_type, exc_value, traceback):
        sys.stderr.close()
        sys.stderr = self.prevstream
        if self.prevfd:
            os.dup2(self.prevfd, sys.stderr.fileno())
            os.close(self.prevfd)


def get_paths_cfg(
        sys_file='pythran.cfg',
        platform_file='pythran-{}.cfg'.format(sys.platform),
        user_file='.pythranrc'
):
    sys_config_dir = os.path.dirname(__file__)
    sys_config_path = os.path.join(sys_config_dir, sys_file)

    platform_config_path = os.path.join(sys_config_dir, platform_file)
    if not os.path.exists(platform_config_path):
        platform_config_path = os.path.join(sys_config_dir,
                                            "pythran-default.cfg")

    user_config_path = os.environ.get('PYTHRANRC', None)
    if user_config_path is None:
        user_config_dir = os.environ.get('XDG_CONFIG_HOME', None)
        if not user_config_dir:
            user_config_dir = os.environ.get('HOME', None)
        if not user_config_dir:
            user_config_dir = '~'
        user_config_path = os.path.expanduser(os.path.join(user_config_dir,
                                                           user_file))

    paths = {"sys": sys_config_path,
             "platform": platform_config_path}

    if user_config_path:
        paths["user"] = user_config_path

    return paths


def init_cfg(sys_file, platform_file, user_file, config_args=None):
    paths = get_paths_cfg(sys_file, platform_file, user_file)

    sys_config_path = paths["sys"]
    platform_config_path = paths["platform"]
    user_config_path = paths.get("user")

    cfgp = ConfigParser()
    for required in (sys_config_path, platform_config_path):
        cfgp.read([required])

    if user_config_path:
        cfgp.read([user_config_path])

    if config_args is not None:
        update_cfg(cfgp, config_args)

    return cfgp


def update_cfg(cfgp, config_args):
    # Override the config options with those provided on the command line
    # e.g. compiler.blas=pythran-openblas.
    for arg in config_args:
        try:
            lhs, rhs = arg.split('=', maxsplit=1)
            section, item = lhs.split('.')
            if not cfgp.has_section(section):
                cfgp.add_section(section)
            cfgp.set(section, item, rhs)
        except Exception:
            pass


def lint_cfg(cfgp, **paths):
    if not paths:
        paths = get_paths_cfg()

    # Use configuration from sys and platform as "reference"
    cfgp_ref = ConfigParser()
    cfgp_ref.read([paths["sys"], paths["platform"]])

    # Check if pythran configuration files exists
    for loc, path in paths.items():
        exists = os.path.exists(path)

        msg = " ".join([
            "{} file".format(loc).rjust(13),
            "exists:" if exists else "does not exist:",
            path
        ])
        logger.info(msg) if exists else logger.warning(msg)

    for section in cfgp.sections():
        # Check if section in the current configuration exists in the
        # reference configuration
        if cfgp_ref.has_section(section):
            options = set(cfgp.options(section))
            options_ref = set(cfgp_ref.options(section))

            # Check if the options in the section are supported by the
            # reference configuration
            if options.issubset(options_ref):
                logger.info(
                    (
                        "pythranrc section [{}] is valid and options are "
                        "correct"
                    ).format(section)
                )
            else:
                logger.warning(
                    (
                        "pythranrc section [{}] is valid but options {} "
                        "are incorrect!"
                    ).format(section, options.difference(options_ref))
                )
        else:
            logger.warning("pythranrc section [{}] is invalid!"
                           .format(section))


def make_extension(python, **extra):
    # load platform specific configuration then user configuration
    cfg = init_cfg('pythran.cfg',
                   'pythran-{}.cfg'.format(sys.platform),
                   '.pythranrc',
                   extra.pop('config', None))


    def parse_define(define):
        index = define.find('=')
        if index < 0:
            return (define, None)
        else:
            return define[:index], define[index + 1:]

    extension = {
        "language": "c++",
        # forcing str conversion to handle Unicode case (the default on MS)
        "define_macros": [str(x) for x in
                          shsplit(cfg.get('compiler', 'defines'))],
        "undef_macros": [str(x) for x in
                         shsplit(cfg.get('compiler', 'undefs'))],
        "include_dirs": [str(x) for x in
                         shsplit(cfg.get('compiler', 'include_dirs'))],
        "library_dirs": [str(x) for x in
                         shsplit(cfg.get('compiler', 'library_dirs'))],
        "libraries": [str(x) for x in
                      shsplit(cfg.get('compiler', 'libs'))],
        "extra_compile_args": [str(x) for x in
                               shsplit(cfg.get('compiler', 'cflags'))],
        "extra_link_args": [str(x) for x in
                            shsplit(cfg.get('compiler', 'ldflags'))],
        "extra_objects": []
    }

    # In case the extension doesn't include any optimization level, make sure we
    # have a decent default. Later options have prcedence so this is still
    # customizable by users.
    extension['extra_compile_args'].insert(0, "-O2")

    if python:
        extension['define_macros'].append('ENABLE_PYTHON_MODULE')

    pythonic_dir = get_include()

    extension["include_dirs"].append(pythonic_dir)

    extra.pop('language', None)  # forced to c++ anyway
    cxx = extra.pop('cxx', None)
    cc = extra.pop('cc', None)

    if cxx is None:
        cxx = compiler()
    if cxx is not None:
        extension['cxx'] = cxx
        extension['cc'] = cc or cxx

    # Honor CXXFLAGS (note: Pythran calls this `cflags` everywhere, however the
    # standard environment variable is `CXXFLAGS` not `CFLAGS`).
    cflags = os.environ.get('CXXFLAGS', None)
    if cflags is not None:
        extension['extra_compile_args'].extend(shsplit(cflags))

    # Honor LDFLAGS
    ldflags = os.environ.get('LDFLAGS', None)
    if ldflags is not None:
        extension['extra_link_args'].extend(shsplit(ldflags))

    major = minor = None
    try:
        limited_api = cfg.getboolean("backend", "limited_api")
        if limited_api:
            if sys.version_info < (3, 10):
                logger.warning("Setting 'backend.limited_api' is only supported starting with Python 3.10, ignoring.")
            else:
                major, minor = 3, 7
    except ValueError:
        limited_api = cfg.get("backend", "limited_api")
        if re.match(r'^3\.[0-9]+$', limited_api):
            major, minor  = map(int, limited_api.split('.'))
            if (major, minor) < (3, 7):
                logger.warning("Invalid 'backend.limited_api' entry, expecting at least 3.7, ignoring.")
                major = minor = None
        else:
            logger.warning("Invalid 'backend.limited_api' entry, expected boolean or '3.7+', ignoring.")

    if major is not None:
        extension["define_macros"].append(f'Py_LIMITED_API={hex(major << 24 | minor << 16)}')
        extension["py_limited_api"] = True

    for k, w in extra.items():
        extension[k].extend(w)
    if cfg.getboolean('pythran', 'complex_hook'):
        # the patch is *not* portable
        extension["include_dirs"].append(pythonic_dir + '/pythonic/patch')

    # numpy specific
    if python:
        extension['include_dirs'].append(numpy.get_include())

    # blas dependency
    reserved_blas_entries = 'scipy-openblas', 'pythran-openblas', 'none'
    user_blas = cfg.get('compiler', 'blas')
    if user_blas == 'pythran-openblas':
        try:
            import pythran_openblas as openblas
            # required to cope with atlas missing extern "C"
            extension['define_macros'].append('PYTHRAN_BLAS_OPENBLAS')
            extension['include_dirs'].extend(openblas.include_dirs)
            extension['extra_objects'].append(
                os.path.join(openblas.library_dir, openblas.static_library)
            )
        except ImportError:
            logger.warning("Failed to find 'pythran-openblas' package. "
                           "Please install it or change the compiler.blas "
                           "setting. Defaulting to 'none'")
            user_blas = 'none'
    elif user_blas == 'scipy-openblas':
        try:
            import scipy_openblas64 as openblas
            # required to cope with atlas missing extern "C"
            extension['define_macros'].append('PYTHRAN_BLAS_SCIPY_OPENBLAS')
            extension['include_dirs'].append(openblas.get_include_dir())
            extension['library_dirs'].append(openblas.get_lib_dir())
            extension['libraries'].append(openblas.get_library())
            extension['extra_link_args'].append("-Wl,-rpath=" + openblas.get_lib_dir())
        except ImportError:
            logger.warning("Failed to find 'scipy-openblas64' package. "
                           "Please install it or change the compiler.blas "
                           "setting. Defaulting to 'none'")
            user_blas = 'none'


    if user_blas == 'none':
        extension['define_macros'].append('PYTHRAN_BLAS_NONE')

    if user_blas not in reserved_blas_entries:
        extension["libraries"].append(user_blas)
        extension['define_macros'].append('PYTHRAN_BLAS_{}'.format(user_blas.upper()))


    # final macro normalization
    extension["define_macros"] = [
        dm if isinstance(dm, tuple) else parse_define(dm)
        for dm in extension["define_macros"]]
    return extension


def compiler():
    """Get compiler to use for C++ to binary process. The precedence for
    choosing the compiler is as follows::

      1. `CXX` environment variable
      2. User configuration (~/.pythranrc)

    Returns None if none is set or if it's set to the empty string

    """
    cfg_cxx = str(cfg.get('compiler', 'CXX'))
    if not cfg_cxx:
        cfg_cxx = None
    return os.environ.get('CXX', cfg_cxx) or None


# load platform specific configuration then user configuration
cfg = init_cfg('pythran.cfg',
               'pythran-{}.cfg'.format(sys.platform),
               '.pythranrc')


def run():
    '''
    Dump on stdout the config flags required to compile pythran-generated code.
    '''
    import argparse
    import pythran

    parser = argparse.ArgumentParser(
        prog='pythran-config',
        description='output build options for pythran-generated code',
        epilog="It's a megablast!"
    )

    parser.add_argument('--compiler', action='store_true',
                        help='print default compiler')

    group = parser.add_mutually_exclusive_group()
    group.add_argument('--cflags', action='store_true',
                        help='print compilation flags to compile extension '
                             'modules directly')

    # The with-BLAS variant could be added if there's user demand. Unclear
    # if there are packages that need this, and also integrate with a build
    # system for Python and NumPy flags. SciPy and scikit-image need the
    # no-BLAS variant.
    group.add_argument('--cflags-pythran-only', action='store_true',
                       help='print compilation flags for usage by a build '
                            'system (doesn\'t include Python, NumPy or BLAS '
                            'flags).'
                       )

    parser.add_argument('--include-dir', action='store_true',
                        help=(
                            'print Pythran include directory '
                            '(matches `pythran.get_include()`).')
                        )

    parser.add_argument('--libs', action='store_true',
                        help='print linker flags')

    parser.add_argument('--no-python', action='store_true',
                        help='do not include Python-related flags')

    parser.add_argument('-V', '--version',
                        action='version',
                        version=pythran.version.__version__)

    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help=(
                            'verbose mode: [-v] prints warnings if pythranrc '
                            'has an invalid configuration; use '
                            '[-vv] for more information')
                        )

    args = parser.parse_args(sys.argv[1:])

    # This should not rely on distutils/setuptools, or anything else not in the
    # stdlib. Please don't add imports higher up.
    if args.cflags_pythran_only:
        compile_flags = [
                "-DENABLE_PYTHON_MODULE",
                "-DPYTHRAN_BLAS_NONE",
                ]
        print(" ".join(compile_flags), f"-I{get_include()}")
        return None


    import distutils.ccompiler
    import distutils.sysconfig
    import numpy

    args.python = not args.no_python

    output = []

    extension = pythran.config.make_extension(python=args.python)

    if args.verbose >= 1:
        if args.verbose == 1:
            logger.setLevel(logging.WARNING)
        else:
            logger.setLevel(logging.INFO)

        lint_cfg(cfg)

    if args.compiler or args.verbose >= 2:
        cxx = compiler() or 'c++'
        logger.info('CXX = '.rjust(10) + cxx)
        if args.compiler:
            output.append(cxx)

    compiler_obj = distutils.ccompiler.new_compiler()
    distutils.sysconfig.customize_compiler(compiler_obj)

    if args.cflags or args.verbose >= 2:
        def fmt_define(define):
            name, value = define
            if value is None:
                return '-D' + name
            else:
                return '-D' + name + '=' + value

        cflags = []
        cflags.extend(fmt_define(define)
                      for define in extension['define_macros'])
        cflags.extend(('-I' + include)
                      for include in extension['include_dirs'])
        if args.python:
            cflags.append('-I' + numpy.get_include())
            cflags.append('-I' + distutils.sysconfig.get_python_inc())

        logger.info('CXXFLAGS = '.rjust(10) + ' '.join(cflags))
        if args.cflags:
            output.extend(cflags)

    if args.libs or args.verbose >= 2:
        ldflags = []
        ldflags.extend(('-L' + compiler_obj.library_dir_option(include))
                       for include in extension['library_dirs'])
        ldflags.extend((compiler_obj.library_option(include))
                       for include in extension['libraries'])

        if args.python:
            libs = distutils.sysconfig.get_config_var('LIBS')
            if libs:
                ldflags.extend(shsplit(libs))
            libdir = distutils.sysconfig.get_config_var('LIBDIR')
            pylib = distutils.sysconfig.get_config_var('LDLIBRARY')
            if libdir and pylib:
                ldflags.append(os.path.join(libdir, pylib))
            else:
                ldflags.append(compiler_obj.library_option('python')
                               + distutils.sysconfig.get_config_var('VERSION'))

        logger.info('LDFLAGS = '.rjust(10) + ' '.join(ldflags))
        if args.libs:
            output.extend(ldflags)

    if args.include_dir:
        output.append(get_include())

    if output:
        print(' '.join(output))


if __name__ == '__main__':
    run()
