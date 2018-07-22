"""
Utility methods to print system info for debugging

adapted from :func:`pandas.show_versions`
"""
# License: BSD 3 clause

import platform
import sys
import importlib


def _get_sys_info():
    """System information

    Return
    ------
    sys_info : dict
        system and Python version information

    """
    python = sys.version.replace('\n', ' ')

    blob = [
        ("python", python),
        ('executable', sys.executable),
        ("machine", platform.platform()),
    ]

    return dict(blob)


def _get_deps_info():
    """Overview of the installed version of main dependencies

    Returns
    -------
    deps_info: dict
        version information on relevant Python libraries

    """
    deps = [
        # (MODULE_NAME, f(mod) -> mod version)
        ("pip", lambda mod: mod.__version__),
        ("setuptools", lambda mod: mod.__version__),
        ("sklearn", lambda mod: mod.__version__),
        ("numpy", lambda mod: mod.__version__),
        ("scipy", lambda mod: mod.version.version),
        ("Cython", lambda mod: mod.__version__),
        ("pandas", lambda mod: mod.__version__),
    ]

    deps_blob = []

    for modname, ver_func in deps:
        try:
            if modname in sys.modules:
                mod = sys.modules[modname]
            else:
                mod = importlib.import_module(modname)
            ver = ver_func(mod)
            deps_blob.append((modname, ver))
        except ImportError:
            deps_blob.append((modname, None))

    return dict(deps_blob)


def _get_blas_info():
    """Information on system BLAS

    Uses the `scikit-learn` builtin method
    :func:`sklearn._build_utils.get_blas_info` which may fail from time to time

    Returns
    -------
    blas_info: dict
        system BLAS information

    """
    from .._build_utils import get_blas_info

    cblas_libs, blas_dict = get_blas_info()

    macros = ['{key}={val}'.format(key=a, val=b)
              for (a, b) in blas_dict.get('define_macros', [])]

    blas_blob = [
        ('macros', ', '.join(macros)),
        ('lib_dirs', ':'.join(blas_dict.get('library_dirs', ''))),
        ('cblas_libs', ', '.join(cblas_libs)),
    ]

    return dict(blas_blob)


def show_versions():
    "Print useful debugging information to the terminal"

    sys_info = _get_sys_info()
    deps_info = _get_deps_info()
    blas_info = _get_blas_info()

    print('\nSystem')
    print('------')
    for k, stat in sys_info.items():
        print("{k:>10}: {stat}".format(k=k, stat=stat))

    print('\nBLAS')
    print('----')
    for k, stat in blas_info.items():
        print("{k:>10}: {stat}".format(k=k, stat=stat))

    print('\nPython deps')
    print('-----------')
    for k, stat in deps_info.items():
        print("{k:>10}: {stat}".format(k=k, stat=stat))
