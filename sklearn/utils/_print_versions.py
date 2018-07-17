"""
Utility methods to print system info for debugging

adapted from :func:`pandas.show_versions`
"""
# License: BSD 3 clause

import os
import platform
import sys
import struct
import subprocess
import codecs
import locale
import importlib

from .._build_utils import get_blas_info


def get_sys_info():
    """System information

    Return
    ------
    sys_info : dict
        system and Python version information

    """

    blob = []

    # get full commit hash
    commit = None
    if os.path.isdir(".git") and os.path.isdir("sklearn"):
        try:
            pipe = subprocess.Popen('git log --format="%H" -n 1'.split(" "),
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            so, serr = pipe.communicate()
        except:
            pass
        else:
            if pipe.returncode == 0:
                commit = so
                try:
                    commit = so.decode('utf-8')
                except ValueError:
                    pass
                commit = commit.strip().strip('"')

    blob.append(('commit', commit))

    try:
        (sysname, nodename, release,
         version, machine, processor) = platform.uname()
        blob.extend([
            ("python", '.'.join(map(str, sys.version_info))),
            ("python-bits", struct.calcsize("P") * 8),
            ("OS", "{sysname}".format(sysname=sysname)),
            ("OS-release", "{release}".format(release=release)),
            ("machine", "{machine}".format(machine=machine)),
            ("processor", "{processor}".format(processor=processor)),
            ("byteorder", "{byteorder}".format(byteorder=sys.byteorder)),
            ("LC_ALL", "{lc}".format(lc=os.environ.get('LC_ALL', "None"))),
            ("LANG", "{lang}".format(lang=os.environ.get('LANG', "None"))),
            ("LOCALE", '.'.join(map(str, locale.getlocale()))),
        ])
    except:
        pass

    return dict(blob)


def get_deps_info():
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
        ("numpy", lambda mod: mod.version.version),
        ("scipy", lambda mod: mod.version.version),
        ("Cython", lambda mod: mod.__version__),
        ("pandas", lambda mod: mod.__version__),
        ("matplotlib", lambda mod: mod.__version__),
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


def show_versions(as_json=False, with_blas=False):
    sys_info = get_sys_info()
    deps_info = get_deps_info()

    if with_blas:
        try:
            cblas_libs, blas_info = get_blas_info()
            blas_info['cblas_libs'] = cblas_libs
        except:
            blas_info = {'error': "Could not retrieve BLAS information"}

    if (as_json):
        try:
            import json
        except ImportError:
            import simplejson as json

        j = dict(system=sys_info,
                 dependencies=deps_info)

        if with_blas:
            j['blas'] = blas_info

        if as_json is True:
            print(j)
        else:
            with codecs.open(as_json, "wb", encoding='utf8') as f:
                json.dump(j, f, indent=2)
    else:
        print("")
        print('System info')
        print('-----------')
        for k, stat in sys_info.items():
            print("{k}: {stat}".format(k=k, stat=stat))

        if with_blas:
            print("")
            print('BLAS info')
            print('---------')
            for k, stat in blas_info.items():
                print("{k}: {stat}".format(k=k, stat=stat))

        print("")
        print('Python libs info')
        print('----------------')
        for k, stat in deps_info.items():
            print("{k}: {stat}".format(k=k, stat=stat))
