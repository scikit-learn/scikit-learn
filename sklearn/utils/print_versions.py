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
    "Returns system information as a dict"

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

    return blob


def show_versions(as_json=False, with_blas=False):
    sys_info = get_sys_info()

    cblas_libs, blas_info = get_blas_info()

    deps = [
        # (MODULE_NAME, f(mod) -> mod version)
        ("pip", lambda mod: mod.__version__),
        ("setuptools", lambda mod: mod.__version__),
        ("numpy", lambda mod: mod.version.version),
        ("scipy", lambda mod: mod.version.version),
        ("Cython", lambda mod: mod.__version__),
        ("pandas", lambda mod: mod.__version__),
        ("matplotlib", lambda mod: mod.__version__),
        # ("sphinx", lambda mod: mod.__version__),
        # ("pytest", lambda mod: mod.__version__),
        # ("patsy", lambda mod: mod.__version__),
        # ("pyarrow", lambda mod: mod.__version__),
        # ("xarray", lambda mod: mod.__version__),
        # ("IPython", lambda mod: mod.__version__),
        # ("dateutil", lambda mod: mod.__version__),
        # ("pytz", lambda mod: mod.VERSION),
        # ("blosc", lambda mod: mod.__version__),
        # ("bottleneck", lambda mod: mod.__version__),
        # ("tables", lambda mod: mod.__version__),
        # ("numexpr", lambda mod: mod.__version__),
        # ("feather", lambda mod: mod.__version__),
        # ("openpyxl", lambda mod: mod.__version__),
        # ("xlrd", lambda mod: mod.__VERSION__),
        # ("xlwt", lambda mod: mod.__VERSION__),
        # ("xlsxwriter", lambda mod: mod.__version__),
        # ("lxml", lambda mod: mod.etree.__version__),
        # ("bs4", lambda mod: mod.__version__),
        # ("html5lib", lambda mod: mod.__version__),
        # ("sqlalchemy", lambda mod: mod.__version__),
        # ("pymysql", lambda mod: mod.__version__),
        # ("psycopg2", lambda mod: mod.__version__),
        # ("jinja2", lambda mod: mod.__version__),
        # ("s3fs", lambda mod: mod.__version__),
        # ("fastparquet", lambda mod: mod.__version__),
        # ("pandas_gbq", lambda mod: mod.__version__),
        # ("pandas_datareader", lambda mod: mod.__version__),
        # ("gcsfs", lambda mod: mod.__version__),
    ]

    deps_blob = list()
    for (modname, ver_f) in deps:
        try:
            if modname in sys.modules:
                mod = sys.modules[modname]
            else:
                mod = importlib.import_module(modname)
            ver = ver_f(mod)
            deps_blob.append((modname, ver))
        except ImportError:
            deps_blob.append((modname, None))

    if (as_json):
        try:
            import json
        except ImportError:
            import simplejson as json

        j = dict(system=dict(sys_info),
                 dependencies=dict(deps_blob))

        if with_blas:
            blas_info['cblas_libs'] = cblas_libs
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
        for k, stat in sys_info:
            print("{k}: {stat}".format(k=k, stat=stat))

        if with_blas:
            print("")
            print('BLAS info')
            print('---------')
            for k, stat in blas_info.items():
                print("{k}: {stat}".format(k=k, stat=stat))
            print('CBLAS libs: {libs}'.format(libs=cblas_libs))

        print("")
        print('Python libs info')
        print('----------------')
        for k, stat in deps_blob:
            print("{k}: {stat}".format(k=k, stat=stat))
