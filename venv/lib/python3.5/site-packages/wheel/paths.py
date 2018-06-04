"""
Installation paths.

Map the .data/ subdirectory names to install paths.
"""

import distutils.command.install as install
import distutils.dist as dist
import os.path
import sys


def get_install_command(name):
    # late binding due to potential monkeypatching
    d = dist.Distribution({'name': name})
    i = install.install(d)
    i.finalize_options()
    return i


def get_install_paths(name):
    """
    Return the (distutils) install paths for the named dist.

    A dict with ('purelib', 'platlib', 'headers', 'scripts', 'data') keys.
    """
    paths = {}

    i = get_install_command(name)

    for key in install.SCHEME_KEYS:
        paths[key] = getattr(i, 'install_' + key)

    # pip uses a similar path as an alternative to the system's (read-only)
    # include directory:
    if hasattr(sys, 'real_prefix'):  # virtualenv
        paths['headers'] = os.path.join(sys.prefix,
                                        'include',
                                        'site',
                                        'python' + sys.version[:3],
                                        name)

    return paths
