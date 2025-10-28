"""
Sphinx Gallery
==============

"""

import os

# dev versions should have "dev" in them, stable should not.
# doc/conf.py makes use of this to set the version drop-down.
__version__ = "0.19.0"


def glr_path_static():
    """Returns path to packaged static files"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), "_static"))
