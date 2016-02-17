"""
==============
Sphinx Gallery
==============

"""
import os
__version__ = '0.1.2'


def glr_path_static():
    """Returns path to packaged static files"""
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '_static'))
