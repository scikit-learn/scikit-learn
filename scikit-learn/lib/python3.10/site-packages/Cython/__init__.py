from .Shadow import __version__

# Void cython.* directives (for case insensitive operating systems).
from .Shadow import *


def load_ipython_extension(ip):
    """Load the extension in IPython."""
    from .Build.IpythonMagic import CythonMagics  # pylint: disable=cyclic-import
    ip.register_magics(CythonMagics)
