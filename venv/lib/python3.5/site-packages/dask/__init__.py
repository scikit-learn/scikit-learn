from __future__ import absolute_import, division, print_function

from .core import istask
from .context import set_options
from .local import get_sync as get
try:
    from .delayed import delayed
except ImportError:
    pass
try:
    from .base import visualize, compute, persist, optimize, is_dask_collection
except ImportError:
    pass

from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
