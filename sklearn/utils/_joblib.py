# We need the absolute_import to avoid the local joblib to override the
# site one
from __future__ import absolute_import
import os as _os

# An environment variable to use the site joblib
if _os.environ.get('SKLEARN_SITE_JOBLIB', False):
    from joblib import __all__
    from joblib import *  # noqa
    from joblib import __version__
    from joblib import logger
else:
    from ..externals.joblib import __all__   # noqa
    from ..externals.joblib import *  # noqa
    from ..externals.joblib import __version__  # noqa
    from ..externals.joblib import logger  # noqa
