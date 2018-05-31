# We need the absolute_import to avoid the local joblib to override the
# site one
from __future__ import absolute_import
import os as _os

# An environment variable to use the site joblib
if _os.environ.get('SKLEARN_SITE_JOBLIB', False):
    from joblib import *
    from joblib import __version__
    from joblib import logger
else:
    from ._joblib import *
    from ._joblib import __version__
    from ._joblib import logger

