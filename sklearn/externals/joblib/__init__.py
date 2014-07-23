
import os
import warnings
if int(os.environ.get('SKLEARN_SYSTEM_JOBLIB', '0')):
    # Try to import joblib from the system for development purpose
    from joblib import *
    from joblib import __version__
    from .embedded_joblib_init import __version__ as orig_version
    warnings.warn('Using system joblib %s instead of embedded joblib %s'
                  % (__version__, orig_version))
else:
    from .embedded_joblib_init import *

