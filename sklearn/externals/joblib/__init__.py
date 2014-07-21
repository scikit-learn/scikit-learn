
import os
if int(os.environ.get('SKLEARN_SYSTEM_JOBLIB', '0')):
    # Try to import joblib from the system for development purpose
    from joblib import *
else:
    from .embedded_joblib_init import *

