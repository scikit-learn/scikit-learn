# Import necessary to preserve backward compatibility of pickles
import sys
import warnings

from joblib import *


msg = ("sklearn.externals.joblib is deprecated in 0.21 and will be removed "
       "in 0.23. Please import this functionality directly from joblib, "
       "which can be installed with: pip install joblib. If this warning is "
       "raised when loading pickled models, you may need to re-serialize "
       "those models with scikit-learn 0.21+.")

if not hasattr(sys, "_is_pytest_session"):
    warnings.warn(msg, category=DeprecationWarning)
