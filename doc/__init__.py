"""Enable experimental modules so that experimental estimators can be discovered
properly by sphinx"""

from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.experimental import enable_iterative_imputer  # noqa

# This doesn't consistently prevent the ImportError when sphinx-gallery is build in
# parallel. Can be deleted later.
