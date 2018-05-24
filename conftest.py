# This file is here so that when running from the root folder
# ./sklearn is added to sys.path by pytest.
# See https://docs.pytest.org/en/latest/pythonpath.html for more details.
# For example, this allows to build extensions in place and run pytest
# doc/modules/clustering.rst and use sklearn from the local folder
# rather than the one from site-packages.

# Set numpy array str/repr to legacy behaviour on numpy > 1.13 to make
# the doctests pass
import numpy as np
try:
    np.set_printoptions(legacy='1.13')
except TypeError:
    pass
