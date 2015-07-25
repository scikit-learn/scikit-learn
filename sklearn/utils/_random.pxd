# Authors: Arnaud Joly
#
# Licence: BSD 3 clause


import numpy as np
cimport numpy as np


cpdef sample_without_replacement(np.npy_intp n_population,
                                 np.npy_intp n_samples,
                                 method=*,
                                 random_state=*)
