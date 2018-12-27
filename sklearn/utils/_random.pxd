# Authors: Arnaud Joly
#
# License: BSD 3 clause
#
# cython: language_level=3


import numpy as np
cimport numpy as np


cpdef sample_without_replacement(np.int_t n_population,
                                 np.int_t n_samples,
                                 method=*,
                                 random_state=*)

