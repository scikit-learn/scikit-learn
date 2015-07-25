# Authors: Arnaud Joly
#
# Licence: BSD 3 clause


import numpy as np
cimport numpy as np


cpdef sample_without_replacement(np.int_t n_population,
                                 np.int_t n_samples,
                                 method=*,
                                 random_state=*)

