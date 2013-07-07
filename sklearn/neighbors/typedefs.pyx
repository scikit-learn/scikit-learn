#!python

import numpy as np
cimport numpy as np
from libc.math cimport sqrt

# use a hack to determine the associated numpy data types
cdef ITYPE_t idummy
cdef ITYPE_t[:] idummy_view = <ITYPE_t[:1]> &idummy
ITYPE = np.asarray(idummy_view).dtype

cdef DTYPE_t ddummy
cdef DTYPE_t[:] ddummy_view = <DTYPE_t[:1]> &ddummy
DTYPE = np.asarray(ddummy_view).dtype

# some handy constants
cdef DTYPE_t INF = np.inf
cdef DTYPE_t PI = np.pi
cdef DTYPE_t ROOT_2PI = sqrt(2 * PI)
