#!python

import numpy as np
cimport numpy as np
from libc.math cimport sqrt

# use a hack to determine the associated numpy data types
# NOTE: the following requires the buffer interface, only available in
#       numpy 1.5+.  We'll choose the DTYPE by hand instead.
#cdef ITYPE_t idummy
#cdef ITYPE_t[:] idummy_view = <ITYPE_t[:1]> &idummy
#ITYPE = np.asarray(idummy_view).dtype
ITYPE = np.intp  # WARNING: this should match ITYPE_t in typedefs.pxd

#cdef DTYPE_t ddummy
#cdef DTYPE_t[:] ddummy_view = <DTYPE_t[:1]> &ddummy
#DTYPE = np.asarray(ddummy_view).dtype
DTYPE = np.float64  # WARNING: this should match DTYPE_t in typedefs.pxd

# some handy constants
cdef DTYPE_t INF = np.inf
cdef DTYPE_t PI = np.pi
cdef DTYPE_t ROOT_2PI = sqrt(2 * PI)
