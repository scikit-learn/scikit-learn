"""
Small collection of auxiliary functions that operate on arrays

"""

from cython cimport floating
from libc.math cimport fabs
from libc.float cimport DBL_MAX, FLT_MAX

from ._cython_blas cimport _copy, _rotg, _rot


def min_pos(const floating[:] X):
    """Find the minimum value of an array over positive values

    Returns the maximum representable value of the input dtype if none of the
    values are positive.
    """
    cdef Py_ssize_t i
    cdef floating min_val = FLT_MAX if floating is float else DBL_MAX
    for i in range(X.size):
        if 0. < X[i] < min_val:
            min_val = X[i]
    return min_val


# General Cholesky Delete.
# Remove an element from the cholesky factorization
# m = columns
# n = rows
#
# TODO: put transpose as an option
def cholesky_delete(const floating[:, :] L, int go_out):
    cdef:
        int n = L.shape[0]
        int m = L.strides[0]
        floating c, s
        floating *L1
        int i

    if floating is float:
        m /= sizeof(float)
    else:
        m /= sizeof(double)

    # delete row go_out
    L1 = &L[0, 0] + (go_out * m)
    for i in range(go_out, n-1):
        _copy(i + 2, L1 + m, 1, L1, 1)
        L1 += m

    L1 = &L[0, 0] + (go_out * m)
    for i in range(go_out, n-1):
        _rotg(L1 + i, L1 + i + 1, &c, &s)
        if L1[i] < 0:
            # Diagonals cannot be negative
            L1[i] = fabs(L1[i])
            c = -c
            s = -s

        L1[i + 1] = 0.  # just for cleanup
        L1 += m

        _rot(n - i - 2, L1 + i, m, L1 + i + 1, m, c, s)
