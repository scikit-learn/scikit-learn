from libc.math cimport lgamma


cdef double lgamma(double x):
    msg = ("sklearn.utils.lgamma is deprecated in v0.22 and will "
           "be removed in v0.24. Use libc.math.lgamma instead.")
    warnings.warn(msg, category=DeprecationWarning)
    if x <= 0:
        raise ValueError("x must be strictly positive, got %f" % x)
    return lgamma(x)
