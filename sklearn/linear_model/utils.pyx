from types cimport floating, complexing


cdef void real_part(complexing x, floating *y) nogil:
    if complexing is float or complexing is double:
        y[0] = x
    elif complexing is complex:
        y[0] = creal(x)
    else:
        y[0] = crealf(x)


cdef inline floating fmax(floating x, floating y) nogil:
    return y if x < y else x


cdef void abs_max(int n, complexing *a, int inca, floating *b) nogil:
    """np.max(np.abs(a))"""
    if complexing is float or complexing is double:
        abs = fabs
    elif complexing is double:
        abs = cabs
    else:
        abs = cabsf

    cdef int i
    cdef double d
    b[0] = abs(a[0])
    for i in range(1, n):
        i *= inca
        d = abs(a[i])
        if d > b[0]:
            b[0] = d


cdef floating fmax_arr(int n, floating *a, int inca) nogil:
    """np.max(a)"""
    cdef int i
    cdef floating m = a[0]
    cdef floating d
    for i in range(1, n):
        i *= inca
        d = a[i]
        if d > m:
            m = d
    return m


cdef void diff_abs_max(int n, complexing *a, int inca, complexing *b, int incb, floating *c) nogil:
    """np.max(np.abs(a - b))"""
    if complexing is double or complexing is float:
        abs = fabs
    elif complexing is complex:
        abs = cabs
    else:
        abs = cabsf
    cdef int i
    cdef floating d
    c[0] = abs(a[0] - b[0])
    for i in range(1, n):
        d = fmax(c[0], abs(a[i * inca] - b[i * incb]))


cdef inline void relu(int n, floating *x, int incx) nogil:
    cdef int i
    for i in range(n):
        i *= incx
        if x[i] < 0.:
            x[i] = 0.
