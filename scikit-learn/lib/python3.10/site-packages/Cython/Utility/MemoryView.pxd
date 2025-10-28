# from cpython cimport ...

cimport cython
cdef extern from "Python.h":
    ctypedef struct PyObject
    ctypedef void *PyThread_type_lock

cdef extern from *:
    ctypedef struct {{memviewslice_name}}:
        pass

    ctypedef struct __pyx_buffer "Py_buffer":
        PyObject *obj

    ctypedef struct __Pyx_TypeInfo:
        pass

    cdef struct __pyx_memoryview "__pyx_memoryview_obj":
        Py_buffer view
        PyObject *obj
        const __Pyx_TypeInfo *typeinfo

    ctypedef int __pyx_atomic_int_type

    {{memviewslice_name}} slice_copy_contig "__pyx_memoryview_copy_new_contig"(
                                 {{memviewslice_name}} *from_mvs,
                                 char *mode, int ndim,
                                 size_t sizeof_dtype, int contig_flag,
                                 bint dtype_is_object) except * nogil
    bint slice_is_contig "__pyx_memviewslice_is_contig" (
                            {{memviewslice_name}} mvs, char order, int ndim) nogil
    bint slices_overlap "__pyx_slices_overlap" ({{memviewslice_name}} *slice1,
                                                {{memviewslice_name}} *slice2,
                                                int ndim, size_t itemsize) nogil

@cname("__pyx_array")
cdef class array:

    cdef:
        char *data
        Py_ssize_t len
        char *format
        int ndim
        Py_ssize_t *_shape
        Py_ssize_t *_strides
        Py_ssize_t itemsize
        unicode mode  # FIXME: this should have been a simple 'char'
        bytes _format
        void (*callback_free_data)(void *data) noexcept
        # cdef object _memview
        cdef bint free_data
        cdef bint dtype_is_object

    @cname('get_memview')
    cdef get_memview(self)

@cname("__pyx_array_allocate_buffer")
cdef int _allocate_buffer(array self) except -1

@cname("__pyx_array_new")
cdef array array_cwrapper(tuple shape, Py_ssize_t itemsize, char *format, const char *c_mode, char *buf)

@cname('__pyx_memoryview')
cdef class memoryview:

    cdef object obj
    cdef object _size
    cdef object _array_interface
    cdef PyThread_type_lock lock
    cdef __pyx_atomic_int_type acquisition_count
    cdef Py_buffer view
    cdef int flags
    cdef bint dtype_is_object
    cdef const __Pyx_TypeInfo *typeinfo

    cdef char *get_item_pointer(memoryview self, object index) except NULL
    cdef is_slice(self, obj)
    cdef setitem_slice_assignment(self, dst, src)
    cdef setitem_slice_assign_scalar(self, memoryview dst, value)
    cdef setitem_indexed(self, index, value)
    cdef convert_item_to_object(self, char *itemp)
    cdef assign_item_from_object(self, char *itemp, object value)
    cdef _get_base(self)

@cname('__pyx_memoryview_new')
cdef memoryview_cwrapper(object o, int flags, bint dtype_is_object, const __Pyx_TypeInfo *typeinfo)

@cname('__pyx_memoryview_check')
cdef inline bint memoryview_check(object o) noexcept:
    return isinstance(o, memoryview)

@cname('__pyx_memview_slice')
cdef memoryview memview_slice(memoryview memview, object indices)

@cname('__pyx_memoryview_slice_memviewslice')
cdef int slice_memviewslice(
        {{memviewslice_name}} *dst,
        Py_ssize_t shape, Py_ssize_t stride, Py_ssize_t suboffset,
        int dim, int new_ndim, int *suboffset_dim,
        Py_ssize_t start, Py_ssize_t stop, Py_ssize_t step,
        int have_start, int have_stop, int have_step,
        bint is_slice) except -1 nogil

@cname('__pyx_pybuffer_index')
cdef char *pybuffer_index(Py_buffer *view, char *bufp, Py_ssize_t index,
                          Py_ssize_t dim) except NULL

@cname('__pyx_memslice_transpose')
cdef int transpose_memslice({{memviewslice_name}} *memslice) except -1 nogil

@cname('__pyx_memoryview_fromslice')
cdef memoryview_fromslice({{memviewslice_name}} memviewslice,
                          int ndim,
                          object (*to_object_func)(char *),
                          int (*to_dtype_func)(char *, object) except 0,
                          bint dtype_is_object)

@cname('__pyx_memoryview_get_slice_from_memoryview')
cdef {{memviewslice_name}} *get_slice_from_memview(memoryview memview,
                                                {{memviewslice_name}} *mslice) except NULL

@cname('__pyx_memoryview_slice_copy')
cdef void slice_copy(memoryview memview, {{memviewslice_name}} *dst) noexcept

@cname('__pyx_memoryview_copy_object')
cdef memoryview_copy(memoryview memview)

@cname('__pyx_memoryview_copy_object_from_slice')
cdef memoryview_copy_from_slice(memoryview memview, {{memviewslice_name}} *memviewslice)

@cname('__pyx_get_best_slice_order')
cdef char get_best_order({{memviewslice_name}} *mslice, int ndim) noexcept nogil

@cname('__pyx_memoryview_slice_get_size')
cdef Py_ssize_t slice_get_size({{memviewslice_name}} *src, int ndim) noexcept nogil

@cname('__pyx_fill_contig_strides_array')
cdef Py_ssize_t fill_contig_strides_array(
                Py_ssize_t *shape, Py_ssize_t *strides, Py_ssize_t stride,
                int ndim, char order) noexcept nogil

@cname('__pyx_memoryview_copy_data_to_temp')
cdef void *copy_data_to_temp({{memviewslice_name}} *src,
                             {{memviewslice_name}} *tmpslice,
                             char order,
                             int ndim) except NULL nogil

@cname('__pyx_memoryview_err_extents')
cdef int _err_extents(int i, Py_ssize_t extent1,
                             Py_ssize_t extent2) except -1 with gil

@cname('__pyx_memoryview_err_dim')
cdef int _err_dim(PyObject *error, str msg, int dim) except -1 with gil

@cname('__pyx_memoryview_err')
cdef int _err(PyObject *error, str msg) except -1 with gil

@cname('__pyx_memoryview_err_no_memory')
cdef int _err_no_memory() except -1 with gil

@cname('__pyx_memoryview_copy_contents')
cdef int memoryview_copy_contents({{memviewslice_name}} src,
                                  {{memviewslice_name}} dst,
                                  int src_ndim, int dst_ndim,
                                  bint dtype_is_object) except -1 nogil

@cname('__pyx_memoryview_broadcast_leading')
cdef void broadcast_leading({{memviewslice_name}} *mslice,
                            int ndim,
                            int ndim_other) noexcept nogil

@cname('__pyx_memoryview_refcount_copying')
cdef void refcount_copying({{memviewslice_name}} *dst, bint dtype_is_object, int ndim, bint inc) noexcept nogil

@cname('__pyx_memoryview_refcount_objects_in_slice')
cdef void refcount_objects_in_slice(char *data, Py_ssize_t *shape,
                                    Py_ssize_t *strides, int ndim, bint inc) noexcept

@cname('__pyx_memoryview_slice_assign_scalar')
cdef void slice_assign_scalar({{memviewslice_name}} *dst, int ndim,
                              size_t itemsize, void *item,
                              bint dtype_is_object) noexcept nogil

@cname('__pyx_memoryview__slice_assign_scalar')
cdef void _slice_assign_scalar(char *data, Py_ssize_t *shape,
                              Py_ssize_t *strides, int ndim,
                              size_t itemsize, void *item) noexcept nogil
