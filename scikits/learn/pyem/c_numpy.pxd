# :Author:    Robert Kern
# :Copyright: 2004, Enthought, Inc.
# :License:   BSD Style


cdef extern from "numpy/arrayobject.h":
    ctypedef enum PyArray_TYPES:
        PyArray_BOOL
        PyArray_BYTE
        PyArray_UBYTE
        PyArray_SHORT
        PyArray_USHORT 
        PyArray_INT
        PyArray_UINT 
        PyArray_LONG
        PyArray_ULONG
        PyArray_LONGLONG
        PyArray_ULONGLONG
        PyArray_FLOAT
        PyArray_DOUBLE 
        PyArray_LONGDOUBLE
        PyArray_CFLOAT
        PyArray_CDOUBLE
        PyArray_CLONGDOUBLE
        PyArray_OBJECT
        PyArray_STRING
        PyArray_UNICODE
        PyArray_VOID
        PyArray_NTYPES
        PyArray_NOTYPE

    ctypedef int intp 

    ctypedef extern class numpy.dtype [object PyArray_Descr]:
        cdef int type_num, elsize, alignment
        cdef char type, kind, byteorder, hasobject
        cdef object fields, typeobj

    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef object base
        cdef dtype descr
        cdef int flags

    ndarray PyArray_SimpleNew(int ndims, intp* dims, int item_type)
    int PyArray_Check(object obj)
    ndarray PyArray_ContiguousFromObject(object obj, PyArray_TYPES type, 
        int mindim, int maxdim)
    intp PyArray_SIZE(ndarray arr)
    void *PyArray_DATA(ndarray arr)
    ndarray PyArray_FromAny(object obj, dtype newtype, int mindim, int maxdim,
		    int requirements, object context)
    ndarray PyArray_NewFromDescr(object subtype, dtype newtype, int nd, intp* dims,
		    intp* strides, void* data, int flags, object parent)

    void import_array()
