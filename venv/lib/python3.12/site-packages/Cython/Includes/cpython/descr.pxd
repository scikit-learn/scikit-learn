from .object cimport PyObject, PyTypeObject

cdef extern from "Python.h":
    ctypedef object (*wrapperfunc)(self, args, void* wrapped)
    ctypedef object (*wrapperfunc_kwds)(self, args, void* wrapped, kwds)

    struct wrapperbase:
        char* name
        int offset
        void* function
        wrapperfunc wrapper
        char* doc
        int flags
        PyObject* name_strobj

    int PyWrapperFlag_KEYWORDS

    ctypedef class __builtin__.wrapper_descriptor [object PyWrapperDescrObject]:
        cdef type d_type
        cdef d_name
        cdef wrapperbase* d_base
        cdef void* d_wrapped

    object PyDescr_NewWrapper(PyTypeObject* cls, wrapperbase* wrapper, void* wrapped)

    int PyDescr_IsData(descr)
