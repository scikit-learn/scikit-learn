# cython: language_level=3, auto_pickle=False

from cpython.ref cimport PyObject, Py_INCREF, Py_CLEAR, Py_XDECREF, Py_XINCREF
from cpython.exc cimport PyErr_Fetch, PyErr_Restore
from cpython.pystate cimport PyThreadState_Get

cimport cython

cdef extern from *:
    """
    #if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING
    #define __Pyx_refnanny_mutex PyMutex
    static CYTHON_INLINE void __Pyx_refnanny_lock_acquire(PyMutex *lock) {
        PyMutex_Lock(lock);
    }

    static CYTHON_INLINE void __Pyx_refnanny_lock_release(PyMutex *lock) {
        PyMutex_Unlock(lock);
    }
    #else
    #define __Pyx_refnanny_mutex void*
    #define __Pyx_refnanny_lock_acquire(lock)
    #define __Pyx_refnanny_lock_release(lock)
    #endif
    """
    ctypedef void *__Pyx_refnanny_mutex
    void __Pyx_refnanny_lock_acquire(__Pyx_refnanny_mutex *lock)
    void __Pyx_refnanny_lock_release(__Pyx_refnanny_mutex *lock)

loglevel = 0
reflog = []

cdef int log(int level, action, obj, lineno) except -1:
    if (<int> loglevel) >= level:
        if reflog is None:
            # can happen during finalisation
            return 0
        reflog.append((lineno, action, id(obj)))
    return 0

LOG_NONE, LOG_ALL = range(2)
cdef int _LOG_NONE = LOG_NONE
cdef int _LOG_ALL = LOG_ALL

cdef object NO_REFS = (0, None)


@cython.final
cdef class Context(object):
    cdef readonly object name, filename
    cdef readonly dict refs
    cdef readonly list errors
    cdef readonly Py_ssize_t start
    cdef __Pyx_refnanny_mutex lock

    def __cinit__(self, name, line=0, filename=None):
        self.name = name
        self.start = line
        self.filename = filename
        self.refs = {} # id -> (count, [lineno])
        self.errors = []

    cdef void acquire_lock(self) noexcept:
        __Pyx_refnanny_lock_acquire(&self.lock)

    cdef void release_lock(self) noexcept:
        __Pyx_refnanny_lock_release(&self.lock)

    cdef int regref(self, obj, Py_ssize_t lineno, bint is_null) except -1:
        log(_LOG_ALL, u'regref', u"<NULL>" if is_null else obj, lineno)
        if is_null:
            self.errors.append(f"NULL argument on line {lineno}")
            return 0
        id_ = id(obj)
        count, linenumbers = self.refs.get(id_, NO_REFS)
        if linenumbers is None:
            linenumbers = []
        self.refs[id_] = (count + 1, linenumbers)
        linenumbers.append(lineno)
        return 0

    cdef bint delref(self, obj, Py_ssize_t lineno, bint is_null) except -1:
        # returns whether it is ok to do the decref operation
        log(_LOG_ALL, u'delref', u"<NULL>" if is_null else obj, lineno)
        if is_null:
            self.errors.append(f"NULL argument on line {lineno}")
            return False
        id_ = id(obj)
        count, linenumbers = self.refs.get(id_, NO_REFS)
        if linenumbers is None:
            linenumbers = []
        if count == 0:
            self.errors.append(f"Too many decrefs on line {lineno}, reference acquired on lines {linenumbers!r}")
            return False
        if count == 1:
            del self.refs[id_]
        else:
            self.refs[id_] = (count - 1, linenumbers)
        return True

    cdef end(self):
        if self.refs:
            msg = u"References leaked:"
            for count, linenos in self.refs.values():
                msg += f"\n  ({count}) acquired on lines: {u', '.join([f'{x}' for x in linenos])}"
            self.errors.append(msg)
        return u"\n".join([f'REFNANNY: {error}' for error in self.errors]) if self.errors else None


cdef void report_unraisable(filename, Py_ssize_t lineno, object e=None) noexcept:
    try:
        if e is None:
            import sys
            e = sys.exc_info()[1]
        print(f"refnanny raised an exception from {filename}:{lineno}: {e}")
    finally:
        return  # We absolutely cannot exit with an exception


# All Python operations must happen after any existing
# exception has been fetched, in case we are called from
# exception-handling code.

cdef PyObject* SetupContext(char* funcname, Py_ssize_t lineno, char* filename) except NULL:
    if Context is None:
        # Context may be None during finalize phase.
        # In that case, we don't want to be doing anything fancy
        # like caching and resetting exceptions.
        return NULL
    cdef (PyObject*) type = NULL, value = NULL, tb = NULL, result = NULL
    PyThreadState_Get()  # Check that we hold the GIL
    PyErr_Fetch(&type, &value, &tb)
    try:
        ctx = Context.__new__(Context, funcname, lineno, filename)
        Py_INCREF(ctx)
        result = <PyObject*>ctx
    except Exception, e:
        report_unraisable(filename, lineno, e)
    PyErr_Restore(type, value, tb)
    return result

cdef void GOTREF(PyObject* _ctx, PyObject* p_obj, Py_ssize_t lineno):
    if _ctx == NULL: return
    cdef (PyObject*) type = NULL, value = NULL, tb = NULL
    cdef Context ctx = <Context> _ctx
    ctx.acquire_lock()
    PyErr_Fetch(&type, &value, &tb)
    try:
        ctx.regref(
            <object>p_obj if p_obj is not NULL else None,
            lineno,
            is_null=p_obj is NULL,
        )
    except:
        report_unraisable(ctx.filename, lineno=ctx.start)
    finally:
        PyErr_Restore(type, value, tb)
        ctx.release_lock()
        return  # swallow any exceptions

cdef bint GIVEREF_and_report(PyObject* _ctx, PyObject* p_obj, Py_ssize_t lineno):
    if _ctx == NULL: return 1
    cdef (PyObject*) type = NULL, value = NULL, tb = NULL
    cdef bint decref_ok = False
    cdef Context ctx = <Context> _ctx
    ctx.acquire_lock()
    PyErr_Fetch(&type, &value, &tb)
    try:
        decref_ok = ctx.delref(
            <object>p_obj if p_obj is not NULL else None,
            lineno,
            is_null=p_obj is NULL,
        )
    except:
        report_unraisable(ctx.filename, lineno=ctx.start)
    finally:
        PyErr_Restore(type, value, tb)
        ctx.release_lock()
        return decref_ok  # swallow any exceptions

cdef void GIVEREF(PyObject* ctx, PyObject* p_obj, Py_ssize_t lineno):
    GIVEREF_and_report(ctx, p_obj, lineno)

cdef void INCREF(PyObject* ctx, PyObject* obj, Py_ssize_t lineno):
    Py_XINCREF(obj)
    PyThreadState_Get()  # Check that we hold the GIL
    GOTREF(ctx, obj, lineno)

cdef void DECREF(PyObject* ctx, PyObject* obj, Py_ssize_t lineno):
    if GIVEREF_and_report(ctx, obj, lineno):
        Py_XDECREF(obj)
    PyThreadState_Get()  # Check that we hold the GIL

cdef void FinishContext(PyObject** ctx):
    if ctx == NULL or ctx[0] == NULL: return
    cdef (PyObject*) type = NULL, value = NULL, tb = NULL
    cdef object errors = None
    cdef Context context
    PyThreadState_Get()  # Check that we hold the GIL
    PyErr_Fetch(&type, &value, &tb)
    try:
        context = <Context>ctx[0]
        errors = context.end()
        if errors:
            print(f"{context.filename.decode('latin1')}: {context.name.decode('latin1')}()")
            print(errors)
        context = None
    except:
        report_unraisable(
            context.filename if context is not None else None,
            lineno=context.start if context is not None else 0,
        )
    finally:
        Py_CLEAR(ctx[0])
        PyErr_Restore(type, value, tb)
        return  # swallow any exceptions

ctypedef struct RefNannyAPIStruct:
    void (*INCREF)(PyObject*, PyObject*, Py_ssize_t)
    void (*DECREF)(PyObject*, PyObject*, Py_ssize_t)
    void (*GOTREF)(PyObject*, PyObject*, Py_ssize_t)
    void (*GIVEREF)(PyObject*, PyObject*, Py_ssize_t)
    PyObject* (*SetupContext)(char*, Py_ssize_t, char*) except NULL
    void (*FinishContext)(PyObject**)

cdef RefNannyAPIStruct api
api.INCREF = INCREF
api.DECREF =  DECREF
api.GOTREF =  GOTREF
api.GIVEREF = GIVEREF
api.SetupContext = SetupContext
api.FinishContext = FinishContext

cdef extern from "Python.h":
    object PyLong_FromVoidPtr(void*)

RefNannyAPI = PyLong_FromVoidPtr(<void*>&api)
