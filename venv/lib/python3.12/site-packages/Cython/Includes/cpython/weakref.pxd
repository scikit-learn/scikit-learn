from .object cimport PyObject

cdef extern from *:
    # Backport PyWeakref_GetRef to Python < 3.13
    """
    #if PY_VERSION_HEX < 0x030d0000
    static CYTHON_INLINE int __Pyx_PyWeakref_GetRef(PyObject *ref, PyObject **pobj)
    {
        PyObject *obj = PyWeakref_GetObject(ref);
        if (obj == NULL) {
            // SystemError if ref is NULL
            *pobj = NULL;
            return -1;
        }
        if (obj == Py_None) {
            *pobj = NULL;
            return 0;
        }
        Py_INCREF(obj);
        *pobj = obj;
        return 1;
    }
    #else
    #define __Pyx_PyWeakref_GetRef PyWeakref_GetRef
    #endif
    """
    bint PyWeakref_GetRef "__Pyx_PyWeakref_GetRef" (object ref, PyObject** pobj)
    # Get a strong reference to the referenced object from a weak reference,
    # ref, into *pobj.
    # - On success, set *pobj to a new strong reference to the referenced
    #   object and return 1.
    # - If the reference is dead, set *pobj to NULL and return 0.
    # - On error, raise an exception and return -1.

cdef extern from "Python.h":

    bint PyWeakref_Check(object ob)
    # Return true if ob is either a reference or proxy object.

    bint PyWeakref_CheckRef(object ob)
    # Return true if ob is a reference object.

    bint PyWeakref_CheckProxy(ob)
    # Return true if *ob* is a proxy object.

    object PyWeakref_NewRef(object ob, object callback)
    # Return a weak reference object for the object ob.  This will
    # always return a new reference, but is not guaranteed to create a
    # new object; an existing reference object may be returned.  The
    # second parameter, callback, can be a callable object that
    # receives notification when ob is garbage collected; it should
    # accept a single parameter, which will be the weak reference
    # object itself. callback may also be None or NULL.  If ob is not
    # a weakly-referencable object, or if callback is not callable,
    # None, or NULL, this will return NULL and raise TypeError.

    object PyWeakref_NewProxy(object ob, object callback)
    # Return a weak reference proxy object for the object ob.  This
    # will always return a new reference, but is not guaranteed to
    # create a new object; an existing proxy object may be returned.
    # The second parameter, callback, can be a callable object that
    # receives notification when ob is garbage collected; it should
    # accept a single parameter, which will be the weak reference
    # object itself. callback may also be None or NULL.  If ob is not
    # a weakly-referencable object, or if callback is not callable,
    # None, or NULL, this will return NULL and raise TypeError.

    PyObject* PyWeakref_GetObject(object ref) except NULL
    # Return a borrowed reference to the referenced object from a weak
    # reference, ref.  If the referent is no longer live, returns None.
    # Deprecated since Python 3.13, will be removed in version 3.15: Use
    # PyWeakref_GetRef instead.

    PyObject* PyWeakref_GET_OBJECT(object ref)
    # Similar to PyWeakref_GetObject, but implemented as a macro that
    # does no error checking.
    # Deprecated since Python 3.13, will be removed in version 3.15: Use
    # PyWeakref_GetRef instead.
