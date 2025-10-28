from .object cimport PyObject
from .object cimport PyTypeObject, Py_TYPE  # legacy imports for re-export

cdef extern from "Python.h":
    #####################################################################
    # 3. Reference Counts
    #####################################################################
    # The macros in this section are used for managing reference counts of Python objects.
    void Py_INCREF(object o)
    # Increment the reference count for object o. The object must not
    # be NULL; if you aren't sure that it isn't NULL, use
    # Py_XINCREF().

    void Py_XINCREF(PyObject* o)
    # Increment the reference count for object o. The object may be NULL, in which case the macro has no effect.

    void Py_DECREF(object o)
    # Decrement the reference count for object o. The object must not
    # be NULL; if you aren't sure that it isn't NULL, use
    # Py_XDECREF(). If the reference count reaches zero, the object's
    # type's deallocation function (which must not be NULL) is
    # invoked.

    # Warning: The deallocation function can cause arbitrary Python
    # code to be invoked (e.g. when a class instance with a __del__()
    # method is deallocated). While exceptions in such code are not
    # propagated, the executed code has free access to all Python
    # global variables. This means that any object that is reachable
    # from a global variable should be in a consistent state before
    # Py_DECREF() is invoked. For example, code to delete an object
    # from a list should copy a reference to the deleted object in a
    # temporary variable, update the list data structure, and then
    # call Py_DECREF() for the temporary variable.

    void Py_XDECREF(PyObject* o)
    # Decrement the reference count for object o. The object may be
    # NULL, in which case the macro has no effect; otherwise the
    # effect is the same as for Py_DECREF(), and the same warning
    # applies.

    void Py_CLEAR(PyObject* o)
    # Decrement the reference count for object o. The object may be
    # NULL, in which case the macro has no effect; otherwise the
    # effect is the same as for Py_DECREF(), except that the argument
    # is also set to NULL. The warning for Py_DECREF() does not apply
    # with respect to the object passed because the macro carefully
    # uses a temporary variable and sets the argument to NULL before
    # decrementing its reference count.
    # It is a good idea to use this macro whenever decrementing the
    # value of a variable that might be traversed during garbage
    # collection.

    Py_ssize_t Py_REFCNT(object o)
    # Get the reference count of the Python object o.

    # Note that the returned value may not actually reflect how many
    # references to the object are actually held. For example, some
    # objects are “immortal” and have a very high refcount that does not
    # reflect the actual number of references. Consequently, do not rely
    # on the returned value to be accurate, other than a value of 0 or
    # 1.

    Py_ssize_t _Py_REFCNT "Py_REFCNT" (PyObject *ptr)
    # Get the reference count for the PyObject pointer ptr.

    # This is useful when it would be awkward to create an owned reference just
    # to get the reference count. See the note for Py_REFCNT above about the
    # accuracy of reference counts.

    int PyUnstable_Object_EnableDeferredRefcount(object o)
    # Enable deferred reference counting on obj, if supported by the runtime.
    # In the free-threaded build, this allows the interpreter to avoid reference count
    # adjustments to obj, which may improve multi-threaded performance.
    # The tradeoff is that obj will only be deallocated by the tracing garbage collector,
    # and not when the interpreter no longer has any references to it.
    #
    # This function returns 1 if deferred reference counting is enabled on obj,
    # and 0 if deferred reference counting is not supported or if the hint was ignored
    # by the interpreter, such as when deferred reference counting is already enabled on obj.
    # This function is thread-safe, and cannot fail.
    #
    # This function does nothing on builds with the GIL enabled, which do not support
    # deferred reference counting. This also does nothing if obj is not an object tracked
    # by the garbage collector (see gc.is_tracked() and PyObject_GC_IsTracked()).
    #
    # This function is intended to be used soon after obj is created, by the code
    # that creates it, such as in the object’s tp_new slot.
    #
    # Added in CPython 3.14.

    int PyUnstable_Object_IsUniqueReferencedTemporary(object o)
    # Check if obj is a unique temporary object.
    # Returns 1 if obj is known to be a unique temporary object, and 0 otherwise.
    # This function cannot fail, but the check is conservative, and may return 0
    # in some cases even if obj is a unique temporary object.
    #
    # If an object is a unique temporary, it is guaranteed that the current code
    # has the only reference to the object. For arguments to C functions, this should
    # be used instead of checking if the reference count is 1.
    # Starting with Python 3.14, the interpreter internally avoids some reference count
    # modifications when loading objects onto the operands stack by borrowing references
    # when possible, which means that a reference count of 1 by itself does not guarantee
    # that a function argument uniquely referenced.
    #
    # Added in CPython 3.14.

    int PyUnstable_IsImmortal(object o)
    # This function returns non-zero if obj is immortal, and zero otherwise.
    # This function cannot fail.
    #
    # Added in CPython 3.14.

    void PyUnstable_EnableTryIncRef(object o)
    # Enables subsequent uses of PyUnstable_TryIncRef() on obj.
    # The caller must hold a strong reference to obj when calling this.
    #
    # Added in CPython 3.14.

    bint PyUnstable_TryIncRef(PyObject *o)
    # Increments the reference count of obj if it is not zero.
    # Returns 1 if the object’s reference count was successfully incremented.
    # Otherwise, this function returns 0.
    #
    # PyUnstable_EnableTryIncRef() must have been called earlier on obj
    # or this function may spuriously return 0 in the free threading build.
    #
    # Added in CPython 3.14.

    int PyUnstable_Object_IsUniquelyReferenced(object o)
    # Determine if op only has one reference.
    #
    # On GIL-enabled builds, this function is equivalent to Py_REFCNT(op) == 1.
    #
    # On a free threaded build, this checks if op’s reference count is equal
    # to one and additionally checks if op is only used by this thread.
    # Py_REFCNT(op) == 1 is not thread-safe on free threaded builds; prefer this function.
    #
    # The caller must hold an attached thread state, despite the fact that this function
    # doesn’t call into the Python interpreter. This function cannot fail.
    #
    # Added in CPython 3.14.
