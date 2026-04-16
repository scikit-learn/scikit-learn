from .object cimport PyObject

cdef extern from *:
    # Backport PyList_GetItemRef, PyList_Extend and PyList_Clear
    """
    #if __PYX_LIMITED_VERSION_HEX < 0x030d0000
    static CYTHON_INLINE PyObject *
    __Pyx_CAPI_PyList_GetItemRef(PyObject *list, Py_ssize_t index)
    {
        PyObject *item = PyList_GetItem(list, index);
        Py_XINCREF(item);
        return item;
    }
    #else
    #define __Pyx_CAPI_PyList_GetItemRef PyList_GetItemRef
    #endif

    #if CYTHON_COMPILING_IN_LIMITED_API || PY_VERSION_HEX < 0x030d0000
    static CYTHON_INLINE int
    __Pyx_CAPI_PyList_Extend(PyObject *list, PyObject *iterable)
    {
        return PyList_SetSlice(list, PY_SSIZE_T_MAX, PY_SSIZE_T_MAX, iterable);
    }

    static CYTHON_INLINE int
    __Pyx_CAPI_PyList_Clear(PyObject *list)
    {
        return PyList_SetSlice(list, 0, PY_SSIZE_T_MAX, NULL);
    }
    #else
    #define __Pyx_CAPI_PyList_Extend PyList_Extend
    #define __Pyx_CAPI_PyList_Clear PyList_Clear
    #endif
    """

cdef extern from "Python.h":

    ############################################################################
    # Lists
    ############################################################################
    list PyList_New(Py_ssize_t len)
    # Return a new list of length len on success, or NULL on failure.
    #
    # Note: If length is greater than zero, the returned list object's
    # items are set to NULL. Thus you cannot use abstract API
    # functions such as PySequence_SetItem() or expose the object to
    # Python code before setting all items to a real object with
    # PyList_SetItem().

    bint PyList_Check(object p)
    # Return true if p is a list object or an instance of a subtype of
    # the list type.

    bint PyList_CheckExact(object p)
    # Return true if p is a list object, but not an instance of a
    # subtype of the list type.

    Py_ssize_t PyList_Size(object list) except -1
    # Return the length of the list object in list; this is equivalent
    # to "len(list)" on a list object.

    Py_ssize_t PyList_GET_SIZE(object list)
    # Macro form of PyList_Size() without error checking.

    object PyList_GetItemRef "__Pyx_CAPI_PyList_GetItemRef" (object list, Py_ssize_t index)
    # Return value: New reference.
    # Return the object at position index in the list pointed to by list.
    # The position must be non-negative; indexing from the end of the
    # list is not supported. If index is out of bounds (<0 or >=len(list)),
    # return NULL and set an IndexError exception.

    PyObject* PyList_GetItem(object list, Py_ssize_t index) except NULL
    # Return value: Borrowed reference.
    # Like PyList_GetItemRef(), but returns a borrowed reference instead of
    # a strong reference.

    PyObject* PyList_GET_ITEM(object list, Py_ssize_t i)
    # Return value: Borrowed reference.
    # Macro form of PyList_GetItem() without error checking.

    int PyList_SetItem(object list, Py_ssize_t index, object item) except -1
    # Set the item at index index in list to item. Return 0 on success
    # or -1 on failure.
    #
    # WARNING: This function _steals_ a reference to item and discards a
    # reference to an item already in the list at the affected position.

    void PyList_SET_ITEM(object list, Py_ssize_t i, object o)
    # Macro form of PyList_SetItem() without error checking. This is
    # normally only used to fill in new lists where there is no
    # previous content.
    #
    # WARNING: This function _steals_ a reference to item, and, unlike
    # PyList_SetItem(), does not discard a reference to any item that
    # it being replaced; any reference in list at position i will be *leaked*.

    int PyList_Insert(object list, Py_ssize_t index, object item) except -1
    # Insert the item item into list list in front of index
    # index. Return 0 if successful; return -1 and set an exception if
    # unsuccessful. Analogous to list.insert(index, item).

    int PyList_Append(object list, object item) except -1
    # Append the object item at the end of list list. Return 0 if
    # successful; return -1 and set an exception if
    # unsuccessful. Analogous to list.append(item).

    list PyList_GetSlice(object list, Py_ssize_t low, Py_ssize_t high)
    # Return value: New reference.
    # Return a list of the objects in list containing the objects
    # between low and high. Return NULL and set an exception if
    # unsuccessful. Analogous to list[low:high].

    int PyList_SetSlice(object list, Py_ssize_t low, Py_ssize_t high, object itemlist) except -1
    # Set the slice of list between low and high to the contents of
    # itemlist. Analogous to list[low:high] = itemlist. The itemlist
    # may be NULL, indicating the assignment of an empty list (slice
    # deletion). Return 0 on success, -1 on failure.

    int PyList_Extend "__Pyx_CAPI_PyList_Extend" (object list, object iterable) except -1
    # Extend list with the contents of iterable. This is the same as
    # PyList_SetSlice(list, PY_SSIZE_T_MAX, PY_SSIZE_T_MAX, iterable)
    # and analogous to list.extend(iterable) or list += iterable.
    # Raise an exception and return -1 if list is not a list object.
    # Return 0 on success.

    int PyList_Clear "__Pyx_CAPI_PyList_Clear" (object list) except -1
    # Remove all items from list. This is the same as
    # PyList_SetSlice(list, 0, PY_SSIZE_T_MAX, NULL) and analogous
    # to list.clear() or del list[:].
    # Raise an exception and return -1 if list is not a list object.
    # Return 0 on success.

    int PyList_Sort(object list) except -1
    # Sort the items of list in place. Return 0 on success, -1 on
    # failure. This is equivalent to "list.sort()".

    int PyList_Reverse(object list) except -1
    # Reverse the items of list in place. Return 0 on success, -1 on
    # failure. This is the equivalent of "list.reverse()".

    tuple PyList_AsTuple(object list)
    # Return value: New reference.
    # Return a new tuple object containing the contents of list;
    # equivalent to "tuple(list)".
