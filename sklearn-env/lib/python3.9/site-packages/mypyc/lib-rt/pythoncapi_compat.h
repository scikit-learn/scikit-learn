// Header file providing new functions of the Python C API to old Python
// versions.
//
// File distributed under the MIT license.
// Copyright Contributors to the pythoncapi_compat project.
//
// Homepage:
// https://github.com/pythoncapi/pythoncapi_compat
//
// Latest version:
// https://raw.githubusercontent.com/pythoncapi/pythoncapi_compat/master/pythoncapi_compat.h
//
// SPDX-License-Identifier: MIT

#ifndef PYTHONCAPI_COMPAT
#define PYTHONCAPI_COMPAT

#ifdef __cplusplus
extern "C" {
#endif

#include <Python.h>
#include "frameobject.h"          // PyFrameObject, PyFrame_GetBack()


// Compatibility with Visual Studio 2013 and older which don't support
// the inline keyword in C (only in C++): use __inline instead.
#if (defined(_MSC_VER) && _MSC_VER < 1900 \
     && !defined(__cplusplus) && !defined(inline))
#  define inline __inline
#  define PYTHONCAPI_COMPAT_MSC_INLINE
   // These two macros are undefined at the end of this file
#endif


// Cast argument to PyObject* type.
#ifndef _PyObject_CAST
#  define _PyObject_CAST(op) ((PyObject*)(op))
#endif
#ifndef _PyObject_CAST_CONST
#  define _PyObject_CAST_CONST(op) ((const PyObject*)(op))
#endif


// bpo-42262 added Py_NewRef() to Python 3.10.0a3
#if PY_VERSION_HEX < 0x030A00A3 && !defined(Py_NewRef)
static inline PyObject* _Py_NewRef(PyObject *obj)
{
    Py_INCREF(obj);
    return obj;
}
#define Py_NewRef(obj) _Py_NewRef(_PyObject_CAST(obj))
#endif


// bpo-42262 added Py_XNewRef() to Python 3.10.0a3
#if PY_VERSION_HEX < 0x030A00A3 && !defined(Py_XNewRef)
static inline PyObject* _Py_XNewRef(PyObject *obj)
{
    Py_XINCREF(obj);
    return obj;
}
#define Py_XNewRef(obj) _Py_XNewRef(_PyObject_CAST(obj))
#endif


// See https://bugs.python.org/issue42522
#if !defined(_Py_StealRef)
static inline PyObject* __Py_StealRef(PyObject *obj)
{
    Py_DECREF(obj);
    return obj;
}
#define _Py_StealRef(obj) __Py_StealRef(_PyObject_CAST(obj))
#endif


// See https://bugs.python.org/issue42522
#if !defined(_Py_XStealRef)
static inline PyObject* __Py_XStealRef(PyObject *obj)
{
    Py_XDECREF(obj);
    return obj;
}
#define _Py_XStealRef(obj) __Py_XStealRef(_PyObject_CAST(obj))
#endif


// bpo-39573 added Py_SET_REFCNT() to Python 3.9.0a4
#if PY_VERSION_HEX < 0x030900A4 && !defined(Py_SET_REFCNT)
static inline void _Py_SET_REFCNT(PyObject *ob, Py_ssize_t refcnt)
{
    ob->ob_refcnt = refcnt;
}
#define Py_SET_REFCNT(ob, refcnt) _Py_SET_REFCNT(_PyObject_CAST(ob), refcnt)
#endif


// Py_SETREF() and Py_XSETREF() were added to Python 3.5.2.
// It is excluded from the limited C API.
#if (PY_VERSION_HEX < 0x03050200 && !defined(Py_SETREF)) && !defined(Py_LIMITED_API)
#define Py_SETREF(op, op2)                      \
    do {                                        \
        PyObject *_py_tmp = _PyObject_CAST(op); \
        (op) = (op2);                           \
        Py_DECREF(_py_tmp);                     \
    } while (0)

#define Py_XSETREF(op, op2)                     \
    do {                                        \
        PyObject *_py_tmp = _PyObject_CAST(op); \
        (op) = (op2);                           \
        Py_XDECREF(_py_tmp);                    \
    } while (0)
#endif


// bpo-43753 added Py_Is(), Py_IsNone(), Py_IsTrue() and Py_IsFalse()
// to Python 3.10.0b1.
#if PY_VERSION_HEX < 0x030A00B1 && !defined(Py_Is)
#  define Py_Is(x, y) ((x) == (y))
#endif
#if PY_VERSION_HEX < 0x030A00B1 && !defined(Py_IsNone)
#  define Py_IsNone(x) Py_Is(x, Py_None)
#endif
#if PY_VERSION_HEX < 0x030A00B1 && !defined(Py_IsTrue)
#  define Py_IsTrue(x) Py_Is(x, Py_True)
#endif
#if PY_VERSION_HEX < 0x030A00B1 && !defined(Py_IsFalse)
#  define Py_IsFalse(x) Py_Is(x, Py_False)
#endif


// bpo-39573 added Py_SET_TYPE() to Python 3.9.0a4
#if PY_VERSION_HEX < 0x030900A4 && !defined(Py_SET_TYPE)
static inline void
_Py_SET_TYPE(PyObject *ob, PyTypeObject *type)
{
    ob->ob_type = type;
}
#define Py_SET_TYPE(ob, type) _Py_SET_TYPE(_PyObject_CAST(ob), type)
#endif


// bpo-39573 added Py_SET_SIZE() to Python 3.9.0a4
#if PY_VERSION_HEX < 0x030900A4 && !defined(Py_SET_SIZE)
static inline void
_Py_SET_SIZE(PyVarObject *ob, Py_ssize_t size)
{
    ob->ob_size = size;
}
#define Py_SET_SIZE(ob, size) _Py_SET_SIZE((PyVarObject*)(ob), size)
#endif


// bpo-40421 added PyFrame_GetCode() to Python 3.9.0b1
#if PY_VERSION_HEX < 0x030900B1
static inline PyCodeObject*
PyFrame_GetCode(PyFrameObject *frame)
{
    assert(frame != NULL);
    assert(frame->f_code != NULL);
    return (PyCodeObject*)Py_NewRef(frame->f_code);
}
#endif

static inline PyCodeObject*
_PyFrame_GetCodeBorrow(PyFrameObject *frame)
{
    return (PyCodeObject *)_Py_StealRef(PyFrame_GetCode(frame));
}


// bpo-40421 added PyFrame_GetCode() to Python 3.9.0b1
#if PY_VERSION_HEX < 0x030900B1 && !defined(PYPY_VERSION)
static inline PyFrameObject*
PyFrame_GetBack(PyFrameObject *frame)
{
    assert(frame != NULL);
    return (PyFrameObject*)Py_XNewRef(frame->f_back);
}
#endif

#if !defined(PYPY_VERSION)
static inline PyFrameObject*
_PyFrame_GetBackBorrow(PyFrameObject *frame)
{
    return (PyFrameObject *)_Py_XStealRef(PyFrame_GetBack(frame));
}
#endif


// bpo-39947 added PyThreadState_GetInterpreter() to Python 3.9.0a5
#if PY_VERSION_HEX < 0x030900A5
static inline PyInterpreterState *
PyThreadState_GetInterpreter(PyThreadState *tstate)
{
    assert(tstate != NULL);
    return tstate->interp;
}
#endif


// bpo-40429 added PyThreadState_GetFrame() to Python 3.9.0b1
#if PY_VERSION_HEX < 0x030900B1 && !defined(PYPY_VERSION)
static inline PyFrameObject*
PyThreadState_GetFrame(PyThreadState *tstate)
{
    assert(tstate != NULL);
    return (PyFrameObject *)Py_XNewRef(tstate->frame);
}
#endif

#if !defined(PYPY_VERSION)
static inline PyFrameObject*
_PyThreadState_GetFrameBorrow(PyThreadState *tstate)
{
    return (PyFrameObject *)_Py_XStealRef(PyThreadState_GetFrame(tstate));
}
#endif


// bpo-39947 added PyInterpreterState_Get() to Python 3.9.0a5
#if PY_VERSION_HEX < 0x030900A5
static inline PyInterpreterState *
PyInterpreterState_Get(void)
{
    PyThreadState *tstate;
    PyInterpreterState *interp;

    tstate = PyThreadState_GET();
    if (tstate == NULL) {
        Py_FatalError("GIL released (tstate is NULL)");
    }
    interp = tstate->interp;
    if (interp == NULL) {
        Py_FatalError("no current interpreter");
    }
    return interp;
}
#endif


// bpo-39947 added PyInterpreterState_Get() to Python 3.9.0a6
#if 0x030700A1 <= PY_VERSION_HEX && PY_VERSION_HEX < 0x030900A6 && !defined(PYPY_VERSION)
static inline uint64_t
PyThreadState_GetID(PyThreadState *tstate)
{
    assert(tstate != NULL);
    return tstate->id;
}
#endif

// bpo-43760 added PyThreadState_EnterTracing() to Python 3.11.0a2
#if PY_VERSION_HEX < 0x030B00A2 && !defined(PYPY_VERSION)
static inline void PyThreadState_EnterTracing(PyThreadState *tstate)
{
    tstate->tracing++;
#if PY_VERSION_HEX >= 0x030A00A1
    tstate->cframe->use_tracing = 0;
#else
    tstate->use_tracing = 0;
#endif
}
#endif

// bpo-43760 added PyThreadState_LeaveTracing() to Python 3.11.0a2
#if PY_VERSION_HEX < 0x030B00A2 && !defined(PYPY_VERSION)
static inline void PyThreadState_LeaveTracing(PyThreadState *tstate)
{
    tstate->tracing--;
    int use_tracing = (tstate->c_tracefunc != NULL
                       || tstate->c_profilefunc != NULL);
#if PY_VERSION_HEX >= 0x030A00A1
    tstate->cframe->use_tracing = use_tracing;
#else
    tstate->use_tracing = use_tracing;
#endif
}
#endif


// bpo-37194 added PyObject_CallNoArgs() to Python 3.9.0a1
#if PY_VERSION_HEX < 0x030900A1
static inline PyObject*
PyObject_CallNoArgs(PyObject *func)
{
    return PyObject_CallFunctionObjArgs(func, NULL);
}
#endif


// bpo-39245 made PyObject_CallOneArg() public (previously called
// _PyObject_CallOneArg) in Python 3.9.0a4
#if PY_VERSION_HEX < 0x030900A4
static inline PyObject*
PyObject_CallOneArg(PyObject *func, PyObject *arg)
{
    return PyObject_CallFunctionObjArgs(func, arg, NULL);
}
#endif


// bpo-1635741 added PyModule_AddObjectRef() to Python 3.10.0a3
#if PY_VERSION_HEX < 0x030A00A3
static inline int
PyModule_AddObjectRef(PyObject *module, const char *name, PyObject *value)
{
    int res;
    Py_XINCREF(value);
    res = PyModule_AddObject(module, name, value);
    if (res < 0) {
        Py_XDECREF(value);
    }
    return res;
}
#endif


// bpo-40024 added PyModule_AddType() to Python 3.9.0a5
#if PY_VERSION_HEX < 0x030900A5
static inline int
PyModule_AddType(PyObject *module, PyTypeObject *type)
{
    const char *name, *dot;

    if (PyType_Ready(type) < 0) {
        return -1;
    }

    // inline _PyType_Name()
    name = type->tp_name;
    assert(name != NULL);
    dot = strrchr(name, '.');
    if (dot != NULL) {
        name = dot + 1;
    }

    return PyModule_AddObjectRef(module, name, (PyObject *)type);
}
#endif


// bpo-40241 added PyObject_GC_IsTracked() to Python 3.9.0a6.
// bpo-4688 added _PyObject_GC_IS_TRACKED() to Python 2.7.0a2.
#if PY_VERSION_HEX < 0x030900A6 && !defined(PYPY_VERSION)
static inline int
PyObject_GC_IsTracked(PyObject* obj)
{
    return (PyObject_IS_GC(obj) && _PyObject_GC_IS_TRACKED(obj));
}
#endif

// bpo-40241 added PyObject_GC_IsFinalized() to Python 3.9.0a6.
// bpo-18112 added _PyGCHead_FINALIZED() to Python 3.4.0 final.
#if PY_VERSION_HEX < 0x030900A6 && PY_VERSION_HEX >= 0x030400F0 && !defined(PYPY_VERSION)
static inline int
PyObject_GC_IsFinalized(PyObject *obj)
{
    return (PyObject_IS_GC(obj) && _PyGCHead_FINALIZED((PyGC_Head *)(obj)-1));
}
#endif


// bpo-39573 added Py_IS_TYPE() to Python 3.9.0a4
#if PY_VERSION_HEX < 0x030900A4 && !defined(Py_IS_TYPE)
static inline int
_Py_IS_TYPE(const PyObject *ob, const PyTypeObject *type) {
    return ob->ob_type == type;
}
#define Py_IS_TYPE(ob, type) _Py_IS_TYPE(_PyObject_CAST_CONST(ob), type)
#endif


// Py_UNUSED() was added to Python 3.4.0b2.
#if PY_VERSION_HEX < 0x030400B2 && !defined(Py_UNUSED)
#  if defined(__GNUC__) || defined(__clang__)
#    define Py_UNUSED(name) _unused_ ## name __attribute__((unused))
#  else
#    define Py_UNUSED(name) _unused_ ## name
#  endif
#endif


#ifdef PYTHONCAPI_COMPAT_MSC_INLINE
#  undef inline
#  undef PYTHONCAPI_COMPAT_MSC_INLINE
#endif

#ifdef __cplusplus
}
#endif
#endif  // PYTHONCAPI_COMPAT
