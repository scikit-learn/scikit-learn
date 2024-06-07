#ifndef MYPYC_UTIL_H
#define MYPYC_UTIL_H

#include <Python.h>
#include <frameobject.h>
#include <assert.h>

#if defined(__clang__) || defined(__GNUC__)
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)
#define CPy_Unreachable() __builtin_unreachable()
#else
#define likely(x)       (x)
#define unlikely(x)     (x)
#define CPy_Unreachable() abort()
#endif

#if defined(__clang__) || defined(__GNUC__)
#define CPy_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define CPy_NOINLINE __declspec(noinline)
#else
#define CPy_NOINLINE
#endif

// INCREF and DECREF that assert the pointer is not NULL.
// asserts are disabled in release builds so there shouldn't be a perf hit.
// I'm honestly kind of surprised that this isn't done by default.
#define CPy_INCREF(p) do { assert(p); Py_INCREF(p); } while (0)
#define CPy_DECREF(p) do { assert(p); Py_DECREF(p); } while (0)
// Here just for consistency
#define CPy_XDECREF(p) Py_XDECREF(p)

// Tagged integer -- our representation of Python 'int' objects.
// Small enough integers are represented as unboxed integers (shifted
// left by 1); larger integers (larger than 63 bits on a 64-bit
// platform) are stored as a tagged pointer (PyObject *)
// representing a Python int object, with the lowest bit set.
// Tagged integers are always normalized. A small integer *must not*
// have the tag bit set.
typedef size_t CPyTagged;

typedef size_t CPyPtr;

#define CPY_INT_BITS (CHAR_BIT * sizeof(CPyTagged))

#define CPY_TAGGED_MAX (((Py_ssize_t)1 << (CPY_INT_BITS - 2)) - 1)
#define CPY_TAGGED_MIN (-((Py_ssize_t)1 << (CPY_INT_BITS - 2)))
#define CPY_TAGGED_ABS_MIN (0-(size_t)CPY_TAGGED_MIN)

typedef PyObject CPyModule;

// Tag bit used for long integers
#define CPY_INT_TAG 1

// Error value for signed fixed-width (low-level) integers
#define CPY_LL_INT_ERROR -113

// Error value for unsigned fixed-width (low-level) integers
#define CPY_LL_UINT_ERROR 239

// Error value for floats
#define CPY_FLOAT_ERROR -113.0

typedef void (*CPyVTableItem)(void);

static inline CPyTagged CPyTagged_ShortFromInt(int x) {
    return x << 1;
}

static inline CPyTagged CPyTagged_ShortFromSsize_t(Py_ssize_t x) {
    return x << 1;
}

// Are we targeting Python 3.12 or newer?
#define CPY_3_12_FEATURES (PY_VERSION_HEX >= 0x030c0000)

#if CPY_3_12_FEATURES

// Same as macros in CPython internal/pycore_long.h, but with a CPY_ prefix
#define CPY_NON_SIZE_BITS 3
#define CPY_SIGN_ZERO 1
#define CPY_SIGN_NEGATIVE 2
#define CPY_SIGN_MASK 3

#define CPY_LONG_DIGIT(o, n) ((o)->long_value.ob_digit[n])

// Only available on Python 3.12 and later
#define CPY_LONG_TAG(o) ((o)->long_value.lv_tag)
#define CPY_LONG_IS_NEGATIVE(o) (((o)->long_value.lv_tag & CPY_SIGN_MASK) == CPY_SIGN_NEGATIVE)
// Only available on Python 3.12 and later
#define CPY_LONG_SIZE(o) ((o)->long_value.lv_tag >> CPY_NON_SIZE_BITS)
// Number of digits; negative for negative ints
#define CPY_LONG_SIZE_SIGNED(o) (CPY_LONG_IS_NEGATIVE(o) ? -CPY_LONG_SIZE(o) : CPY_LONG_SIZE(o))
// Number of digits, assuming int is non-negative
#define CPY_LONG_SIZE_UNSIGNED(o) CPY_LONG_SIZE(o)

static inline void CPyLong_SetUnsignedSize(PyLongObject *o, Py_ssize_t n) {
    if (n == 0)
        o->long_value.lv_tag = CPY_SIGN_ZERO;
    else
        o->long_value.lv_tag = n << CPY_NON_SIZE_BITS;
}

#else

#define CPY_LONG_DIGIT(o, n) ((o)->ob_digit[n])
#define CPY_LONG_IS_NEGATIVE(o) (((o)->ob_base.ob_size < 0)
#define CPY_LONG_SIZE_SIGNED(o) ((o)->ob_base.ob_size)
#define CPY_LONG_SIZE_UNSIGNED(o) ((o)->ob_base.ob_size)

static inline void CPyLong_SetUnsignedSize(PyLongObject *o, Py_ssize_t n) {
    o->ob_base.ob_size = n;
}

#endif

#endif
