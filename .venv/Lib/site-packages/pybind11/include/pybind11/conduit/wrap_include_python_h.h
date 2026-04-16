#pragma once

// Copyright (c) 2024 The pybind Community.

// STRONG REQUIREMENT:
//   This header is a wrapper around `#include <Python.h>`, therefore it
//   MUST BE INCLUDED BEFORE ANY STANDARD HEADERS are included.
// See also:
//   https://docs.python.org/3/c-api/intro.html#include-files
// Quoting from there:
//   Note: Since Python may define some pre-processor definitions which affect
//   the standard headers on some systems, you must include Python.h before
//   any standard headers are included.

// To maximize reusability:
// DO NOT ADD CODE THAT REQUIRES C++ EXCEPTION HANDLING.

// Disable linking to pythonX_d.lib on Windows in debug mode.
#if defined(_MSC_VER) && defined(_DEBUG) && !defined(Py_DEBUG)
// Workaround for a VS 2022 issue.
// See https://github.com/pybind/pybind11/pull/3497 for full context.
// NOTE: This workaround knowingly violates the Python.h include order
//       requirement (see above).
#    include <yvals.h>
#    if _MSVC_STL_VERSION >= 143
#        include <crtdefs.h>
#    endif
#    define PYBIND11_DEBUG_MARKER
#    undef _DEBUG
#endif

// Don't let Python.h #define (v)snprintf as macro because they are implemented
// properly in Visual Studio since 2015.
#if defined(_MSC_VER)
#    define HAVE_SNPRINTF 1
#endif

#if defined(_MSC_VER)
#    pragma warning(push)
#    pragma warning(disable : 4505)
// C4505: 'PySlice_GetIndicesEx': unreferenced local function has been removed
#endif

#include <Python.h>
#include <frameobject.h>
#include <pythread.h>

#if defined(_MSC_VER)
#    pragma warning(pop)
#endif

#if defined(PYBIND11_DEBUG_MARKER)
#    define _DEBUG 1
#    undef PYBIND11_DEBUG_MARKER
#endif

// Python #defines overrides on all sorts of core functions, which
// tends to wreak havok in C++ codebases that expect these to work
// like regular functions (potentially with several overloads).
#if defined(isalnum)
#    undef isalnum
#    undef isalpha
#    undef islower
#    undef isspace
#    undef isupper
#    undef tolower
#    undef toupper
#endif

#if defined(copysign)
#    undef copysign
#endif
