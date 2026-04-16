// Copyright (c) 2016-2025 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

// PLEASE DO NOT ADD ANY INCLUDES HERE

// Define some generic pybind11 helper macros for warning management.
//
// Note that compiler-specific push/pop pairs are baked into the
// PYBIND11_NAMESPACE_BEGIN/PYBIND11_NAMESPACE_END pair of macros. Therefore manual
// PYBIND11_WARNING_PUSH/PYBIND11_WARNING_POP are usually only needed in `#include` sections.
//
// If you find you need to suppress a warning, please try to make the suppression as local as
// possible using these macros. Please also be sure to push/pop with the pybind11 macros. Please
// only use compiler specifics if you need to check specific versions, e.g. Apple Clang vs. vanilla
// Clang.
#if defined(__INTEL_COMPILER)
#    define PYBIND11_COMPILER_INTEL
#    define PYBIND11_PRAGMA(...) _Pragma(#__VA_ARGS__)
#    define PYBIND11_WARNING_PUSH PYBIND11_PRAGMA(warning push)
#    define PYBIND11_WARNING_POP PYBIND11_PRAGMA(warning pop)
#elif defined(__clang__)
#    define PYBIND11_COMPILER_CLANG
#    define PYBIND11_PRAGMA(...) _Pragma(#__VA_ARGS__)
#    define PYBIND11_WARNING_PUSH PYBIND11_PRAGMA(clang diagnostic push)
#    define PYBIND11_WARNING_POP PYBIND11_PRAGMA(clang diagnostic pop)
#elif defined(__GNUC__)
#    define PYBIND11_COMPILER_GCC
#    define PYBIND11_PRAGMA(...) _Pragma(#__VA_ARGS__)
#    define PYBIND11_WARNING_PUSH PYBIND11_PRAGMA(GCC diagnostic push)
#    define PYBIND11_WARNING_POP PYBIND11_PRAGMA(GCC diagnostic pop)
#elif defined(_MSC_VER) // Must be after the clang branch because clang-cl also defines _MSC_VER
#    define PYBIND11_COMPILER_MSVC
#    define PYBIND11_PRAGMA(...) __pragma(__VA_ARGS__)
#    define PYBIND11_WARNING_PUSH PYBIND11_PRAGMA(warning(push))
#    define PYBIND11_WARNING_POP PYBIND11_PRAGMA(warning(pop))
#endif

#ifdef PYBIND11_COMPILER_MSVC
#    define PYBIND11_WARNING_DISABLE_MSVC(name) PYBIND11_PRAGMA(warning(disable : name))
#else
#    define PYBIND11_WARNING_DISABLE_MSVC(name)
#endif

#ifdef PYBIND11_COMPILER_CLANG
#    define PYBIND11_WARNING_DISABLE_CLANG(name) PYBIND11_PRAGMA(clang diagnostic ignored name)
#else
#    define PYBIND11_WARNING_DISABLE_CLANG(name)
#endif

#ifdef PYBIND11_COMPILER_GCC
#    define PYBIND11_WARNING_DISABLE_GCC(name) PYBIND11_PRAGMA(GCC diagnostic ignored name)
#else
#    define PYBIND11_WARNING_DISABLE_GCC(name)
#endif

#ifdef PYBIND11_COMPILER_INTEL
#    define PYBIND11_WARNING_DISABLE_INTEL(name) PYBIND11_PRAGMA(warning disable name)
#else
#    define PYBIND11_WARNING_DISABLE_INTEL(name)
#endif

#define PYBIND11_NAMESPACE_BEGIN(name)                                                            \
    namespace name {                                                                              \
    PYBIND11_WARNING_PUSH

#define PYBIND11_NAMESPACE_END(name)                                                              \
    PYBIND11_WARNING_POP                                                                          \
    }

// Robust support for some features and loading modules compiled against different pybind versions
// requires forcing hidden visibility on pybind code, so we enforce this by setting the attribute
// on the main `pybind11` namespace.
#if !defined(PYBIND11_NAMESPACE)
#    if defined(__GNUG__) && !defined(_WIN32)
#        define PYBIND11_NAMESPACE pybind11 __attribute__((visibility("hidden")))
#    else
#        define PYBIND11_NAMESPACE pybind11
#    endif
#endif
