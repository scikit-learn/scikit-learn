#pragma once

// Copyright (c) 2024 The pybind Community.

// To maximize reusability:
// DO NOT ADD CODE THAT REQUIRES C++ EXCEPTION HANDLING.

#include "wrap_include_python_h.h"

// Implementation details. DO NOT USE ELSEWHERE. (Unfortunately we cannot #undef them.)
// This is duplicated here to maximize portability.
#define PYBIND11_PLATFORM_ABI_ID_STRINGIFY(x) #x
#define PYBIND11_PLATFORM_ABI_ID_TOSTRING(x) PYBIND11_PLATFORM_ABI_ID_STRINGIFY(x)

#ifdef PYBIND11_COMPILER_TYPE
//   // To maintain backward compatibility (see PR #5439).
#    define PYBIND11_COMPILER_TYPE_LEADING_UNDERSCORE ""
#else
#    define PYBIND11_COMPILER_TYPE_LEADING_UNDERSCORE "_"
#    if defined(__MINGW32__)
#        define PYBIND11_COMPILER_TYPE "mingw"
#    elif defined(__CYGWIN__)
#        define PYBIND11_COMPILER_TYPE "gcc_cygwin"
#    elif defined(_MSC_VER)
#        define PYBIND11_COMPILER_TYPE "msvc"
#    elif defined(__clang__) || defined(__GNUC__)
#        define PYBIND11_COMPILER_TYPE "system" // Assumed compatible with system compiler.
#    else
#        error "Unknown PYBIND11_COMPILER_TYPE: PLEASE REVISE THIS CODE."
#    endif
#endif

// PR #5439 made this macro obsolete. However, there are many manipulations of this macro in the
// wild. Therefore, to maintain backward compatibility, it is kept around.
#ifndef PYBIND11_STDLIB
#    define PYBIND11_STDLIB ""
#endif

#ifndef PYBIND11_BUILD_ABI
#    if defined(_MSC_VER)                 // See PR #4953.
#        if defined(_MT) && defined(_DLL) // Corresponding to CL command line options /MD or /MDd.
#            if (_MSC_VER) / 100 == 19
#                define PYBIND11_BUILD_ABI "_md_mscver19"
#            else
#                error "Unknown major version for MSC_VER: PLEASE REVISE THIS CODE."
#            endif
#        elif defined(_MT) // Corresponding to CL command line options /MT or /MTd.
#            define PYBIND11_BUILD_ABI "_mt_mscver" PYBIND11_PLATFORM_ABI_ID_TOSTRING(_MSC_VER)
#        else
#            if (_MSC_VER) / 100 == 19
#                define PYBIND11_BUILD_ABI "_none_mscver19"
#            else
#                error "Unknown major version for MSC_VER: PLEASE REVISE THIS CODE."
#            endif
#        endif
#    elif defined(_LIBCPP_ABI_VERSION) // https://libcxx.llvm.org/DesignDocs/ABIVersioning.html
#        define PYBIND11_BUILD_ABI                                                                \
            "_libcpp_abi" PYBIND11_PLATFORM_ABI_ID_TOSTRING(_LIBCPP_ABI_VERSION)
#    elif defined(_GLIBCXX_USE_CXX11_ABI) // See PR #5439.
#        if defined(__NVCOMPILER)
//           // Assume that NVHPC is in the 1xxx ABI family.
//           // THIS ASSUMPTION IS NOT FUTURE PROOF but apparently the best we can do.
//           // Please let us know if there is a way to validate the assumption here.
#        elif !defined(__GXX_ABI_VERSION)
#            error                                                                                \
                "Unknown platform or compiler (_GLIBCXX_USE_CXX11_ABI): PLEASE REVISE THIS CODE."
#        endif
#        if defined(__GXX_ABI_VERSION) && __GXX_ABI_VERSION < 1002 || __GXX_ABI_VERSION >= 2000
#            error "Unknown platform or compiler (__GXX_ABI_VERSION): PLEASE REVISE THIS CODE."
#        endif
#        define PYBIND11_BUILD_ABI                                                                \
            "_libstdcpp_gxx_abi_1xxx_use_cxx11_abi_" PYBIND11_PLATFORM_ABI_ID_TOSTRING(           \
                _GLIBCXX_USE_CXX11_ABI)
#    else
#        error "Unknown platform or compiler: PLEASE REVISE THIS CODE."
#    endif
#endif

// On MSVC, debug and release builds are not ABI-compatible!
#if defined(_MSC_VER) && defined(_DEBUG)
#    define PYBIND11_BUILD_TYPE "_debug"
#else
#    define PYBIND11_BUILD_TYPE ""
#endif

#define PYBIND11_PLATFORM_ABI_ID                                                                  \
    PYBIND11_COMPILER_TYPE PYBIND11_STDLIB PYBIND11_BUILD_ABI PYBIND11_BUILD_TYPE
