// Copyright (c) 2021 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include <pybind11/cast.h>
#include <pybind11/detail/common.h>
#include <pybind11/detail/descr.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>

#include <string>

#if defined(PYBIND11_HAS_FILESYSTEM)
#    include <filesystem>
#elif defined(PYBIND11_HAS_EXPERIMENTAL_FILESYSTEM)
#    include <experimental/filesystem>
#else
#    error "Neither #include <filesystem> nor #include <experimental/filesystem> is available."
#endif

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

#ifdef PYPY_VERSION
#    define PYBIND11_REINTERPRET_CAST_VOID_PTR_IF_NOT_PYPY(...) (__VA_ARGS__)
#else
#    define PYBIND11_REINTERPRET_CAST_VOID_PTR_IF_NOT_PYPY(...)                                   \
        (reinterpret_cast<void *>(__VA_ARGS__))
#endif

#if defined(PYBIND11_HAS_FILESYSTEM) || defined(PYBIND11_HAS_EXPERIMENTAL_FILESYSTEM)
template <typename T>
struct path_caster {

private:
    static PyObject *unicode_from_fs_native(const std::string &w) {
#    if !defined(PYPY_VERSION)
        return PyUnicode_DecodeFSDefaultAndSize(w.c_str(), ssize_t(w.size()));
#    else
        // PyPy mistakenly declares the first parameter as non-const.
        return PyUnicode_DecodeFSDefaultAndSize(const_cast<char *>(w.c_str()), ssize_t(w.size()));
#    endif
    }

    static PyObject *unicode_from_fs_native(const std::wstring &w) {
        return PyUnicode_FromWideChar(w.c_str(), ssize_t(w.size()));
    }

public:
    static handle cast(const T &path, return_value_policy, handle) {
        if (auto py_str = unicode_from_fs_native(path.native())) {
            return module_::import("pathlib")
                .attr("Path")(reinterpret_steal<object>(py_str))
                .release();
        }
        return nullptr;
    }

    bool load(handle handle, bool) {
        // PyUnicode_FSConverter and PyUnicode_FSDecoder normally take care of
        // calling PyOS_FSPath themselves, but that's broken on PyPy (PyPy
        // issue #3168) so we do it ourselves instead.
        PyObject *buf = PyOS_FSPath(handle.ptr());
        if (!buf) {
            PyErr_Clear();
            return false;
        }
        PyObject *native = nullptr;
        if constexpr (std::is_same_v<typename T::value_type, char>) {
            if (PyUnicode_FSConverter(buf, PYBIND11_REINTERPRET_CAST_VOID_PTR_IF_NOT_PYPY(&native))
                != 0) {
                if (auto *c_str = PyBytes_AsString(native)) {
                    // AsString returns a pointer to the internal buffer, which
                    // must not be free'd.
                    value = c_str;
                }
            }
        } else if constexpr (std::is_same_v<typename T::value_type, wchar_t>) {
            if (PyUnicode_FSDecoder(buf, PYBIND11_REINTERPRET_CAST_VOID_PTR_IF_NOT_PYPY(&native))
                != 0) {
                if (auto *c_str = PyUnicode_AsWideCharString(native, nullptr)) {
                    // AsWideCharString returns a new string that must be free'd.
                    value = c_str; // Copies the string.
                    PyMem_Free(c_str);
                }
            }
        }
        Py_XDECREF(native);
        Py_DECREF(buf);
        if (PyErr_Occurred()) {
            PyErr_Clear();
            return false;
        }
        return true;
    }

    PYBIND11_TYPE_CASTER(T, io_name("os.PathLike | str | bytes", "pathlib.Path"));
};

#endif // PYBIND11_HAS_FILESYSTEM || defined(PYBIND11_HAS_EXPERIMENTAL_FILESYSTEM)

#if defined(PYBIND11_HAS_FILESYSTEM)
template <>
struct type_caster<std::filesystem::path> : public path_caster<std::filesystem::path> {};
#elif defined(PYBIND11_HAS_EXPERIMENTAL_FILESYSTEM)
template <>
struct type_caster<std::experimental::filesystem::path>
    : public path_caster<std::experimental::filesystem::path> {};
#endif

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
