/*
    pybind11/warnings.h: Python warnings wrappers.

    Copyright (c) 2024 Jan Iwaszkiewicz <jiwaszkiewicz6@gmail.com>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "pybind11.h"
#include "detail/common.h"

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

PYBIND11_NAMESPACE_BEGIN(detail)

inline bool PyWarning_Check(PyObject *obj) {
    int result = PyObject_IsSubclass(obj, PyExc_Warning);
    if (result == 1) {
        return true;
    }
    if (result == -1) {
        raise_from(PyExc_SystemError,
                   "pybind11::detail::PyWarning_Check(): PyObject_IsSubclass() call failed.");
        throw error_already_set();
    }
    return false;
}

PYBIND11_NAMESPACE_END(detail)

PYBIND11_NAMESPACE_BEGIN(warnings)

inline object
new_warning_type(handle scope, const char *name, handle base = PyExc_RuntimeWarning) {
    if (!detail::PyWarning_Check(base.ptr())) {
        pybind11_fail("pybind11::warnings::new_warning_type(): cannot create custom warning, base "
                      "must be a subclass of "
                      "PyExc_Warning!");
    }
    if (hasattr(scope, name)) {
        pybind11_fail("pybind11::warnings::new_warning_type(): an attribute with name \""
                      + std::string(name) + "\" exists already.");
    }
    std::string full_name = scope.attr("__name__").cast<std::string>() + std::string(".") + name;
    handle h(PyErr_NewException(full_name.c_str(), base.ptr(), nullptr));
    if (!h) {
        raise_from(PyExc_SystemError,
                   "pybind11::warnings::new_warning_type(): PyErr_NewException() call failed.");
        throw error_already_set();
    }
    auto obj = reinterpret_steal<object>(h);
    scope.attr(name) = obj;
    return obj;
}

// Similar to Python `warnings.warn()`
inline void
warn(const char *message, handle category = PyExc_RuntimeWarning, int stack_level = 2) {
    if (!detail::PyWarning_Check(category.ptr())) {
        pybind11_fail(
            "pybind11::warnings::warn(): cannot raise warning, category must be a subclass of "
            "PyExc_Warning!");
    }

    if (PyErr_WarnEx(category.ptr(), message, stack_level) == -1) {
        throw error_already_set();
    }
}

PYBIND11_NAMESPACE_END(warnings)

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
