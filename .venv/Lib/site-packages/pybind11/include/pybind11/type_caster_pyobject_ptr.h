// Copyright (c) 2023 The pybind Community.

#pragma once

#include "detail/common.h"
#include "detail/descr.h"
#include "cast.h"
#include "pytypes.h"

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

template <>
class type_caster<PyObject> {
public:
    static constexpr auto name = const_name("object"); // See discussion under PR #4601.

    // This overload is purely to guard against accidents.
    template <typename T,
              detail::enable_if_t<!is_same_ignoring_cvref<T, PyObject *>::value, int> = 0>
    static handle cast(T &&, return_value_policy, handle /*parent*/) {
        static_assert(is_same_ignoring_cvref<T, PyObject *>::value,
                      "Invalid C++ type T for to-Python conversion (type_caster<PyObject>).");
        return nullptr; // Unreachable.
    }

    static handle cast(PyObject *src, return_value_policy policy, handle /*parent*/) {
        if (src == nullptr) {
            throw error_already_set();
        }
        if (PyErr_Occurred()) {
            raise_from(PyExc_SystemError, "src != nullptr but PyErr_Occurred()");
            throw error_already_set();
        }
        if (policy == return_value_policy::take_ownership) {
            return src;
        }
        if (policy == return_value_policy::reference
            || policy == return_value_policy::automatic_reference) {
            return handle(src).inc_ref();
        }
        pybind11_fail("type_caster<PyObject>::cast(): unsupported return_value_policy: "
                      + std::to_string(static_cast<int>(policy)));
    }

    bool load(handle src, bool) {
        value = reinterpret_borrow<object>(src);
        return true;
    }

    template <typename T>
    using cast_op_type = PyObject *;

    explicit operator PyObject *() { return value.ptr(); }

private:
    object value;
};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
