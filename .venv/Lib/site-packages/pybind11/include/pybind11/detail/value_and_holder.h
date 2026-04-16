// Copyright (c) 2016-2024 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include "common.h"

#include <cstddef>
#include <cstdint>
#include <typeinfo>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

struct value_and_holder {
    instance *inst = nullptr;
    size_t index = 0u;
    const detail::type_info *type = nullptr;
    void **vh = nullptr;

    // Main constructor for a found value/holder:
    value_and_holder(instance *i, const detail::type_info *type, size_t vpos, size_t index)
        : inst{i}, index{index}, type{type},
          vh{inst->simple_layout ? inst->simple_value_holder
                                 : &inst->nonsimple.values_and_holders[vpos]} {}

    // Default constructor (used to signal a value-and-holder not found by get_value_and_holder())
    value_and_holder() = default;

    // Used for past-the-end iterator
    explicit value_and_holder(size_t index) : index{index} {}

    template <typename V = void>
    V *&value_ptr() const {
        return reinterpret_cast<V *&>(vh[0]);
    }
    // True if this `value_and_holder` has a non-null value pointer
    explicit operator bool() const { return value_ptr() != nullptr; }

    template <typename H>
    H &holder() const {
        return reinterpret_cast<H &>(vh[1]);
    }
    bool holder_constructed() const {
        return inst->simple_layout
                   ? inst->simple_holder_constructed
                   : (inst->nonsimple.status[index] & instance::status_holder_constructed) != 0u;
    }
    // NOLINTNEXTLINE(readability-make-member-function-const)
    void set_holder_constructed(bool v = true) {
        if (inst->simple_layout) {
            inst->simple_holder_constructed = v;
        } else if (v) {
            inst->nonsimple.status[index] |= instance::status_holder_constructed;
        } else {
            inst->nonsimple.status[index]
                &= static_cast<std::uint8_t>(~instance::status_holder_constructed);
        }
    }
    bool instance_registered() const {
        return inst->simple_layout
                   ? inst->simple_instance_registered
                   : ((inst->nonsimple.status[index] & instance::status_instance_registered) != 0);
    }
    // NOLINTNEXTLINE(readability-make-member-function-const)
    void set_instance_registered(bool v = true) {
        if (inst->simple_layout) {
            inst->simple_instance_registered = v;
        } else if (v) {
            inst->nonsimple.status[index] |= instance::status_instance_registered;
        } else {
            inst->nonsimple.status[index]
                &= static_cast<std::uint8_t>(~instance::status_instance_registered);
        }
    }
};

// This is a semi-public API to check if the corresponding instance has been constructed with a
// holder. That is, if the instance has been constructed with a holder, the `__init__` method is
// called and the C++ object is valid. Otherwise, the C++ object might only be allocated, but not
// initialized. This will lead to **SEGMENTATION FAULTS** if the C++ object is used in any way.
// Example usage: https://pybind11.readthedocs.io/en/stable/advanced/classes.html#custom-type-setup
//                for `tp_traverse` and `tp_clear` implementations.
// WARNING: The caller is responsible for ensuring that the `reinterpret_cast` is valid.
inline bool is_holder_constructed(PyObject *obj) {
    auto *const instance = reinterpret_cast<pybind11::detail::instance *>(obj);
    return instance->get_value_and_holder().holder_constructed();
}

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
