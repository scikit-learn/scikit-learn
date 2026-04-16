// Copyright (c) 2022-2025 The pybind Community.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include "detail/common.h"
#include "detail/native_enum_data.h"
#include "detail/type_caster_base.h"
#include "cast.h"

#include <cassert>
#include <limits>
#include <type_traits>
#include <typeindex>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

/// Conversions between Python's native (stdlib) enum types and C++ enums.
template <typename EnumType>
class native_enum : public detail::native_enum_data {
public:
    using Underlying = typename std::underlying_type<EnumType>::type;

    native_enum(handle parent_scope,
                const char *name,
                const char *native_type_name,
                const char *class_doc = "")
        : detail::native_enum_data(
              parent_scope, name, native_type_name, class_doc, make_record()) {
        if (detail::get_local_type_info(typeid(EnumType)) != nullptr
            || detail::get_global_type_info(typeid(EnumType)) != nullptr) {
            pybind11_fail(
                "pybind11::native_enum<...>(\"" + enum_name_encoded
                + "\") is already registered as a `pybind11::enum_` or `pybind11::class_`!");
        }
        if (detail::global_internals_native_enum_type_map_contains(enum_type_index)) {
            pybind11_fail("pybind11::native_enum<...>(\"" + enum_name_encoded
                          + "\") is already registered!");
        }
        arm_finalize_check();
    }

    /// Export enumeration entries into the parent scope
    native_enum &export_values() {
        assert(!export_values_flag); // Catch redundant calls.
        export_values_flag = true;
        return *this;
    }

    /// Add an enumeration entry
    native_enum &value(char const *name, EnumType value, const char *doc = nullptr) {
        // Disarm for the case that the native_enum_data dtor runs during exception unwinding.
        disarm_finalize_check("value after finalize");
        members.append(make_tuple(name, static_cast<Underlying>(value)));
        if (doc) {
            member_docs.append(make_tuple(name, doc));
        }
        arm_finalize_check(); // There was no exception.
        return *this;
    }

    native_enum(const native_enum &) = delete;
    native_enum &operator=(const native_enum &) = delete;

private:
    static detail::native_enum_record make_record() {
        detail::native_enum_record ret;
        ret.cpptype = &typeid(EnumType);
        ret.size_bytes = sizeof(EnumType);
        ret.is_signed = std::is_signed<Underlying>::value;
        return ret;
    }
};

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
