// Copyright (c) 2022-2025 The pybind Community.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include "../pytypes.h"
#include "common.h"
#include "internals.h"

#include <cassert>
#include <sstream>
#include <string>
#include <typeindex>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

// This is a separate function only to enable easy unit testing.
inline std::string
native_enum_missing_finalize_error_message(const std::string &enum_name_encoded) {
    return "pybind11::native_enum<...>(\"" + enum_name_encoded + "\", ...): MISSING .finalize()";
}

// Internals for pybind11::native_enum; one native_enum_data object exists
// inside each pybind11::native_enum and lives only for the duration of the
// native_enum binding statement.
class native_enum_data {
public:
    native_enum_data(handle parent_scope_,
                     const char *enum_name,
                     const char *native_type_name,
                     const char *class_doc,
                     const native_enum_record &enum_record_)
        : enum_name_encoded{enum_name}, native_type_name_encoded{native_type_name},
          enum_type_index{*enum_record_.cpptype},
          parent_scope(reinterpret_borrow<object>(parent_scope_)), enum_name{enum_name},
          native_type_name{native_type_name}, class_doc(class_doc), export_values_flag{false},
          finalize_needed{false} {
        // Create the enum record capsule. It will be installed on the enum
        // type object during finalize(). Its destructor removes the enum
        // mapping from our internals, so that we won't try to convert to an
        // enum type that's been destroyed.
        enum_record = capsule(
            new native_enum_record{enum_record_},
            native_enum_record::attribute_name(),
            +[](void *record_) {
                auto *record = static_cast<native_enum_record *>(record_);
                with_internals([&](internals &internals) {
                    internals.native_enum_type_map.erase(*record->cpptype);
                });
                delete record;
            });
    }

    void finalize();

    native_enum_data(const native_enum_data &) = delete;
    native_enum_data &operator=(const native_enum_data &) = delete;

#if !defined(NDEBUG)
    // This dtor cannot easily be unit tested because it terminates the process.
    ~native_enum_data() {
        if (finalize_needed) {
            pybind11_fail(native_enum_missing_finalize_error_message(enum_name_encoded));
        }
    }
#endif

protected:
    void disarm_finalize_check(const char *error_context) {
        if (!finalize_needed) {
            pybind11_fail("pybind11::native_enum<...>(\"" + enum_name_encoded
                          + "\"): " + error_context);
        }
        finalize_needed = false;
    }

    void arm_finalize_check() {
        assert(!finalize_needed); // Catch redundant calls.
        finalize_needed = true;
    }

    std::string enum_name_encoded;
    std::string native_type_name_encoded;
    std::type_index enum_type_index;

private:
    object parent_scope;
    str enum_name;
    str native_type_name;
    std::string class_doc;
    capsule enum_record;

protected:
    list members;
    list member_docs;
    bool export_values_flag : 1; // Attention: It is best to keep the bools together.

private:
    bool finalize_needed : 1;
};

inline handle
global_internals_native_enum_type_map_get_item(const std::type_index &enum_type_index) {
    return with_internals([&](internals &internals) {
        auto found = internals.native_enum_type_map.find(enum_type_index);
        if (found != internals.native_enum_type_map.end()) {
            return handle(found->second);
        }
        return handle();
    });
}

inline bool
global_internals_native_enum_type_map_contains(const std::type_index &enum_type_index) {
    return with_internals([&](internals &internals) {
        return internals.native_enum_type_map.count(enum_type_index) != 0;
    });
}

inline object import_or_getattr(const std::string &fully_qualified_name,
                                const std::string &append_to_exception_message) {
    std::istringstream stream(fully_qualified_name);
    std::string part;

    if (!std::getline(stream, part, '.') || part.empty()) {
        std::string msg = "Invalid fully-qualified name `";
        msg += fully_qualified_name;
        msg += "`";
        msg += append_to_exception_message;
        throw value_error(msg);
    }

    auto curr_scope = reinterpret_steal<object>(PyImport_ImportModule(part.c_str()));
    if (!curr_scope) {
        std::string msg = "Failed to import top-level module `";
        msg += part;
        msg += "`";
        msg += append_to_exception_message;
        raise_from(PyExc_ImportError, msg.c_str());
        throw error_already_set();
    }

    // Now recursively getattr or import remaining parts
    std::string curr_path = part;
    while (std::getline(stream, part, '.')) {
        if (part.empty()) {
            std::string msg = "Invalid fully-qualified name `";
            msg += fully_qualified_name;
            msg += "`";
            msg += append_to_exception_message;
            throw value_error(msg);
        }
        std::string next_path = curr_path;
        next_path += ".";
        next_path += part;
        auto next_scope
            = reinterpret_steal<object>(PyObject_GetAttrString(curr_scope.ptr(), part.c_str()));
        if (!next_scope) {
            error_fetch_and_normalize stored_getattr_error("getattr");
            // Try importing the next level
            next_scope = reinterpret_steal<object>(PyImport_ImportModule(next_path.c_str()));
            if (!next_scope) {
                error_fetch_and_normalize stored_import_error("import");
                std::string msg = "Failed to import or getattr `";
                msg += part;
                msg += "` from `";
                msg += curr_path;
                msg += "`";
                msg += append_to_exception_message;
                msg += "\n-------- getattr exception --------\n";
                msg += stored_getattr_error.error_string();
                msg += "\n-------- import exception --------\n";
                msg += stored_import_error.error_string();
                throw import_error(msg.c_str());
            }
        }
        curr_scope = next_scope;
        curr_path = next_path;
    }
    return curr_scope;
}

inline void native_enum_data::finalize() {
    disarm_finalize_check("DOUBLE finalize");
    if (hasattr(parent_scope, enum_name)) {
        pybind11_fail("pybind11::native_enum<...>(\"" + enum_name_encoded
                      + "\"): an object with that name is already defined");
    }
    auto py_enum_type = import_or_getattr(native_type_name, " (native_type_name)");
    auto py_enum = py_enum_type(enum_name, members);
    object module_name = get_module_name_if_available(parent_scope);
    if (module_name) {
        py_enum.attr("__module__") = module_name;
    }
    if (hasattr(parent_scope, "__qualname__")) {
        const auto parent_qualname = parent_scope.attr("__qualname__").cast<std::string>();
        py_enum.attr("__qualname__") = str(parent_qualname + "." + enum_name.cast<std::string>());
    }
    parent_scope.attr(enum_name) = py_enum;
    if (export_values_flag) {
        for (auto member : members) {
            auto member_name = member[int_(0)];
            if (hasattr(parent_scope, member_name)) {
                pybind11_fail("pybind11::native_enum<...>(\"" + enum_name_encoded + "\").value(\""
                              + member_name.cast<std::string>()
                              + "\"): an object with that name is already defined");
            }
            parent_scope.attr(member_name) = py_enum[member_name];
        }
    }
    if (!class_doc.empty()) {
        py_enum.attr("__doc__") = class_doc.c_str();
    }
    for (auto doc : member_docs) {
        py_enum[doc[int_(0)]].attr("__doc__") = doc[int_(1)];
    }

    py_enum.attr(native_enum_record::attribute_name()) = enum_record;
    with_internals([&](internals &internals) {
        internals.native_enum_type_map[enum_type_index] = py_enum.ptr();
    });
}

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
