/*
    pybind11/detail/holder_caster_foreign_helpers.h: Logic to implement
    set_foreign_holder() in copyable_ and movable_holder_caster.

    Copyright (c) 2025 Hudson River Trading LLC <opensource@hudson-trading.com>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <pybind11/gil.h>

#include "common.h"

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

struct holder_caster_foreign_helpers {
    struct py_deleter {
        void operator()(const void *) const noexcept {
            // Don't run the deleter if the interpreter has been shut down
            if (Py_IsInitialized() == 0) {
                return;
            }
            gil_scoped_acquire guard;
            Py_DECREF(o);
        }

        PyObject *o;
    };

    // Downcast shared_ptr from the enable_shared_from_this base to the target type.
    // SFINAE probe: use static_pointer_cast when the static downcast is valid (common case),
    // fall back to dynamic_pointer_cast when it isn't (virtual inheritance — issue #5989).
    // We can't use dynamic_pointer_cast unconditionally because it requires polymorphic types;
    // we can't use is_polymorphic to choose because that's orthogonal to virtual inheritance.
    // (The implementation uses the "tag dispatch via overload priority" trick.)
    template <typename type, typename esft_base>
    static auto esft_downcast(const std::shared_ptr<esft_base> &existing, int /*preferred*/)
        -> decltype(static_cast<type *>(std::declval<esft_base *>()), std::shared_ptr<type>()) {
        return std::static_pointer_cast<type>(existing);
    }

    template <typename type, typename esft_base>
    static std::shared_ptr<type> esft_downcast(const std::shared_ptr<esft_base> &existing,
                                               ... /*fallback*/) {
        return std::dynamic_pointer_cast<type>(existing);
    }

    template <typename type>
    static auto set_via_shared_from_this(type *value, std::shared_ptr<type> *holder_out)
        -> decltype(value->shared_from_this(), bool()) {
        // object derives from enable_shared_from_this;
        // try to reuse an existing shared_ptr if one is known
        if (auto existing = try_get_shared_from_this(value)) {
            *holder_out = esft_downcast<type>(existing, 0);
            return true;
        }
        return false;
    }

    template <typename type>
    static bool set_via_shared_from_this(void *, std::shared_ptr<type> *) {
        return false;
    }

    template <typename type>
    static bool set_foreign_holder(handle src, type *value, std::shared_ptr<type> *holder_out) {
        // We only support using std::shared_ptr<T> for foreign T, and
        // it's done by creating a new shared_ptr control block that
        // owns a reference to the original Python object.
        if (value == nullptr) {
            *holder_out = {};
            return true;
        }
        if (set_via_shared_from_this(value, holder_out)) {
            return true;
        }
        *holder_out = std::shared_ptr<type>(value, py_deleter{src.inc_ref().ptr()});
        return true;
    }

    template <typename type>
    static bool
    set_foreign_holder(handle src, const type *value, std::shared_ptr<const type> *holder_out) {
        std::shared_ptr<type> holder_mut;
        if (set_foreign_holder(src, const_cast<type *>(value), &holder_mut)) {
            *holder_out = holder_mut;
            return true;
        }
        return false;
    }

    template <typename type>
    static bool set_foreign_holder(handle, type *, ...) {
        throw cast_error("Unable to cast foreign type to held instance -- "
                         "only std::shared_ptr<T> is supported in this case");
    }
};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
