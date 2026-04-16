// Copyright (c) 2020-2024 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

/* Proof-of-Concept for smart pointer interoperability.

High-level aspects:

* Support all `unique_ptr`, `shared_ptr` interops that are feasible.

* Cleanly and clearly report all interops that are infeasible.

* Meant to fit into a `PyObject`, as a holder for C++ objects.

* Support a system design that makes it impossible to trigger
  C++ Undefined Behavior, especially from Python.

* Support a system design with clean runtime inheritance casting. From this
  it follows that the `smart_holder` needs to be type-erased (`void*`).

* Handling of RTTI for the type-erased held pointer is NOT implemented here.
  It is the responsibility of the caller to ensure that `static_cast<T *>`
  is well-formed when calling `as_*` member functions. Inheritance casting
  needs to be handled in a different layer (similar to the code organization
  in boost/python/object/inheritance.hpp).

Details:

* The "root holder" chosen here is a `shared_ptr<void>` (named `vptr` in this
  implementation). This choice is practically inevitable because `shared_ptr`
  has only very limited support for inspecting and accessing its deleter.

* If created from a raw pointer, or a `unique_ptr` without a custom deleter,
  `vptr` always uses a custom deleter, to support `unique_ptr`-like disowning.
  The custom deleters could be extended to included life-time management for
  external objects (e.g. `PyObject`).

* If created from an external `shared_ptr`, or a `unique_ptr` with a custom
  deleter, including life-time management for external objects is infeasible.

* By choice, the smart_holder is movable but not copyable, to keep the design
  simple, and to guard against accidental copying overhead.

* The `void_cast_raw_ptr` option is needed to make the `smart_holder` `vptr`
  member invisible to the `shared_from_this` mechanism, in case the lifetime
  of a `PyObject` is tied to the pointee.
*/

#pragma once

#include "pybind11_namespace_macros.h"

#include <cstring>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <utility>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(memory)

// Default fallback.
static constexpr bool type_has_shared_from_this(...) { return false; }

// This overload uses SFINAE to skip enable_shared_from_this checks when the
// base is inaccessible (e.g. private inheritance).
template <typename T>
static auto type_has_shared_from_this(const T *ptr)
    -> decltype(static_cast<const std::enable_shared_from_this<T> *>(ptr), true) {
    return true;
}

// Inaccessible base → substitution failure → fallback overload selected
template <typename T>
static constexpr bool type_has_shared_from_this(const void *) {
    return false;
}

struct guarded_delete {
    // NOTE: PYBIND11_INTERNALS_VERSION needs to be bumped if changes are made to this struct.
    std::weak_ptr<void> released_ptr;    // Trick to keep the smart_holder memory footprint small.
    std::function<void(void *)> del_fun; // Rare case.
    void (*del_ptr)(void *);             // Common case.
    bool use_del_fun;
    bool armed_flag;
    guarded_delete(std::function<void(void *)> &&del_fun, bool armed_flag)
        : del_fun{std::move(del_fun)}, del_ptr{nullptr}, use_del_fun{true},
          armed_flag{armed_flag} {}
    guarded_delete(void (*del_ptr)(void *), bool armed_flag)
        : del_ptr{del_ptr}, use_del_fun{false}, armed_flag{armed_flag} {}
    void operator()(void *raw_ptr) const {
        if (armed_flag) {
            if (use_del_fun) {
                del_fun(raw_ptr);
            } else {
                del_ptr(raw_ptr);
            }
        }
    }
};

inline guarded_delete *get_guarded_delete(const std::shared_ptr<void> &ptr) {
    return std::get_deleter<guarded_delete>(ptr);
}

using get_guarded_delete_fn = guarded_delete *(*) (const std::shared_ptr<void> &);

template <typename T, typename std::enable_if<std::is_destructible<T>::value, int>::type = 0>
inline void std_default_delete_if_destructible(void *raw_ptr) {
    std::default_delete<T>{}(static_cast<T *>(raw_ptr));
}

template <typename T, typename std::enable_if<!std::is_destructible<T>::value, int>::type = 0>
inline void std_default_delete_if_destructible(void *) {
    // This noop operator is needed to avoid a compilation error (for `delete raw_ptr;`), but
    // throwing an exception from a destructor will std::terminate the process. Therefore the
    // runtime check for lifetime-management correctness is implemented elsewhere (in
    // ensure_pointee_is_destructible()).
}

template <typename T>
guarded_delete make_guarded_std_default_delete(bool armed_flag) {
    return guarded_delete(std_default_delete_if_destructible<T>, armed_flag);
}

template <typename T, typename D>
struct custom_deleter {
    // NOTE: PYBIND11_INTERNALS_VERSION needs to be bumped if changes are made to this struct.
    D deleter;
    explicit custom_deleter(D &&deleter) : deleter{std::forward<D>(deleter)} {}
    void operator()(void *raw_ptr) { deleter(static_cast<T *>(raw_ptr)); }
};

template <typename T, typename D>
guarded_delete make_guarded_custom_deleter(D &&uqp_del, bool armed_flag) {
    return guarded_delete(
        std::function<void(void *)>(custom_deleter<T, D>(std::forward<D>(uqp_del))), armed_flag);
}

template <typename T, typename D>
constexpr bool uqp_del_is_std_default_delete() {
    return std::is_same<D, std::default_delete<T>>::value
           || std::is_same<D, std::default_delete<T const>>::value;
}

inline bool type_info_equal_across_dso_boundaries(const std::type_info &a,
                                                  const std::type_info &b) {
    // RTTI pointer comparison may fail across DSOs (e.g., macOS libc++).
    // Fallback to name comparison, which is generally safe and ABI-stable enough for our use.
    return a == b || std::strcmp(a.name(), b.name()) == 0;
}

struct smart_holder {
    // NOTE: PYBIND11_INTERNALS_VERSION needs to be bumped if changes are made to this struct.
    const std::type_info *rtti_uqp_del = nullptr;
    std::shared_ptr<void> vptr;
    bool vptr_is_using_noop_deleter : 1;
    bool vptr_is_using_std_default_delete : 1;
    bool vptr_is_external_shared_ptr : 1;
    bool is_populated : 1;
    bool is_disowned : 1;

    // Design choice: smart_holder is movable but not copyable.
    smart_holder(smart_holder &&) = default;
    smart_holder(const smart_holder &) = delete;
    smart_holder &operator=(smart_holder &&) = delete;
    smart_holder &operator=(const smart_holder &) = delete;

    smart_holder()
        : vptr_is_using_noop_deleter{false}, vptr_is_using_std_default_delete{false},
          vptr_is_external_shared_ptr{false}, is_populated{false}, is_disowned{false} {}

    bool has_pointee() const { return vptr != nullptr; }

    template <typename T>
    static void ensure_pointee_is_destructible(const char *context) {
        if (!std::is_destructible<T>::value) {
            throw std::invalid_argument(std::string("Pointee is not destructible (") + context
                                        + ").");
        }
    }

    void ensure_is_populated(const char *context) const {
        if (!is_populated) {
            throw std::runtime_error(std::string("Unpopulated holder (") + context + ").");
        }
    }
    void ensure_is_not_disowned(const char *context) const {
        if (is_disowned) {
            throw std::runtime_error(std::string("Holder was disowned already (") + context
                                     + ").");
        }
    }

    void ensure_vptr_is_using_std_default_delete(const char *context) const {
        if (vptr_is_external_shared_ptr) {
            throw std::invalid_argument(std::string("Cannot disown external shared_ptr (")
                                        + context + ").");
        }
        if (vptr_is_using_noop_deleter) {
            throw std::invalid_argument(std::string("Cannot disown non-owning holder (") + context
                                        + ").");
        }
        if (!vptr_is_using_std_default_delete) {
            throw std::invalid_argument(std::string("Cannot disown custom deleter (") + context
                                        + ").");
        }
    }

    template <typename T, typename D>
    void ensure_compatible_uqp_del(const char *context) const {
        if (!rtti_uqp_del) {
            if (!uqp_del_is_std_default_delete<T, D>()) {
                throw std::invalid_argument(std::string("Missing unique_ptr deleter (") + context
                                            + ").");
            }
            ensure_vptr_is_using_std_default_delete(context);
            return;
        }
        if (uqp_del_is_std_default_delete<T, D>() && vptr_is_using_std_default_delete) {
            return;
        }
        if (!type_info_equal_across_dso_boundaries(typeid(D), *rtti_uqp_del)) {
            throw std::invalid_argument(std::string("Incompatible unique_ptr deleter (") + context
                                        + ").");
        }
    }

    void ensure_has_pointee(const char *context) const {
        if (!has_pointee()) {
            throw std::invalid_argument(std::string("Disowned holder (") + context + ").");
        }
    }

    void ensure_use_count_1(const char *context) const {
        if (vptr == nullptr) {
            throw std::invalid_argument(std::string("Cannot disown nullptr (") + context + ").");
        }
        // In multithreaded environments accessing use_count can lead to
        // race conditions, but in the context of Python it is a bug (elsewhere)
        // if the Global Interpreter Lock (GIL) is not being held when this code
        // is reached.
        // PYBIND11:REMINDER: This may need to be protected by a mutex in free-threaded Python.
        if (vptr.use_count() != 1) {
            throw std::invalid_argument(std::string("Cannot disown use_count != 1 (") + context
                                        + ").");
        }
    }

    void reset_vptr_deleter_armed_flag(const get_guarded_delete_fn ggd_fn, bool armed_flag) const {
        auto *gd = ggd_fn(vptr);
        if (gd == nullptr) {
            throw std::runtime_error(
                "smart_holder::reset_vptr_deleter_armed_flag() called in an invalid context.");
        }
        gd->armed_flag = armed_flag;
    }

    // Caller is responsible for precondition: ensure_compatible_uqp_del<T, D>() must succeed.
    template <typename T, typename D>
    std::unique_ptr<D> extract_deleter(const char *context,
                                       const get_guarded_delete_fn ggd_fn) const {
        auto *gd = ggd_fn(vptr);
        if (gd && gd->use_del_fun) {
            const auto &custom_deleter_ptr = gd->del_fun.template target<custom_deleter<T, D>>();
            if (custom_deleter_ptr == nullptr) {
                throw std::runtime_error(
                    std::string("smart_holder::extract_deleter() precondition failure (") + context
                    + ").");
            }
            static_assert(std::is_copy_constructible<D>::value,
                          "Required for compatibility with smart_holder functionality.");
            return std::unique_ptr<D>(new D(custom_deleter_ptr->deleter));
        }
        return nullptr;
    }

    static smart_holder from_raw_ptr_unowned(void *raw_ptr) {
        smart_holder hld;
        hld.vptr.reset(raw_ptr, [](void *) {});
        hld.vptr_is_using_noop_deleter = true;
        hld.is_populated = true;
        return hld;
    }

    template <typename T>
    T *as_raw_ptr_unowned() const {
        return static_cast<T *>(vptr.get());
    }

    template <typename T>
    static smart_holder from_raw_ptr_take_ownership(T *raw_ptr, bool void_cast_raw_ptr = false) {
        ensure_pointee_is_destructible<T>("from_raw_ptr_take_ownership");
        smart_holder hld;
        auto gd = make_guarded_std_default_delete<T>(true);
        if (void_cast_raw_ptr) {
            hld.vptr.reset(static_cast<void *>(raw_ptr), std::move(gd));
        } else {
            hld.vptr.reset(raw_ptr, std::move(gd));
        }
        hld.vptr_is_using_std_default_delete = true;
        hld.is_populated = true;
        return hld;
    }

    // Caller is responsible for ensuring the complex preconditions
    // (see `smart_holder_type_caster_support::load_helper`).
    void disown(const get_guarded_delete_fn ggd_fn) {
        reset_vptr_deleter_armed_flag(ggd_fn, false);
        is_disowned = true;
    }

    // Caller is responsible for ensuring the complex preconditions
    // (see `smart_holder_type_caster_support::load_helper`).
    void reclaim_disowned(const get_guarded_delete_fn ggd_fn) {
        reset_vptr_deleter_armed_flag(ggd_fn, true);
        is_disowned = false;
    }

    // Caller is responsible for ensuring the complex preconditions
    // (see `smart_holder_type_caster_support::load_helper`).
    void release_disowned() { vptr.reset(); }

    void ensure_can_release_ownership(const char *context = "ensure_can_release_ownership") const {
        ensure_is_not_disowned(context);
        ensure_vptr_is_using_std_default_delete(context);
        ensure_use_count_1(context);
    }

    // Caller is responsible for ensuring the complex preconditions
    // (see `smart_holder_type_caster_support::load_helper`).
    void release_ownership(const get_guarded_delete_fn ggd_fn) {
        reset_vptr_deleter_armed_flag(ggd_fn, false);
        release_disowned();
    }

    template <typename T, typename D>
    static smart_holder from_unique_ptr(std::unique_ptr<T, D> &&unq_ptr,
                                        void *mi_subobject_ptr = nullptr) {
        smart_holder hld;
        hld.rtti_uqp_del = &typeid(D);
        hld.vptr_is_using_std_default_delete = uqp_del_is_std_default_delete<T, D>();

        // Build the owning control block on the *real object start* (T*).
        guarded_delete gd
            = hld.vptr_is_using_std_default_delete
                  ? make_guarded_std_default_delete<T>(true)
                  : make_guarded_custom_deleter<T, D>(std::move(unq_ptr.get_deleter()), true);
        // Critical: construct owner with pointer we intend to delete
        std::shared_ptr<T> owner(unq_ptr.get(), std::move(gd));
        // Relinquish ownership only after successful construction of owner
        (void) unq_ptr.release();

        // Publish either the MI/VI subobject pointer (if provided) or the full object.
        // Why this is needed:
        //   * The `owner` shared_ptr must always manage the true object start (T*).
        //     That ensures the deleter is invoked on a valid object header, so the
        //     virtual destructor can dispatch safely (critical on MSVC with virtual
        //     inheritance, where base subobjects are not at offset 0).
        //   * However, pybind11 needs to *register* and expose the subobject pointer
        //     appropriate for the type being bound.
        //     This pointer may differ from the T* object start under multiple/virtual
        //     inheritance.
        // This is achieved by using an aliasing shared_ptr<void>:
        //   - `owner` retains lifetime of the actual T* object start for deletion.
        //   - `vptr` points at the adjusted subobject (mi_subobject_ptr), giving
        //     Python the correct identity/registration address.
        // If no subobject pointer is passed, we simply publish the full object.
        if (mi_subobject_ptr) {
            hld.vptr = std::shared_ptr<void>(owner, mi_subobject_ptr);
        } else {
            hld.vptr = std::static_pointer_cast<void>(owner);
        }

        hld.is_populated = true;
        return hld;
    }

    template <typename T>
    static smart_holder from_shared_ptr(const std::shared_ptr<T> &shd_ptr) {
        smart_holder hld;
        hld.vptr = std::static_pointer_cast<void>(shd_ptr);
        hld.vptr_is_external_shared_ptr = true;
        hld.is_populated = true;
        return hld;
    }

    template <typename T>
    std::shared_ptr<T> as_shared_ptr() const {
        return std::static_pointer_cast<T>(vptr);
    }
};

PYBIND11_NAMESPACE_END(memory)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
