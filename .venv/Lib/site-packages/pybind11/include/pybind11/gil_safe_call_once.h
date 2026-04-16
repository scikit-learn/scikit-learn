// Copyright (c) 2023 The pybind Community.

#pragma once

#include "detail/common.h"
#include "detail/internals.h"
#include "gil.h"

#include <cassert>
#include <mutex>

#if defined(Py_GIL_DISABLED) || defined(PYBIND11_HAS_SUBINTERPRETER_SUPPORT)
#    include <atomic>
#endif
#ifdef PYBIND11_HAS_SUBINTERPRETER_SUPPORT
#    include <cstdint>
#    include <memory>
#    include <string>
#endif

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

PYBIND11_NAMESPACE_BEGIN(detail)
#if defined(Py_GIL_DISABLED) || defined(PYBIND11_HAS_SUBINTERPRETER_SUPPORT)
using atomic_bool = std::atomic_bool;
#else
using atomic_bool = bool;
#endif
PYBIND11_NAMESPACE_END(detail)

// Use the `gil_safe_call_once_and_store` class below instead of the naive
//
//   static auto imported_obj = py::module_::import("module_name"); // BAD, DO NOT USE!
//
// which has two serious issues:
//
//     1. Py_DECREF() calls potentially after the Python interpreter was finalized already, and
//     2. deadlocks in multi-threaded processes (because of missing lock ordering).
//
// The following alternative avoids both problems:
//
//   PYBIND11_CONSTINIT static py::gil_safe_call_once_and_store<py::object> storage;
//   auto &imported_obj = storage // Do NOT make this `static`!
//       .call_once_and_store_result([]() {
//           return py::module_::import("module_name");
//       })
//       .get_stored();
//
// The parameter of `call_once_and_store_result()` must be callable. It can make
// CPython API calls, and in particular, it can temporarily release the GIL.
//
// `T` can be any C++ type, it does not have to involve CPython API types.
//
// The behavior with regard to signals, e.g. `SIGINT` (`KeyboardInterrupt`),
// is not ideal. If the main thread is the one to actually run the `Callable`,
// then a `KeyboardInterrupt` will interrupt it if it is running normal Python
// code. The situation is different if a non-main thread runs the
// `Callable`, and then the main thread starts waiting for it to complete:
// a `KeyboardInterrupt` will not interrupt the non-main thread, but it will
// get processed only when it is the main thread's turn again and it is running
// normal Python code. However, this will be unnoticeable for quick call-once
// functions, which is usually the case.
//
// For in-depth background, see docs/advanced/deadlock.md
#ifndef PYBIND11_HAS_SUBINTERPRETER_SUPPORT
// Subinterpreter support is disabled.
// In this case, we can store the result globally, because there is only a single interpreter.
//
// The life span of the stored result is the entire process lifetime. It is leaked on process
// termination to avoid destructor calls after the Python interpreter was finalized.
template <typename T>
class gil_safe_call_once_and_store {
public:
    // PRECONDITION: The GIL must be held when `call_once_and_store_result()` is called.
    //
    // NOTE: The second parameter (finalize callback) is intentionally unused when subinterpreter
    // support is disabled. In that case, storage is process-global and intentionally leaked to
    // avoid calling destructors after the Python interpreter has been finalized.
    template <typename Callable>
    gil_safe_call_once_and_store &call_once_and_store_result(Callable &&fn,
                                                             void (*)(T &) /*unused*/ = nullptr) {
        if (!is_initialized_) { // This read is guarded by the GIL.
            // Multiple threads may enter here, because the GIL is released in the next line and
            // CPython API calls in the `fn()` call below may release and reacquire the GIL.
            gil_scoped_release gil_rel; // Needed to establish lock ordering.
            std::call_once(once_flag_, [&] {
                // Only one thread will ever enter here.
                gil_scoped_acquire gil_acq;
                ::new (storage_) T(fn()); // fn may release, but will reacquire, the GIL.
                is_initialized_ = true;   // This write is guarded by the GIL.
            });
            // All threads will observe `is_initialized_` as true here.
        }
        // Intentionally not returning `T &` to ensure the calling code is self-documenting.
        return *this;
    }

    // This must only be called after `call_once_and_store_result()` was called.
    T &get_stored() {
        assert(is_initialized_);
        PYBIND11_WARNING_PUSH
#    if !defined(__clang__) && defined(__GNUC__) && __GNUC__ < 5
        // Needed for gcc 4.8.5
        PYBIND11_WARNING_DISABLE_GCC("-Wstrict-aliasing")
#    endif
        return *reinterpret_cast<T *>(storage_);
        PYBIND11_WARNING_POP
    }

    constexpr gil_safe_call_once_and_store() = default;
    // The instance is a global static, so its destructor runs when the process
    // is terminating. Therefore, do nothing here because the Python interpreter
    // may have been finalized already.
    PYBIND11_DTOR_CONSTEXPR ~gil_safe_call_once_and_store() = default;

    // Disable copy and move operations.
    gil_safe_call_once_and_store(const gil_safe_call_once_and_store &) = delete;
    gil_safe_call_once_and_store(gil_safe_call_once_and_store &&) = delete;
    gil_safe_call_once_and_store &operator=(const gil_safe_call_once_and_store &) = delete;
    gil_safe_call_once_and_store &operator=(gil_safe_call_once_and_store &&) = delete;

private:
    // The global static storage (per-process) when subinterpreter support is disabled.
    alignas(T) char storage_[sizeof(T)] = {};
    std::once_flag once_flag_;

    // The `is_initialized_`-`storage_` pair is very similar to `std::optional`,
    // but the latter does not have the triviality properties of former,
    // therefore `std::optional` is not a viable alternative here.
    detail::atomic_bool is_initialized_{false};
};
#else
// Subinterpreter support is enabled.
// In this case, we should store the result per-interpreter instead of globally, because each
// subinterpreter has its own separate state. The cached result may not shareable across
// interpreters (e.g., imported modules and their members).

PYBIND11_NAMESPACE_BEGIN(detail)

template <typename T>
struct call_once_storage {
    alignas(T) char storage[sizeof(T)] = {};
    std::once_flag once_flag;
    void (*finalize)(T &) = nullptr;
    std::atomic_bool is_initialized{false};

    call_once_storage() = default;
    ~call_once_storage() {
        if (is_initialized) {
            if (finalize != nullptr) {
                finalize(*reinterpret_cast<T *>(storage));
            } else {
                reinterpret_cast<T *>(storage)->~T();
            }
        }
    }
    call_once_storage(const call_once_storage &) = delete;
    call_once_storage(call_once_storage &&) = delete;
    call_once_storage &operator=(const call_once_storage &) = delete;
    call_once_storage &operator=(call_once_storage &&) = delete;
};

PYBIND11_NAMESPACE_END(detail)

// Prefix for storage keys in the interpreter state dict.
#    define PYBIND11_CALL_ONCE_STORAGE_KEY_PREFIX PYBIND11_INTERNALS_ID "_call_once_storage__"

// The life span of the stored result is the entire interpreter lifetime. An additional
// `finalize_fn` can be provided to clean up the stored result when the interpreter is destroyed.
template <typename T>
class gil_safe_call_once_and_store {
public:
    // PRECONDITION: The GIL must be held when `call_once_and_store_result()` is called.
    template <typename Callable>
    gil_safe_call_once_and_store &call_once_and_store_result(Callable &&fn,
                                                             void (*finalize_fn)(T &) = nullptr) {
        if (!is_last_storage_valid()) {
            // Multiple threads may enter here, because the GIL is released in the next line and
            // CPython API calls in the `fn()` call below may release and reacquire the GIL.
            gil_scoped_release gil_rel; // Needed to establish lock ordering.
            // There can be multiple threads going through here.
            storage_type *value = nullptr;
            {
                gil_scoped_acquire gil_acq; // Restore lock ordering.
                // This function is thread-safe under free-threading.
                value = get_or_create_storage_in_state_dict();
            }
            assert(value != nullptr);
            std::call_once(value->once_flag, [&] {
                // Only one thread will ever enter here.
                gil_scoped_acquire gil_acq;
                // fn may release, but will reacquire, the GIL.
                ::new (value->storage) T(fn());
                value->finalize = finalize_fn;
                value->is_initialized = true;
                last_storage_ptr_ = reinterpret_cast<T *>(value->storage);
                is_initialized_by_at_least_one_interpreter_ = true;
            });
            // All threads will observe `is_initialized_by_at_least_one_interpreter_` as true here.
        }
        // Intentionally not returning `T &` to ensure the calling code is self-documenting.
        return *this;
    }

    // This must only be called after `call_once_and_store_result()` was called.
    T &get_stored() {
        T *result = last_storage_ptr_;
        if (!is_last_storage_valid()) {
            gil_scoped_acquire gil_acq;
            auto *value = get_or_create_storage_in_state_dict();
            result = last_storage_ptr_ = reinterpret_cast<T *>(value->storage);
        }
        assert(result != nullptr);
        return *result;
    }

    constexpr gil_safe_call_once_and_store() = default;
    // The instance is a global static, so its destructor runs when the process
    // is terminating. Therefore, do nothing here because the Python interpreter
    // may have been finalized already.
    PYBIND11_DTOR_CONSTEXPR ~gil_safe_call_once_and_store() = default;

    // Disable copy and move operations because the memory address is used as key.
    gil_safe_call_once_and_store(const gil_safe_call_once_and_store &) = delete;
    gil_safe_call_once_and_store(gil_safe_call_once_and_store &&) = delete;
    gil_safe_call_once_and_store &operator=(const gil_safe_call_once_and_store &) = delete;
    gil_safe_call_once_and_store &operator=(gil_safe_call_once_and_store &&) = delete;

private:
    using storage_type = detail::call_once_storage<T>;

    // Indicator of fast path for single-interpreter case.
    bool is_last_storage_valid() const {
        return is_initialized_by_at_least_one_interpreter_
               && !detail::has_seen_non_main_interpreter();
    }

    // Get the unique key for this storage instance in the interpreter's state dict.
    // The return type should not be `py::str` because PyObject is interpreter-dependent.
    std::string get_storage_key() const {
        // The instance is expected to be global static, so using its address as unique identifier.
        // The typical usage is like:
        //
        //   PYBIND11_CONSTINIT static gil_safe_call_once_and_store<T> storage;
        //
        return PYBIND11_CALL_ONCE_STORAGE_KEY_PREFIX
               + std::to_string(reinterpret_cast<std::uintptr_t>(this));
    }

    // Get or create per-storage capsule in the current interpreter's state dict.
    // The storage is interpreter-dependent and will not be shared across interpreters.
    storage_type *get_or_create_storage_in_state_dict() {
        return detail::atomic_get_or_create_in_state_dict<storage_type>(get_storage_key().c_str())
            .first;
    }

    // No storage needed when subinterpreter support is enabled.
    // The actual storage is stored in the per-interpreter state dict via
    // `get_or_create_storage_in_state_dict()`.

    // Fast local cache to avoid repeated lookups when there are no multiple interpreters.
    // This is only valid if there is a single interpreter. Otherwise, it is not used.
    // WARNING: We cannot use thread local cache similar to `internals_pp_manager::internals_p_tls`
    //          because the thread local storage cannot be explicitly invalidated when interpreters
    //          are destroyed (unlike `internals_pp_manager` which has explicit hooks for that).
    T *last_storage_ptr_ = nullptr;
    // This flag is true if the value has been initialized by any interpreter (may not be the
    // current one).
    detail::atomic_bool is_initialized_by_at_least_one_interpreter_{false};
};
#endif

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
