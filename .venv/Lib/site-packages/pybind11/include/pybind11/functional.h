/*
    pybind11/functional.h: std::function<> support

    Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "pybind11.h"

#include <functional>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)
PYBIND11_NAMESPACE_BEGIN(type_caster_std_function_specializations)

// ensure GIL is held during functor destruction
struct func_handle {
    function f;
#if !(defined(_MSC_VER) && _MSC_VER == 1916 && defined(PYBIND11_CPP17))
    // This triggers a syntax error under very special conditions (very weird indeed).
    explicit
#endif
        func_handle(function &&f_) noexcept
        : f(std::move(f_)) {
    }
    func_handle(const func_handle &f_) { operator=(f_); }
    func_handle &operator=(const func_handle &f_) {
        gil_scoped_acquire acq;
        f = f_.f;
        return *this;
    }
    ~func_handle() {
        gil_scoped_acquire acq;
        function kill_f(std::move(f));
    }
};

// to emulate 'move initialization capture' in C++11
struct func_wrapper_base {
    func_handle hfunc;
    explicit func_wrapper_base(func_handle &&hf) noexcept : hfunc(hf) {}
};

template <typename Return, typename... Args>
struct func_wrapper : func_wrapper_base {
    using func_wrapper_base::func_wrapper_base;
    Return operator()(Args... args) const { // NOLINT(performance-unnecessary-value-param)
        gil_scoped_acquire acq;
        // casts the returned object as a rvalue to the return type
        return hfunc.f(std::forward<Args>(args)...).template cast<Return>();
    }
};

PYBIND11_NAMESPACE_END(type_caster_std_function_specializations)

template <typename Return, typename... Args>
struct type_caster<std::function<Return(Args...)>> {
    using type = std::function<Return(Args...)>;
    using retval_type = conditional_t<std::is_same<Return, void>::value, void_type, Return>;
    using function_type = Return (*)(Args...);

public:
    bool load(handle src, bool convert) {
        if (src.is_none()) {
            // Defer accepting None to other overloads (if we aren't in convert mode):
            if (!convert) {
                return false;
            }
            return true;
        }

        if (!isinstance<function>(src)) {
            return false;
        }

        auto func = reinterpret_borrow<function>(src);

        /*
           When passing a C++ function as an argument to another C++
           function via Python, every function call would normally involve
           a full C++ -> Python -> C++ roundtrip, which can be prohibitive.
           Here, we try to at least detect the case where the function is
           stateless (i.e. function pointer or lambda function without
           captured variables), in which case the roundtrip can be avoided.
         */
        if (auto cfunc = func.cpp_function()) {
            auto *cfunc_self = PyCFunction_GET_SELF(cfunc.ptr());
            if (cfunc_self == nullptr) {
                PyErr_Clear();
            } else {
                function_record *rec = function_record_ptr_from_PyObject(cfunc_self);
                while (rec != nullptr) {
                    if (rec->is_stateless
                        && same_type(typeid(function_type),
                                     *reinterpret_cast<const std::type_info *>(rec->data[1]))) {
                        struct capture {
                            function_type f;

                            static capture *from_data(void **data) {
                                return PYBIND11_STD_LAUNDER(reinterpret_cast<capture *>(data));
                            }
                        };
                        PYBIND11_ENSURE_PRECONDITION_FOR_FUNCTIONAL_H_PERFORMANCE_OPTIMIZATIONS(
                            std::is_standard_layout<capture>::value);
                        value = capture::from_data(rec->data)->f;
                        return true;
                    }
                    rec = rec->next;
                }
            }
            // PYPY segfaults here when passing builtin function like sum.
            // Raising an fail exception here works to prevent the segfault, but only on gcc.
            // See PR #1413 for full details
        }

        value = type_caster_std_function_specializations::func_wrapper<Return, Args...>(
            type_caster_std_function_specializations::func_handle(std::move(func)));
        return true;
    }

    template <typename Func>
    static handle cast(Func &&f_, return_value_policy policy, handle /* parent */) {
        if (!f_) {
            return none().release();
        }

        auto result = f_.template target<function_type>();
        if (result) {
            return cpp_function(*result, policy).release();
        }
        return cpp_function(std::forward<Func>(f_), policy).release();
    }

    PYBIND11_TYPE_CASTER(
        type,
        const_name("collections.abc.Callable[[")
            + ::pybind11::detail::concat(::pybind11::detail::arg_descr(make_caster<Args>::name)...)
            + const_name("], ") + ::pybind11::detail::return_descr(make_caster<retval_type>::name)
            + const_name("]"));
};

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
