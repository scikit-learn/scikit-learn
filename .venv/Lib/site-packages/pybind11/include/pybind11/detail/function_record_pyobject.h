// Copyright (c) 2024-2025 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

// For background see the description of PR google/pybind11clif#30099.

#pragma once

#include <pybind11/attr.h>
#include <pybind11/conduit/pybind11_platform_abi_id.h>
#include <pybind11/pytypes.h>

#include "common.h"

#include <cstring>
#include <utility>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

struct function_record_PyObject {
    PyObject_HEAD
    function_record *cpp_func_rec;
};

PYBIND11_NAMESPACE_BEGIN(function_record_PyTypeObject_methods)

PyObject *tp_new_impl(PyTypeObject *type, PyObject *args, PyObject *kwds);
PyObject *tp_alloc_impl(PyTypeObject *type, Py_ssize_t nitems);
int tp_init_impl(PyObject *self, PyObject *args, PyObject *kwds);
void tp_dealloc_impl(PyObject *self);
void tp_free_impl(void *self);

static PyObject *reduce_ex_impl(PyObject *self, PyObject *, PyObject *);

static PyMethodDef tp_methods_impl[]
    = {{"__reduce_ex__",
        // reduce_ex_impl is a PyCFunctionWithKeywords, but PyMethodDef
        // requires a PyCFunction. The cast through void* is safe and
        // idiomatic with METH_KEYWORDS, and it successfully sidesteps
        // unhelpful compiler warnings.
        // NOLINTNEXTLINE(bugprone-casting-through-void)
        reinterpret_cast<PyCFunction>(reinterpret_cast<void *>(reduce_ex_impl)),
        METH_VARARGS | METH_KEYWORDS,
        nullptr},
       {nullptr, nullptr, 0, nullptr}};

// Python 3.12+ emits a DeprecationWarning for heap types whose tp_name does
// not contain a dot ('.') and that lack a __module__ attribute. For pybind11's
// internal function_record type, we do not have an actual module object to
// attach, so we cannot use PyType_FromModuleAndSpec (introduced in Python 3.9)
// to set __module__ automatically.
//
// As a workaround, we define a "qualified" type name that includes a dummy
// module name (PYBIND11_DUMMY_MODULE_NAME). This is non‑idiomatic but avoids
// the deprecation warning, and results in reprs like
//
//     <class 'pybind11_builtins.pybind11_detail_function_record_...'>
//
// even though no real pybind11_builtins module exists. If pybind11 gains an
// actual module object in the future, this code should switch to
// PyType_FromModuleAndSpec for Python 3.9+ and drop the dummy module
// workaround.
//
// Note that this name is versioned.
#define PYBIND11_DETAIL_FUNCTION_RECORD_TP_PLAINNAME                                              \
    "pybind11_detail_function_record_" PYBIND11_DETAIL_FUNCTION_RECORD_ABI_ID                     \
    "_" PYBIND11_PLATFORM_ABI_ID
constexpr char tp_plainname_impl[] = PYBIND11_DETAIL_FUNCTION_RECORD_TP_PLAINNAME;
constexpr char tp_qualname_impl[]
    = PYBIND11_DUMMY_MODULE_NAME "." PYBIND11_DETAIL_FUNCTION_RECORD_TP_PLAINNAME;

PYBIND11_NAMESPACE_END(function_record_PyTypeObject_methods)

static PyType_Slot function_record_PyType_Slots[] = {
    {Py_tp_dealloc,
     reinterpret_cast<void *>(function_record_PyTypeObject_methods::tp_dealloc_impl)},
    {Py_tp_methods,
     reinterpret_cast<void *>(function_record_PyTypeObject_methods::tp_methods_impl)},
    {Py_tp_init, reinterpret_cast<void *>(function_record_PyTypeObject_methods::tp_init_impl)},
    {Py_tp_alloc, reinterpret_cast<void *>(function_record_PyTypeObject_methods::tp_alloc_impl)},
    {Py_tp_new, reinterpret_cast<void *>(function_record_PyTypeObject_methods::tp_new_impl)},
    {Py_tp_free, reinterpret_cast<void *>(function_record_PyTypeObject_methods::tp_free_impl)},
    {0, nullptr}};

static PyType_Spec function_record_PyType_Spec
    = {function_record_PyTypeObject_methods::tp_qualname_impl,
       sizeof(function_record_PyObject),
       0,
       Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HEAPTYPE,
       function_record_PyType_Slots};

inline PyTypeObject *get_function_record_PyTypeObject() {
    PYBIND11_LOCK_INTERNALS(get_internals());
    PyTypeObject *&py_type_obj = detail::get_local_internals().function_record_py_type;
    if (!py_type_obj) {
        PyObject *py_obj = PyType_FromSpec(&function_record_PyType_Spec);
        if (py_obj == nullptr) {
            throw error_already_set();
        }
        py_type_obj = reinterpret_cast<PyTypeObject *>(py_obj);
    }
    return py_type_obj;
}

inline bool is_function_record_PyObject(PyObject *obj) {
    if (PyType_Check(obj) != 0) {
        return false;
    }
    PyTypeObject *obj_type = Py_TYPE(obj);

    PyTypeObject *frtype = get_function_record_PyTypeObject();

    // Fast path (pointer comparison).
    if (obj_type == frtype) {
        return true;
    }
    // This works across extension modules. Note that tp_name is versioned.
    if (strcmp(obj_type->tp_name, function_record_PyTypeObject_methods::tp_qualname_impl) == 0
        || strcmp(obj_type->tp_name, function_record_PyTypeObject_methods::tp_plainname_impl)
               == 0) {
        return true;
    }
    return false;
}

inline function_record *function_record_ptr_from_PyObject(PyObject *obj) {
    if (is_function_record_PyObject(obj)) {
        return (reinterpret_cast<detail::function_record_PyObject *>(obj))->cpp_func_rec;
    }
    return nullptr;
}

inline object function_record_PyObject_New() {
    auto *py_func_rec = PyObject_New(function_record_PyObject, get_function_record_PyTypeObject());
    if (py_func_rec == nullptr) {
        throw error_already_set();
    }
    py_func_rec->cpp_func_rec = nullptr; // For clarity/purity. Redundant in practice.
    return reinterpret_steal<object>(reinterpret_cast<PyObject *>(py_func_rec));
}

PYBIND11_NAMESPACE_BEGIN(function_record_PyTypeObject_methods)

// Guard against accidents & oversights, in particular when porting to future Python versions.
inline PyObject *tp_new_impl(PyTypeObject *, PyObject *, PyObject *) {
    pybind11_fail("UNEXPECTED CALL OF function_record_PyTypeObject_methods::tp_new_impl");
    // return nullptr; // Unreachable.
}

inline PyObject *tp_alloc_impl(PyTypeObject *, Py_ssize_t) {
    pybind11_fail("UNEXPECTED CALL OF function_record_PyTypeObject_methods::tp_alloc_impl");
    // return nullptr; // Unreachable.
}

inline int tp_init_impl(PyObject *, PyObject *, PyObject *) {
    pybind11_fail("UNEXPECTED CALL OF function_record_PyTypeObject_methods::tp_init_impl");
    // return -1; // Unreachable.
}

inline void tp_free_impl(void *) {
    pybind11_fail("UNEXPECTED CALL OF function_record_PyTypeObject_methods::tp_free_impl");
}

inline PyObject *reduce_ex_impl(PyObject *self, PyObject *, PyObject *) {
    // Deliberately ignoring the arguments for simplicity (expected is `protocol: int`).
    const function_record *rec = function_record_ptr_from_PyObject(self);
    if (rec == nullptr) {
        pybind11_fail(
            "FATAL: function_record_PyTypeObject reduce_ex_impl(): cannot obtain cpp_func_rec.");
    }
    if (rec->name != nullptr && rec->name[0] != '\0' && rec->scope
        && PyModule_Check(rec->scope.ptr()) != 0) {
        object scope_module = get_scope_module(rec->scope);
        if (scope_module) {
            auto builtins = reinterpret_borrow<dict>(PyEval_GetBuiltins());
            auto builtins_eval = builtins["eval"];
            auto reconstruct_args = make_tuple(str("__import__('importlib').import_module('")
                                               + scope_module + str("')"));
            return make_tuple(std::move(builtins_eval), std::move(reconstruct_args))
                .release()
                .ptr();
        }
    }
    set_error(PyExc_RuntimeError, repr(self) + str(" is not pickleable."));
    return nullptr;
}

PYBIND11_NAMESPACE_END(function_record_PyTypeObject_methods)

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
