// Copyright (c) 2024 The pybind Community.

/* The pybind11_conduit_v1 feature enables type-safe interoperability between

* different independent Python/C++ bindings systems,

* including pybind11 versions with different PYBIND11_INTERNALS_VERSION's.

    * NOTE: The conduit feature
            only covers    from-Python-to-C++ conversions, it
            does not cover from-C++-to-Python conversions.
            (For the latter, a different feature would have to be added.)

The naming of the feature is a bit misleading:

* The feature is in no way tied to pybind11 internals.

* It just happens to originate from pybind11 and currently still lives there.

* The only external dependency is <Python.h>.

The implementation is a VERY light-weight dependency. It is designed to be
compatible with any ISO C++11 (or higher) compiler, and does NOT require
C++ Exception Handling to be enabled.

Please see https://github.com/pybind/pybind11/pull/5296 for more background.

The implementation involves a

def _pybind11_conduit_v1_(
    self,
    pybind11_platform_abi_id: bytes,
    cpp_type_info_capsule: capsule,
    pointer_kind: bytes) -> capsule

method that is meant to be added to Python objects wrapping C++ objects
(e.g. pybind11::class_-wrapped types).

The design of the _pybind11_conduit_v1_ feature provides two layers of
protection against C++ ABI mismatches:

* The first and most important layer is that the pybind11_platform_abi_id's
  must match between extensions. — This will never be perfect, but is the same
  pragmatic approach used in pybind11 since 2017
  (https://github.com/pybind/pybind11/commit/96997a4b9d4ec3d389a570604394af5d5eee2557,
  PYBIND11_INTERNALS_ID).

* The second layer is that the typeid(std::type_info).name()'s must match
  between extensions.

The implementation below (which is shorter than this comment!), serves as a
battle-tested specification. The main API is this one function:

auto *cpp_pointer = pybind11_conduit_v1::get_type_pointer_ephemeral<YourType>(py_obj);

It is meant to be a minimalistic reference implementation, intentionally
without comprehensive error reporting. It is expected that major bindings
systems will roll their own, compatible implementations, potentially with
system-specific error reporting. The essential specifications all bindings
systems need to agree on are merely:

* PYBIND11_PLATFORM_ABI_ID (const char* literal).

* The cpp_type_info capsule (see below: a void *ptr and a const char *name).

* The cpp_conduit capsule (see below: a void *ptr and a const char *name).

* "raw_pointer_ephemeral" means: the lifetime of the pointer is the lifetime
  of the py_obj.

*/

// THIS MUST STAY AT THE TOP!
#include "pybind11_platform_abi_id.h"

#include <Python.h>
#include <typeinfo>

namespace pybind11_conduit_v1 {

inline void *get_raw_pointer_ephemeral(PyObject *py_obj, const std::type_info *cpp_type_info) {
    PyObject *cpp_type_info_capsule
        = PyCapsule_New(const_cast<void *>(static_cast<const void *>(cpp_type_info)),
                        typeid(std::type_info).name(),
                        nullptr);
    if (cpp_type_info_capsule == nullptr) {
        return nullptr;
    }
    PyObject *cpp_conduit = PyObject_CallMethod(py_obj,
                                                "_pybind11_conduit_v1_",
                                                "yOy",
                                                PYBIND11_PLATFORM_ABI_ID,
                                                cpp_type_info_capsule,
                                                "raw_pointer_ephemeral");
    Py_DECREF(cpp_type_info_capsule);
    if (cpp_conduit == nullptr) {
        return nullptr;
    }
    void *raw_ptr = PyCapsule_GetPointer(cpp_conduit, cpp_type_info->name());
    Py_DECREF(cpp_conduit);
    if (PyErr_Occurred()) {
        return nullptr;
    }
    return raw_ptr;
}

template <typename T>
T *get_type_pointer_ephemeral(PyObject *py_obj) {
    void *raw_ptr = get_raw_pointer_ephemeral(py_obj, &typeid(T));
    if (raw_ptr == nullptr) {
        return nullptr;
    }
    return static_cast<T *>(raw_ptr);
}

} // namespace pybind11_conduit_v1
