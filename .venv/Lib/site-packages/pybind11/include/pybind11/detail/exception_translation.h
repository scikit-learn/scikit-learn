/*
    pybind11/detail/exception_translation.h: means to translate C++ exceptions to Python exceptions

    Copyright (c) 2024 The Pybind Development Team.

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include "common.h"
#include "internals.h"

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
PYBIND11_NAMESPACE_BEGIN(detail)

// Apply all the extensions translators from a list
// Return true if one of the translators completed without raising an exception
// itself. Return of false indicates that if there are other translators
// available, they should be tried.
inline bool apply_exception_translators(std::forward_list<ExceptionTranslator> &translators) {
    auto last_exception = std::current_exception();

    for (auto &translator : translators) {
        try {
            translator(last_exception);
            return true;
        } catch (...) {
            last_exception = std::current_exception();
        }
    }
    return false;
}

inline void try_translate_exceptions() {
    /* When an exception is caught, give each registered exception
        translator a chance to translate it to a Python exception. First
        all module-local translators will be tried in reverse order of
        registration. If none of the module-locale translators handle
        the exception (or there are no module-locale translators) then
        the global translators will be tried, also in reverse order of
        registration.

        A translator may choose to do one of the following:

        - catch the exception and call py::set_error()
            to set a standard (or custom) Python exception, or
        - do nothing and let the exception fall through to the next translator, or
        - delegate translation to the next translator by throwing a new type of exception.
        */

    bool handled = with_exception_translators(
        [&](std::forward_list<ExceptionTranslator> &exception_translators,
            std::forward_list<ExceptionTranslator> &local_exception_translators) {
            if (detail::apply_exception_translators(local_exception_translators)) {
                return true;
            }
            if (detail::apply_exception_translators(exception_translators)) {
                return true;
            }
            return false;
        });

    if (!handled) {
        set_error(PyExc_SystemError, "Exception escaped from default exception translator!");
    }
}

PYBIND11_NAMESPACE_END(detail)
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
