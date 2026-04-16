// Copyright (c) 2016-2025 The Pybind Development Team.
// All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

#pragma once

#include "detail/common.h"

#include <cassert>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

class gil_scoped_acquire_simple {
    PyGILState_STATE state;

public:
    gil_scoped_acquire_simple() : state{PyGILState_Ensure()} {}
    gil_scoped_acquire_simple(const gil_scoped_acquire_simple &) = delete;
    gil_scoped_acquire_simple &operator=(const gil_scoped_acquire_simple &) = delete;
    ~gil_scoped_acquire_simple() { PyGILState_Release(state); }
};

class gil_scoped_release_simple {
    PyThreadState *state;

public:
    // PRECONDITION: The GIL must be held when this constructor is called.
    gil_scoped_release_simple() {
        assert(PyGILState_Check());
        state = PyEval_SaveThread();
    }
    gil_scoped_release_simple(const gil_scoped_release_simple &) = delete;
    gil_scoped_release_simple &operator=(const gil_scoped_release_simple &) = delete;
    ~gil_scoped_release_simple() { PyEval_RestoreThread(state); }
};

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)
