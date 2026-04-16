/*
    pybind11/gil.h: RAII helpers for managing the GIL

    Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#if defined(PYBIND11_SIMPLE_GIL_MANAGEMENT)

#    include "detail/common.h"
#    include "gil_simple.h"

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

using gil_scoped_acquire = gil_scoped_acquire_simple;
using gil_scoped_release = gil_scoped_release_simple;

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

#else

#    include "detail/common.h"
#    include "detail/internals.h"

#    include <cassert>

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)

PYBIND11_NAMESPACE_BEGIN(detail)

PYBIND11_WARNING_PUSH
PYBIND11_WARNING_DISABLE_GCC("-Wredundant-decls")

// forward declarations
PyThreadState *get_thread_state_unchecked();

PYBIND11_WARNING_POP

PYBIND11_NAMESPACE_END(detail)

/* The functions below essentially reproduce the PyGILState_* API using a RAII
 * pattern, but there are a few important differences:
 *
 * 1. When acquiring the GIL from an non-main thread during the finalization
 *    phase, the GILState API blindly terminates the calling thread, which
 *    is often not what is wanted. This API does not do this.
 *
 * 2. The gil_scoped_release function can optionally cut the relationship
 *    of a PyThreadState and its associated thread, which allows moving it to
 *    another thread (this is a fairly rare/advanced use case).
 *
 * 3. The reference count of an acquired thread state can be controlled. This
 *    can be handy to prevent cases where callbacks issued from an external
 *    thread would otherwise constantly construct and destroy thread state data
 *    structures.
 *
 * See the Python bindings of NanoGUI (http://github.com/wjakob/nanogui) for an
 * example which uses features 2 and 3 to migrate the Python thread of
 * execution to another thread (to run the event loop on the original thread,
 * in this case).
 */

class gil_scoped_acquire {
public:
    PYBIND11_NOINLINE gil_scoped_acquire() {
        auto &internals = detail::get_internals();
        tstate = internals.tstate.get();

        if (!tstate) {
            /* Check if the GIL was acquired using the PyGILState_* API instead (e.g. if
               calling from a Python thread). Since we use a different key, this ensures
               we don't create a new thread state and deadlock in PyEval_AcquireThread
               below. Note we don't save this state with internals.tstate, since we don't
               create it we would fail to clear it (its reference count should be > 0). */
            tstate = PyGILState_GetThisThreadState();
        }

        if (!tstate) {
            tstate = PyThreadState_New(internals.istate);
#    if defined(PYBIND11_DETAILED_ERROR_MESSAGES)
            if (!tstate) {
                pybind11_fail("scoped_acquire: could not create thread state!");
            }
#    endif
            tstate->gilstate_counter = 0;
            internals.tstate = tstate;
        } else {
            release = detail::get_thread_state_unchecked() != tstate;
        }

        if (release) {
            PyEval_AcquireThread(tstate);
        }

        inc_ref();
    }

    gil_scoped_acquire(const gil_scoped_acquire &) = delete;
    gil_scoped_acquire &operator=(const gil_scoped_acquire &) = delete;

    void inc_ref() { ++tstate->gilstate_counter; }

    PYBIND11_NOINLINE void dec_ref() {
        --tstate->gilstate_counter;
#    if defined(PYBIND11_DETAILED_ERROR_MESSAGES)
        if (detail::get_thread_state_unchecked() != tstate) {
            pybind11_fail("scoped_acquire::dec_ref(): thread state must be current!");
        }
        if (tstate->gilstate_counter < 0) {
            pybind11_fail("scoped_acquire::dec_ref(): reference count underflow!");
        }
#    endif
        if (tstate->gilstate_counter == 0) {
#    if defined(PYBIND11_DETAILED_ERROR_MESSAGES)
            if (!release) {
                pybind11_fail("scoped_acquire::dec_ref(): internal error!");
            }
#    endif
            // Make sure that PyThreadState_Clear is not recursively called by finalizers.
            // See issue #5827
            ++tstate->gilstate_counter;
            PyThreadState_Clear(tstate);
            --tstate->gilstate_counter;
            if (active) {
                PyThreadState_DeleteCurrent();
            }
            detail::get_internals().tstate.reset();
            release = false;
        }
    }

    /// This method will disable the PyThreadState_DeleteCurrent call and the
    /// GIL won't be released. This method should be used if the interpreter
    /// could be shutting down when this is called, as thread deletion is not
    /// allowed during shutdown. Check _Py_IsFinalizing() on Python 3.7+, and
    /// protect subsequent code.
    PYBIND11_NOINLINE void disarm() { active = false; }

    PYBIND11_NOINLINE ~gil_scoped_acquire() {
        dec_ref();
        if (release) {
            PyEval_SaveThread();
        }
    }

private:
    PyThreadState *tstate = nullptr;
    bool release = true;
    bool active = true;
};

class gil_scoped_release {
public:
    // PRECONDITION: The GIL must be held when this constructor is called.
    explicit gil_scoped_release(bool disassoc = false) : disassoc(disassoc) {
        assert(PyGILState_Check());
        // `get_internals()` must be called here unconditionally in order to initialize
        // `internals.tstate` for subsequent `gil_scoped_acquire` calls. Otherwise, an
        // initialization race could occur as multiple threads try `gil_scoped_acquire`.
        auto &internals = detail::get_internals();
        // NOLINTNEXTLINE(cppcoreguidelines-prefer-member-initializer)
        tstate = PyEval_SaveThread();
        if (disassoc) {
            internals.tstate.reset();
        }
    }

    gil_scoped_release(const gil_scoped_release &) = delete;
    gil_scoped_release &operator=(const gil_scoped_release &) = delete;

    /// This method will disable the PyThreadState_DeleteCurrent call and the
    /// GIL won't be acquired. This method should be used if the interpreter
    /// could be shutting down when this is called, as thread deletion is not
    /// allowed during shutdown. Check _Py_IsFinalizing() on Python 3.7+, and
    /// protect subsequent code.
    PYBIND11_NOINLINE void disarm() { active = false; }

    ~gil_scoped_release() {
        if (!tstate) {
            return;
        }
        // `PyEval_RestoreThread()` should not be called if runtime is finalizing
        if (active) {
            PyEval_RestoreThread(tstate);
        }
        if (disassoc) {
            detail::get_internals().tstate = tstate;
        }
    }

private:
    PyThreadState *tstate;
    bool disassoc;
    bool active = true;
};

PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

#endif // !PYBIND11_SIMPLE_GIL_MANAGEMENT
