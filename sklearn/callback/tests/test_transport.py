# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import threading
import time
from multiprocessing.connection import Client

import pytest

from sklearn.callback._transport import close_listener, open_listener


def _assert_no_leftover_threads(n_expected, timeout=1):
    """Assert that the background threads of the transport have exited.

    They exit asynchronously, hence the polling.
    """
    while (n_threads := threading.active_count()) > n_expected:
        timeout -= 0.05
        assert timeout > 0, f"{n_threads - n_expected} leftover threads"
        time.sleep(0.05)


@pytest.mark.thread_unsafe  # counts the threads of the whole process
def test_close_listener_stops_accept_thread():
    """Check that closing a listener doesn't leave its accept thread behind.

    An accept thread that is left behind also holds on to the last connection it
    accepted, because its frame keeps a reference to it. Both ends of that socket then
    stay open even though nothing reads from them anymore.
    """
    n_threads = threading.active_count()
    listener_handle = open_listener(lambda message: None)

    close_listener(listener_handle)
    close_listener(listener_handle)  # closing twice must be harmless.

    _assert_no_leftover_threads(n_threads)


@pytest.mark.thread_unsafe  # counts the threads of the whole process
def test_close_listener_whose_accept_thread_is_already_gone():
    """Check that closing a listener that failed on its own is harmless.

    An accept thread that stops because `accept` failed closes its listener on the way
    out, and the address stays registered until someone closes it, so `close_listener`
    then has an address it cannot connect to.
    """
    n_threads = threading.active_count()
    listener_handle = open_listener(lambda message: None)

    # Connecting without the authkey makes `accept` fail, which is how the accept thread
    # is stopped, here without unregistering the listener first.
    Client(listener_handle.address).close()
    _assert_no_leftover_threads(n_threads)

    close_listener(listener_handle)
