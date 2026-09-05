# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import gc
import pickle
import threading
import time

import pytest

from sklearn.callback._transport import Channel


def _assert_no_leftover_threads(n_expected, timeout=1):
    """Assert that the background threads of the transport have exited.

    They exit asynchronously, hence the polling.
    """
    while (n_threads := threading.active_count()) > n_expected:
        timeout -= 0.05
        assert timeout > 0, f"{n_threads - n_expected} leftover threads"
        time.sleep(0.05)


@pytest.mark.thread_unsafe  # counts the threads of the whole process
def test_no_leftover_connection_threads_after_worker_copy_is_dropped():
    """Check that dropping a worker copy of a channel frees its connection.

    The copy a worker unpickles opens a connection on first `send`. That connection is
    held by the copy, so collecting the copy closes it, which makes the handler thread
    on the listening side exit.
    """
    channel = Channel(lambda message: None)

    # Count after initializing the channel so that the number already includes the
    # listener's accept thread.
    n_threads = threading.active_count()

    copy = pickle.loads(pickle.dumps(channel))
    copy.send(None)

    assert threading.active_count() > n_threads

    del copy
    gc.collect()

    _assert_no_leftover_threads(n_threads)

    channel.close()
