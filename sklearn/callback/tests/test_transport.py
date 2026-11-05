# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import gc
import pickle
import threading
import time

import pytest

from sklearn.callback._transport import Channel


def _assert_no_leftover_threads(baseline, timeout=1):
    """Assert that the background threads of the transport have exited.

    They exit asynchronously, hence the polling.
    """
    while extra := set(threading.enumerate()) - baseline:
        timeout -= 0.05
        assert timeout > 0, f"{len(extra)} leftover threads"
        time.sleep(0.05)


@pytest.mark.thread_unsafe  # inspects the threads of the whole process
def test_no_leftover_connection_threads_after_worker_copy_is_dropped():
    """Check that dropping a worker copy of a channel frees its connection.

    The copy a worker unpickles opens a connection on first `send`. That connection is
    held by the copy, so collecting the copy closes it, which makes the handler thread
    on the listening side exit.
    """
    channel = Channel(lambda message: None)

    # Snapshot after initializing the channel so the listener's accept thread is
    # part of the baseline.
    baseline = set(threading.enumerate())

    copy = pickle.loads(pickle.dumps(channel))
    copy.send(None)

    assert set(threading.enumerate()) - baseline

    del copy
    gc.collect()

    _assert_no_leftover_threads(baseline)

    channel.close()
