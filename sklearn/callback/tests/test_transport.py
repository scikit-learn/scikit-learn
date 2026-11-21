# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import gc
import pickle
import threading
import time
from collections import Counter
from multiprocessing.connection import Client

import pytest

from sklearn.callback._transport import (
    _TRANSPORT_THREAD_NAME_PREFIX,
    Channel,
)


def _transport_threads():
    return [
        thread
        for thread in threading.enumerate()
        if thread.name.startswith(_TRANSPORT_THREAD_NAME_PREFIX)
    ]


def _snapshot_transport_threads():
    # Thread objects, not their names or counts. A listener finalized on another
    # thread can exit between the snapshot and the check; comparing counts then
    # hides the thread this test just started, or reports it as a leak.
    return set(_transport_threads())


def _new_transport_threads(baseline):
    return [thread for thread in _transport_threads() if thread not in baseline]


def _assert_accept_thread_running(baseline, timeout=1):
    """This test's accept thread is alive and carries the transport name.

    Threads already in `baseline` are ignored. Without this check,
    `_assert_no_leftover_threads` would pass vacuously if the threads were no
    longer named with the prefix.
    """
    accept = f"{_TRANSPORT_THREAD_NAME_PREFIX}accept"
    while True:
        new = _new_transport_threads(baseline)
        if [thread.name for thread in new] == [accept] and new[0].is_alive():
            return
        timeout -= 0.05
        assert timeout > 0, [thread.name for thread in new]
        time.sleep(0.05)


def _assert_no_leftover_threads(baseline, timeout=1):
    """Assert that the background threads started by this test have exited.

    They exit asynchronously, hence the polling. Only threads whose names start
    with the transport prefix are considered, and only the ones that were not
    already in `baseline`, so unrelated threads and threads left by earlier
    tests cannot hide a leak or fail the test.
    """
    while new := _new_transport_threads(baseline):
        timeout -= 0.05
        assert timeout > 0, f"leftover threads: {[thread.name for thread in new]}"
        time.sleep(0.05)


@pytest.mark.thread_unsafe  # module-level listener registries
def test_no_leftover_connection_threads_after_worker_copy_is_dropped():
    """Check that dropping a worker copy of a channel frees its connection.

    The copy a worker unpickles opens a connection on first `send`. That connection is
    held by the copy, so collecting the copy closes it, which makes the handler thread
    on the listening side exit.
    """
    channel = Channel(lambda message: None)

    # Snapshot after initializing the channel so the listener's accept thread is
    # part of the baseline.
    baseline = _snapshot_transport_threads()

    copy = pickle.loads(pickle.dumps(channel))
    copy.send(None)

    assert _new_transport_threads(baseline)

    del copy
    gc.collect()

    _assert_no_leftover_threads(baseline)

    channel.close()


@pytest.mark.thread_unsafe  # module-level listener registries
def test_disconnect_frees_connection_while_copy_is_alive():
    """Check that `disconnect` closes the connection before the copy is collected.

    Simulates a non-propagated callback (e.g. `ScoringMonitor`) running in a loky
    worker: `teardown` calls `disconnect` to free the socket and the main-process
    handler thread right away, without relying on the unpickled copy being collected.

    The callback context currently doesn't create reference cycles that keep callbacks
    alive, but this must not be a requirement for the transport to free its resources:
    the copy is made part of a reference cycle here, so that it can only be reclaimed
    by the cyclic GC.
    """
    channel = Channel(lambda message: None)
    # Snapshot after initializing the channel so the listener's accept thread is
    # part of the baseline.
    baseline = _snapshot_transport_threads()

    copy = pickle.loads(pickle.dumps(channel))
    copy.self_ref = copy
    copy.send(None)
    assert _new_transport_threads(baseline)

    copy.disconnect()
    _assert_no_leftover_threads(baseline)
    assert copy._connection is None

    channel.close()


@pytest.mark.thread_unsafe  # module-level listener registries
def test_close_listener_stops_accept_thread():
    """Check that closing a channel doesn't leave its accept thread behind.

    An accept thread that is left behind also holds on to the last connection it
    accepted, because its frame keeps a reference to it. Both ends of that socket then
    stay open even though nothing reads from them anymore.
    """
    baseline = _snapshot_transport_threads()
    channel = Channel(lambda message: None)
    _assert_accept_thread_running(baseline)

    channel.close()
    channel.close()  # closing twice must be harmless.

    _assert_no_leftover_threads(baseline)


@pytest.mark.thread_unsafe  # module-level listener registries
def test_close_listener_whose_accept_thread_is_already_gone():
    """Check that closing a channel whose accept thread already exited is harmless.

    An accept thread that stops because `accept` failed closes its listener on the way
    out, and the address stays registered until someone closes it, so `close` then has
    an address it cannot connect to.
    """
    baseline = _snapshot_transport_threads()
    channel = Channel(lambda message: None)
    _assert_accept_thread_running(baseline)

    # Connecting without the authkey makes `accept` fail, which is how the accept thread
    # is stopped, here without unregistering the listener first.
    Client(channel._listener_address).close()
    _assert_no_leftover_threads(baseline)

    channel.close()


@pytest.mark.thread_unsafe  # module-level listener registries
def test_handle_thread_is_named_and_exits():
    """The per-connection handler thread uses the transport name and then exits.

    The consumer blocks until the name has been observed, otherwise the handler
    can return before `threading.enumerate` sees it.
    """
    started = threading.Event()
    release = threading.Event()

    def consumer(_message):
        started.set()
        release.wait()

    baseline = _snapshot_transport_threads()
    channel = Channel(consumer)

    def send_one():
        client = Client(channel._listener_address, authkey=channel._authkey)
        try:
            client.send(None)
            client.recv()
        finally:
            client.close()

    sender = threading.Thread(target=send_one, daemon=True)
    sender.start()
    try:
        assert started.wait(timeout=2)
        new = _new_transport_threads(baseline)
        assert Counter(thread.name for thread in new) == Counter(
            {
                f"{_TRANSPORT_THREAD_NAME_PREFIX}accept": 1,
                f"{_TRANSPORT_THREAD_NAME_PREFIX}handle": 1,
            }
        ), [thread.name for thread in new]
    finally:
        release.set()
        sender.join(timeout=2)
        channel.close()

    assert not sender.is_alive()
    _assert_no_leftover_threads(baseline)


@pytest.mark.thread_unsafe  # module-level listener registries
def test_rejected_handshake_keeps_accept_thread():
    """A wrong auth digest does not stop the accept thread.

    `AuthenticationError` is not an `OSError`. `close` wakes the accept thread by
    closing without answering, which raises `EOFError`. A peer that answers and is
    rejected must leave the listener open.
    """
    received = []
    baseline = _snapshot_transport_threads()
    channel = Channel(received.append)
    try:
        _assert_accept_thread_running(baseline)

        client = Client(channel._listener_address)  # no authkey: handshake is manual
        try:
            client.recv_bytes(256)  # the challenge
            client.send_bytes(b"not-a-valid-mac")
            try:
                client.recv_bytes(256)  # the server's rejection, if it is still open
            except (EOFError, OSError):
                pass
        finally:
            client.close()

        good = Client(channel._listener_address, authkey=channel._authkey)
        try:
            good.send("ok")
            good.recv()
        finally:
            good.close()

        assert received == ["ok"]
        # The handler for "ok" exits after the client closes, so wait it out.
        # The accept thread itself must stay.
        accept = f"{_TRANSPORT_THREAD_NAME_PREFIX}accept"
        timeout = 1
        while Counter(thread.name for thread in _new_transport_threads(baseline)) != (
            Counter({accept: 1})
        ):
            names = [thread.name for thread in _new_transport_threads(baseline)]
            assert names.count(accept) >= 1, names
            timeout -= 0.05
            assert timeout > 0, names
            time.sleep(0.05)
    finally:
        channel.close()

    _assert_no_leftover_threads(baseline)
