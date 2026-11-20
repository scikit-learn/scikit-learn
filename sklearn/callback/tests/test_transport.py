# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import threading
import time
from multiprocessing.connection import Client

import pytest

from sklearn.callback._transport import (
    _TRANSPORT_THREAD_NAME_PREFIX,
    close_listener,
    open_listener,
)


def _transport_threads():
    return [
        thread
        for thread in threading.enumerate()
        if thread.name.startswith(_TRANSPORT_THREAD_NAME_PREFIX)
    ]


def _assert_accept_thread_running():
    """The listener's accept thread is alive and carries the transport name.

    Without this check, `_assert_no_leftover_threads` would pass vacuously if the
    threads were no longer named with the prefix.
    """
    names = [thread.name for thread in _transport_threads()]
    assert names == [f"{_TRANSPORT_THREAD_NAME_PREFIX}accept"], names


def _assert_no_leftover_threads(timeout=1):
    """Assert that the background threads of the transport have exited.

    They exit asynchronously, hence the polling. Only threads whose names start
    with the transport prefix are considered, so unrelated threads spawned or
    stopped in the background cannot hide a leak or fail the test.
    """
    while leftover := _transport_threads():
        timeout -= 0.05
        names = [thread.name for thread in leftover]
        assert timeout > 0, f"{len(leftover)} leftover threads: {names}"
        time.sleep(0.05)


@pytest.mark.thread_unsafe  # module-level listener registries
def test_close_listener_stops_accept_thread():
    """Check that closing a listener doesn't leave its accept thread behind.

    An accept thread that is left behind also holds on to the last connection it
    accepted, because its frame keeps a reference to it. Both ends of that socket then
    stay open even though nothing reads from them anymore.
    """
    listener_handle = open_listener(lambda message: None)
    _assert_accept_thread_running()

    close_listener(listener_handle)
    close_listener(listener_handle)  # closing twice must be harmless.

    _assert_no_leftover_threads()


@pytest.mark.thread_unsafe  # module-level listener registries
def test_close_listener_whose_accept_thread_is_already_gone():
    """Check that closing a listener that failed on its own is harmless.

    An accept thread that stops because `accept` failed closes its listener on the way
    out, and the address stays registered until someone closes it, so `close_listener`
    then has an address it cannot connect to.
    """
    listener_handle = open_listener(lambda message: None)
    _assert_accept_thread_running()

    # Connecting without the authkey makes `accept` fail, which is how the accept thread
    # is stopped, here without unregistering the listener first.
    Client(listener_handle.address).close()
    _assert_no_leftover_threads()

    close_listener(listener_handle)


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

    listener_handle = open_listener(consumer)

    def send_one():
        client = Client(listener_handle.address, authkey=listener_handle.authkey)
        try:
            client.send(None)
            client.recv()
        finally:
            client.close()

    sender = threading.Thread(target=send_one, daemon=True)
    sender.start()
    try:
        assert started.wait(timeout=2)
        assert sorted(thread.name for thread in _transport_threads()) == [
            f"{_TRANSPORT_THREAD_NAME_PREFIX}accept",
            f"{_TRANSPORT_THREAD_NAME_PREFIX}handle",
        ]
    finally:
        release.set()
        sender.join(timeout=2)
        close_listener(listener_handle)

    assert not sender.is_alive()
    _assert_no_leftover_threads()


@pytest.mark.thread_unsafe  # module-level listener registries
def test_rejected_handshake_keeps_accept_thread():
    """A wrong auth digest does not stop the accept thread.

    `AuthenticationError` is not an `OSError`. `close_listener` wakes the accept
    thread by closing without answering, which raises `EOFError`. A peer that
    answers and is rejected must leave the listener open.
    """
    received = []
    listener_handle = open_listener(received.append)
    try:
        _assert_accept_thread_running()

        client = Client(listener_handle.address)  # no authkey: handshake is manual
        try:
            client.recv_bytes(256)  # the challenge
            client.send_bytes(b"not-a-valid-mac")
            try:
                client.recv_bytes(256)  # the server's rejection, if it is still open
            except (EOFError, OSError):
                pass
        finally:
            client.close()

        good = Client(listener_handle.address, authkey=listener_handle.authkey)
        try:
            good.send("ok")
            good.recv()
        finally:
            good.close()

        assert received == ["ok"]
        # The handler for "ok" exits after the client closes, so wait it out.
        # The accept thread itself must stay.
        timeout = 1
        accept_name = f"{_TRANSPORT_THREAD_NAME_PREFIX}accept"
        while (names := [thread.name for thread in _transport_threads()]) != [
            accept_name
        ]:
            assert accept_name in names, names
            timeout -= 0.05
            assert timeout > 0, names
            time.sleep(0.05)
    finally:
        close_listener(listener_handle)

    _assert_no_leftover_threads()
