# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

"""Cross-process transport for callback messages.

When a callback is registered on an estimator that uses multiple worker processes, every
worker ends up with its own copy of the callback (sent there by pickling). However, the
user-visible state (e.g. logs being filled in or progress bars advancing) lives on the
main process. This module provides a way for the worker copies to ship their messages
back to the main process over a local endpoint (a UNIX socket on Unix, a Windows named
pipe on Windows).

Remark: we don't use a `multiprocessing.Manager` because its proxy objects become
unusable once the Manager subprocess is gone, e.g. after unpickling in a fresh
interpreter. The only ways to work around that either rely on multiprocessing private
API or give up cross-process capabilities.
"""

import os
from multiprocessing.connection import Client, Connection, Listener
from threading import Thread
from typing import Callable, NamedTuple


class ListenerHandle(NamedTuple):
    """A picklable reference to a main-process listener.

    Attributes
    ----------
    address : str
        Address of the local endpoint the listener is bound to.
        Workers use this to connect, and the main process also uses it as a key into
        this module's registries to find the live listener.

    authkey : bytes
        Shared secret used to authenticate connections to the listener.
    """

    address: str
    authkey: bytes


# Process-local registries of live listeners and their message consumers, keyed by the
# listener's address:
#
# - _listeners: live listeners. It contains one listener per call to `open_listener`.
#               - ScoringMonitor opens one listener per callback instance.
#               - ProgressBar opens one listener per fit.
#
# - _message_consumers: callables that handle messages from the listeners, e.g. putting
#                       in the progress bar queue or appending to the scoring monitor
#                       log. There's one per listener.
#
# These are intentionally module-level: each process (main or worker) has its
# own copy when the module is imported, and a worker that receives a pickled
# `ListenerHandle` from the main process will find these dicts empty, which is how
# `send` decides to go over the socket instead of directly calling the message consumer.
_listeners: dict[str, Listener] = {}
_message_consumers: dict[str, Callable] = {}

# Worker-side cache of client connections back to the main-process listeners, keyed by
# listener address. The first cross-process send opens a client. Subsequent sends to the
# same listener reuse the cached connection.
_worker_connections: dict[str, Connection] = {}


def open_listener(message_consumer):
    """Create a listener for incoming messages on the main process.

    Also registers the listener and its message consumer in the module-level dicts.

    Parameters
    ----------
    message_consumer : callable
        A one-argument function, `message_consumer(message)`, that processes incoming
        message to update the callback's state.

    Returns
    -------
    listener_handle : ListenerHandle
        A reference to the listener.
    """
    authkey = os.urandom(32)
    # `backlog` is the kernel's accept queue size. The stdlib default of 1 is too
    # small: while the accept thread is busy with the authentication handshake of
    # an in-flight connection, any concurrent worker that calls `Client(...)` on a
    # full queue gets `ConnectionRefusedError` (macOS in particular enforces this
    # strictly). With the worker-side connection cache (see `_worker_connections`)
    # only one `connect()` per worker per listener is ever needed, so 128 is comfortably
    # more than enough.
    listener = Listener(authkey=authkey, backlog=128)
    listener_handle = ListenerHandle(address=listener.address, authkey=authkey)

    _listeners[listener_handle.address] = listener
    _message_consumers[listener_handle.address] = message_consumer

    def _handle(conn):
        # Read messages until the worker disconnects. After processing each message,
        # send a one-byte acknowledgement so that the worker-side `send` only returns
        # once the message has actually been consumed here.
        try:
            while True:
                message_consumer(conn.recv())
                conn.send(None)
        except (EOFError, OSError):
            return

    def _accept():
        while True:
            try:
                conn = listener.accept()
            except OSError:
                return
            Thread(target=_handle, args=(conn,), daemon=True).start()

    Thread(target=_accept, daemon=True).start()
    return listener_handle


def close_listener(listener_handle):
    """Stop listening for `listener_handle` and free its background threads."""
    _message_consumers.pop(listener_handle.address)
    _listeners.pop(listener_handle.address).close()


def can_reuse_listener(listener_handle):
    """Whether the listener at `listener_handle` is usable from this process.

    Helper for callbacks that open their listener eagerly (e.g. in `__init__`) and
    therefore have to decide, on unpickling, whether to keep the inherited handle
    or open a fresh listener. The listener is not reusable when:

    - We are the process that originally opened the listener. Reusing the handle
      would route messages through the in-process fast path of `send`, into the
      original instance's message consumer instead of the unpickled instance's.

    - The listener is no longer reachable, e.g. unpickling in a fresh interpreter,
      or on a host that cannot reach the original listener.
    """
    if listener_handle.address in _listeners:
        return False
    try:
        Client(listener_handle.address, authkey=listener_handle.authkey).close()
    except OSError:
        return False
    return True


def send(listener_handle, message):
    """Deliver `message` to whoever is listening on `listener_handle`.

    There are two possible delivery paths:

    - In-process fast path: `send` is called in the same process that called
      `open_listener` for this listener handle. The message consumer can directly be
      called without any serialization overhead.

    - Cross-process path: `send` is called in a different process. The worker opens
      a `Client` connection to the main-process listener on first use and caches it
      in `_worker_connections`, so all subsequent messages reuse the same socket.
      `send` then waits for an acknowledgement from the main process so that, by
      the time it returns, the message has actually been processed by the consumer.
    """
    message_consumer = _message_consumers.get(listener_handle.address)
    if message_consumer is not None:  # fast path
        message_consumer(message)
        return

    address = listener_handle.address
    connection = _worker_connections.get(address)
    if connection is None:
        connection = Client(address, authkey=listener_handle.authkey)
        _worker_connections[address] = connection
    connection.send(message)
    connection.recv()
