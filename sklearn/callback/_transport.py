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
from multiprocessing.connection import Client, Listener
from threading import Thread
from typing import NamedTuple


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
_listeners = {}
_message_consumers = {}


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
    listener = Listener(authkey=authkey)
    listener_handle = ListenerHandle(address=listener.address, authkey=authkey)

    _listeners[listener_handle.address] = listener
    _message_consumers[listener_handle.address] = message_consumer

    def _handle(conn):
        try:
            while True:
                message_consumer(conn.recv())
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


def send(listener_handle, message):
    """Deliver `message` to whoever is listening on `listener_handle`.

    There are two possible delivery paths:

    - In-process fast path: `send` is called in the same process that called
      `open_listener` for this listener handle. The message consumer can directly be
      called without any serialization overhead.

    - Cross-process path: `send` is called in a different process. A fresh socket
      connection is opened to the main-process listener and the message is sent over it.
    """
    message_consumer = _message_consumers.get(listener_handle.address)
    if message_consumer is not None:  # fast path
        message_consumer(message)
        return

    try:
        connection = Client(listener_handle.address, authkey=listener_handle.authkey)
        connection.send(message)
    finally:
        connection.close()
