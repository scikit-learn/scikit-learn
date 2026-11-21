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
import weakref
from multiprocessing import AuthenticationError
from multiprocessing.connection import Client, Listener
from threading import Lock, Thread

# Process-local registry of live listeners, keyed by address. There is one listener
# per call to `open_listener` (one per ScoringMonitor instance, one per ProgressBar
# fit). It is intentionally module-level so that `is_reusable` can tell, from an
# unpickled copy that does not hold the listener, whether this process is the one that
# opened it.
_listeners: dict[str, Listener] = {}

# Names need not be unique. The shared prefix lets tests and debuggers tell these
# threads apart from unrelated ones (joblib, pytest plugins, ...).
_TRANSPORT_THREAD_NAME_PREFIX = "sklearn_callback_transport_"


class Channel:
    """A two-way endpoint for the messages of a callback.

    The process that creates a channel opens a listener that receives the messages and
    hands them to `message_consumer`.

    The copy of the channel that a worker unpickles opens its own connection on the
    first `send` and reuses it for the following ones.

    Parameters
    ----------
    message_consumer : callable
        A one-argument function, `message_consumer(message)`, that processes an
        incoming message to update the callback's state. This callable may be called
        from multiple different threads, and must therefore behave in a thread-safe
        manner. In particular, it must not be a bound method of the callback.
    """

    # Class attributes so that a copy unpickled in a worker starts out without
    # connection and consumer, and opens its own connection on first send.
    _connection = None
    _message_consumer = None

    def __init__(self, message_consumer):
        self._message_consumer = message_consumer
        self._lock = Lock()
        self._listener_address, self._authkey = open_listener(message_consumer)
        weakref.finalize(self, close_listener, self._listener_address)

    def __getstate__(self):
        # exclude the connection and lock: a worker that is itself pickling this channel
        # should not share its live client, and the copied channel opens its own
        # connection and lock when unpickled.
        return {"_listener_address": self._listener_address, "_authkey": self._authkey}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._lock = Lock()

    def send(self, message):
        """Deliver `message` to whoever is listening on the other end.

        There are two possible delivery paths:

        - In-process fast path: `send` is called on the instance that created the
          channel. The message consumer can directly be called without any serialization
          overhead.

        - Cross-process path: `send` is called from an unpickled copy of the channel.
          It connects to the listener on first use and reuses that connection for the
          following messages. `send` then waits for an acknowledgement from the main
          process so that, by the time it returns, the message has actually been
          processed by the consumer. Access to the connection is serialized with a lock
          because nested thread parallelism inside a worker process shares the same
          unpickled channel across threads.
        """
        if self._message_consumer is not None:  # fast path
            self._message_consumer(message)
            return

        with self._lock:
            if self._connection is None:
                self._connection = Client(self._listener_address, authkey=self._authkey)
            try:
                self._connection.send(message)
                self._connection.recv()
            except BaseException:
                connection = self._connection
                self._connection = None
                connection.close()
                raise

    def disconnect(self):
        """Release the connection in this process, leaving the listener open.

        A later `send` reconnects, so this only frees resources. Non-propagated
        callbacks, i.e. callbacks that are meant to be torn down in the workers, should
        call this in their `teardown` hook to free resources.
        """
        # Close explicitly instead of relying on the Client being collected: the
        # handler thread in the listening process stays blocked in `recv` until this end
        # closes, and this channel may stay reachable (e.g. cached by a long-lived
        # worker) well after the fit it was used for.
        with self._lock:
            connection = self._connection
            self._connection = None
            if connection is not None:
                connection.close()

    def close(self):
        """Close the channel, i.e. its listener and its connection, if it has one.

        Auto-propagated callbacks, i.e. callbacks that create one channel per fit,
        should call this in their `teardown` hook because the channel is not meant to be
        reused afterwards.
        """
        self.disconnect()
        close_listener(self._listener_address)

    def is_reusable(self):
        """Whether this channel, freshly unpickled, can be used as is.

        Helper for callbacks that create their channel eagerly (e.g. in `__init__`) and
        therefore have to decide, on unpickling, whether to keep the inherited channel
        or create a fresh one. The channel is not reusable when:

        - We are the process that originally opened the listener. Reusing the channel
          would deliver messages to the original instance's consumer, through that
          listener, instead of the unpickled instance's.

        - The listener is no longer reachable, e.g. unpickling in a fresh interpreter,
          or on a host that cannot reach the original listener.
        """
        if self._listener_address in _listeners:
            return False
        try:
            Client(self._listener_address, authkey=self._authkey).close()
        except OSError:
            return False
        return True


def open_listener(message_consumer):
    """Create a listener for incoming messages on the main process.

    Also registers the listener in the module-level dict.

    Parameters
    ----------
    message_consumer : callable
        Called with each incoming message.

    Returns
    -------
    address : str
        Address of the local endpoint the listener is bound to.

    authkey : bytes
        Shared secret used to authenticate connections to the listener.
    """
    authkey = os.urandom(32)
    # `backlog` is the kernel's accept queue size. The stdlib default of 1 is too
    # small: while the accept thread is busy with the authentication handshake of
    # an in-flight connection, any concurrent worker that calls `Client(...)` on a
    # full queue gets `ConnectionRefusedError` (macOS in particular enforces this
    # strictly). At most one connection is active at a time per process-local copy, so
    # 128 is comfortably more than enough.
    listener = Listener(authkey=authkey, backlog=128)
    _listeners[listener.address] = listener

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
        finally:
            conn.close()

    def _accept():
        try:
            while listener.address in _listeners:
                try:
                    conn = listener.accept()
                except AuthenticationError:
                    # Rejected handshake. The listening socket is still open.
                    continue
                except (OSError, EOFError):
                    break  # woken up by `close_listener`, or the listener failed
                Thread(
                    target=_handle,
                    args=(conn,),
                    name=f"{_TRANSPORT_THREAD_NAME_PREFIX}handle",
                    daemon=True,
                ).start()
                del conn
        finally:
            listener.close()

    Thread(
        target=_accept, name=f"{_TRANSPORT_THREAD_NAME_PREFIX}accept", daemon=True
    ).start()
    return listener.address, authkey


def close_listener(address):
    """Stop listening at `address` and free its background threads."""
    if _listeners.pop(address, None) is None:
        return  # already closed

    try:
        # The accept thread is waiting in `accept`, and closing the listener would not
        # interrupt that call. Connecting to it does: we connect without the authkey so
        # that the authentication handshake fails, and that failure is what makes the
        # accept thread break out of its loop.
        Client(address).close()
    except OSError:
        pass  # the accept thread is already gone, e.g. its listener failed
