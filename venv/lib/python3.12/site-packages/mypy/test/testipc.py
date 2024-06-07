from __future__ import annotations

import sys
import time
from multiprocessing import Queue, get_context
from unittest import TestCase, main

import pytest

from mypy.ipc import IPCClient, IPCServer

CONNECTION_NAME = "dmypy-test-ipc"


def server(msg: str, q: Queue[str]) -> None:
    server = IPCServer(CONNECTION_NAME)
    q.put(server.connection_name)
    data = ""
    while not data:
        with server:
            server.write(msg)
            data = server.read()
    server.cleanup()


def server_multi_message_echo(q: Queue[str]) -> None:
    server = IPCServer(CONNECTION_NAME)
    q.put(server.connection_name)
    data = ""
    with server:
        while data != "quit":
            data = server.read()
            server.write(data)
    server.cleanup()


class IPCTests(TestCase):
    def setUp(self) -> None:
        if sys.platform == "linux":
            # The default "fork" start method is potentially unsafe
            self.ctx = get_context("forkserver")
        else:
            self.ctx = get_context("spawn")

    def test_transaction_large(self) -> None:
        queue: Queue[str] = self.ctx.Queue()
        msg = "t" * 200000  # longer than the max read size of 100_000
        p = self.ctx.Process(target=server, args=(msg, queue), daemon=True)
        p.start()
        connection_name = queue.get()
        with IPCClient(connection_name, timeout=1) as client:
            assert client.read() == msg
            client.write("test")
        queue.close()
        queue.join_thread()
        p.join()

    def test_connect_twice(self) -> None:
        queue: Queue[str] = self.ctx.Queue()
        msg = "this is a test message"
        p = self.ctx.Process(target=server, args=(msg, queue), daemon=True)
        p.start()
        connection_name = queue.get()
        with IPCClient(connection_name, timeout=1) as client:
            assert client.read() == msg
            client.write(
                ""
            )  # don't let the server hang up yet, we want to connect again.

        with IPCClient(connection_name, timeout=1) as client:
            assert client.read() == msg
            client.write("test")
        queue.close()
        queue.join_thread()
        p.join()
        assert p.exitcode == 0

    def test_multiple_messages(self) -> None:
        queue: Queue[str] = self.ctx.Queue()
        p = self.ctx.Process(
            target=server_multi_message_echo, args=(queue,), daemon=True
        )
        p.start()
        connection_name = queue.get()
        with IPCClient(connection_name, timeout=1) as client:
            # "foo bar" with extra accents on letters.
            # In UTF-8 encoding so we don't confuse editors opening this file.
            fancy_text = b"f\xcc\xb6o\xcc\xb2\xf0\x9d\x91\x9c \xd0\xb2\xe2\xb7\xa1a\xcc\xb6r\xcc\x93\xcd\x98\xcd\x8c"
            client.write(fancy_text.decode("utf-8"))
            assert client.read() == fancy_text.decode("utf-8")

            client.write("Test with spaces")
            client.write("Test write before reading previous")
            time.sleep(
                0
            )  # yield to the server to force reading of all messages by server.
            assert client.read() == "Test with spaces"
            assert client.read() == "Test write before reading previous"

            client.write("quit")
            assert client.read() == "quit"
        queue.close()
        queue.join_thread()
        p.join()
        assert p.exitcode == 0

    # Run test_connect_twice a lot, in the hopes of finding issues.
    # This is really slow, so it is skipped, but can be enabled if
    # needed to debug IPC issues.
    @pytest.mark.skip
    def test_connect_alot(self) -> None:
        t0 = time.time()
        for i in range(1000):
            try:
                print(i, "start")
                self.test_connect_twice()
            finally:
                t1 = time.time()
                print(i, t1 - t0)
                sys.stdout.flush()
                t0 = t1


if __name__ == "__main__":
    main()
