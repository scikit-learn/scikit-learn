"""
Telnet server.

Example usage::

    class MyTelnetApplication(TelnetApplication):
        def client_connected(self, telnet_connection):
            # Set CLI with simple prompt.
            telnet_connection.set_application(
                telnet_connection.create_prompt_application(...))

        def handle_command(self, telnet_connection, document):
            # When the client enters a command, just reply.
            telnet_connection.send('You said: %r\n\n' % document.text)

        ...

    a = MyTelnetApplication()
    TelnetServer(application=a, host='127.0.0.1', port=23).run()
"""
from __future__ import unicode_literals

import socket
import select

import threading
import os
import fcntl

from six import int2byte, text_type, binary_type
from codecs import getincrementaldecoder

from prompt_toolkit.enums import DEFAULT_BUFFER
from prompt_toolkit.eventloop.base import EventLoop
from prompt_toolkit.interface import CommandLineInterface, Application
from prompt_toolkit.layout.screen import Size
from prompt_toolkit.shortcuts import create_prompt_application
from prompt_toolkit.terminal.vt100_input import InputStream
from prompt_toolkit.terminal.vt100_output import Vt100_Output

from .log import logger
from .protocol import IAC, DO, LINEMODE, SB, MODE, SE, WILL, ECHO, NAWS, SUPPRESS_GO_AHEAD
from .protocol import TelnetProtocolParser
from .application import TelnetApplication

__all__ = (
    'TelnetServer',
)


def _initialize_telnet(connection):
    logger.info('Initializing telnet connection')

    # Iac Do Linemode
    connection.send(IAC + DO + LINEMODE)

    # Suppress Go Ahead. (This seems important for Putty to do correct echoing.)
    # This will allow bi-directional operation.
    connection.send(IAC + WILL + SUPPRESS_GO_AHEAD)

    # Iac sb
    connection.send(IAC + SB + LINEMODE + MODE + int2byte(0) + IAC + SE)

    # IAC Will Echo
    connection.send(IAC + WILL + ECHO)

    # Negotiate window size
    connection.send(IAC + DO + NAWS)


class _ConnectionStdout(object):
    """
    Wrapper around socket which provides `write` and `flush` methods for the
    Vt100_Output output.
    """
    def __init__(self, connection, encoding):
        self._encoding = encoding
        self._connection = connection
        self._buffer = []

    def write(self, data):
        assert isinstance(data, text_type)
        self._buffer.append(data.encode(self._encoding))
        self.flush()

    def flush(self):
        try:
            self._connection.send(b''.join(self._buffer))
        except socket.error as e:
            logger.error("Couldn't send data over socket: %s" % e)

        self._buffer = []


class TelnetConnection(object):
    """
    Class that represents one Telnet connection.
    """
    def __init__(self, conn, addr, application, server, encoding):
        assert isinstance(addr, tuple)  # (addr, port) tuple
        assert isinstance(application, TelnetApplication)
        assert isinstance(server, TelnetServer)
        assert isinstance(encoding, text_type)  # e.g. 'utf-8'

        self.conn = conn
        self.addr = addr
        self.application = application
        self.closed = False
        self.handling_command = True
        self.server = server
        self.encoding = encoding
        self.callback = None  # Function that handles the CLI result.

        # Create "Output" object.
        self.size = Size(rows=40, columns=79)

        # Initialize.
        _initialize_telnet(conn)

        # Create output.
        def get_size():
            return self.size
        self.stdout = _ConnectionStdout(conn, encoding=encoding)
        self.vt100_output = Vt100_Output(self.stdout, get_size, write_binary=False)

        # Create an eventloop (adaptor) for the CommandLineInterface.
        self.eventloop = _TelnetEventLoopInterface(server)

        # Set default CommandLineInterface.
        self.set_application(create_prompt_application())

        # Call client_connected
        application.client_connected(self)

        # Draw for the first time.
        self.handling_command = False
        self.cli._redraw()

    def set_application(self, app, callback=None):
        """
        Set ``CommandLineInterface`` instance for this connection.
        (This can be replaced any time.)

        :param cli: CommandLineInterface instance.
        :param callback: Callable that takes the result of the CLI.
        """
        assert isinstance(app, Application)
        assert callback is None or callable(callback)

        self.cli = CommandLineInterface(
            application=app,
            eventloop=self.eventloop,
            output=self.vt100_output)
        self.callback = callback

        # Create a parser, and parser callbacks.
        cb = self.cli.create_eventloop_callbacks()
        inputstream = InputStream(cb.feed_key)

        # Input decoder for stdin. (Required when working with multibyte
        # characters, like chinese input.)
        stdin_decoder_cls = getincrementaldecoder(self.encoding)
        stdin_decoder = [stdin_decoder_cls()]  # nonlocal

        # Tell the CLI that it's running. We don't start it through the run()
        # call, but will still want _redraw() to work.
        self.cli._is_running = True

        def data_received(data):
            """ TelnetProtocolParser 'data_received' callback """
            assert isinstance(data, binary_type)

            try:
                result = stdin_decoder[0].decode(data)
                inputstream.feed(result)
            except UnicodeDecodeError:
                stdin_decoder[0] = stdin_decoder_cls()
                return ''

        def size_received(rows, columns):
            """ TelnetProtocolParser 'size_received' callback """
            self.size = Size(rows=rows, columns=columns)
            cb.terminal_size_changed()

        self.parser = TelnetProtocolParser(data_received, size_received)

    def feed(self, data):
        """
        Handler for incoming data. (Called by TelnetServer.)
        """
        assert isinstance(data, binary_type)

        self.parser.feed(data)

        # Render again.
        self.cli._redraw()

        # When a return value has been set (enter was pressed), handle command.
        if self.cli.is_returning:
            try:
                return_value = self.cli.return_value()
            except (EOFError, KeyboardInterrupt) as e:
                # Control-D or Control-C was pressed.
                logger.info('%s, closing connection.', type(e).__name__)
                self.close()
                return

            # Handle CLI command
            self._handle_command(return_value)

    def _handle_command(self, command):
        """
        Handle command. This will run in a separate thread, in order not
        to block the event loop.
        """
        logger.info('Handle command %r', command)

        def in_executor():
            self.handling_command = True
            try:
                if self.callback is not None:
                    self.callback(self, command)
            finally:
                self.server.call_from_executor(done)

        def done():
            self.handling_command = False

            # Reset state and draw again. (If the connection is still open --
            # the application could have called TelnetConnection.close()
            if not self.closed:
                self.cli.reset()
                self.cli.buffers[DEFAULT_BUFFER].reset()
                self.cli.renderer.request_absolute_cursor_position()
                self.vt100_output.flush()
                self.cli._redraw()

        self.server.run_in_executor(in_executor)

    def erase_screen(self):
        """
        Erase output screen.
        """
        self.vt100_output.erase_screen()
        self.vt100_output.cursor_goto(0, 0)
        self.vt100_output.flush()

    def send(self, data):
        """
        Send text to the client.
        """
        assert isinstance(data, text_type)

        # When data is send back to the client, we should replace the line
        # endings. (We didn't allocate a real pseudo terminal, and the telnet
        # connection is raw, so we are responsible for inserting \r.)
        self.stdout.write(data.replace('\n', '\r\n'))
        self.stdout.flush()

    def close(self):
        """
        Close the connection.
        """
        self.application.client_leaving(self)

        self.conn.close()
        self.closed = True


class _TelnetEventLoopInterface(EventLoop):
    """
    Eventloop object to be assigned to `CommandLineInterface`.
    """
    def __init__(self, server):
        self._server = server

    def close(self):
        " Ignore. "

    def stop(self):
        " Ignore. "

    def run_in_executor(self, callback):
        self._server.run_in_executor(callback)

    def call_from_executor(self, callback, _max_postpone_until=None):
        self._server.call_from_executor(callback)

    def add_reader(self, fd, callback):
        raise NotImplementedError

    def remove_reader(self, fd):
        raise NotImplementedError


class TelnetServer(object):
    """
    Telnet server implementation.
    """
    def __init__(self, host='127.0.0.1', port=23, application=None, encoding='utf-8'):
        assert isinstance(host, text_type)
        assert isinstance(port, int)
        assert isinstance(application, TelnetApplication)
        assert isinstance(encoding, text_type)

        self.host = host
        self.port = port
        self.application = application
        self.encoding = encoding

        self.connections = set()

        self._calls_from_executor = []

        # Create a pipe for inter thread communication.
        self._schedule_pipe = os.pipe()
        fcntl.fcntl(self._schedule_pipe[0], fcntl.F_SETFL, os.O_NONBLOCK)

    @classmethod
    def create_socket(cls, host, port):
        # Create and bind socket
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        s.bind((host, port))

        s.listen(4)
        return s

    def run_in_executor(self, callback):
        threading.Thread(target=callback).start()

    def call_from_executor(self, callback):
        self._calls_from_executor.append(callback)

        if self._schedule_pipe:
            os.write(self._schedule_pipe[1], b'x')

    def _process_callbacks(self):
        """
        Process callbacks from `call_from_executor` in eventloop.
        """
        # Flush all the pipe content.
        os.read(self._schedule_pipe[0], 1024)

        # Process calls from executor.
        calls_from_executor, self._calls_from_executor = self._calls_from_executor, []
        for c in calls_from_executor:
            c()

    def run(self):
        """
        Run the eventloop for the telnet server.
        """
        listen_socket = self.create_socket(self.host, self.port)
        logger.info('Listening for telnet connections on %s port %r', self.host, self.port)

        try:
            while True:
                # Removed closed connections.
                self.connections = set([c for c in self.connections if not c.closed])

                # Ignore connections handling commands.
                connections = set([c for c in self.connections if not c.handling_command])

                # Wait for next event.
                read_list = (
                    [listen_socket, self._schedule_pipe[0]] +
                    [c.conn for c in connections])

                read, _, _ = select.select(read_list, [], [])

                for s in read:
                    # When the socket itself is ready, accept a new connection.
                    if s == listen_socket:
                        self._accept(listen_socket)

                    # If we receive something on our "call_from_executor" pipe, process
                    # these callbacks in a thread safe way.
                    elif s == self._schedule_pipe[0]:
                        self._process_callbacks()

                    # Handle incoming data on socket.
                    else:
                        self._handle_incoming_data(s)
        finally:
            listen_socket.close()

    def _accept(self, listen_socket):
        """
        Accept new incoming connection.
        """
        conn, addr = listen_socket.accept()
        connection = TelnetConnection(conn, addr, self.application, self, encoding=self.encoding)
        self.connections.add(connection)

        logger.info('New connection %r %r', *addr)

    def _handle_incoming_data(self, conn):
        """
        Handle incoming data on socket.
        """
        connection = [c for c in self.connections if c.conn == conn][0]
        data = conn.recv(1024)
        if data:
            connection.feed(data)
        else:
            self.connections.remove(connection)
