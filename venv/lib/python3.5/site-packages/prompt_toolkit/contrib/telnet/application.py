"""
Interface for Telnet applications.
"""
from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from six import with_metaclass

__all__ = (
    'TelnetApplication',
)


class TelnetApplication(with_metaclass(ABCMeta, object)):
    """
    The interface which has to be implemented for any telnet application.
    An instance of this class has to be passed to `TelnetServer`.
    """
    @abstractmethod
    def client_connected(self, telnet_connection):
        """
        Called when a new client was connected.

        Probably you want to call `telnet_connection.set_cli` here to set a
        the CommandLineInterface instance to be used.
        Hint: Use the following shortcut: `prompt_toolkit.shortcuts.create_cli`
        """

    @abstractmethod
    def client_leaving(self, telnet_connection):
        """
        Called when a client quits.
        """
