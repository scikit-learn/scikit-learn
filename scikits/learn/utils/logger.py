"""
Logging utility.

Author: Feth Arezki <feth A+ tuttu.info>

This module exposes 'log()' in an attempt to provide a drop in replacement
to this use of print statements throughout the code.


LOGGERS is a dict with logging domain as keys and loggers as values.
Every logger created or accessed through logtools() is referenced in it.

Simple uses
~~~~~~~~~~~

log() will trigger the general log system (to stdout):
>>> log("message with %(placeholder)s but no replacement")
message with %(placeholder)s but no replacement
>>> log("message with %s placeholder", 1)
message with 1 placeholder
>>> log("message with %s + %s placeholder", (1, 1))
message with 1 + 1 placeholder
>>> log("message with %(placeholder)s %(quantifier)s",
... {'placeholder': 'placeholder replaced', 'quantifier': 'twice'})
message with placeholder replaced twice


Criticity
~~~~~~~~~

However, log() takes an optional 'verbosity' int argument that must be
higher than the logger's threshold for the message to be displayed.
set_log_threshold() lets you set this value.

You want to use the generic log levels:
    10 DEBUG=10
    20 INFO=20
    30 WARNING
    40 ERROR
    50 CRITICAL

>>> from logging import DEBUG
>>> log("You should not be reading this", verbosity=DEBUG)
>>> from logging import INFO
>>> log("This must have been displayed", verbosity=INFO)
This must have been displayed
>>> set_log_threshold(DEBUG)
>>> log("You should now be reading this", verbosity=DEBUG)
You should now be reading this


Later
~~~~~

We might use
* a nice formatter to the logs,
* a hierarchy of loggers, so as to adjust loglevel for a particular context.
"""


#LICENSE:BSD 3 clause


__all__ = ['log', 'logtools', 'set_log_threshold', 'LOGGERS']


#Importing sys and not sys.stdout because doctest is going to replace
#sys.stdout and we need to access the same fd as doctest
import sys

from logging import CRITICAL, INFO, getLogger, StreamHandler
from operator import isSequenceType


def _non_str_sequence(item):
    """
    Identifies tuples or lists from strings.
    Strings indeed have a special meaning for message formatting :-)
    """
    return isSequenceType(item) and not isinstance(item, (unicode, str))


DEFAULT_THRESHOLD = INFO
LOGGERS = {}


def _getdefaulthandler():
    """
    The default handler is writing to sys.stdout.
    """
    if not hasattr(_getdefaulthandler, 'handler'):
        _getdefaulthandler.handler = StreamHandler(sys.stdout)
    return _getdefaulthandler.handler


def _getlogger(domain):
    """
    gets or makes a logger with default threshold.
    logger.handlers[] will still be empty.
    """
    if domain in LOGGERS:
        return LOGGERS[domain]
    logger = getLogger(domain)
    logger.setLevel(DEFAULT_THRESHOLD)
    LOGGERS[domain] = logger
    return logger


def _log(logger, message, substitutions, verbosity):
    """
    domain: str
    other parameters: see the (generated) log() function

    Performs a late initialization of a log handler if there was none.
    """
    if not logger.handlers:
        logger.addHandler(_getdefaulthandler())

    if verbosity is None:
        verbosity = CRITICAL

    if substitutions is None:
        logger.log(verbosity, message)
        return
    if isinstance(substitutions, dict):
        logger.log(verbosity, message, substitutions)
        return
    if _non_str_sequence(substitutions):
        #pylint: disable=W0142
        #It's not me! It's logger.log signature that wants * magic!
        logger.log(verbosity, message, *substitutions)
        #pylint: enable=W0142
        return
    logger.log(verbosity, message, substitutions)


def logtools(domain):
    """
    Retrieves or creates a logger for the specified domain.

    Parameters:
    domain: str

    returns (logger, log_helper(), threshold_helper())
    """
    logger = _getlogger(domain)

    def log_helper(message, substitutions=None, verbosity=None):
        """
        Outputs a message.

        Parameters:
        -----------
        message: str
            with or without placeholders
        substitutions: filling for placeholders, optional
        verbosity: int, optional (defaults None)
            See this module's doc about log level.
            If unset, the message is critical: always displayed.
        """
        _log(logger, message, substitutions, verbosity)

    def threshold_helper(value):
        """
        Sets global logging threshold for the logger of domain 'scikits.learn'.
        Messages with a verbosity under that threshold will be blocked.
        See 'Criticity' in this module's doc.

        Parameters:
        -----------
        value: int
        """
        logger.setLevel(value)

    return logger, log_helper, threshold_helper


BASE_LOGGER_NAME = "scikits.learn"
#pylint: disable=C0103
#log and set_log_threshold are functions, that's why the function naming
#convention applies here
DEFAULT_LOGGER, log, set_log_threshold = logtools(BASE_LOGGER_NAME)
#pylint: enable=C0103
