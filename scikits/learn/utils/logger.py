"""
Logging utility.

Author: Feth Arezki <feth A+ tuttu.info>

This module exposes 'log()' in an attempt to provide a drop in replacement
to this use of print statements throughout the code.

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
BASE_LOGGER_NAME = "scikits.learn"


def _defaultlogger():
    """
    Gets or makes a default logger for scikits.learn to use.
    """
    if not hasattr(_defaultlogger, 'logger'):
        logger = getLogger(BASE_LOGGER_NAME)
        logger.addHandler(StreamHandler(sys.stdout))
        logger.setLevel(DEFAULT_THRESHOLD)
        _defaultlogger.logger = logger
    return _defaultlogger.logger


def set_log_threshold(value):
    """
    Sets global logging threshold for the logger of domain 'scikits.learn'.
    Messages with a verbosity under that threshold will be blocked.
    """
    _defaultlogger().setLevel(value)


def log(message, substitutions=None, verbosity=None):
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
    logger = _defaultlogger()
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
