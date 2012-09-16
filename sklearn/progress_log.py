"""
A logger with a nested control of verbosity for progress messages.
"""

import sys
import inspect
import logging
import pickle
import copy_reg
# To make users' life easier, import useful symbols
from logging import DEBUG, INFO, ERROR, WARNING, WARN, NOTSET

# XXX: must make logger objects picklable

def get_logger(verbosity=0, name=None, caller_name=None):
    if isinstance(verbosity, ProgressLog):
        # Should we clone and set the caller_name?
        return verbosity
    if name is None:
        # Retrieve the name of the calling module
        frame = inspect.currentframe()
        name = frame.f_back.f_globals['__name__'].split('.')[0]
    logger = ProgressLog(name=name, verbosity=verbosity,
                         caller_name=caller_name)
    # Retrieve the module-level logger and set it as a parent of our
    # class-specific logger to benefit from its settings
    module_log = logging.getLogger(name)
    logger.parent = module_log
    return logger


class HasLog(object):
    """ A class with a getter for a logger.
    """

    def _get_logger(self):
        verbosity = getattr(self, 'verbose', 0)
        name = self.__class__.__module__.split('.')[0]
        caller_name = self.__class__.__name__
        return get_logger(verbosity=verbosity, name=name,
                          caller_name=caller_name)


class ProgressLog(logging.Logger):
    """ A logger with a nested control of verbosity for progress
        messages.
    """

    def __init__(self, name, verbosity=0, level=logging.NOTSET,
                 caller_name=None):
        self.verbosity = verbosity
        self.caller_name = caller_name
        # Standard logging functionality
        super(ProgressLog, self).__init__(name=name, level=level)

    def progress(self, message, msg_vars=(), short_message=None,
                 verbosity_offset=0):
        verbosity_offset += self.verbosity
        if verbosity_offset <= 0:
            return
        caller_name = self.caller_name
        caller_frame = inspect.currentframe().f_back
        if caller_name is None:
            # XXX: following code may be fragile -> try/except?
            if 'self' in caller_frame.f_locals:
                caller_name = caller_frame.f_locals['self'].__class__.__name__
                if verbosity_offset > 10:
                    caller_name = '%s.%s' % (caller_name,
                                            caller_frame.f_code.co_name)
            else:
                caller_name = caller_frame.f_code.co_name
        if self.isEnabledFor(logging.DEBUG):
            caller_name = "%s %s:%i" % (caller_name,
                                        caller_frame.f_code.co_filename,
                                        caller_frame.f_lineno)
        message = '[%s] %s' % (caller_name, message)
        if self.isEnabledFor(logging.DEBUG):
            # High verbosity logging. We pass everything on to the
            # standard logger
            return self.debug(message, *msg_vars)
        # XXX: need to deal with the short_message
        self.info(message, *msg_vars)

    def clone(self, verbosity_offset=-1):
        logger = ProgressLog(name=self.name,
                             verbosity=self.verbosity + verbosity_offset,
                             level=self.level)
        # No need to pass in the level: it is set through the parent
        logger.parent = self
        return logger

    def warning(self, message, msg_vars=()):
        # Does some logic so that warnings.warn is called in addition to
        # logging.
        # XXX: should call warnings.warn and be clever so that the
        # message does not pop up twice
        return super(ProgressLog, self).warning(message, msg_vars)

    #def progress_context(self, verbosity_offset=-1):
    #    pass

    def __repr__(self):
        return '%s(verbosity=%s)' % (self.__class__.__name__, self.verbosity)

    def __reduce__(self):
        # Need to pop the manager at reduce and during reconstruction
        # grab the one from the root logger
        # XXX: what do we do with 'parent'?
        pass


###############################################################################
# Code to make sure that our logging framework always pickles cleanly

class _StdErr(object):
    " No-op tracer for pickling"
    pass


class _StdOut(object):
    " No-op tracer for pickling"
    pass


class FileStreamHandler(logging.StreamHandler):
    """ StreamHandler subclass to deal with pickling problems."""

    def __init__(self, stream=None):
        self._original_stream = stream
        if isinstance(stream, basestring):
            stream = open(stream, 'ab')
        elif (not stream is sys.__stderr__
                        or stream is sys.__stdout__):
            # Just try pickling to raise an error if the stream doesn't
            # pickle
            pickle.dumps(stream)
        if stream is None:
            stream = sys.stdout
        super(FileStreamHandler, self).__init__(stream)

    def __reduce__(self):
        return _file_stream_handler_unpickle, (self._original_stream, )


def _file_stream_handler_unpickle(stream):
    if stream is _StdOut:
        stream = sys.stdout
    if stream is _StdErr:
        stream = sys.stderr
    return FileStreamHandler(stream)

# XXX: the 'logger.manager.logDict' dictionary doesn't pickle right
# The reason is that we don't control what goes in it, and at some point,
# there will be a logger with a handler that doesn't pickle

class PickeableLogger(logging.Logger):
    """ A logger that pickles safely by not saving its
        manager.

        The manager contains a reference to all loggers, and some of
        these may have open file handles that do not pickle.
    """
    def __reduce__(self):
        reconstructor, object_init, state = \
            super(PickeableLogger, self).__reduce__
        state.pop('manager', None)
        return _picklable_logger_unpickle, object_init, state


def _picklable_logger_unpickle(cls, base, state):
    obj = copy_reg._reconstructor(cls, base, state)
    # need to register the logger
    manager = logging.Logger.manager
    logging._acquireLock()
    try:
        manager.loggerDict[obj.name] = obj
    finally:
        manager._releaseLock()

    return obj

# XXX need to be able to set the logger

def get_pickleable_logger(name):
    """ Replacement of logging.getLogger to control the class used in
        order to plug in our class, that is picklable.
    """
    manager = logging.Logger.manager
    logging._acquireLock()
    try:
        if name in manager.loggerDict:
            logger = manager.loggerDict[name]
        else:
            # XXX: Should I use the ProgressLog here or another smaller
            # subclass?
            logger = PickeableLogger(name)
            logger.manager = manager
            manager.loggerDict[name] = logger
            manager._fixupParents(logger)
    finally:
        manager._releaseLock()
    return logger

###############################################################################

def setup_logger(name, level=logging.INFO, log_file=None, dots=True,
                 display_name=False, time_stamp=False,
                 clear_previous_handlers=True):
    """
    Parameters
    ----------

    log_file: string or open file, optional, default: sys.stdout
        we'll log messages and progress there.
        if a string is supplied, it is assumed to be the path where to log
        The file will be created or appended to.

    dots: boolean
        do you want dots printed in this log_file?

    time_stamp: boolean, defaults to None
        do you want logs to be prefixed by a timestamp?
        if unset (None), the value set at object
        initialization (in __init__) is reused
    """
    logger = get_pickleable_logger(name)
    logger.setLevel(level)

    log_format = []
    if time_stamp:
        log_format.append("[%(asctime)s]")
    if display_name:
        log_format.append("[%(name)s]")
    log_format.append("%(message)s")
    log_format = "".join(log_format)

    formatter = logging.Formatter(fmt=log_format)

    if clear_previous_handlers:
        for handler in logger.handlers:
            logger.removeHandler(handler)

    handler = FileStreamHandler(log_file)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


setup_logger('__main__')
setup_logger('sklearn')
