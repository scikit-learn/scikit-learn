"""
A logger with a nested control of verbosity for progress messages.
"""

import sys
import inspect
import logging
# To make users' life easier, import useful symbols
from logging import DEBUG, INFO, ERROR, WARNING, WARN, NOTSET

# XXX: timer, for use e.g. in kmeans?
def get_logger(verbosity=0, name=None, caller_name=None):
    if isinstance(verbosity, ProgressLogger):
        # Should we clone and set the caller_name?
        return verbosity
    if name is None:
        # Retrieve the name of the calling module
        frame = inspect.currentframe()
        name = frame.f_back.f_globals['__name__'].split('.')[0]
    logger = ProgressLogger(logger=name, verbosity=verbosity,
                         caller_name=caller_name)
    # Retrieve the module-level logger and set it as a parent of our
    # class-specific logger to benefit from its settings
    module_log = logging.getLogger(name)
    logger.parent = module_log
    return logger


class HasLogger(object):
    """ A class with a getter for a logger.
    """

    def _get_logger(self):
        verbosity = getattr(self, 'verbose', 0)
        name = self.__class__.__module__.split('.')[0]
        caller_name = self.__class__.__name__
        return get_logger(verbosity=verbosity, name=name,
                          caller_name=caller_name)


class ProgressLogger(object):
    """ A logger with a nested control of verbosity for progress
        messages.
    """

    def __init__(self, logger, verbosity=0, caller_name=None):
        self.verbosity = verbosity
        self.caller_name = caller_name
        # Standard logging functionality
        if isinstance(logger, basestring):
            logger = logging.getLogger(name=logger)
        self.logger = logger
        self.debug = logger.debug
        self.info = logger.info
        self.error = logger.error

    def progress(self, message, *msg_vars, **options):
        verbosity_offset = options.pop('verbosity_offset', 0)
        short_message = options.pop('short_message', None)
        if options:
            raise TypeError('Unknown arguments %s' % options.keys())
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
        if self.logger.isEnabledFor(logging.DEBUG):
            caller_name = "%s %s:%i" % (caller_name,
                                        caller_frame.f_code.co_filename,
                                        caller_frame.f_lineno)
        message = '[%s] %s' % (caller_name, message)
        if self.logger.isEnabledFor(logging.DEBUG):
            # High verbosity logging. We pass everything on to the
            # standard logger
            return self.debug(message, *msg_vars)
        # XXX: need to deal with the short_message
        self.info(message, *msg_vars)

    def clone(self, verbosity_offset=-1):
        logger = self.__class__(logger=self.logger,
                             verbosity=self.verbosity + verbosity_offset)
        return logger

    def warning(self, message, *msg_vars):
        # Does some logic so that warnings.warn is called in addition to
        # logging.
        # XXX: should call warnings.warn and be clever so that the
        # message does not pop up twice
        return self.logger.warning(message, *msg_vars)

    #def progress_context(self, verbosity_offset=-1):
    #    pass

    def __repr__(self):
        return '%s(verbosity=%s)' % (self.__class__.__name__, self.verbosity)

    def __reduce__(self):
        """ We don't store the logger when pickling, just its name, as
            loggers don't pickle.
        """
        return (self.__class__, (self.logger.name, self.verbosity,
                self.caller_name))

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
    logger = logging.getLogger(name)
    logger.setLevel(level)

    log_format = []
    if time_stamp:
        log_format.append("[%(asctime)s]")
    if display_name:
        log_format.append("[%(name)s]")
    log_format.append("%(message)s")
    log_format = "".join(log_format)

    formatter = logging.Formatter(fmt=log_format)

    if log_file is None:
        log_file = sys.stdout
    elif isinstance(log_file, basestring):
        log_file = file(log_file, 'ab')

    if clear_previous_handlers:
        for handler in logger.handlers:
            logger.removeHandler(handler)

    handler = logging.StreamHandler(log_file)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


# Setup loggers with sane defaults
setup_logger('sklearn')
if not '__main__' in logging.Logger.manager.loggerDict:
    # Setup the __main__ logger, because if people define classes
    # in a script, this is what is used
    setup_logger('__main__')

