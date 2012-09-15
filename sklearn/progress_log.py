
# The verbosity offset logic is enforced in each call to the logger

# We want to try to use the 'parent' of the logger as much as
# possible.
# Our _get_logger method should have a name that defaults to the top
# module in the class module path (ie sklearn for
# sklearn.lda.LDA)

# We need to remember to define logging rules somewhere for each module

import sys
import inspect
import logging
# To make users' life easier, import useful symbols
from logging import DEBUG, INFO, ERROR, WARNING, WARN, NOTSET

class HasLog(object):

    def _get_logger(self):
        if isinstance(self.verbose, ProgressLog):
            return self.verbose
        name = self.__class__.__module__.split('.')[0]
        log = ProgressLog(name=name,
                          verbosity_offset=self.verbose)
        # This is the module-level logger. We set it as a parent of our
        # class-specific logger to benefit from its settings
        module_log = logging.getLogger(name)
        log.parent = module_log
        return log


class ProgressLog(logging.Logger):

    def __init__(self, name, verbosity_offset=0,
                 level=logging.NOTSET):
        self.verbosity_offset = verbosity_offset
        # Standard logging functionality
        super(ProgressLog, self).__init__(name=name, level=level)

    def progress(self, message, msg_vars=(), short_message=None,
                 verbosity_offset=0):
        verbosity_offset += self.verbosity_offset
        if verbosity_offset <= 0:
            return
        # XXX: following code may be fragile -> try/except?
        caller_frame = inspect.currentframe().f_back
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
                verbosity_offset=self.verbosity_offset + verbosity_offset)
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


def setup_logger(name, level=logging.INFO, log_file=None, dots=True,
                 display_name=False, timestamp=False,
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

    timestamp: boolean, defaults to None
        do you want logs to be prefixed by a timestamp?
        if unset (None), the value set at object
        initialization (in __init__) is reused
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    if log_file is None:
        log_file = sys.stdout
    elif isinstance(log_file, basestring):
        log_file = open(log_file, 'ab')

    log_format = []
    if timestamp:
        log_format.append("[%(asctime)s]")
    if display_name:
        log_format.append("[%(name)s]")
    log_format.append("%(message)s")
    log_format = "".join(log_format)

    formatter = logging.Formatter(fmt=log_format)

    if clear_previous_handlers:
        for handler in logger.handlers:
            logger.removeHandler(handler)

    handler = logging.StreamHandler(log_file)
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    return logger


setup_logger('__main__')
