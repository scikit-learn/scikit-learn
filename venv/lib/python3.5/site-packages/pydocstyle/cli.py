"""Command line interface for pydocstyle."""
import logging
import sys

from .utils import log
from .violations import Error
from .config import ConfigurationParser, IllegalConfiguration
from .checker import check


__all__ = ('main', )


class ReturnCode(object):
    no_violations_found = 0
    violations_found = 1
    invalid_options = 2


def run_pydocstyle():
    log.setLevel(logging.DEBUG)
    conf = ConfigurationParser()
    setup_stream_handlers(conf.get_default_run_configuration())

    try:
        conf.parse()
    except IllegalConfiguration:
        return ReturnCode.invalid_options

    run_conf = conf.get_user_run_configuration()

    # Reset the logger according to the command line arguments
    setup_stream_handlers(run_conf)

    log.debug("starting in debug mode.")

    Error.explain = run_conf.explain
    Error.source = run_conf.source

    errors = []
    try:
        for filename, checked_codes, ignore_decorators in \
                conf.get_files_to_check():
            errors.extend(check((filename,), select=checked_codes,
                                ignore_decorators=ignore_decorators))
    except IllegalConfiguration as error:
        # An illegal configuration file was found during file generation.
        log.error(error.args[0])
        return ReturnCode.invalid_options

    count = 0
    for error in errors:
        if hasattr(error, 'code'):
            sys.stdout.write('%s\n' % error)
        count += 1
    if count == 0:
        exit_code = ReturnCode.no_violations_found
    else:
        exit_code = ReturnCode.violations_found
    if run_conf.count:
        print(count)
    return exit_code


def main():
    """Run pydocstyle as a script."""
    try:
        sys.exit(run_pydocstyle())
    except KeyboardInterrupt:
        pass


def setup_stream_handlers(conf):
    """Setup logging stream handlers according to the options."""
    class StdoutFilter(logging.Filter):
        def filter(self, record):
            return record.levelno in (logging.DEBUG, logging.INFO)

    log.handlers = []

    stdout_handler = logging.StreamHandler(sys.stdout)
    stdout_handler.setLevel(logging.WARNING)
    stdout_handler.addFilter(StdoutFilter())
    if conf.debug:
        stdout_handler.setLevel(logging.DEBUG)
    elif conf.verbose:
        stdout_handler.setLevel(logging.INFO)
    else:
        stdout_handler.setLevel(logging.WARNING)
    log.addHandler(stdout_handler)

    stderr_handler = logging.StreamHandler(sys.stderr)
    msg_format = "%(levelname)s: %(message)s"
    stderr_handler.setFormatter(logging.Formatter(fmt=msg_format))
    stderr_handler.setLevel(logging.WARNING)
    log.addHandler(stderr_handler)
