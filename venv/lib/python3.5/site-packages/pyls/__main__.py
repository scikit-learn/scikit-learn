# Copyright 2017 Palantir Technologies, Inc.
import argparse
import json
import logging
import logging.config
import sys
from .python_ls import start_io_lang_server, start_tcp_lang_server, PythonLanguageServer

LOG_FORMAT = "%(asctime)s UTC - %(levelname)s - %(name)s - %(message)s"


def add_arguments(parser):
    parser.description = "Python Language Server"

    parser.add_argument(
        "--tcp", action="store_true",
        help="Use TCP server instead of stdio"
    )
    parser.add_argument(
        "--host", default="127.0.0.1",
        help="Bind to this address"
    )
    parser.add_argument(
        "--port", type=int, default=2087,
        help="Bind to this port"
    )

    log_group = parser.add_mutually_exclusive_group()
    log_group.add_argument(
        "--log-config",
        help="Path to a JSON file containing Python logging config."
    )
    log_group.add_argument(
        "--log-file",
        help="Redirect logs to the given file instead of writing to stderr."
        "Has no effect if used with --log-config."
    )

    parser.add_argument(
        '-v', '--verbose', action='count', default=0,
        help="Increase verbosity of log output, overrides log config file"
    )


def main():
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args = parser.parse_args()
    _configure_logger(args.verbose, args.log_config, args.log_file)

    if args.tcp:
        start_tcp_lang_server(args.host, args.port, PythonLanguageServer)
    else:
        stdin, stdout = _binary_stdio()
        start_io_lang_server(stdin, stdout, PythonLanguageServer)


def _binary_stdio():
    """Construct binary stdio streams (not text mode).

    This seems to be different for Window/Unix Python2/3, so going by:
        https://stackoverflow.com/questions/2850893/reading-binary-data-from-stdin
    """
    PY3K = sys.version_info >= (3, 0)

    if PY3K:
        # pylint: disable=no-member
        stdin, stdout = sys.stdin.buffer, sys.stdout.buffer
    else:
        # Python 2 on Windows opens sys.stdin in text mode, and
        # binary data that read from it becomes corrupted on \r\n
        if sys.platform == "win32":
            # set sys.stdin to binary mode
            # pylint: disable=no-member,import-error
            import os
            import msvcrt
            msvcrt.setmode(sys.stdin.fileno(), os.O_BINARY)
            msvcrt.setmode(sys.stdout.fileno(), os.O_BINARY)
        stdin, stdout = sys.stdin, sys.stdout

    return stdin, stdout


def _configure_logger(verbose=0, log_config=None, log_file=None):
    root_logger = logging.root

    if log_config:
        with open(log_config, 'r') as f:
            logging.config.dictConfig(json.load(f))
    else:
        formatter = logging.Formatter(LOG_FORMAT)
        if log_file:
            log_handler = logging.handlers.RotatingFileHandler(
                log_file, mode='a', maxBytes=50*1024*1024,
                backupCount=10, encoding=None, delay=0
            )
        else:
            log_handler = logging.StreamHandler()
        log_handler.setFormatter(formatter)
        root_logger.addHandler(log_handler)

    if verbose == 0:
        level = logging.WARNING
    elif verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG

    root_logger.setLevel(level)


if __name__ == '__main__':
    main()
