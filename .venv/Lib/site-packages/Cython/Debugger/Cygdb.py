#!/usr/bin/env python3

"""
The Cython debugger

The current directory should contain a directory named 'cython_debug', or a
path to the cython project directory should be given (the parent directory of
cython_debug).

Additional gdb args can be provided only if a path to the project directory is
given.
"""

import os
import sys
import glob
import tempfile
import textwrap
import subprocess
import argparse
import logging

logger = logging.getLogger(__name__)


def make_command_file(path_to_debug_info, prefix_code='',
                      no_import=False, skip_interpreter=False):
    if not no_import:
        pattern = os.path.join(path_to_debug_info,
                               'cython_debug',
                               'cython_debug_info_*')
        debug_files = glob.glob(pattern)

        if not debug_files:
            sys.exit('No cython debug files were found in %s. Aborting.' % (
                                   os.path.abspath(path_to_debug_info)))

    fd, tempfilename = tempfile.mkstemp()
    f = os.fdopen(fd, 'w')
    try:
        f.write(prefix_code)
        f.write(textwrap.dedent('''\
            # This is a gdb command file
            # See https://sourceware.org/gdb/onlinedocs/gdb/Command-Files.html

            set breakpoint pending on
            set print pretty on

            python
            try:
                # Activate virtualenv, if we were launched from one
                import os
                import sys
                virtualenv = os.getenv('VIRTUAL_ENV')
                if virtualenv:
                    scripts_dir = 'Scripts' if sys.platform == "win32" else 'bin'
                    path_to_activate_this_py = os.path.join(virtualenv, scripts_dir, 'activate_this.py')
                    print("gdb command file: Activating virtualenv: %s; path_to_activate_this_py: %s" % (
                        virtualenv, path_to_activate_this_py))
                    with open(path_to_activate_this_py) as f:
                        exec(f.read(), dict(__file__=path_to_activate_this_py))
                from Cython.Debugger import libcython, libpython
            except Exception as ex:
                from traceback import print_exc
                print("There was an error in Python code originating from the file " + ''' + repr(__file__) + ''')
                print("It used the Python interpreter " + str(sys.executable))
                print_exc()
                exit(1)
            end
            '''))

        if no_import:
            # don't do this, this overrides file command in .gdbinit
            # f.write("file %s\n" % sys.executable)
            pass
        else:
            if not skip_interpreter:
                # Point Cygdb to the interpreter that was used to generate
                # the debugging information.
                path = os.path.join(path_to_debug_info, "cython_debug", "interpreter")
                interpreter_file = open(path)
                try:
                    interpreter = interpreter_file.read()
                finally:
                    interpreter_file.close()
                f.write("file %s\n" % interpreter)

            f.write('\n'.join('cy import %s\n' % fn for fn in debug_files))

            if not skip_interpreter:
                f.write(textwrap.dedent('''\
                    python
                    import sys
                    # Check if the Python executable provides a symbol table.
                    if not hasattr(gdb.selected_inferior().progspace, "symbol_file"):
                        sys.stderr.write(
                            "''' + interpreter + ''' was not compiled with debug symbols (or it was "
                            "stripped). Some functionality may not work (properly).\\n")
                    end
                '''))

            f.write("source .cygdbinit\n")
    finally:
        f.close()

    return tempfilename


def main():
    """
    Start the Cython debugger. This tells gdb to import the Cython and Python
    extensions (libcython.py and libpython.py) and it enables gdb's pending
    breakpoints.
    """
    parser = argparse.ArgumentParser(
        prog="cygdb",
        description="Cython debugger",
    )
    parser.add_argument("gdb_argv", nargs="*",
                        help="Arguments to forward to gdb; specified after --")
    parser.add_argument("--build-dir", dest="build_dir", default=None,
                        help="Directory containing cython_build/ files")
    parser.add_argument("--gdb-executable",
        dest="gdb", default='gdb',
        help="gdb executable to use [default: gdb]")
    parser.add_argument("--verbose", "-v",
        dest="verbosity", action="count", default=0,
        help="Verbose mode. Multiple -v options increase the verbosity")
    parser.add_argument("--skip-interpreter",
                      dest="skip_interpreter", default=False, action="store_true",
                      help="Do not automatically point GDB to the same interpreter "
                           "used to generate debugging information")

    options = parser.parse_args()
    path_to_debug_info = options.build_dir
    gdb_argv = options.gdb_argv
    no_import = path_to_debug_info is None

    if options.build_dir is None and gdb_argv and os.path.isdir(gdb_argv[0]):
        import warnings
        gdb_argv = options.gdb_argv[1:]
        path_to_debug_info = options.gdb_argv[0]
        warnings.warn(f'Using deprecated positional parameter to find build directory. Use "--build-dir {path_to_debug_info}" argument instead.')

    logging_level = logging.WARN
    if options.verbosity == 1:
        logging_level = logging.INFO
    if options.verbosity >= 2:
        logging_level = logging.DEBUG
    logging.basicConfig(level=logging_level)

    skip_interpreter = options.skip_interpreter

    logger.debug("options = %r", options)
    tempfilename = make_command_file(path_to_debug_info,
                                     no_import=no_import,
                                     skip_interpreter=skip_interpreter)
    logger.info("Launching %s with command file: %s and gdb_argv: %s",
        options.gdb, tempfilename, gdb_argv)
    with open(tempfilename) as tempfile:
        logger.debug('Command file (%s) contains: """\n%s"""', tempfilename, tempfile.read())
        logger.info("Spawning %s...", options.gdb)
        p = subprocess.Popen([options.gdb, '-command', tempfilename] + gdb_argv)
        logger.info("Spawned %s (pid %d)", options.gdb, p.pid)
        while True:
            try:
                logger.debug("Waiting for gdb (pid %d) to exit...", p.pid)
                ret = p.wait()
                logger.debug("Wait for gdb (pid %d) to exit is done. Returned: %r", p.pid, ret)
            except KeyboardInterrupt:
                pass
            else:
                break
        logger.debug("Closing temp command file with fd: %s", tempfile.fileno())
    logger.debug("Removing temp command file: %s", tempfilename)
    os.remove(tempfilename)
    logger.debug("Removed temp command file: %s", tempfilename)
