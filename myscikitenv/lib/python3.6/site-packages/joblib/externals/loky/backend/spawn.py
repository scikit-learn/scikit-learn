###############################################################################
# Prepares and processes the data to setup the new process environment
#
# author: Thomas Moreau and Olivier Grisel
#
# adapted from multiprocessing/spawn.py (17/02/2017)
#  * Improve logging data
#
import os
import sys
import runpy
import types
from multiprocessing import process, util

from joblib.externals.loky.backend import context


if sys.platform != 'win32':
    WINEXE = False
    WINSERVICE = False
else:
    WINEXE = (sys.platform == 'win32' and getattr(sys, 'frozen', False))
    WINSERVICE = sys.executable.lower().endswith("pythonservice.exe")

if WINSERVICE:
    _python_exe = os.path.join(sys.exec_prefix, 'python.exe')
else:
    _python_exe = sys.executable


def get_executable():
    return _python_exe


def _check_not_importing_main():
    if getattr(process.current_process(), '_inheriting', False):
        raise RuntimeError('''
        An attempt has been made to start a new process before the
        current process has finished its bootstrapping phase.

        This probably means that you are not using fork to start your
        child processes and you have forgotten to use the proper idiom
        in the main module:

            if __name__ == '__main__':
                freeze_support()
                ...

        The "freeze_support()" line can be omitted if the program
        is not going to be frozen to produce an executable.''')


def get_preparation_data(name, init_main_module=True):
    '''
    Return info about parent needed by child to unpickle process object
    '''
    _check_not_importing_main()
    d = dict(
        log_to_stderr=util._log_to_stderr,
        authkey=bytes(process.current_process().authkey),
    )

    if util._logger is not None:
        d['log_level'] = util._logger.getEffectiveLevel()
        if len(util._logger.handlers) > 0:
            h = util._logger.handlers[0]
            d['log_fmt'] = h.formatter._fmt

    sys_path = [p for p in sys.path]
    try:
        i = sys_path.index('')
    except ValueError:
        pass
    else:
        sys_path[i] = process.ORIGINAL_DIR

    d.update(
        name=name,
        sys_path=sys_path,
        sys_argv=sys.argv,
        orig_dir=process.ORIGINAL_DIR,
        dir=os.getcwd()
    )

    if sys.platform != "win32":
        # Pass the semaphore_tracker pid to avoid re-spawning it in every child
        from . import semaphore_tracker
        semaphore_tracker.ensure_running()
        d['tracker_pid'] = semaphore_tracker._semaphore_tracker._pid

    # Figure out whether to initialise main in the subprocess as a module
    # or through direct execution (or to leave it alone entirely)
    if init_main_module:
        main_module = sys.modules['__main__']
        try:
            main_mod_name = getattr(main_module.__spec__, "name", None)
        except BaseException:
            main_mod_name = None
        if main_mod_name is not None:
            d['init_main_from_name'] = main_mod_name
        elif sys.platform != 'win32' or (not WINEXE and not WINSERVICE):
            main_path = getattr(main_module, '__file__', None)
            if main_path is not None:
                if (not os.path.isabs(main_path) and
                        process.ORIGINAL_DIR is not None):
                    main_path = os.path.join(process.ORIGINAL_DIR, main_path)
                d['init_main_from_path'] = os.path.normpath(main_path)
                # Compat for python2.7
                d['main_path'] = d['init_main_from_path']

    return d


#
# Prepare current process
#
old_main_modules = []


def prepare(data):
    '''
    Try to get current process ready to unpickle process object
    '''
    if 'name' in data:
        process.current_process().name = data['name']

    if 'authkey' in data:
        process.current_process().authkey = data['authkey']

    if 'log_to_stderr' in data and data['log_to_stderr']:
        util.log_to_stderr()

    if 'log_level' in data:
        util.get_logger().setLevel(data['log_level'])

    if 'log_fmt' in data:
        import logging
        util.get_logger().handlers[0].setFormatter(
            logging.Formatter(data['log_fmt'])
        )

    if 'sys_path' in data:
        sys.path = data['sys_path']

    if 'sys_argv' in data:
        sys.argv = data['sys_argv']

    if 'dir' in data:
        os.chdir(data['dir'])

    if 'orig_dir' in data:
        process.ORIGINAL_DIR = data['orig_dir']

    if 'tacker_pid' in data:
        from . import semaphore_tracker
        semaphore_tracker._semaphore_tracker._pid = data["tracker_pid"]

    if 'init_main_from_name' in data:
        _fixup_main_from_name(data['init_main_from_name'])
    elif 'init_main_from_path' in data:
        _fixup_main_from_path(data['init_main_from_path'])


# Multiprocessing module helpers to fix up the main module in
# spawned subprocesses
def _fixup_main_from_name(mod_name):
    # __main__.py files for packages, directories, zip archives, etc, run
    # their "main only" code unconditionally, so we don't even try to
    # populate anything in __main__, nor do we make any changes to
    # __main__ attributes
    current_main = sys.modules['__main__']
    if mod_name == "__main__" or mod_name.endswith(".__main__"):
        return

    # If this process was forked, __main__ may already be populated
    if getattr(current_main.__spec__, "name", None) == mod_name:
        return

    # Otherwise, __main__ may contain some non-main code where we need to
    # support unpickling it properly. We rerun it as __mp_main__ and make
    # the normal __main__ an alias to that
    old_main_modules.append(current_main)
    main_module = types.ModuleType("__mp_main__")
    main_content = runpy.run_module(mod_name,
                                    run_name="__mp_main__",
                                    alter_sys=True)
    main_module.__dict__.update(main_content)
    sys.modules['__main__'] = sys.modules['__mp_main__'] = main_module


def _fixup_main_from_path(main_path):
    # If this process was forked, __main__ may already be populated
    current_main = sys.modules['__main__']

    # Unfortunately, the main ipython launch script historically had no
    # "if __name__ == '__main__'" guard, so we work around that
    # by treating it like a __main__.py file
    # See https://github.com/ipython/ipython/issues/4698
    main_name = os.path.splitext(os.path.basename(main_path))[0]
    if main_name == 'ipython':
        return

    # Otherwise, if __file__ already has the setting we expect,
    # there's nothing more to do
    if getattr(current_main, '__file__', None) == main_path:
        return

    # If the parent process has sent a path through rather than a module
    # name we assume it is an executable script that may contain
    # non-main code that needs to be executed
    old_main_modules.append(current_main)
    main_module = types.ModuleType("__mp_main__")
    main_content = runpy.run_path(main_path,
                                  run_name="__mp_main__")
    main_module.__dict__.update(main_content)
    sys.modules['__main__'] = sys.modules['__mp_main__'] = main_module


def import_main_path(main_path):
    '''
    Set sys.modules['__main__'] to module at main_path
    '''
    _fixup_main_from_path(main_path)
