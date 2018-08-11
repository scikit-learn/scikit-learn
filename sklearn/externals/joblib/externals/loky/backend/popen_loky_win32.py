import os
import sys

from .context import get_spawning_popen, set_spawning_popen
from . import spawn
from . import reduction
from multiprocessing import util

if sys.platform == "win32":
    # Avoid import error by code introspection tools such as test runners
    # trying to import this module while running on non-Windows systems.
    import msvcrt
    from .compat_win32 import _winapi
    from .compat_win32 import Popen as _Popen
else:
    _Popen = object

if sys.version_info[:2] < (3, 3):
    from os import fdopen as open

__all__ = ['Popen']

#
#
#

TERMINATE = 0x10000
WINEXE = (sys.platform == 'win32' and getattr(sys, 'frozen', False))
WINSERVICE = sys.executable.lower().endswith("pythonservice.exe")


#
# We define a Popen class similar to the one from subprocess, but
# whose constructor takes a process object as its argument.
#

class Popen(_Popen):
    '''
    Start a subprocess to run the code of a process object
    '''
    method = 'loky'

    def __init__(self, process_obj):
        prep_data = spawn.get_preparation_data(
            process_obj._name, process_obj.init_main_module)

        # read end of pipe will be "stolen" by the child process
        # -- see spawn_main() in spawn.py.
        rhandle, wfd = _winapi.CreatePipe(None, 0)
        if sys.version_info[:2] > (3, 3):
            wfd = msvcrt.open_osfhandle(wfd, 0)

        cmd = spawn.get_command_line(parent_pid=os.getpid(),
                                     pipe_handle=rhandle)
        cmd = ' '.join('"%s"' % x for x in cmd)

        try:
            with open(wfd, 'wb') as to_child:
                # start process
                try:
                    inherit = sys.version_info[:2] < (3, 4)
                    hp, ht, pid, tid = _winapi.CreateProcess(
                        spawn.get_executable(), cmd,
                        None, None, inherit, 0,
                        None, None, None)
                    _winapi.CloseHandle(ht)
                except:
                    _winapi.CloseHandle(rhandle)
                    raise

                # set attributes of self
                self.pid = pid
                self.returncode = None
                self._handle = hp
                self.sentinel = int(hp)
                util.Finalize(self, _winapi.CloseHandle, (self.sentinel,))

                # send information to child
                set_spawning_popen(self)
                if sys.version_info[:2] < (3, 4):
                    Popen._tls.process_handle = int(hp)
                try:
                    reduction.dump(prep_data, to_child)
                    reduction.dump(process_obj, to_child)
                finally:
                    set_spawning_popen(None)
                    if sys.version_info[:2] < (3, 4):
                        del Popen._tls.process_handle
        except IOError as exc:
            # IOError 22 happens when the launched subprocess terminated before
            # wfd.close is called. Thus we can safely ignore it.
            if exc.errno != 22:
                raise
            util.debug("While starting {}, ignored a IOError 22"
                       .format(process_obj._name))

    def duplicate_for_child(self, handle):
        assert self is get_spawning_popen()
        return reduction.duplicate(handle, self.sentinel)
