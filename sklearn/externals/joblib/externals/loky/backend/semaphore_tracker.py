###############################################################################
# Server process to keep track of unlinked semaphores and clean them.
#
# author: Thomas Moreau
#
# adapted from multiprocessing/semaphore_tracker.py  (17/02/2017)
#  * include custom spawnv_passfds to start the process
#  * use custom unlink from our own SemLock implementation
#  * add some VERBOSE logging
#

#
# On Unix we run a server process which keeps track of unlinked
# semaphores. The server ignores SIGINT and SIGTERM and reads from a
# pipe.  Every other process of the program has a copy of the writable
# end of the pipe, so we get EOF when all other processes have exited.
# Then the server process unlinks any remaining semaphore names.
#
# This is important because the system only supports a limited number
# of named semaphores, and they will not be automatically removed till
# the next reboot.  Without this semaphore tracker process, "killall
# python" would probably leave unlinked semaphores.
#

import os
import signal
import sys
import threading
import warnings

from . import spawn
from multiprocessing import util

try:
    from _multiprocessing import sem_unlink
except ImportError:
    from .semlock import sem_unlink

__all__ = ['ensure_running', 'register', 'unregister']

VERBOSE = False


class SemaphoreTracker(object):

    def __init__(self):
        self._lock = threading.Lock()
        self._fd = None

    def getfd(self):
        self.ensure_running()
        return self._fd

    def ensure_running(self):
        '''Make sure that semaphore tracker process is running.

        This can be run from any process.  Usually a child process will use
        the semaphore created by its parent.'''
        with self._lock:
            if self._fd is not None:
                return
            fds_to_pass = []
            try:
                fds_to_pass.append(sys.stderr.fileno())
            except Exception:
                pass

            cmd = 'from {} import main; main(%d)'.format(main.__module__)
            r, w = os.pipe()
            try:
                fds_to_pass.append(r)
                # process will out live us, so no need to wait on pid
                exe = spawn.get_executable()
                args = [exe] + util._args_from_interpreter_flags()
                # In python 3.3, there is a bug which put `-RRRRR..` instead of
                # `-R` in args. Replace it to get the correct flags.
                # See https://github.com/python/cpython/blob/3.3/Lib/subprocess.py#L488
                if sys.version_info[:2] <= (3, 3):
                    import re
                    for i in range(1, len(args)):
                        args[i] = re.sub("-R+", "-R", args[i])
                args += ['-c', cmd % r]
                util.debug("launching Semaphore tracker: {}".format(args))
                spawnv_passfds(exe, args, fds_to_pass)
            except:
                os.close(w)
                raise
            else:
                self._fd = w
            finally:
                os.close(r)

    def register(self, name):
        '''Register name of semaphore with semaphore tracker.'''
        self._send('REGISTER', name)

    def unregister(self, name):
        '''Unregister name of semaphore with semaphore tracker.'''
        self._send('UNREGISTER', name)

    def _send(self, cmd, name):
        self.ensure_running()
        msg = '{0}:{1}\n'.format(cmd, name).encode('ascii')
        if len(name) > 512:
            # posix guarantees that writes to a pipe of less than PIPE_BUF
            # bytes are atomic, and that PIPE_BUF >= 512
            raise ValueError('name too long')
        nbytes = os.write(self._fd, msg)
        assert nbytes == len(msg)


_semaphore_tracker = SemaphoreTracker()
ensure_running = _semaphore_tracker.ensure_running
register = _semaphore_tracker.register
unregister = _semaphore_tracker.unregister
getfd = _semaphore_tracker.getfd


def main(fd):
    '''Run semaphore tracker.'''
    # protect the process from ^C and "killall python" etc
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_IGN)

    for f in (sys.stdin, sys.stdout):
        try:
            f.close()
        except Exception:
            pass

    if VERBOSE:  # pragma: no cover
        sys.stderr.write("Main semaphore tracker is running\n")
        sys.stderr.flush()

    cache = set()
    try:
        # keep track of registered/unregistered semaphores
        with os.fdopen(fd, 'rb') as f:
            for line in f:
                try:
                    cmd, name = line.strip().split(b':')
                    if cmd == b'REGISTER':
                        cache.add(name)
                        if VERBOSE:  # pragma: no cover
                            name = name.decode('ascii')
                            sys.stderr.write("[SemaphoreTracker] register {}\n"
                                             .format(name))
                            sys.stderr.flush()
                    elif cmd == b'UNREGISTER':
                        cache.remove(name)
                        if VERBOSE:  # pragma: no cover
                            name = name.decode('ascii')
                            sys.stderr.write("[SemaphoreTracker] unregister {}"
                                             ": cache({})\n"
                                             .format(name, len(cache)))
                            sys.stderr.flush()
                    else:
                        raise RuntimeError('unrecognized command %r' % cmd)
                except Exception:
                    try:
                        sys.excepthook(*sys.exc_info())
                    except:
                        pass
    finally:
        # all processes have terminated; cleanup any remaining semaphores
        if cache:
            try:
                warnings.warn('semaphore_tracker: There appear to be %d '
                              'leaked semaphores to clean up at shutdown' %
                              len(cache))
            except Exception:
                pass
        for name in cache:
            # For some reason the process which created and registered this
            # semaphore has failed to unregister it. Presumably it has died.
            # We therefore unlink it.
            try:
                try:
                    sem_unlink(name)
                    if VERBOSE:  # pragma: no cover
                        name = name.decode('ascii')
                        sys.stderr.write("[SemaphoreTracker] unregister {}\n"
                                         .format(name))
                        sys.stderr.flush()
                except Exception as e:
                    warnings.warn('semaphore_tracker: %r: %s' % (name, e))
            finally:
                pass

    if VERBOSE:  # pragma: no cover
        sys.stderr.write("semaphore tracker shut down\n")
        sys.stderr.flush()


#
# Start a program with only specified fds kept open
#

def spawnv_passfds(path, args, passfds):
    passfds = sorted(passfds)
    errpipe_read, errpipe_write = os.pipe()
    try:
        if sys.version_info >= (3, 3):
            import _posixsubprocess
            return _posixsubprocess.fork_exec(
                args, [os.fsencode(path)], True, passfds, None, None,
                -1, -1, -1, -1, -1, -1, errpipe_read, errpipe_write,
                False, False, None)
        else:
            from .reduction import _mk_inheritable
            _pass = []
            for fd in passfds:
                _pass += [_mk_inheritable(fd)]
            from .fork_exec import fork_exec
            fork_exec(args, _pass)
    finally:
        os.close(errpipe_read)
        os.close(errpipe_write)
