# Copyright (c) 2005-2022, NumPy Developers.
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:

#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.

#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.

#     * Neither the name of the NumPy Developers nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
import os
import sys
import subprocess
import re
from distutils.errors import DistutilsExecError

from numpy.distutils import log


def is_sequence(seq):
    if isinstance(seq, str):
        return False
    try:
        len(seq)
    except Exception:
        return False
    return True


def forward_bytes_to_stdout(val):
    """
    Forward bytes from a subprocess call to the console, without attempting to
    decode them.

    The assumption is that the subprocess call already returned bytes in
    a suitable encoding.
    """
    if hasattr(sys.stdout, "buffer"):
        # use the underlying binary output if there is one
        sys.stdout.buffer.write(val)
    elif hasattr(sys.stdout, "encoding"):
        # round-trip the encoding if necessary
        sys.stdout.write(val.decode(sys.stdout.encoding))
    else:
        # make a best-guess at the encoding
        sys.stdout.write(val.decode("utf8", errors="replace"))


def CCompiler_spawn(self, cmd, display=None, env=None):
    """
    Execute a command in a sub-process.

    Parameters
    ----------
    cmd : str
        The command to execute.
    display : str or sequence of str, optional
        The text to add to the log file kept by `numpy.distutils`.
        If not given, `display` is equal to `cmd`.
    env: a dictionary for environment variables, optional

    Returns
    -------
    None

    Raises
    ------
    DistutilsExecError
        If the command failed, i.e. the exit status was not 0.

    """
    env = env if env is not None else dict(os.environ)
    if display is None:
        display = cmd
        if is_sequence(display):
            display = " ".join(list(display))
    log.info(display)
    try:
        if self.verbose:
            subprocess.check_output(cmd, env=env)
        else:
            subprocess.check_output(cmd, stderr=subprocess.STDOUT, env=env)
    except subprocess.CalledProcessError as exc:
        o = exc.output
        s = exc.returncode
    except OSError as e:
        # OSError doesn't have the same hooks for the exception
        # output, but exec_command() historically would use an
        # empty string for EnvironmentError (base class for
        # OSError)
        # o = b''
        # still that would make the end-user lost in translation!
        o = f"\n\n{e}\n\n\n"
        try:
            o = o.encode(sys.stdout.encoding)
        except AttributeError:
            o = o.encode("utf8")
        # status previously used by exec_command() for parent
        # of OSError
        s = 127
    else:
        # use a convenience return here so that any kind of
        # caught exception will execute the default code after the
        # try / except block, which handles various exceptions
        return None

    if is_sequence(cmd):
        cmd = " ".join(list(cmd))

    if self.verbose:
        forward_bytes_to_stdout(o)

    if re.search(b"Too many open files", o):
        msg = "\nTry rerunning setup command until build succeeds."
    else:
        msg = ""
    raise DistutilsExecError(
        'Command "%s" failed with exit status %d%s' % (cmd, s, msg)
    )
