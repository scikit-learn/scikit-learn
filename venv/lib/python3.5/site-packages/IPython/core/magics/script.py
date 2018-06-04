"""Magic functions for running cells in various scripts."""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import errno
import os
import sys
import signal
import time
from subprocess import Popen, PIPE
import atexit

from IPython.core import magic_arguments
from IPython.core.magic import  (
    Magics, magics_class, line_magic, cell_magic
)
from IPython.lib.backgroundjobs import BackgroundJobManager
from IPython.utils import py3compat
from IPython.utils.process import arg_split
from traitlets import List, Dict, default

#-----------------------------------------------------------------------------
# Magic implementation classes
#-----------------------------------------------------------------------------

def script_args(f):
    """single decorator for adding script args"""
    args = [
        magic_arguments.argument(
            '--out', type=str,
            help="""The variable in which to store stdout from the script.
            If the script is backgrounded, this will be the stdout *pipe*,
            instead of the stderr text itself.
            """
        ),
        magic_arguments.argument(
            '--err', type=str,
            help="""The variable in which to store stderr from the script.
            If the script is backgrounded, this will be the stderr *pipe*,
            instead of the stderr text itself.
            """
        ),
        magic_arguments.argument(
            '--bg', action="store_true",
            help="""Whether to run the script in the background.
            If given, the only way to see the output of the command is
            with --out/err.
            """
        ),
        magic_arguments.argument(
            '--proc', type=str,
            help="""The variable in which to store Popen instance.
            This is used only when --bg option is given.
            """
        ),
    ]
    for arg in args:
        f = arg(f)
    return f

@magics_class
class ScriptMagics(Magics):
    """Magics for talking to scripts
    
    This defines a base `%%script` cell magic for running a cell
    with a program in a subprocess, and registers a few top-level
    magics that call %%script with common interpreters.
    """
    script_magics = List(
        help="""Extra script cell magics to define
        
        This generates simple wrappers of `%%script foo` as `%%foo`.
        
        If you want to add script magics that aren't on your path,
        specify them in script_paths
        """,
    ).tag(config=True)
    @default('script_magics')
    def _script_magics_default(self):
        """default to a common list of programs"""
        
        defaults = [
            'sh',
            'bash',
            'perl',
            'ruby',
            'python',
            'python2',
            'python3',
            'pypy',
        ]
        if os.name == 'nt':
            defaults.extend([
                'cmd',
            ])
        
        return defaults
    
    script_paths = Dict(
        help="""Dict mapping short 'ruby' names to full paths, such as '/opt/secret/bin/ruby'
        
        Only necessary for items in script_magics where the default path will not
        find the right interpreter.
        """
    ).tag(config=True)
    
    def __init__(self, shell=None):
        super(ScriptMagics, self).__init__(shell=shell)
        self._generate_script_magics()
        self.job_manager = BackgroundJobManager()
        self.bg_processes = []
        atexit.register(self.kill_bg_processes)

    def __del__(self):
        self.kill_bg_processes()
    
    def _generate_script_magics(self):
        cell_magics = self.magics['cell']
        for name in self.script_magics:
            cell_magics[name] = self._make_script_magic(name)
    
    def _make_script_magic(self, name):
        """make a named magic, that calls %%script with a particular program"""
        # expand to explicit path if necessary:
        script = self.script_paths.get(name, name)
        
        @magic_arguments.magic_arguments()
        @script_args
        def named_script_magic(line, cell):
            # if line, add it as cl-flags
            if line:
                 line = "%s %s" % (script, line)
            else:
                line = script
            return self.shebang(line, cell)
        
        # write a basic docstring:
        named_script_magic.__doc__ = \
        """%%{name} script magic
        
        Run cells with {script} in a subprocess.
        
        This is a shortcut for `%%script {script}`
        """.format(**locals())
        
        return named_script_magic
    
    @magic_arguments.magic_arguments()
    @script_args
    @cell_magic("script")
    def shebang(self, line, cell):
        """Run a cell via a shell command
        
        The `%%script` line is like the #! line of script,
        specifying a program (bash, perl, ruby, etc.) with which to run.
        
        The rest of the cell is run by that program.
        
        Examples
        --------
        ::
        
            In [1]: %%script bash
               ...: for i in 1 2 3; do
               ...:   echo $i
               ...: done
            1
            2
            3
        """
        argv = arg_split(line, posix = not sys.platform.startswith('win'))
        args, cmd = self.shebang.parser.parse_known_args(argv)
        
        try:
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, stdin=PIPE)
        except OSError as e:
            if e.errno == errno.ENOENT:
                print("Couldn't find program: %r" % cmd[0])
                return
            else:
                raise
        
        if not cell.endswith('\n'):
            cell += '\n'
        cell = cell.encode('utf8', 'replace')
        if args.bg:
            self.bg_processes.append(p)
            self._gc_bg_processes()
            if args.out:
                self.shell.user_ns[args.out] = p.stdout
            if args.err:
                self.shell.user_ns[args.err] = p.stderr
            self.job_manager.new(self._run_script, p, cell, daemon=True)
            if args.proc:
                self.shell.user_ns[args.proc] = p
            return
        
        try:
            out, err = p.communicate(cell)
        except KeyboardInterrupt:
            try:
                p.send_signal(signal.SIGINT)
                time.sleep(0.1)
                if p.poll() is not None:
                    print("Process is interrupted.")
                    return
                p.terminate()
                time.sleep(0.1)
                if p.poll() is not None:
                    print("Process is terminated.")
                    return
                p.kill()
                print("Process is killed.")
            except OSError:
                pass
            except Exception as e:
                print("Error while terminating subprocess (pid=%i): %s" \
                    % (p.pid, e))
            return
        out = py3compat.decode(out)
        err = py3compat.decode(err)
        if args.out:
            self.shell.user_ns[args.out] = out
        else:
            sys.stdout.write(out)
            sys.stdout.flush()
        if args.err:
            self.shell.user_ns[args.err] = err
        else:
            sys.stderr.write(err)
            sys.stderr.flush()
    
    def _run_script(self, p, cell):
        """callback for running the script in the background"""
        p.stdin.write(cell)
        p.stdin.close()
        p.wait()

    @line_magic("killbgscripts")
    def killbgscripts(self, _nouse_=''):
        """Kill all BG processes started by %%script and its family."""
        self.kill_bg_processes()
        print("All background processes were killed.")

    def kill_bg_processes(self):
        """Kill all BG processes which are still running."""
        if not self.bg_processes:
            return
        for p in self.bg_processes:
            if p.poll() is None:
                try:
                    p.send_signal(signal.SIGINT)
                except:
                    pass
        time.sleep(0.1)
        self._gc_bg_processes()
        if not self.bg_processes:
            return
        for p in self.bg_processes:
            if p.poll() is None:
                try:
                    p.terminate()
                except:
                    pass
        time.sleep(0.1)
        self._gc_bg_processes()
        if not self.bg_processes:
            return
        for p in self.bg_processes:
            if p.poll() is None:
                try:
                    p.kill()
                except:
                    pass
        self._gc_bg_processes()

    def _gc_bg_processes(self):
        self.bg_processes = [p for p in self.bg_processes if p.poll() is None]
