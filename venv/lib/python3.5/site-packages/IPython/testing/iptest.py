# -*- coding: utf-8 -*-
"""IPython Test Suite Runner.

This module provides a main entry point to a user script to test IPython
itself from the command line. There are two ways of running this script:

1. With the syntax `iptest all`.  This runs our entire test suite by
   calling this script (with different arguments) recursively.  This
   causes modules and package to be tested in different processes, using nose
   or trial where appropriate.
2. With the regular nose syntax, like `iptest -vvs IPython`.  In this form
   the script simply calls nose, but with special command line flags and
   plugins loaded.

"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.


import glob
from io import BytesIO
import os
import os.path as path
import sys
from threading import Thread, Lock, Event
import warnings

import nose.plugins.builtin
from nose.plugins.xunit import Xunit
from nose import SkipTest
from nose.core import TestProgram
from nose.plugins import Plugin
from nose.util import safe_str

from IPython import version_info
from IPython.utils.py3compat import decode
from IPython.utils.importstring import import_item
from IPython.testing.plugin.ipdoctest import IPythonDoctest
from IPython.external.decorators import KnownFailure, knownfailureif

pjoin = path.join


# Enable printing all warnings raise by IPython's modules
warnings.filterwarnings('ignore', message='.*Matplotlib is building the font cache.*', category=UserWarning, module='.*')
warnings.filterwarnings('error', message='.*', category=ResourceWarning, module='.*')
warnings.filterwarnings('error', message=".*{'config': True}.*", category=DeprecationWarning, module='IPy.*')
warnings.filterwarnings('default', message='.*', category=Warning, module='IPy.*')

warnings.filterwarnings('error', message='.*apply_wrapper.*', category=DeprecationWarning, module='.*')
warnings.filterwarnings('error', message='.*make_label_dec', category=DeprecationWarning, module='.*')
warnings.filterwarnings('error', message='.*decorated_dummy.*', category=DeprecationWarning, module='.*')
warnings.filterwarnings('error', message='.*skip_file_no_x11.*', category=DeprecationWarning, module='.*')
warnings.filterwarnings('error', message='.*onlyif_any_cmd_exists.*', category=DeprecationWarning, module='.*')

warnings.filterwarnings('error', message='.*disable_gui.*', category=DeprecationWarning, module='.*')

warnings.filterwarnings('error', message='.*ExceptionColors global is deprecated.*', category=DeprecationWarning, module='.*')

# Jedi older versions
warnings.filterwarnings(
    'error', message='.*elementwise != comparison failed and.*', category=FutureWarning, module='.*')

if version_info < (6,):
    # nose.tools renames all things from `camelCase` to `snake_case` which raise an
    # warning with the runner they also import from standard import library. (as of Dec 2015)
    # Ignore, let's revisit that in a couple of years for IPython 6.
    warnings.filterwarnings(
        'ignore', message='.*Please use assertEqual instead', category=Warning, module='IPython.*')

if version_info < (7,):
    warnings.filterwarnings('ignore', message='.*Completer.complete.*',
                            category=PendingDeprecationWarning, module='.*')
else:
    warnings.warn(
        'Completer.complete was pending deprecation and should be changed to Deprecated', FutureWarning)



# ------------------------------------------------------------------------------
# Monkeypatch Xunit to count known failures as skipped.
# ------------------------------------------------------------------------------
def monkeypatch_xunit():
    try:
        knownfailureif(True)(lambda: None)()
    except Exception as e:
        KnownFailureTest = type(e)

    def addError(self, test, err, capt=None):
        if issubclass(err[0], KnownFailureTest):
            err = (SkipTest,) + err[1:]
        return self.orig_addError(test, err, capt)

    Xunit.orig_addError = Xunit.addError
    Xunit.addError = addError

#-----------------------------------------------------------------------------
# Check which dependencies are installed and greater than minimum version.
#-----------------------------------------------------------------------------
def extract_version(mod):
    return mod.__version__

def test_for(item, min_version=None, callback=extract_version):
    """Test to see if item is importable, and optionally check against a minimum
    version.

    If min_version is given, the default behavior is to check against the
    `__version__` attribute of the item, but specifying `callback` allows you to
    extract the value you are interested in. e.g::

        In [1]: import sys

        In [2]: from IPython.testing.iptest import test_for

        In [3]: test_for('sys', (2,6), callback=lambda sys: sys.version_info)
        Out[3]: True

    """
    try:
        check = import_item(item)
    except (ImportError, RuntimeError):
        # GTK reports Runtime error if it can't be initialized even if it's
        # importable.
        return False
    else:
        if min_version:
            if callback:
                # extra processing step to get version to compare
                check = callback(check)

            return check >= min_version
        else:
            return True

# Global dict where we can store information on what we have and what we don't
# have available at test run time
have = {'matplotlib': test_for('matplotlib'),
        'pygments': test_for('pygments'),
        'sqlite3': test_for('sqlite3')}

#-----------------------------------------------------------------------------
# Test suite definitions
#-----------------------------------------------------------------------------

test_group_names = ['core',
                    'extensions', 'lib', 'terminal', 'testing', 'utils',
                   ]

class TestSection(object):
    def __init__(self, name, includes):
        self.name = name
        self.includes = includes
        self.excludes = []
        self.dependencies = []
        self.enabled = True
    
    def exclude(self, module):
        if not module.startswith('IPython'):
            module = self.includes[0] + "." + module
        self.excludes.append(module.replace('.', os.sep))
    
    def requires(self, *packages):
        self.dependencies.extend(packages)
    
    @property
    def will_run(self):
        return self.enabled and all(have[p] for p in self.dependencies)

# Name -> (include, exclude, dependencies_met)
test_sections = {n:TestSection(n, ['IPython.%s' % n]) for n in test_group_names}


# Exclusions and dependencies
# ---------------------------

# core:
sec = test_sections['core']
if not have['sqlite3']:
    sec.exclude('tests.test_history')
    sec.exclude('history')
if not have['matplotlib']:
    sec.exclude('pylabtools'),
    sec.exclude('tests.test_pylabtools')

# lib:
sec = test_sections['lib']
sec.exclude('kernel')
if not have['pygments']:
    sec.exclude('tests.test_lexers')
# We do this unconditionally, so that the test suite doesn't import
# gtk, changing the default encoding and masking some unicode bugs.
sec.exclude('inputhookgtk')
# We also do this unconditionally, because wx can interfere with Unix signals.
# There are currently no tests for it anyway.
sec.exclude('inputhookwx')
# Testing inputhook will need a lot of thought, to figure out
# how to have tests that don't lock up with the gui event
# loops in the picture
sec.exclude('inputhook')

# testing:
sec = test_sections['testing']
# These have to be skipped on win32 because they use echo, rm, cd, etc.
# See ticket https://github.com/ipython/ipython/issues/87
if sys.platform == 'win32':
    sec.exclude('plugin.test_exampleip')
    sec.exclude('plugin.dtexample')

# don't run jupyter_console tests found via shim
test_sections['terminal'].exclude('console')

# extensions:
sec = test_sections['extensions']
# This is deprecated in favour of rpy2
sec.exclude('rmagic')
# autoreload does some strange stuff, so move it to its own test section
sec.exclude('autoreload')
sec.exclude('tests.test_autoreload')
test_sections['autoreload'] = TestSection('autoreload',
        ['IPython.extensions.autoreload', 'IPython.extensions.tests.test_autoreload'])
test_group_names.append('autoreload')


#-----------------------------------------------------------------------------
# Functions and classes
#-----------------------------------------------------------------------------

def check_exclusions_exist():
    from IPython.paths import get_ipython_package_dir
    from warnings import warn
    parent = os.path.dirname(get_ipython_package_dir())
    for sec in test_sections:
        for pattern in sec.exclusions:
            fullpath = pjoin(parent, pattern)
            if not os.path.exists(fullpath) and not glob.glob(fullpath + '.*'):
                warn("Excluding nonexistent file: %r" % pattern)


class ExclusionPlugin(Plugin):
    """A nose plugin to effect our exclusions of files and directories.
    """
    name = 'exclusions'
    score = 3000  # Should come before any other plugins
    
    def __init__(self, exclude_patterns=None):
        """
        Parameters
        ----------

        exclude_patterns : sequence of strings, optional
          Filenames containing these patterns (as raw strings, not as regular
          expressions) are excluded from the tests.
        """
        self.exclude_patterns = exclude_patterns or []
        super(ExclusionPlugin, self).__init__()

    def options(self, parser, env=os.environ):
        Plugin.options(self, parser, env)
    
    def configure(self, options, config):
        Plugin.configure(self, options, config)
        # Override nose trying to disable plugin.
        self.enabled = True
        
    def wantFile(self, filename):
        """Return whether the given filename should be scanned for tests.
        """
        if any(pat in filename for pat in self.exclude_patterns):
            return False
        return None

    def wantDirectory(self, directory):
        """Return whether the given directory should be scanned for tests.
        """
        if any(pat in directory for pat in self.exclude_patterns):
            return False
        return None


class StreamCapturer(Thread):
    daemon = True  # Don't hang if main thread crashes
    started = False
    def __init__(self, echo=False):
        super(StreamCapturer, self).__init__()
        self.echo = echo
        self.streams = []
        self.buffer = BytesIO()
        self.readfd, self.writefd = os.pipe()
        self.buffer_lock = Lock()
        self.stop = Event()

    def run(self):
        self.started = True

        while not self.stop.is_set():
            chunk = os.read(self.readfd, 1024)

            with self.buffer_lock:
                self.buffer.write(chunk)
            if self.echo:
                sys.stdout.write(decode(chunk))

        os.close(self.readfd)
        os.close(self.writefd)

    def reset_buffer(self):
        with self.buffer_lock:
            self.buffer.truncate(0)
            self.buffer.seek(0)

    def get_buffer(self):
        with self.buffer_lock:
            return self.buffer.getvalue()

    def ensure_started(self):
        if not self.started:
            self.start()

    def halt(self):
        """Safely stop the thread."""
        if not self.started:
            return

        self.stop.set()
        os.write(self.writefd, b'\0')  # Ensure we're not locked in a read()
        self.join()

class SubprocessStreamCapturePlugin(Plugin):
    name='subprocstreams'
    def __init__(self):
        Plugin.__init__(self)
        self.stream_capturer = StreamCapturer()
        self.destination = os.environ.get('IPTEST_SUBPROC_STREAMS', 'capture')
        # This is ugly, but distant parts of the test machinery need to be able
        # to redirect streams, so we make the object globally accessible.
        nose.iptest_stdstreams_fileno = self.get_write_fileno

    def get_write_fileno(self):
        if self.destination == 'capture':
            self.stream_capturer.ensure_started()
            return self.stream_capturer.writefd
        elif self.destination == 'discard':
            return os.open(os.devnull, os.O_WRONLY)
        else:
            return sys.__stdout__.fileno()
    
    def configure(self, options, config):
        Plugin.configure(self, options, config)
        # Override nose trying to disable plugin.
        if self.destination == 'capture':
            self.enabled = True
    
    def startTest(self, test):
        # Reset log capture
        self.stream_capturer.reset_buffer()
    
    def formatFailure(self, test, err):
        # Show output
        ec, ev, tb = err
        captured = self.stream_capturer.get_buffer().decode('utf-8', 'replace')
        if captured.strip():
            ev = safe_str(ev)
            out = [ev, '>> begin captured subprocess output <<',
                    captured,
                    '>> end captured subprocess output <<']
            return ec, '\n'.join(out), tb

        return err
    
    formatError = formatFailure
    
    def finalize(self, result):
        self.stream_capturer.halt()


def run_iptest():
    """Run the IPython test suite using nose.

    This function is called when this script is **not** called with the form
    `iptest all`.  It simply calls nose with appropriate command line flags
    and accepts all of the standard nose arguments.
    """
    # Apply our monkeypatch to Xunit
    if '--with-xunit' in sys.argv and not hasattr(Xunit, 'orig_addError'):
        monkeypatch_xunit()

    arg1 = sys.argv[1]
    if arg1 in test_sections:
        section = test_sections[arg1]
        sys.argv[1:2] = section.includes
    elif arg1.startswith('IPython.') and arg1[8:] in test_sections:
        section = test_sections[arg1[8:]]
        sys.argv[1:2] = section.includes
    else:
        section = TestSection(arg1, includes=[arg1])
        

    argv = sys.argv + [ '--detailed-errors',  # extra info in tracebacks
                        # We add --exe because of setuptools' imbecility (it
                        # blindly does chmod +x on ALL files).  Nose does the
                        # right thing and it tries to avoid executables,
                        # setuptools unfortunately forces our hand here.  This
                        # has been discussed on the distutils list and the
                        # setuptools devs refuse to fix this problem!
                        '--exe',
                        ]
    if '-a' not in argv and '-A' not in argv:
        argv = argv + ['-a', '!crash']

    if nose.__version__ >= '0.11':
        # I don't fully understand why we need this one, but depending on what
        # directory the test suite is run from, if we don't give it, 0 tests
        # get run.  Specifically, if the test suite is run from the source dir
        # with an argument (like 'iptest.py IPython.core', 0 tests are run,
        # even if the same call done in this directory works fine).  It appears
        # that if the requested package is in the current dir, nose bails early
        # by default.  Since it's otherwise harmless, leave it in by default
        # for nose >= 0.11, though unfortunately nose 0.10 doesn't support it.
        argv.append('--traverse-namespace')

    plugins = [ ExclusionPlugin(section.excludes), KnownFailure(),
               SubprocessStreamCapturePlugin() ]
    
    # we still have some vestigial doctests in core
    if (section.name.startswith(('core', 'IPython.core', 'IPython.utils'))):
        plugins.append(IPythonDoctest())
        argv.extend([
            '--with-ipdoctest',
            '--ipdoctest-tests',
            '--ipdoctest-extension=txt',
        ])

    
    # Use working directory set by parent process (see iptestcontroller)
    if 'IPTEST_WORKING_DIR' in os.environ:
        os.chdir(os.environ['IPTEST_WORKING_DIR'])
    
    # We need a global ipython running in this process, but the special
    # in-process group spawns its own IPython kernels, so for *that* group we
    # must avoid also opening the global one (otherwise there's a conflict of
    # singletons).  Ultimately the solution to this problem is to refactor our
    # assumptions about what needs to be a singleton and what doesn't (app
    # objects should, individual shells shouldn't).  But for now, this
    # workaround allows the test suite for the inprocess module to complete.
    if 'kernel.inprocess' not in section.name:
        from IPython.testing import globalipapp
        globalipapp.start_ipython()

    # Now nose can run
    TestProgram(argv=argv, addplugins=plugins)

if __name__ == '__main__':
    run_iptest()
