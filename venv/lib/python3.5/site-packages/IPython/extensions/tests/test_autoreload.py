"""Tests for autoreload extension.
"""
#-----------------------------------------------------------------------------
#  Copyright (c) 2012 IPython Development Team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os
import sys
import tempfile
import textwrap
import shutil
import random
import time
from io import StringIO

import nose.tools as nt
import IPython.testing.tools as tt

from IPython.testing.decorators import skipif

from IPython.extensions.autoreload import AutoreloadMagics
from IPython.core.events import EventManager, pre_run_cell

#-----------------------------------------------------------------------------
# Test fixture
#-----------------------------------------------------------------------------

noop = lambda *a, **kw: None

class FakeShell(object):

    def __init__(self):
        self.ns = {}
        self.events = EventManager(self, {'pre_run_cell', pre_run_cell})
        self.auto_magics = AutoreloadMagics(shell=self)
        self.events.register('pre_run_cell', self.auto_magics.pre_run_cell)

    register_magics = set_hook = noop

    def run_code(self, code):
        self.events.trigger('pre_run_cell')
        exec(code, self.ns)
        self.auto_magics.post_execute_hook()

    def push(self, items):
        self.ns.update(items)

    def magic_autoreload(self, parameter):
        self.auto_magics.autoreload(parameter)

    def magic_aimport(self, parameter, stream=None):
        self.auto_magics.aimport(parameter, stream=stream)
        self.auto_magics.post_execute_hook()


class Fixture(object):
    """Fixture for creating test module files"""

    test_dir = None
    old_sys_path = None
    filename_chars = "abcdefghijklmopqrstuvwxyz0123456789"

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.old_sys_path = list(sys.path)
        sys.path.insert(0, self.test_dir)
        self.shell = FakeShell()

    def tearDown(self):
        shutil.rmtree(self.test_dir)
        sys.path = self.old_sys_path

        self.test_dir = None
        self.old_sys_path = None
        self.shell = None

    def get_module(self):
        module_name = "tmpmod_" + "".join(random.sample(self.filename_chars,20))
        if module_name in sys.modules:
            del sys.modules[module_name]
        file_name = os.path.join(self.test_dir, module_name + ".py")
        return module_name, file_name

    def write_file(self, filename, content):
        """
        Write a file, and force a timestamp difference of at least one second

        Notes
        -----
        Python's .pyc files record the timestamp of their compilation
        with a time resolution of one second.

        Therefore, we need to force a timestamp difference between .py
        and .pyc, without having the .py file be timestamped in the
        future, and without changing the timestamp of the .pyc file
        (because that is stored in the file).  The only reliable way
        to achieve this seems to be to sleep.
        """

        # Sleep one second + eps
        time.sleep(1.05)

        # Write
        f = open(filename, 'w')
        try:
            f.write(content)
        finally:
            f.close()

    def new_module(self, code):
        mod_name, mod_fn = self.get_module()
        f = open(mod_fn, 'w')
        try:
            f.write(code)
        finally:
            f.close()
        return mod_name, mod_fn

#-----------------------------------------------------------------------------
# Test automatic reloading
#-----------------------------------------------------------------------------

class TestAutoreload(Fixture):

    @skipif(sys.version_info < (3, 6))
    def test_reload_enums(self):
        import enum
        mod_name, mod_fn = self.new_module(textwrap.dedent("""
                                from enum import Enum
                                class MyEnum(Enum):
                                    A = 'A'
                                    B = 'B'
                            """))
        self.shell.magic_autoreload("2")
        self.shell.magic_aimport(mod_name)
        self.write_file(mod_fn, textwrap.dedent("""
                                from enum import Enum
                                class MyEnum(Enum):
                                    A = 'A'
                                    B = 'B'
                                    C = 'C'
                            """))
        with tt.AssertNotPrints(('[autoreload of %s failed:' % mod_name), channel='stderr'):
            self.shell.run_code("pass")  # trigger another reload


    def _check_smoketest(self, use_aimport=True):
        """
        Functional test for the automatic reloader using either
        '%autoreload 1' or '%autoreload 2'
        """

        mod_name, mod_fn = self.new_module("""
x = 9

z = 123  # this item will be deleted

def foo(y):
    return y + 3

class Baz(object):
    def __init__(self, x):
        self.x = x
    def bar(self, y):
        return self.x + y
    @property
    def quux(self):
        return 42
    def zzz(self):
        '''This method will be deleted below'''
        return 99

class Bar:    # old-style class: weakref doesn't work for it on Python < 2.7
    def foo(self):
        return 1
""")

        #
        # Import module, and mark for reloading
        #
        if use_aimport:
            self.shell.magic_autoreload("1")
            self.shell.magic_aimport(mod_name)
            stream = StringIO()
            self.shell.magic_aimport("", stream=stream)
            nt.assert_in(("Modules to reload:\n%s" % mod_name), stream.getvalue())

            with nt.assert_raises(ImportError):
                self.shell.magic_aimport("tmpmod_as318989e89ds")
        else:
            self.shell.magic_autoreload("2")
            self.shell.run_code("import %s" % mod_name)
            stream = StringIO()
            self.shell.magic_aimport("", stream=stream)
            nt.assert_true("Modules to reload:\nall-except-skipped" in
                           stream.getvalue())
        nt.assert_in(mod_name, self.shell.ns)

        mod = sys.modules[mod_name]

        #
        # Test module contents
        #
        old_foo = mod.foo
        old_obj = mod.Baz(9)
        old_obj2 = mod.Bar()

        def check_module_contents():
            nt.assert_equal(mod.x, 9)
            nt.assert_equal(mod.z, 123)

            nt.assert_equal(old_foo(0), 3)
            nt.assert_equal(mod.foo(0), 3)

            obj = mod.Baz(9)
            nt.assert_equal(old_obj.bar(1), 10)
            nt.assert_equal(obj.bar(1), 10)
            nt.assert_equal(obj.quux, 42)
            nt.assert_equal(obj.zzz(), 99)

            obj2 = mod.Bar()
            nt.assert_equal(old_obj2.foo(), 1)
            nt.assert_equal(obj2.foo(), 1)

        check_module_contents()

        #
        # Simulate a failed reload: no reload should occur and exactly
        # one error message should be printed
        #
        self.write_file(mod_fn, """
a syntax error
""")

        with tt.AssertPrints(('[autoreload of %s failed:' % mod_name), channel='stderr'):
            self.shell.run_code("pass") # trigger reload
        with tt.AssertNotPrints(('[autoreload of %s failed:' % mod_name), channel='stderr'):
            self.shell.run_code("pass") # trigger another reload
        check_module_contents()

        #
        # Rewrite module (this time reload should succeed)
        #
        self.write_file(mod_fn, """
x = 10

def foo(y):
    return y + 4

class Baz(object):
    def __init__(self, x):
        self.x = x
    def bar(self, y):
        return self.x + y + 1
    @property
    def quux(self):
        return 43

class Bar:    # old-style class
    def foo(self):
        return 2
""")

        def check_module_contents():
            nt.assert_equal(mod.x, 10)
            nt.assert_false(hasattr(mod, 'z'))

            nt.assert_equal(old_foo(0), 4) # superreload magic!
            nt.assert_equal(mod.foo(0), 4)

            obj = mod.Baz(9)
            nt.assert_equal(old_obj.bar(1), 11) # superreload magic!
            nt.assert_equal(obj.bar(1), 11)

            nt.assert_equal(old_obj.quux, 43)
            nt.assert_equal(obj.quux, 43)

            nt.assert_false(hasattr(old_obj, 'zzz'))
            nt.assert_false(hasattr(obj, 'zzz'))

            obj2 = mod.Bar()
            nt.assert_equal(old_obj2.foo(), 2)
            nt.assert_equal(obj2.foo(), 2)

        self.shell.run_code("pass") # trigger reload
        check_module_contents()

        #
        # Another failure case: deleted file (shouldn't reload)
        #
        os.unlink(mod_fn)

        self.shell.run_code("pass") # trigger reload
        check_module_contents()

        #
        # Disable autoreload and rewrite module: no reload should occur
        #
        if use_aimport:
            self.shell.magic_aimport("-" + mod_name)
            stream = StringIO()
            self.shell.magic_aimport("", stream=stream)
            nt.assert_true(("Modules to skip:\n%s" % mod_name) in
                           stream.getvalue())

            # This should succeed, although no such module exists
            self.shell.magic_aimport("-tmpmod_as318989e89ds")
        else:
            self.shell.magic_autoreload("0")

        self.write_file(mod_fn, """
x = -99
""")

        self.shell.run_code("pass") # trigger reload
        self.shell.run_code("pass")
        check_module_contents()

        #
        # Re-enable autoreload: reload should now occur
        #
        if use_aimport:
            self.shell.magic_aimport(mod_name)
        else:
            self.shell.magic_autoreload("")

        self.shell.run_code("pass") # trigger reload
        nt.assert_equal(mod.x, -99)

    def test_smoketest_aimport(self):
        self._check_smoketest(use_aimport=True)

    def test_smoketest_autoreload(self):
        self._check_smoketest(use_aimport=False)
