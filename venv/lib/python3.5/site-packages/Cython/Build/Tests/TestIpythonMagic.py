# -*- coding: utf-8 -*-
# tag: ipython

"""Tests for the Cython magics extension."""

from __future__ import absolute_import

import os
import sys
from contextlib import contextmanager
from Cython.Build import IpythonMagic
from Cython.TestUtils import CythonTest

try:
    import IPython.testing.globalipapp
    from IPython.utils import py3compat
except ImportError:
    # Disable tests and fake helpers for initialisation below.
    class _py3compat(object):
        def str_to_unicode(self, s):
            return s

    py3compat = _py3compat()

    def skip_if_not_installed(_):
        return None
else:
    def skip_if_not_installed(c):
        return c

try:
    # disable IPython history thread before it gets started to avoid having to clean it up
    from IPython.core.history import HistoryManager
    HistoryManager.enabled = False
except ImportError:
    pass

code = py3compat.str_to_unicode("""\
def f(x):
    return 2*x
""")

cython3_code = py3compat.str_to_unicode("""\
def f(int x):
    return 2 / x

def call(x):
    return f(*(x,))
""")

pgo_cython3_code = cython3_code + py3compat.str_to_unicode("""\
def main():
    for _ in range(100): call(5)
main()
""")


if sys.platform == 'win32':
    # not using IPython's decorators here because they depend on "nose"
    try:
        from unittest import skip as skip_win32
    except ImportError:
        # poor dev's silent @unittest.skip()
        def skip_win32(dummy):
            def _skip_win32(func):
                return None
            return _skip_win32
else:
    def skip_win32(dummy):
        def _skip_win32(func):
            def wrapper(*args, **kwargs):
                func(*args, **kwargs)
            return wrapper
        return _skip_win32


@skip_if_not_installed
class TestIPythonMagic(CythonTest):

    @classmethod
    def setUpClass(cls):
        CythonTest.setUpClass()
        cls._ip = IPython.testing.globalipapp.get_ipython()

    def setUp(self):
        CythonTest.setUp(self)
        self._ip.extension_manager.load_extension('cython')

    def test_cython_inline(self):
        ip = self._ip
        ip.ex('a=10; b=20')
        result = ip.run_cell_magic('cython_inline', '', 'return a+b')
        self.assertEqual(result, 30)

    @skip_win32('Skip on Windows')
    def test_cython_pyximport(self):
        ip = self._ip
        module_name = '_test_cython_pyximport'
        ip.run_cell_magic('cython_pyximport', module_name, code)
        ip.ex('g = f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)
        ip.run_cell_magic('cython_pyximport', module_name, code)
        ip.ex('h = f(-10)')
        self.assertEqual(ip.user_ns['h'], -20.0)
        try:
            os.remove(module_name + '.pyx')
        except OSError:
            pass

    def test_cython(self):
        ip = self._ip
        ip.run_cell_magic('cython', '', code)
        ip.ex('g = f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)

    def test_cython_name(self):
        # The Cython module named 'mymodule' defines the function f.
        ip = self._ip
        ip.run_cell_magic('cython', '--name=mymodule', code)
        # This module can now be imported in the interactive namespace.
        ip.ex('import mymodule; g = mymodule.f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)

    def test_cython_language_level(self):
        # The Cython cell defines the functions f() and call().
        ip = self._ip
        ip.run_cell_magic('cython', '', cython3_code)
        ip.ex('g = f(10); h = call(10)')
        if sys.version_info[0] < 3:
            self.assertEqual(ip.user_ns['g'], 2 // 10)
            self.assertEqual(ip.user_ns['h'], 2 // 10)
        else:
            self.assertEqual(ip.user_ns['g'], 2.0 / 10.0)
            self.assertEqual(ip.user_ns['h'], 2.0 / 10.0)

    def test_cython3(self):
        # The Cython cell defines the functions f() and call().
        ip = self._ip
        ip.run_cell_magic('cython', '-3', cython3_code)
        ip.ex('g = f(10); h = call(10)')
        self.assertEqual(ip.user_ns['g'], 2.0 / 10.0)
        self.assertEqual(ip.user_ns['h'], 2.0 / 10.0)

    def test_cython2(self):
        # The Cython cell defines the functions f() and call().
        ip = self._ip
        ip.run_cell_magic('cython', '-2', cython3_code)
        ip.ex('g = f(10); h = call(10)')
        self.assertEqual(ip.user_ns['g'], 2 // 10)
        self.assertEqual(ip.user_ns['h'], 2 // 10)

    @skip_win32('Skip on Windows')
    def test_cython3_pgo(self):
        # The Cython cell defines the functions f() and call().
        ip = self._ip
        ip.run_cell_magic('cython', '-3 --pgo', pgo_cython3_code)
        ip.ex('g = f(10); h = call(10); main()')
        self.assertEqual(ip.user_ns['g'], 2.0 / 10.0)
        self.assertEqual(ip.user_ns['h'], 2.0 / 10.0)

    @skip_win32('Skip on Windows')
    def test_extlibs(self):
        ip = self._ip
        code = py3compat.str_to_unicode("""
from libc.math cimport sin
x = sin(0.0)
        """)
        ip.user_ns['x'] = 1
        ip.run_cell_magic('cython', '-l m', code)
        self.assertEqual(ip.user_ns['x'], 0)


    def test_cython_verbose(self):
        ip = self._ip
        ip.run_cell_magic('cython', '--verbose', code)
        ip.ex('g = f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)

    def test_cython_verbose_thresholds(self):
        @contextmanager
        def mock_distutils():
            class MockLog:
                DEBUG = 1
                INFO = 2
                thresholds = [INFO]

                def set_threshold(self, val):
                    self.thresholds.append(val)
                    return self.thresholds[-2]


            new_log = MockLog()
            old_log = IpythonMagic.distutils.log
            try:
                IpythonMagic.distutils.log = new_log
                yield new_log
            finally:
                IpythonMagic.distutils.log = old_log

        ip = self._ip
        with mock_distutils() as verbose_log:
            ip.run_cell_magic('cython', '--verbose', code)
            ip.ex('g = f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)
        self.assertEquals([verbose_log.INFO, verbose_log.DEBUG, verbose_log.INFO],
                          verbose_log.thresholds)

        with mock_distutils() as normal_log:
            ip.run_cell_magic('cython', '', code)
            ip.ex('g = f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)
        self.assertEquals([normal_log.INFO], normal_log.thresholds)
