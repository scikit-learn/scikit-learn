# tag: ipython

"""Tests for the Cython magics extension."""


import os
import io
import sys
from contextlib import contextmanager
from unittest import skipIf

from Cython.Build import IpythonMagic
from Cython.TestUtils import CythonTest
from Cython.Compiler.Annotate import AnnotationCCodeWriter

try:
    import IPython.testing.globalipapp
except ImportError:
    # Disable tests and fake helpers for initialisation below.
    def skip_if_not_installed(_):
        return None
else:
    def skip_if_not_installed(c):
        return c

# not using IPython's decorators here because they depend on "nose"
skip_win32 = skipIf(sys.platform == 'win32', "Skip on Windows")

try:
    # disable IPython history thread before it gets started to avoid having to clean it up
    from IPython.core.history import HistoryManager
    HistoryManager.enabled = False
except ImportError:
    pass


@contextmanager
def capture_output():
    backup = sys.stdout, sys.stderr
    try:
        replacement = [
             io.TextIOWrapper(io.BytesIO(), encoding=sys.stdout.encoding),
             io.TextIOWrapper(io.BytesIO(), encoding=sys.stderr.encoding),
        ]
        sys.stdout, sys.stderr = replacement
        output = []
        yield output
    finally:
        sys.stdout, sys.stderr = backup
        for wrapper in replacement:
            wrapper.seek(0)  # rewind
            output.append(wrapper.read())
            wrapper.close()


code = """\
def f(x):
    return 2*x
"""

cython3_code = """\
def f(int x):
    return 2 / x

def call(x):
    return f(*(x,))
"""

pgo_cython3_code = cython3_code + """\
def main():
    for _ in range(100): call(5)
main()
"""

compile_error_code = '''\
cdef extern from *:
    """
    xxx a=1;
    """
    int a;
def doit():
    return a
'''

compile_warning_code = '''\
cdef extern from *:
    """
    #pragma message ( "CWarning" )
    int a = 42;
    """
    int a;
def doit():
    return a
'''


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

    @skip_win32
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

    def test_cython_compile_error_shown(self):
        ip = self._ip
        with capture_output() as out:
            ip.run_cell_magic('cython', '-3', compile_error_code)
        captured_out, captured_err = out

        # it could be that c-level output is captured by distutil-extension
        # (and not by us) and is printed to stdout:
        captured_all = captured_out + "\n" + captured_err
        self.assertTrue("error" in captured_all, msg="error in " + captured_all)

    def test_cython_link_error_shown(self):
        ip = self._ip
        with capture_output() as out:
            ip.run_cell_magic('cython', '-3 -l=xxxxxxxx', code)
        captured_out, captured_err = out

        # it could be that c-level output is captured by distutil-extension
        # (and not by us) and is printed to stdout:
        captured_all = captured_out + "\n!" + captured_err
        self.assertTrue("error" in captured_all, msg="error in " + captured_all)

    def test_cython_warning_shown(self):
        ip = self._ip
        with capture_output() as out:
            # force rebuild, otherwise no warning as after the first success
            # no build step is performed
            ip.run_cell_magic('cython', '-3 -f', compile_warning_code)
        captured_out, captured_err = out

        # check that warning was printed to stdout even if build hasn't failed
        self.assertTrue("CWarning" in captured_out)

    @skip_win32
    def test_cython3_pgo(self):
        # The Cython cell defines the functions f() and call().
        ip = self._ip
        ip.run_cell_magic('cython', '-3 --pgo', pgo_cython3_code)
        ip.ex('g = f(10); h = call(10); main()')
        self.assertEqual(ip.user_ns['g'], 2.0 / 10.0)
        self.assertEqual(ip.user_ns['h'], 2.0 / 10.0)

    @skip_win32
    def test_extlibs(self):
        ip = self._ip
        code = """
from libc.math cimport sin
x = sin(0.0)
        """
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
        self.assertEqual([verbose_log.INFO, verbose_log.DEBUG, verbose_log.INFO],
                          verbose_log.thresholds)

        with mock_distutils() as normal_log:
            ip.run_cell_magic('cython', '', code)
            ip.ex('g = f(10)')
        self.assertEqual(ip.user_ns['g'], 20.0)
        self.assertEqual([normal_log.INFO], normal_log.thresholds)

    def test_cython_no_annotate(self):
        ip = self._ip
        html = ip.run_cell_magic('cython', '', code)
        self.assertTrue(html is None)

    def test_cython_annotate(self):
        ip = self._ip
        html = ip.run_cell_magic('cython', '--annotate', code)
        # somewhat brittle way to differentiate between annotated htmls
        # with/without complete source code:
        self.assertTrue(AnnotationCCodeWriter.COMPLETE_CODE_TITLE not in html.data)

    def test_cython_annotate_default(self):
        ip = self._ip
        html = ip.run_cell_magic('cython', '-a', code)
        # somewhat brittle way to differentiate between annotated htmls
        # with/without complete source code:
        self.assertTrue(AnnotationCCodeWriter.COMPLETE_CODE_TITLE not in html.data)

    def test_cython_annotate_complete_c_code(self):
        ip = self._ip
        html = ip.run_cell_magic('cython', '--annotate-fullc', code)
        # somewhat brittle way to differentiate between annotated htmls
        # with/without complete source code:
        self.assertTrue(AnnotationCCodeWriter.COMPLETE_CODE_TITLE in html.data)
