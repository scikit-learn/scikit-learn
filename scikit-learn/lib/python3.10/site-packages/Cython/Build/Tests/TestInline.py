import os
import tempfile
import unittest
from Cython.Shadow import inline
from Cython.Build.Inline import safe_type, cymeit
from Cython.TestUtils import CythonTest

try:
    import numpy
    has_numpy = True
except:
    has_numpy = False

test_kwds = dict(force=True, quiet=True)

global_value = 100

class TestInline(CythonTest):
    def setUp(self):
        CythonTest.setUp(self)
        self._call_kwds = dict(test_kwds)
        if os.path.isdir('TEST_TMP'):
            lib_dir = os.path.join('TEST_TMP','inline')
        else:
            lib_dir = tempfile.mkdtemp(prefix='cython_inline_')
        self._call_kwds['lib_dir'] = lib_dir

    def test_simple(self):
        self.assertEqual(inline("return 1+2", **self._call_kwds), 3)

    def test_types(self):
        self.assertEqual(inline("""
            cimport cython
            return cython.typeof(a), cython.typeof(b)
        """, a=1.0, b=[], **self._call_kwds), ('double', 'list object'))

    def test_locals(self):
        a = 1
        b = 2
        self.assertEqual(inline("return a+b", **self._call_kwds), 3)

    def test_globals(self):
        self.assertEqual(inline("return global_value + 1", **self._call_kwds), global_value + 1)

    def test_no_return(self):
        self.assertEqual(inline("""
            a = 1
            cdef double b = 2
            cdef c = []
        """, **self._call_kwds), dict(a=1, b=2.0, c=[]))

    def test_def_node(self):
        foo = inline("def foo(x): return x * x", **self._call_kwds)['foo']
        self.assertEqual(foo(7), 49)

    def test_class_ref(self):
        class Type:
            pass
        tp = inline("Type")['Type']
        self.assertEqual(tp, Type)

    def test_pure(self):
        import cython as cy
        b = inline("""
        b = cy.declare(float, a)
        c = cy.declare(cy.pointer(cy.float), &b)
        return b
        """, a=3, **self._call_kwds)
        self.assertEqual(type(b), float)

    def test_compiler_directives(self):
        self.assertEqual(
            inline('return sum(x)',
                   x=[1, 2, 3],
                   cython_compiler_directives={'boundscheck': False}),
            6
        )

    def test_lang_version(self):
        # GH-3419. Caching for inline code didn't always respect compiler directives.
        inline_divcode = "def f(int a, int b): return a/b"
        self.assertEqual(
            inline(inline_divcode, language_level=2)['f'](5,2),
            2
        )
        self.assertEqual(
            inline(inline_divcode, language_level=3)['f'](5,2),
            2.5
        )
        self.assertEqual(
            inline(inline_divcode, language_level=2)['f'](5,2),
            2
        )

    def test_repeated_use(self):
        inline_mulcode = "def f(int a, int b): return a * b"
        self.assertEqual(inline(inline_mulcode)['f'](5, 2), 10)
        self.assertEqual(inline(inline_mulcode)['f'](5, 3), 15)
        self.assertEqual(inline(inline_mulcode)['f'](6, 2), 12)
        self.assertEqual(inline(inline_mulcode)['f'](5, 2), 10)

        f = inline(inline_mulcode)['f']
        self.assertEqual(f(5, 2), 10)
        self.assertEqual(f(5, 3), 15)

    @unittest.skipIf(not has_numpy, "NumPy is not available")
    def test_numpy(self):
        import numpy
        a = numpy.ndarray((10, 20))
        a[0,0] = 10
        self.assertEqual(safe_type(a), 'numpy.ndarray[numpy.float64_t, ndim=2]')
        self.assertEqual(inline("return a[0,0]", a=a, **self._call_kwds), 10.0)


class TestCymeit(unittest.TestCase):
    def _run(self, code, setup_code=None, **kwargs):
        timings, number = cymeit(code, setup_code=setup_code, **kwargs)

        self.assertGreater(min(timings), 0)

        # Guard that autoscaling leads to reasonable timings.
        # Note: we cannot compare against the expected 0.2 due to large timing variations on CI.
        max_time = max(timing * number for timing in timings)
        if isinstance(max_time, int):
            self.assertGreaterEqual(max_time, 100_000)
        else:
            self.assertGreaterEqual(max_time, 0.0001)
        self.assertGreater(number, 10)  # arbitrary lower bound for our very quick benchmarks

        return timings

    def test_benchmark_simple(self):
        setup_code = "numbers = list(range(0, 1000, 3))"
        self._run("sum([num for num in numbers])", setup_code, repeat=3)

    def test_benchmark_timer(self):
        import time
        setup_code = "numbers = list(range(0, 1000, 3))"
        timings = self._run("sum([num for num in numbers])", setup_code, timer=time.perf_counter, repeat=3)

        for timing in timings:
            self.assertIsInstance(timing, float)

    def test_benchmark_timer_ns(self):
        import time
        setup_code = "numbers = list(range(0, 1000, 3))"
        timings = self._run("sum([num for num in numbers])", setup_code, timer=time.perf_counter_ns, repeat=3)

        for timing in timings:
            self.assertIsInstance(timing, int)

    def test_benchmark_multiline_setup(self):
        setup_code = """
        numbers = list(range(0, 100, 3))

        def csum(numbers):
            result = 0
            for number in numbers:
                result += number
            return result
        """
        self._run("csum(numbers)", setup_code)

    def test_benchmark_multiline_code(self):
        setup_code = "numbers = list(range(0, 100, 3))"
        self._run("""
            sum([
                num
                    for num in numbers
                ])
            """,
            setup_code,
            repeat=3
        )

    def test_benchmark_in_module(self):
        self._run("fsum(range(100))", import_module='math', repeat=2)
