import sys
import unittest
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO  # doesn't accept 'str' in Py2

from Cython.Utils import (
    _CACHE_NAME_PATTERN, _build_cache_name, _find_cache_attributes,
    build_hex_version, cached_method, clear_method_caches, try_finally_contextmanager,
    print_version, normalise_float_repr,
)

METHOD_NAME = "cached_next"
CACHE_NAME = _build_cache_name(METHOD_NAME)
NAMES = CACHE_NAME, METHOD_NAME

class Cached(object):
    @cached_method
    def cached_next(self, x):
        return next(x)


class TestCythonUtils(unittest.TestCase):
    def test_build_hex_version(self):
        self.assertEqual('0x001D00A1', build_hex_version('0.29a1'))
        self.assertEqual('0x001D03C4', build_hex_version('0.29.3rc4'))
        self.assertEqual('0x001D00F0', build_hex_version('0.29'))
        self.assertEqual('0x040000F0', build_hex_version('4.0'))

    ############################## Cached Methods ##############################

    def test_cache_method_name(self):
        method_name = "foo"
        cache_name = _build_cache_name(method_name)
        match = _CACHE_NAME_PATTERN.match(cache_name)

        self.assertIsNot(match, None)
        self.assertEqual(match.group(1), method_name)

    def test_requirements_for_Cached(self):
        obj = Cached()

        self.assertFalse(hasattr(obj, CACHE_NAME))
        self.assertTrue(hasattr(obj, METHOD_NAME))
        self.set_of_names_equal(obj, set())

    def set_of_names_equal(self, obj, value):
        self.assertEqual(set(_find_cache_attributes(obj)), value)

    def test_find_cache_attributes(self):
        obj = Cached()
        method_name = "bar"
        cache_name = _build_cache_name(method_name)

        setattr(obj, CACHE_NAME, {})
        setattr(obj, cache_name, {})

        self.assertFalse(hasattr(obj, method_name))
        self.set_of_names_equal(obj, {NAMES, (cache_name, method_name)})

    def test_cached_method(self):
        obj = Cached()
        value = iter(range(3))  # iter for Py2
        cache = {(value,): 0}

        # cache args
        self.assertEqual(obj.cached_next(value), 0)
        self.set_of_names_equal(obj, {NAMES})
        self.assertEqual(getattr(obj, CACHE_NAME), cache)

        # use cache
        self.assertEqual(obj.cached_next(value), 0)
        self.set_of_names_equal(obj, {NAMES})
        self.assertEqual(getattr(obj, CACHE_NAME), cache)

    def test_clear_method_caches(self):
        obj = Cached()
        value = iter(range(3))  # iter for Py2
        cache = {(value,): 1}

        obj.cached_next(value)  # cache args

        clear_method_caches(obj)
        self.set_of_names_equal(obj, set())

        self.assertEqual(obj.cached_next(value), 1)
        self.set_of_names_equal(obj, {NAMES})
        self.assertEqual(getattr(obj, CACHE_NAME), cache)

    def test_clear_method_caches_with_missing_method(self):
        obj = Cached()
        method_name = "bar"
        cache_name = _build_cache_name(method_name)
        names = cache_name, method_name

        setattr(obj, cache_name, object())

        self.assertFalse(hasattr(obj, method_name))
        self.set_of_names_equal(obj, {names})

        clear_method_caches(obj)
        self.set_of_names_equal(obj, {names})

    def test_try_finally_contextmanager(self):
        states = []
        @try_finally_contextmanager
        def gen(*args, **kwargs):
            states.append("enter")
            yield (args, kwargs)
            states.append("exit")

        with gen(1, 2, 3, x=4) as call_args:
            assert states == ["enter"]
            self.assertEqual(call_args, ((1, 2, 3), {'x': 4}))
        assert states == ["enter", "exit"]

        class MyException(RuntimeError):
            pass

        del states[:]
        with self.assertRaises(MyException):
            with gen(1, 2, y=4) as call_args:
                assert states == ["enter"]
                self.assertEqual(call_args, ((1, 2), {'y': 4}))
                raise MyException("FAIL INSIDE")
            assert states == ["enter", "exit"]

        del states[:]
        with self.assertRaises(StopIteration):
            with gen(1, 2, y=4) as call_args:
                assert states == ["enter"]
                self.assertEqual(call_args, ((1, 2), {'y': 4}))
                raise StopIteration("STOP")
            assert states == ["enter", "exit"]

    def test_print_version(self):
        orig_stderr = sys.stderr
        orig_stdout = sys.stdout
        stderr = sys.stderr = StringIO()
        stdout = sys.stdout = StringIO()
        try:
            print_version()
        finally:
            sys.stdout = orig_stdout
            sys.stderr = orig_stderr

        stdout = stdout.getvalue()
        stderr = stderr.getvalue()

        from .. import __version__ as version
        self.assertIn(version, stdout)
        if stderr:  # Depends on os.fstat(1/2).
            self.assertIn(version, stderr)

    def test_print_version_stdouterr(self):
        orig_stderr = sys.stderr
        orig_stdout = sys.stdout
        stdout = sys.stdout = sys.stderr = StringIO()  # same!
        try:
            print_version()
        finally:
            sys.stdout = orig_stdout
            sys.stderr = orig_stderr

        stdout = stdout.getvalue()

        from .. import __version__ as version
        self.assertIn(version, stdout)
        self.assertEqual(stdout.count(version), 1)

    def test_normalise_float_repr(self):
        examples = [
            ('.0', '.0'),
            ('.000000', '.0'),
            ('.1', '.1'),
            ('1.', '1.'),
            ('1.0', '1.'),
            ('1.000000000000000000000', '1.'),
            ('00000000000000000000001.000000000000000000000', '1.'),
            ('12345.0025', '12345.0025'),
            ('1E5', '100000.'),
            ('.1E-5', '.000001'),
            ('1.1E-5', '.000011'),
            ('12.3E-5', '.000123'),
            ('.1E10', '1000000000.'),
            ('1.1E10', '11000000000.'),
            ('123.4E10', '1234000000000.'),
            ('123.456E0', '123.456'),
            ('123.456E-1', '12.3456'),
            ('123.456E-2', '1.23456'),
            ('123.456E1', '1234.56'),
            ('123.456E2', '12345.6'),
            ('2.1E80', '210000000000000000000000000000000000000000000000000000000000000000000000000000000.'),
        ]

        for float_str, norm_str in examples:
            self.assertEqual(float(float_str), float(norm_str))  # safety check for test data

            result = normalise_float_repr(float_str)
            self.assertEqual(float(float_str), float(result))
            self.assertEqual(
                result, norm_str,
                "normalise_float_repr(%r) == %r != %r  (%.330f)" % (float_str, result, norm_str, float(float_str))
            )
