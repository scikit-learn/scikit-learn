import builtins
import sys
import unittest

from ..Builtin import (
    inferred_method_return_types, find_return_type_of_builtin_method,
    unsafe_compile_time_methods, is_safe_compile_time_method,
    builtin_scope,
)

from ..Code import (
    KNOWN_PYTHON_BUILTINS_VERSION, KNOWN_PYTHON_BUILTINS,
)

class TestBuiltinReturnTypes(unittest.TestCase):
    def test_find_return_type_of_builtin_method(self):
        # It's enough to test the method existence in a recent Python that likely has them.
        look_up_methods = sys.version_info >= (3,10)

        for type_name, methods in inferred_method_return_types.items():
            py_type = getattr(builtins, type_name if type_name != 'unicode' else 'str')

            for method_name, return_type_name in methods.items():
                builtin_type = builtin_scope.lookup(type_name).type
                return_type = find_return_type_of_builtin_method(builtin_type, method_name)

                if return_type.is_builtin_type:
                    if '[' in return_type_name:
                        return_type_name = return_type_name.partition('[')[0]
                    if return_type_name == 'T':
                        return_type_name = type_name
                    self.assertEqual(return_type.name, return_type_name)
                    if look_up_methods:
                        self.assertTrue(hasattr(py_type, method_name), f"{type_name}.{method_name}")
                else:
                    self.assertEqual(return_type.empty_declaration_code(pyrex=True), return_type_name)


class TestBuiltinCompatibility(unittest.TestCase):
    def test_python_builtin_compatibility(self):
        expected_builtins = set(KNOWN_PYTHON_BUILTINS)
        if sys.platform != 'win32':
            expected_builtins.discard("WindowsError")
        runtime_builtins = frozenset(
            name for name in dir(builtins)
            if name not in ('__doc__', '__loader__', '__name__', '__package__', '__spec__'))
        if sys.version_info < KNOWN_PYTHON_BUILTINS_VERSION:
            missing_builtins = expected_builtins - runtime_builtins
            if missing_builtins:
                self.skipTest(f'skipping test, older Python release found. Missing builtins: {", ".join(sorted(missing_builtins))}')
            self.skipTest('skipping test, older Python release found.')
        self.assertSetEqual(runtime_builtins, expected_builtins)

    def test_unsafe_compile_time_methods(self):
        """Validate the table of builtin methods that are not safe for compile time evaluation
        against the table of known builtin methods (and their types).
        """
        for builtin_type_name, unsafe_methods in unsafe_compile_time_methods.items():
            self.assertIsInstance(unsafe_methods, set)

            builtin_type = getattr(builtins, builtin_type_name)  # All named types must exist as builtin types.

            known_methods = sorted(
                inferred_method_return_types[builtin_type_name])  # All types are also in "inferred_method_return_types".

            self.assertFalse(unsafe_methods.difference(known_methods))  # Only known methods are listed.

            for method_name in known_methods:
                builtin_method = getattr(builtin_type, method_name, None)
                if builtin_method is None:
                    self.assertIn(method_name, unsafe_methods)  # Non-portable methods are always unsafe.
                    continue
