import unittest

from Cython import Shadow
from Cython.Compiler import Options, CythonScope, PyrexTypes, Errors

class TestShadow(unittest.TestCase):
    def test_all_directives_in_shadow(self):
        missing_directives = []
        extra_directives = []
        for full_directive in Options.directive_types.keys():
            # Python 2 doesn't support "directive, *rest = full_directive.split('.')"
            split_directive = full_directive.split('.')
            directive, rest = split_directive[0], split_directive[1:]

            scope = Options.directive_scopes.get(full_directive)
            if scope and len(scope) == 1 and scope[0] == "module":
                # module-scoped things can't be used from Cython
                if hasattr(Shadow, directive):
                    extra_directives.append(full_directive)
                continue
            if full_directive == "collection_type":
                # collection_type is current restricted to utility code only
                # so doesn't need to be in Shadow
                continue
            if full_directive == "staticmethod":
                # staticmethod is a weird special-case and not really intended to be
                # used from the cython module
                continue

            if not hasattr(Shadow, directive):
                missing_directives.append(full_directive)
            elif rest:
                directive_value = getattr(Shadow, directive)
                for subdirective in rest:
                    if (hasattr(type(directive_value), '__getattr__') or
                            hasattr(type(directive_value), '__getattribute__')):
                        # skip things like "dataclasses" which override attribute lookup
                        break
        self.assertEqual(missing_directives, [])
        self.assertEqual(extra_directives, [])

    def test_all_types_in_shadow(self):
        cython_scope = CythonScope.create_cython_scope(None)
        # Not doing load_cythonscope at this stage because it requires a proper context and
        # Errors.py to be set up

        missing_types = []
        for key in cython_scope.entries.keys():
            if key.startswith('__') and key.endswith('__'):
                continue
            if key in ('PyTypeObject', 'PyObject_TypeCheck'):
                # These are declared in Shadow.py for reasons that look to
                # be an implementation detail, but it isn't our intention for
                # users to access them from Pure Python mode.
                continue
            if not hasattr(Shadow, key):
                missing_types.append(key)
        self.assertEqual(missing_types, [])

    def test_int_types_in_shadow(self):
        missing_types = []
        for int_name in Shadow.int_types:
            for sign in ['', 'u', 's']:
                name = sign + int_name

                if sign and (
                        int_name in ['Py_UNICODE', 'Py_UCS4', 'Py_ssize_t',
                                     'ssize_t', 'ptrdiff_t', 'Py_hash_t'] or
                        name == "usize_t"):
                    # size_t is special-cased here a little since ssize_t legitimate
                    # but usize_t isn't
                    self.assertNotIn(name, dir(Shadow))
                    self.assertNotIn('p_' + name, dir(Shadow))
                    continue

                if not hasattr(Shadow, name):
                    missing_types.append(name)

                for ptr in range(1, 4):
                    ptr_name = 'p' * ptr + '_' + name
                    if not hasattr(Shadow, ptr_name):
                        missing_types.append(ptr_name)
        self.assertEqual(missing_types, [])

    def test_most_types(self):
        # TODO it's unfortunately hard to get a definite list of types to confirm that they're
        # present (because they're obtained by on-the-fly string parsing in `cython_scope.lookup_type`)

        cython_scope = CythonScope.create_cython_scope(None)
        # Set up just enough of "Context" and "Errors" that CythonScope.lookup_type can fail
        class Context:
            cpp = False
            language_level = 3
            future_directives = []
        cython_scope._context = Context
        Errors.init_thread()

        missing_types = []
        missing_lookups = []
        for (signed, longness, name), type_ in PyrexTypes.modifiers_and_name_to_type.items():
            if name == 'object':
                continue  # This probably shouldn't be in Shadow
            if not hasattr(Shadow, name):
                missing_types.append(name)
            if not cython_scope.lookup_type(name):
                missing_lookups.append(name)
            for ptr in range(1, 4):
                ptr_name = 'p' * ptr + '_' + name
                if not hasattr(Shadow, ptr_name):
                    missing_types.append(ptr_name)
                if not cython_scope.lookup_type(ptr_name):
                    missing_lookups.append(ptr_name)
        self.assertEqual(missing_types, [])
        self.assertEqual(missing_lookups, [])
