from __future__ import absolute_import

import unittest

import Cython.Compiler.PyrexTypes as PT


class TestMethodDispatcherTransform(unittest.TestCase):

    def test_widest_numeric_type(self):
        def assert_widest(type1, type2, widest):
            self.assertEqual(widest, PT.widest_numeric_type(type1, type2))

        assert_widest(PT.c_int_type, PT.c_long_type, PT.c_long_type)
        assert_widest(PT.c_double_type, PT.c_long_type, PT.c_double_type)
        assert_widest(PT.c_longdouble_type, PT.c_long_type, PT.c_longdouble_type)

        cenum = PT.CEnumType("E", "cenum", typedef_flag=False)
        assert_widest(PT.c_int_type, cenum, PT.c_int_type)


class TestTypeIdentifiers(unittest.TestCase):

    TEST_DATA = [
        ("char*", "char__ptr"),
        ("char *", "char__ptr"),
        ("char **", "char__ptr__ptr"),
        ("_typedef", "_typedef"),
        ("__typedef", "__dundertypedef"),
        ("___typedef", "__dunder_typedef"),
        ("____typedef", "__dunder__dundertypedef"),
        ("_____typedef", "__dunder__dunder_typedef"),
        ("const __typedef", "__const___dundertypedef"),
        ("int[42]", "int__lArr42__rArr"),
        ("int[:]", "int__lArr__D__rArr"),
        ("int[:,:]", "int__lArr__D__comma___D__rArr"),
        ("int[:,:,:]", "int__lArr__D__comma___D__comma___D__rArr"),
        ("int[:,:,...]", "int__lArr__D__comma___D__comma___EL__rArr"),
        ("std::vector", "std__in_vector"),
        ("std::vector&&", "std__in_vector__fwref"),
        ("const std::vector", "__const_std__in_vector"),
        ("const std::vector&", "__const_std__in_vector__ref"),
        ("const_std", "const_std"),
    ]

    def test_escape_special_type_characters(self):
        test_func = PT._escape_special_type_characters  # keep test usage visible for IDEs
        function_name = "_escape_special_type_characters"
        self._test_escape(function_name)

    def test_type_identifier_for_declaration(self):
        test_func = PT.type_identifier_from_declaration  # keep test usage visible for IDEs
        function_name = test_func.__name__
        self._test_escape(function_name)

        # differences due to whitespace removal
        test_data = [
            ("const &std::vector", "const__refstd__in_vector"),
            ("const &std::vector<int>", "const__refstd__in_vector__lAngint__rAng"),
            ("const &&std::vector", "const__fwrefstd__in_vector"),
            ("const &&&std::vector", "const__fwref__refstd__in_vector"),
            ("const &&std::vector", "const__fwrefstd__in_vector"),
            ("void (*func)(int x, float y)",
             "975d51__void__lParen__ptrfunc__rParen__lParenint__spac__etc"),
            ("float ** (*func)(int x, int[:] y)",
             "31883a__float__ptr__ptr__lParen__ptrfunc__rParen__lPar__etc"),
        ]
        self._test_escape(function_name, test_data)

    def _test_escape(self, func_name, test_data=TEST_DATA):
        escape = getattr(PT, func_name)
        for declaration, expected in test_data:
            escaped_value = escape(declaration)
            self.assertEqual(escaped_value, expected, "%s('%s') == '%s' != '%s'" % (
                func_name, declaration, escaped_value, expected))
            # test that the length has been successfully capped
            self.assertLessEqual(len(escaped_value), 64)
