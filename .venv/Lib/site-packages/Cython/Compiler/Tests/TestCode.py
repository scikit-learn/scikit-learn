import textwrap
from unittest import TestCase

from .. import Naming
from ..Code import _indent_chunk, UtilityCode, process_utility_ccode

class TestIndent(TestCase):
    def _test_indentations(self, chunk, expected):
        for indentation in range(16):
            expected_indented = textwrap.indent(expected, ' ' * indentation)
            for line in expected_indented.splitlines():
                # Validate before the comparison that empty lines got stripped also by textwrap.indent().
                self.assertTrue(line == '' or line.strip(), repr(line))

            with self.subTest(indentation=indentation):
                result = _indent_chunk(chunk, indentation_length=indentation)
                self.assertEqual(expected_indented, result)

    def test_indent_empty(self):
        self._test_indentations('', '')

    def test_indent_empty_lines(self):
        self._test_indentations('\n', '\n')
        self._test_indentations('\n'*2, '\n'*2)
        self._test_indentations('\n'*3, '\n'*3)
        self._test_indentations(' \n'*2, '\n'*2)
        self._test_indentations('\n  \n \n    \n', '\n'*4)

    def test_indent_one_line(self):
        self._test_indentations('abc', 'abc')

    def test_indent_chunk(self):
        chunk = """
            x = 1
            if x == 2:
                print("False")
            else:
                print("True")
        """
        expected = """
x = 1
if x == 2:
    print("False")
else:
    print("True")
"""
        self._test_indentations(chunk, expected)

    def test_indent_empty_line(self):
        chunk = """
            x = 1

            if x == 2:
                print("False")
            else:
                print("True")
        """
        expected = """
x = 1

if x == 2:
    print("False")
else:
    print("True")
"""
        self._test_indentations(chunk, expected)

    def test_indent_empty_line_unclean(self):
        lines = """
            x = 1

            if x == 2:
                print("False")
            else:
                print("True")
        """.splitlines(keepends=True)
        lines[2] = '            \n'
        chunk = ''.join(lines)
        expected = """
x = 1

if x == 2:
    print("False")
else:
    print("True")
"""
        self._test_indentations(chunk, expected)


class TestUtilityCodeProcessing(TestCase):
    def _process(self, code):
        utility_code = UtilityCode()
        formatted_code, is_module_specific = process_utility_ccode(utility_code, None, code)
        self.assertFalse(is_module_specific)  # cannot currently test this case
        return formatted_code

    def assert_formatted_code(self, code: str, expected: str, dedent=False):
        if dedent:
            expected = textwrap.dedent(expected)
        expected = expected.strip() + '\n\n'
        formatted = self._process(code)
        self.assertEqual(formatted, expected)

    def test_format_cstring(self):
        self.assert_formatted_code('''
        Some Text and a CSTRING("""
        spanning "multiple" 'lines'.
        Really.
        """);   # end of C string
        ''',
        expected=r'''
        Some Text and a "\n"
        "        spanning \042multiple\042 'lines'.\n"
        "        Really.\n"
        "        \n"
        ;   # end of C string
        ''',
        dedent=True)

    def test_cglobal(self):
        self.assert_formatted_code("""
        CGLOBAL(name)
        NAMED_CGLOBAL(empty_tuple)
        """,
        expected=f"""
        {Naming.modulestateglobal_cname}->name
        {Naming.modulestateglobal_cname}->{Naming.empty_tuple}
        """)

    def test_empty_builtin(self):
        self.assert_formatted_code("""
        EMPTY(tuple)EMPTY(bytes)
        EMPTY(tuple);EMPTY(bytes)
        EMPTY(unicode)
        EMPTY(bytes)
        EMPTY(tuple)
        """,
        expected=f"""
        {Naming.modulestateglobal_cname}->{Naming.empty_tuple}{Naming.modulestateglobal_cname}->{Naming.empty_bytes}
        {Naming.modulestateglobal_cname}->{Naming.empty_tuple};{Naming.modulestateglobal_cname}->{Naming.empty_bytes}
        {Naming.modulestateglobal_cname}->{Naming.empty_unicode}
        {Naming.modulestateglobal_cname}->{Naming.empty_bytes}
        {Naming.modulestateglobal_cname}->{Naming.empty_tuple}
        """)
