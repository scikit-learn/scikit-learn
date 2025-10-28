import textwrap
from unittest import TestCase

from ..Code import _indent_chunk

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
