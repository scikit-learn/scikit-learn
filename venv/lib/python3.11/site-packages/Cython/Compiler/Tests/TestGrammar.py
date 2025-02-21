# mode: run
# tag: syntax

"""
Uses TreeFragment to test invalid syntax.
"""

from __future__ import absolute_import

import ast
import textwrap

from ...TestUtils import CythonTest
from .. import ExprNodes
from ..Errors import CompileError

# Copied from CPython's test_grammar.py
VALID_UNDERSCORE_LITERALS = [
    '0_0_0',
    '4_2',
    '1_0000_0000',
    '0b1001_0100',
    '0xffff_ffff',
    '0o5_7_7',
    '1_00_00.5',
    '1_00_00.5j',
    '1_00_00.5e5',
    '1_00_00j',
    '1_00_00e5_1',
    '1e1_0',
    '.1_4',
    '.1_4e1',
    '0b_0',
    '0x_f',
    '0o_5',
    '1_00_00j',
    '1_00_00.5j',
    '1_00_00e5_1j',
    '.1_4j',
    '(1_2.5+3_3j)',
    '(.5_6j)',
]

# Copied from CPython's test_grammar.py
INVALID_UNDERSCORE_LITERALS = [
    # Trailing underscores:
    '0_',
    '42_',
    '1.4j_',
    '0x_',
    '0b1_',
    '0xf_',
    '0o5_',
    '0 if 1_Else 1',
    # Underscores in the base selector:
    '0_b0',
    '0_xf',
    '0_o5',
    # Old-style octal, still disallowed:
    # FIXME: still need to support PY_VERSION_HEX < 3
    '0_7',
    '09_99',
    # Multiple consecutive underscores:
    '4_______2',
    '0.1__4',
    '0.1__4j',
    '0b1001__0100',
    '0xffff__ffff',
    '0x___',
    '0o5__77',
    '1e1__0',
    '1e1__0j',
    # Underscore right before a dot:
    '1_.4',
    '1_.4j',
    # Underscore right after a dot:
    '1._4',
    '1._4j',
    '._5',
    '._5j',
    # Underscore right after a sign:
    '1.0e+_1',
    '1.0e+_1j',
    # Underscore right before j:
    '1.4_j',
    '1.4e5_j',
    # Underscore right before e:
    '1_e1',
    '1.4_e1',
    '1.4_e1j',
    # Underscore right after e:
    '1e_1',
    '1.4e_1',
    '1.4e_1j',
    # Complex cases with parens:
    '(1+1.5_j_)',
    '(1+1.5_j)',
    # Whitespace in literals
    '1_ 2',
    '1 _2',
    '1_2.2_ 1',
    '1_2.2 _1',
    '1_2e _1',
    '1_2e2 _1',
    '1_2e 2_1',
]


INVALID_ELLIPSIS = [
    (". . .", 2, 0),
    (". ..", 2, 0),
    (".. .", 2, 0),
    (". ...", 2, 0),
    (". ... .", 2, 0),
    (".. ... .", 2, 0),
    (". ... ..", 2, 0),
    ("""
    (
        .
        ..
    )
    """, 3, 4),
    ("""
    [
        ..
        .,
        None
    ]
    """, 3, 4),
    ("""
    {
        None,
        .
        .

        .
    }
    """, 4, 4)
]


class TestGrammar(CythonTest):

    def test_invalid_number_literals(self):
        for literal in INVALID_UNDERSCORE_LITERALS:
            for expression in ['%s', '1 + %s', '%s + 1', '2 * %s', '%s * 2']:
                code = 'x = ' + expression % literal
                try:
                    self.fragment(u'''\
                    # cython: language_level=3
                    ''' + code)
                except CompileError as exc:
                    assert code in [s.strip() for s in str(exc).splitlines()], str(exc)
                else:
                    assert False, "Invalid Cython code '%s' failed to raise an exception" % code

    def test_valid_number_literals(self):
        for literal in VALID_UNDERSCORE_LITERALS:
            for i, expression in enumerate(['%s', '1 + %s', '%s + 1', '2 * %s', '%s * 2']):
                code = 'x = ' + expression % literal
                node = self.fragment(u'''\
                    # cython: language_level=3
                    ''' + code).root
                assert node is not None

                literal_node = node.stats[0].rhs  # StatListNode([SingleAssignmentNode('x', expr)])
                if i > 0:
                    # Add/MulNode() -> literal is first or second operand
                    literal_node = literal_node.operand2 if i % 2 else literal_node.operand1
                if 'j' in literal or 'J' in literal:
                    if '+' in literal:
                        # FIXME: tighten this test
                        assert isinstance(literal_node, ExprNodes.AddNode), (literal, literal_node)
                    else:
                        assert isinstance(literal_node, ExprNodes.ImagNode), (literal, literal_node)
                elif '.' in literal or 'e' in literal or 'E' in literal and not ('0x' in literal or '0X' in literal):
                    assert isinstance(literal_node, ExprNodes.FloatNode), (literal, literal_node)
                else:
                    assert isinstance(literal_node, ExprNodes.IntNode), (literal, literal_node)

    def test_invalid_ellipsis(self):
        ERR = ":{0}:{1}: Expected an identifier or literal"
        for code, line, col in INVALID_ELLIPSIS:
            try:
                ast.parse(textwrap.dedent(code))
            except SyntaxError as exc:
                assert True
            else:
                assert False, "Invalid Python code '%s' failed to raise an exception" % code

            try:
                self.fragment(u'''\
                # cython: language_level=3
                ''' + code)
            except CompileError as exc:
                assert ERR.format(line, col) in str(exc), str(exc)
            else:
                assert False, "Invalid Cython code '%s' failed to raise an exception" % code


if __name__ == "__main__":
    import unittest
    unittest.main()
