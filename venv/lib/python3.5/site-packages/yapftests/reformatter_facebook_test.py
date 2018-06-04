# Copyright 2016 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Facebook tests for yapf.reformatter."""

import textwrap
import unittest

from yapf.yapflib import reformatter
from yapf.yapflib import style

from yapftests import yapf_test_helper


class TestsForFacebookStyle(yapf_test_helper.YAPFTest):

  @classmethod
  def setUpClass(cls):
    style.SetGlobalStyle(style.CreateFacebookStyle())

  def testNoNeedForLineBreaks(self):
    unformatted_code = textwrap.dedent("""\
        def overly_long_function_name(
          just_one_arg, **kwargs):
          pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        def overly_long_function_name(just_one_arg, **kwargs):
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testDedentClosingBracket(self):
    unformatted_code = textwrap.dedent("""\
        def overly_long_function_name(
          first_argument_on_the_same_line,
          second_argument_makes_the_line_too_long):
          pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        def overly_long_function_name(
            first_argument_on_the_same_line, second_argument_makes_the_line_too_long
        ):
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testBreakAfterOpeningBracketIfContentsTooBig(self):
    unformatted_code = textwrap.dedent("""\
        def overly_long_function_name(a, b, c, d, e, f, g, h, i, j, k, l, m,
          n, o, p, q, r, s, t, u, v, w, x, y, z):
          pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        def overly_long_function_name(
            a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t, u, \
v, w, x, y, z
        ):
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testDedentClosingBracketWithComments(self):
    unformatted_code = textwrap.dedent("""\
        def overly_long_function_name(
          # comment about the first argument
          first_argument_with_a_very_long_name_or_so,
          # comment about the second argument
          second_argument_makes_the_line_too_long):
          pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        def overly_long_function_name(
            # comment about the first argument
            first_argument_with_a_very_long_name_or_so,
            # comment about the second argument
            second_argument_makes_the_line_too_long
        ):
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testDedentImportAsNames(self):
    code = textwrap.dedent("""\
        from module import (
            internal_function as function,
            SOME_CONSTANT_NUMBER1,
            SOME_CONSTANT_NUMBER2,
            SOME_CONSTANT_NUMBER3,
        )
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self.assertCodeEqual(code, reformatter.Reformat(uwlines))

  def testDedentTestListGexp(self):
    unformatted_code = textwrap.dedent("""\
        try:
            pass
        except (
            IOError, OSError, LookupError, RuntimeError, OverflowError
        ) as exception:
            pass

        try:
            pass
        except (
            IOError, OSError, LookupError, RuntimeError, OverflowError,
        ) as exception:
            pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        try:
            pass
        except (
            IOError, OSError, LookupError, RuntimeError, OverflowError
        ) as exception:
            pass

        try:
            pass
        except (
            IOError,
            OSError,
            LookupError,
            RuntimeError,
            OverflowError,
        ) as exception:
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testBrokenIdempotency(self):
    # TODO(ambv): The following behaviour should be fixed.
    pass0_code = textwrap.dedent("""\
        try:
            pass
        except (IOError, OSError, LookupError, RuntimeError, OverflowError) as exception:
            pass
        """)
    pass1_code = textwrap.dedent("""\
        try:
            pass
        except (
            IOError, OSError, LookupError, RuntimeError, OverflowError
        ) as exception:
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(pass0_code)
    self.assertCodeEqual(pass1_code, reformatter.Reformat(uwlines))

    pass2_code = textwrap.dedent("""\
        try:
            pass
        except (
            IOError, OSError, LookupError, RuntimeError, OverflowError
        ) as exception:
            pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(pass1_code)
    self.assertCodeEqual(pass2_code, reformatter.Reformat(uwlines))

  def testIfExprHangingIndent(self):
    unformatted_code = textwrap.dedent("""\
        if True:
            if True:
                if True:
                    if not self.frobbies and (
                       self.foobars.counters['db.cheeses'] != 1 or
                       self.foobars.counters['db.marshmellow_skins'] != 1):
                        pass
        """)
    expected_formatted_code = textwrap.dedent("""\
        if True:
            if True:
                if True:
                    if not self.frobbies and (
                        self.foobars.counters['db.cheeses'] != 1 or
                        self.foobars.counters['db.marshmellow_skins'] != 1
                    ):
                        pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testSimpleDedenting(self):
    unformatted_code = textwrap.dedent("""\
        if True:
            self.assertEqual(result.reason_not_added, "current preflight is still running")
        """)
    expected_formatted_code = textwrap.dedent("""\
        if True:
            self.assertEqual(
                result.reason_not_added, "current preflight is still running"
            )
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testDedentingWithSubscripts(self):
    unformatted_code = textwrap.dedent("""\
        class Foo:
            class Bar:
                @classmethod
                def baz(cls, clues_list, effect, constraints, constraint_manager):
                    if clues_lists:
                       return cls.single_constraint_not(clues_lists, effect, constraints[0], constraint_manager)

        """)
    expected_formatted_code = textwrap.dedent("""\
        class Foo:
            class Bar:
                @classmethod
                def baz(cls, clues_list, effect, constraints, constraint_manager):
                    if clues_lists:
                        return cls.single_constraint_not(
                            clues_lists, effect, constraints[0], constraint_manager
                        )
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testDedentingCallsWithInnerLists(self):
    code = textwrap.dedent("""\
        class _():
            def _():
                cls.effect_clues = {
                    'effect': Clue((cls.effect_time, 'apache_host'), effect_line, 40)
                }
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self.assertCodeEqual(code, reformatter.Reformat(uwlines))

  def testDedentingListComprehension(self):
    unformatted_code = textwrap.dedent("""\
        class Foo():
            def _pack_results_for_constraint_or():
                self.param_groups = dict(
                    (
                        key + 1, ParamGroup(groups[key], default_converter)
                    ) for key in six.moves.range(len(groups))
                )

                for combination in cls._clues_combinations(clues_lists):
                    if all(
                        cls._verify_constraint(combination, effect, constraint)
                        for constraint in constraints
                    ):
                        pass

                guessed_dict = dict(
                    (
                        key, guessed_pattern_matches[key]
                    ) for key in six.moves.range(len(guessed_pattern_matches))
                )

                content = "".join(
                    itertools.chain(
                        (first_line_fragment, ), lines_between, (last_line_fragment, )
                    )
                )

                rule = Rule(
                    [self.cause1, self.cause2, self.cause1, self.cause2], self.effect, constraints1,
                    Rule.LINKAGE_AND
                )

                assert sorted(log_type.files_to_parse) == [
                    ('localhost', os.path.join(path, 'node_1.log'), super_parser),
                    ('localhost', os.path.join(path, 'node_2.log'), super_parser)
                ]
        """)
    expected_formatted_code = textwrap.dedent("""\
        class Foo():
            def _pack_results_for_constraint_or():
                self.param_groups = dict(
                    (key + 1, ParamGroup(groups[key], default_converter))
                    for key in six.moves.range(len(groups))
                )

                for combination in cls._clues_combinations(clues_lists):
                    if all(
                        cls._verify_constraint(combination, effect, constraint)
                        for constraint in constraints
                    ):
                        pass

                guessed_dict = dict(
                    (key, guessed_pattern_matches[key])
                    for key in six.moves.range(len(guessed_pattern_matches))
                )

                content = "".join(
                    itertools.chain(
                        (first_line_fragment, ), lines_between, (last_line_fragment, )
                    )
                )

                rule = Rule(
                    [self.cause1, self.cause2, self.cause1, self.cause2], self.effect,
                    constraints1, Rule.LINKAGE_AND
                )

                assert sorted(log_type.files_to_parse) == [
                    ('localhost', os.path.join(path, 'node_1.log'), super_parser),
                    ('localhost', os.path.join(path, 'node_2.log'), super_parser)
                ]
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testMustSplitDedenting(self):
    code = textwrap.dedent("""\
        class _():
            def _():
                effect_line = FrontInput(
                    effect_line_offset, line_content,
                    LineSource('localhost', xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx)
                )
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self.assertCodeEqual(code, reformatter.Reformat(uwlines))

  def testDedentIfConditional(self):
    code = textwrap.dedent("""\
        class _():
            def _():
                if True:
                    if not self.frobbies and (
                        self.foobars.counters['db.cheeses'] != 1 or
                        self.foobars.counters['db.marshmellow_skins'] != 1
                    ):
                        pass
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self.assertCodeEqual(code, reformatter.Reformat(uwlines))

  def testDedentSet(self):
    code = textwrap.dedent("""\
        class _():
            def _():
                assert set(self.constraint_links.get_links()) == set(
                    [
                        (2, 10, 100),
                        (2, 10, 200),
                        (2, 20, 100),
                        (2, 20, 200),
                    ]
                )
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    self.assertCodeEqual(code, reformatter.Reformat(uwlines))

  def testDedentingInnerScope(self):
    code = textwrap.dedent("""\
        class Foo():
            @classmethod
            def _pack_results_for_constraint_or(cls, combination, constraints):
                return cls._create_investigation_result(
                    (clue for clue in combination if not clue == Verifier.UNMATCHED),
                    constraints, InvestigationResult.OR
                )
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(code)
    reformatted_code = reformatter.Reformat(uwlines)
    self.assertCodeEqual(code, reformatted_code)

    uwlines = yapf_test_helper.ParseAndUnwrap(reformatted_code)
    reformatted_code = reformatter.Reformat(uwlines)
    self.assertCodeEqual(code, reformatted_code)

  def testCommentWithNewlinesInPrefix(self):
    unformatted_code = textwrap.dedent("""\
        def foo():
            if 0:
                return False
                
                
            #a deadly comment
            elif 1:
                return True


        print(foo())
        """)
    expected_formatted_code = textwrap.dedent("""\
        def foo():
            if 0:
                return False
                
            #a deadly comment
            elif 1:
                return True


        print(foo())
        """)
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))

  def testIfStmtClosingBracket(self):
    unformatted_code = """\
if (isinstance(value  , (StopIteration  , StopAsyncIteration  )) and exc.__cause__ is value_asdfasdfasdfasdfsafsafsafdasfasdfs):
    return False
"""
    expected_formatted_code = """\
if (
    isinstance(value, (StopIteration, StopAsyncIteration)) and
    exc.__cause__ is value_asdfasdfasdfasdfsafsafsafdasfasdfs
):
    return False
"""
    uwlines = yapf_test_helper.ParseAndUnwrap(unformatted_code)
    self.assertCodeEqual(expected_formatted_code, reformatter.Reformat(uwlines))


if __name__ == '__main__':
  unittest.main()
