"""Test cases for IR generation."""

import os.path

from mypy.test.config import test_temp_dir
from mypy.test.data import DataDrivenTestCase
from mypy.errors import CompileError

from mypyc.common import TOP_LEVEL_NAME
from mypyc.ir.pprint import format_func
from mypyc.test.testutil import (
    ICODE_GEN_BUILTINS, use_custom_builtins, MypycDataSuite, build_ir_for_single_file,
    assert_test_output, remove_comment_lines, replace_word_size,
    infer_ir_build_options_from_test_name
)

files = [
    'irbuild-basic.test',
    'irbuild-int.test',
    'irbuild-lists.test',
    'irbuild-tuple.test',
    'irbuild-dict.test',
    'irbuild-set.test',
    'irbuild-str.test',
    'irbuild-bytes.test',
    'irbuild-statements.test',
    'irbuild-nested.test',
    'irbuild-classes.test',
    'irbuild-optional.test',
    'irbuild-any.test',
    'irbuild-generics.test',
    'irbuild-try.test',
    'irbuild-strip-asserts.test',
    'irbuild-vectorcall.test',
    'irbuild-unreachable.test',
    'irbuild-isinstance.test',
    'irbuild-dunders.test',
    'irbuild-singledispatch.test',
    'irbuild-constant-fold.test',
]


class TestGenOps(MypycDataSuite):
    files = files
    base_path = test_temp_dir
    optional_out = True

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        """Perform a runtime checking transformation test case."""
        options = infer_ir_build_options_from_test_name(testcase.name)
        if options is None:
            # Skipped test case
            return
        with use_custom_builtins(os.path.join(self.data_prefix, ICODE_GEN_BUILTINS), testcase):
            expected_output = remove_comment_lines(testcase.output)
            expected_output = replace_word_size(expected_output)
            name = testcase.name
            try:
                ir = build_ir_for_single_file(testcase.input, options)
            except CompileError as e:
                actual = e.messages
            else:
                actual = []
                for fn in ir:
                    if (fn.name == TOP_LEVEL_NAME
                            and not name.endswith('_toplevel')):
                        continue
                    actual.extend(format_func(fn))

            assert_test_output(testcase, actual, 'Invalid source code output',
                               expected_output)
