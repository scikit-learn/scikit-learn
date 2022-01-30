"""Test runner for exception handling transform test cases.

The transform inserts exception handling branch operations to IR.
"""

import os.path

from mypy.test.config import test_temp_dir
from mypy.test.data import DataDrivenTestCase
from mypy.errors import CompileError

from mypyc.common import TOP_LEVEL_NAME
from mypyc.ir.pprint import format_func
from mypyc.transform.uninit import insert_uninit_checks
from mypyc.transform.exceptions import insert_exception_handling
from mypyc.transform.refcount import insert_ref_count_opcodes
from mypyc.test.testutil import (
    ICODE_GEN_BUILTINS, use_custom_builtins, MypycDataSuite, build_ir_for_single_file,
    assert_test_output, remove_comment_lines
)
from mypyc.analysis.blockfreq import frequently_executed_blocks

files = [
    'exceptions.test',
    'exceptions-freq.test',
]


class TestExceptionTransform(MypycDataSuite):
    files = files
    base_path = test_temp_dir

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        """Perform a runtime checking transformation test case."""
        with use_custom_builtins(os.path.join(self.data_prefix, ICODE_GEN_BUILTINS), testcase):
            expected_output = remove_comment_lines(testcase.output)
            try:
                ir = build_ir_for_single_file(testcase.input)
            except CompileError as e:
                actual = e.messages
            else:
                actual = []
                for fn in ir:
                    if (fn.name == TOP_LEVEL_NAME
                            and not testcase.name.endswith('_toplevel')):
                        continue
                    insert_uninit_checks(fn)
                    insert_exception_handling(fn)
                    insert_ref_count_opcodes(fn)
                    actual.extend(format_func(fn))
                    if testcase.name.endswith('_freq'):
                        common = frequently_executed_blocks(fn.blocks[0])
                        actual.append('hot blocks: %s' % sorted(b.label for b in common))

            assert_test_output(testcase, actual, 'Invalid source code output',
                               expected_output)
