"""Tests for mypy incremental error output."""
from typing import List

from mypy import build
from mypy.test.helpers import assert_string_arrays_equal
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.modulefinder import BuildSource
from mypy.errors import CompileError
from mypy.options import Options


class ErrorStreamSuite(DataSuite):
    required_out_section = True
    base_path = '.'
    files = ['errorstream.test']

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        test_error_stream(testcase)


def test_error_stream(testcase: DataDrivenTestCase) -> None:
    """Perform a single error streaming test case.

    The argument contains the description of the test case.
    """
    options = Options()
    options.show_traceback = True

    logged_messages: List[str] = []

    def flush_errors(msgs: List[str], serious: bool) -> None:
        if msgs:
            logged_messages.append('==== Errors flushed ====')
            logged_messages.extend(msgs)

    sources = [BuildSource('main', '__main__', '\n'.join(testcase.input))]
    try:
        build.build(sources=sources,
                    options=options,
                    flush_errors=flush_errors)
    except CompileError as e:
        assert e.messages == []

    assert_string_arrays_equal(testcase.output, logged_messages,
                               'Invalid output ({}, line {})'.format(
                                   testcase.file, testcase.line))
