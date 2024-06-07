"""
A "meta test" which tests the `--update-data` feature for updating .test files.
Updating the expected output, especially when it's in the form of inline (comment) assertions,
can be brittle, which is why we're "meta-testing" here.
"""

from mypy.test.helpers import Suite
from mypy.test.meta._pytest import PytestResult, dedent_docstring, run_pytest_data_suite


def _run_pytest_update_data(data_suite: str) -> PytestResult:
    """
    Runs a suite of data test cases through 'pytest --update-data' until either tests pass
    or until a maximum number of attempts (needed for incremental tests).
    """
    return run_pytest_data_suite(
        data_suite, extra_args=["--update-data"], max_attempts=3
    )


class UpdateDataSuite(Suite):
    def test_update_data(self) -> None:
        # Note: We test multiple testcases rather than 'test case per test case'
        #       so we could also exercise rewriting multiple testcases at once.
        result = _run_pytest_update_data(
            """
            [case testCorrect]
            s: str = 42  # E: Incompatible types in assignment (expression has type "int", variable has type "str")

            [case testWrong]
            s: str = 42  # E: wrong error

            [case testXfail-xfail]
            s: str = 42  # E: wrong error

            [case testWrongMultiline]
            s: str = 42  # E: foo \
                         # N: bar

            [case testMissingMultiline]
            s: str = 42;  i: int = 'foo'

            [case testExtraneous]
            s: str = 'foo'  # E: wrong error

            [case testExtraneousMultiline]
            s: str = 'foo'  # E: foo \
                            # E: bar

            [case testExtraneousMultilineNonError]
            s: str = 'foo'  # W: foo \
                            # N: bar

            [case testOutCorrect]
            s: str = 42
            [out]
            main:1: error: Incompatible types in assignment (expression has type "int", variable has type "str")

            [case testOutWrong]
            s: str = 42
            [out]
            main:1: error: foobar

            [case testOutWrongIncremental]
            s: str = 42
            [out]
            main:1: error: foobar
            [out2]
            main:1: error: foobar

            [case testWrongMultipleFiles]
            import a, b
            s: str = 42  # E: foo
            [file a.py]
            s1: str = 42  # E: bar
            [file b.py]
            s2: str = 43  # E: baz
            [builtins fixtures/list.pyi]
            """
        )

        # Assert
        expected = dedent_docstring(
            """
        [case testCorrect]
        s: str = 42  # E: Incompatible types in assignment (expression has type "int", variable has type "str")

        [case testWrong]
        s: str = 42  # E: Incompatible types in assignment (expression has type "int", variable has type "str")

        [case testXfail-xfail]
        s: str = 42  # E: wrong error

        [case testWrongMultiline]
        s: str = 42  # E: Incompatible types in assignment (expression has type "int", variable has type "str")

        [case testMissingMultiline]
        s: str = 42;  i: int = 'foo'  # E: Incompatible types in assignment (expression has type "int", variable has type "str") \\
                                      # E: Incompatible types in assignment (expression has type "str", variable has type "int")

        [case testExtraneous]
        s: str = 'foo'

        [case testExtraneousMultiline]
        s: str = 'foo'

        [case testExtraneousMultilineNonError]
        s: str = 'foo'

        [case testOutCorrect]
        s: str = 42
        [out]
        main:1: error: Incompatible types in assignment (expression has type "int", variable has type "str")

        [case testOutWrong]
        s: str = 42
        [out]
        main:1: error: Incompatible types in assignment (expression has type "int", variable has type "str")

        [case testOutWrongIncremental]
        s: str = 42
        [out]
        main:1: error: Incompatible types in assignment (expression has type "int", variable has type "str")
        [out2]
        main:1: error: Incompatible types in assignment (expression has type "int", variable has type "str")

        [case testWrongMultipleFiles]
        import a, b
        s: str = 42  # E: Incompatible types in assignment (expression has type "int", variable has type "str")
        [file a.py]
        s1: str = 42  # E: Incompatible types in assignment (expression has type "int", variable has type "str")
        [file b.py]
        s2: str = 43  # E: Incompatible types in assignment (expression has type "int", variable has type "str")
        [builtins fixtures/list.pyi]
        """
        )
        assert result.input_updated == expected
