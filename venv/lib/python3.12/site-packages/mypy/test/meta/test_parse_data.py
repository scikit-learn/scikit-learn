"""
A "meta test" which tests the parsing of .test files. This is not meant to become exhaustive
but to ensure we maintain a basic level of ergonomics for mypy contributors.
"""

from mypy.test.helpers import Suite
from mypy.test.meta._pytest import PytestResult, run_pytest_data_suite


def _run_pytest(data_suite: str) -> PytestResult:
    return run_pytest_data_suite(data_suite, extra_args=[], max_attempts=1)


class ParseTestDataSuite(Suite):
    def test_parse_invalid_case(self) -> None:
        # Act
        result = _run_pytest(
            """
            [case abc]
            s: str
            [case foo-XFAIL]
            s: str
            """
        )

        # Assert
        assert "Invalid testcase id 'foo-XFAIL'" in result.stdout

    def test_parse_invalid_section(self) -> None:
        # Act
        result = _run_pytest(
            """
            [case abc]
            s: str
            [unknownsection]
            abc
            """
        )

        # Assert
        expected_lineno = result.input.splitlines().index("[unknownsection]") + 1
        expected = f".test:{expected_lineno}: Invalid section header [unknownsection] in case 'abc'"
        assert expected in result.stdout

    def test_bad_ge_version_check(self) -> None:
        # Act
        actual = _run_pytest(
            """
            [case abc]
            s: str
            [out version>=3.8]
            abc
            """
        )

        # Assert
        assert (
            "version>=3.8 always true since minimum runtime version is (3, 8)"
            in actual.stdout
        )

    def test_bad_eq_version_check(self) -> None:
        # Act
        actual = _run_pytest(
            """
            [case abc]
            s: str
            [out version==3.7]
            abc
            """
        )

        # Assert
        assert (
            "version==3.7 always false since minimum runtime version is (3, 8)"
            in actual.stdout
        )
