from __future__ import annotations

import os
from unittest import TestCase, mock

from mypy.inspections import parse_location
from mypy.util import _generate_junit_contents, get_terminal_width


class TestGetTerminalSize(TestCase):
    def test_get_terminal_size_in_pty_defaults_to_80(self) -> None:
        # when run using a pty, `os.get_terminal_size()` returns `0, 0`
        ret = os.terminal_size((0, 0))
        mock_environ = os.environ.copy()
        mock_environ.pop("COLUMNS", None)
        with mock.patch.object(os, "get_terminal_size", return_value=ret):
            with mock.patch.dict(os.environ, values=mock_environ, clear=True):
                assert get_terminal_width() == 80

    def test_parse_location_windows(self) -> None:
        assert parse_location(r"C:\test.py:1:1") == (r"C:\test.py", [1, 1])
        assert parse_location(r"C:\test.py:1:1:1:1") == (r"C:\test.py", [1, 1, 1, 1])


class TestWriteJunitXml(TestCase):
    def test_junit_pass(self) -> None:
        serious = False
        messages_by_file: dict[str | None, list[str]] = {}
        expected = """<?xml version="1.0" encoding="utf-8"?>
<testsuite errors="0" failures="0" name="mypy" skips="0" tests="1" time="1.230">
  <testcase classname="mypy" file="mypy" line="1" name="mypy-py3.14-test-plat" time="1.230">
  </testcase>
</testsuite>
"""
        result = _generate_junit_contents(
            dt=1.23,
            serious=serious,
            messages_by_file=messages_by_file,
            version="3.14",
            platform="test-plat",
        )
        assert result == expected

    def test_junit_fail_escape_xml_chars(self) -> None:
        serious = False
        messages_by_file: dict[str | None, list[str]] = {
            "file1.py": ["Test failed", "another line < > &"]
        }
        expected = """<?xml version="1.0" encoding="utf-8"?>
<testsuite errors="0" failures="1" name="mypy" skips="0" tests="1" time="1.230">
  <testcase classname="mypy" file="file1.py" line="1" name="mypy-py3.14-test-plat file1.py" time="1.230">
    <failure message="mypy produced messages">Test failed
another line &lt; &gt; &amp;</failure>
  </testcase>
</testsuite>
"""
        result = _generate_junit_contents(
            dt=1.23,
            serious=serious,
            messages_by_file=messages_by_file,
            version="3.14",
            platform="test-plat",
        )
        assert result == expected

    def test_junit_fail_two_files(self) -> None:
        serious = False
        messages_by_file: dict[str | None, list[str]] = {
            "file1.py": ["Test failed", "another line"],
            "file2.py": ["Another failure", "line 2"],
        }
        expected = """<?xml version="1.0" encoding="utf-8"?>
<testsuite errors="0" failures="2" name="mypy" skips="0" tests="2" time="1.230">
  <testcase classname="mypy" file="file1.py" line="1" name="mypy-py3.14-test-plat file1.py" time="1.230">
    <failure message="mypy produced messages">Test failed
another line</failure>
  </testcase>
  <testcase classname="mypy" file="file2.py" line="1" name="mypy-py3.14-test-plat file2.py" time="1.230">
    <failure message="mypy produced messages">Another failure
line 2</failure>
  </testcase>
</testsuite>
"""
        result = _generate_junit_contents(
            dt=1.23,
            serious=serious,
            messages_by_file=messages_by_file,
            version="3.14",
            platform="test-plat",
        )
        assert result == expected

    def test_serious_error(self) -> None:
        serious = True
        messages_by_file: dict[str | None, list[str]] = {
            None: ["Error line 1", "Error line 2"]
        }
        expected = """<?xml version="1.0" encoding="utf-8"?>
<testsuite errors="1" failures="0" name="mypy" skips="0" tests="1" time="1.230">
  <testcase classname="mypy" file="mypy" line="1" name="mypy-py3.14-test-plat" time="1.230">
    <failure message="mypy produced messages">Error line 1
Error line 2</failure>
  </testcase>
</testsuite>
"""
        result = _generate_junit_contents(
            dt=1.23,
            serious=serious,
            messages_by_file=messages_by_file,
            version="3.14",
            platform="test-plat",
        )
        assert result == expected
