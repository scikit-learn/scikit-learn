"""Test cases for the command line.

To begin we test that "mypy <directory>[/]" always recurses down the
whole tree.
"""

import os
import re
import subprocess
import sys

from typing import List
from typing import Optional

from mypy.test.config import test_temp_dir, PREFIX
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.helpers import (
    assert_string_arrays_equal, normalize_error_messages, check_test_output_files
)

# Path to Python 3 interpreter
python3_path = sys.executable

# Files containing test case descriptions.
cmdline_files = [
    'cmdline.test',
    'cmdline.pyproject.test',
    'reports.test',
    'envvars.test',
]


class PythonCmdlineSuite(DataSuite):
    files = cmdline_files
    native_sep = True

    def run_case(self, testcase: DataDrivenTestCase) -> None:
        for step in [1] + sorted(testcase.output2):
            test_python_cmdline(testcase, step)


def test_python_cmdline(testcase: DataDrivenTestCase, step: int) -> None:
    assert testcase.old_cwd is not None, "test was not properly set up"
    # Write the program to a file.
    program = '_program.py'
    program_path = os.path.join(test_temp_dir, program)
    with open(program_path, 'w', encoding='utf8') as file:
        for s in testcase.input:
            file.write('{}\n'.format(s))
    args = parse_args(testcase.input[0])
    custom_cwd = parse_cwd(testcase.input[1]) if len(testcase.input) > 1 else None
    args.append('--show-traceback')
    if '--error-summary' not in args:
        args.append('--no-error-summary')
    # Type check the program.
    fixed = [python3_path, '-m', 'mypy']
    env = os.environ.copy()
    env.pop('COLUMNS', None)
    env['PYTHONPATH'] = PREFIX
    process = subprocess.Popen(fixed + args,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=os.path.join(
                                   test_temp_dir,
                                   custom_cwd or ""
                                   ),
                               env=env)
    outb, errb = process.communicate()
    result = process.returncode
    # Split output into lines.
    out = [s.rstrip('\n\r') for s in str(outb, 'utf8').splitlines()]
    err = [s.rstrip('\n\r') for s in str(errb, 'utf8').splitlines()]

    if "PYCHARM_HOSTED" in os.environ:
        for pos, line in enumerate(err):
            if line.startswith('pydev debugger: '):
                # Delete the attaching debugger message itself, plus the extra newline added.
                del err[pos:pos + 2]
                break

    # Remove temp file.
    os.remove(program_path)
    # Compare actual output to expected.
    if testcase.output_files:
        # Ignore stdout, but we insist on empty stderr and zero status.
        if err or result:
            raise AssertionError(
                'Expected zero status and empty stderr%s, got %d and\n%s' %
                (' on step %d' % step if testcase.output2 else '',
                 result, '\n'.join(err + out)))
        check_test_output_files(testcase, step)
    else:
        if testcase.normalize_output:
            out = normalize_error_messages(err + out)
        obvious_result = 1 if out else 0
        if obvious_result != result:
            out.append('== Return code: {}'.format(result))
        expected_out = testcase.output if step == 1 else testcase.output2[step]
        # Strip "tmp/" out of the test so that # E: works...
        expected_out = [s.replace("tmp" + os.sep, "") for s in expected_out]
        assert_string_arrays_equal(expected_out, out,
                                   'Invalid output ({}, line {}){}'.format(
                                       testcase.file, testcase.line,
                                       ' on step %d' % step if testcase.output2 else ''))


def parse_args(line: str) -> List[str]:
    """Parse the first line of the program for the command line.

    This should have the form

      # cmd: mypy <options>

    For example:

      # cmd: mypy pkg/
    """
    m = re.match('# cmd: mypy (.*)$', line)
    if not m:
        return []  # No args; mypy will spit out an error.
    return m.group(1).split()


def parse_cwd(line: str) -> Optional[str]:
    """Parse the second line of the program for the command line.

    This should have the form

      # cwd: <directory>

    For example:

      # cwd: main/subdir
    """
    m = re.match('# cwd: (.*)$', line)
    return m.group(1) if m else None
