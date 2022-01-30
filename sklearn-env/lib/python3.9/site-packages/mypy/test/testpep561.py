from contextlib import contextmanager
import os
import pytest
import re
import subprocess
from subprocess import PIPE
import sys
import tempfile
from typing import Tuple, List, Generator

import mypy.api
from mypy.test.config import package_path
from mypy.util import try_find_python2_interpreter
from mypy.test.data import DataDrivenTestCase, DataSuite
from mypy.test.config import test_temp_dir
from mypy.test.helpers import assert_string_arrays_equal, perform_file_operations


# NOTE: options.use_builtins_fixtures should not be set in these
# tests, otherwise mypy will ignore installed third-party packages.


class PEP561Suite(DataSuite):
    files = [
        'pep561.test',
    ]
    base_path = '.'

    def run_case(self, test_case: DataDrivenTestCase) -> None:
        test_pep561(test_case)


@contextmanager
def virtualenv(
                python_executable: str = sys.executable
                ) -> Generator[Tuple[str, str], None, None]:
    """Context manager that creates a virtualenv in a temporary directory

    returns the path to the created Python executable"""
    # Sadly, we need virtualenv, as the Python 3 venv module does not support creating a venv
    # for Python 2, and Python 2 does not have its own venv.
    with tempfile.TemporaryDirectory() as venv_dir:
        proc = subprocess.run([sys.executable,
                               '-m',
                               'virtualenv',
                               '-p{}'.format(python_executable),
                               venv_dir], cwd=os.getcwd(), stdout=PIPE, stderr=PIPE)
        if proc.returncode != 0:
            err = proc.stdout.decode('utf-8') + proc.stderr.decode('utf-8')
            raise Exception("Failed to create venv. Do you have virtualenv installed?\n" + err)
        if sys.platform == 'win32':
            yield venv_dir, os.path.abspath(os.path.join(venv_dir, 'Scripts', 'python'))
        else:
            yield venv_dir, os.path.abspath(os.path.join(venv_dir, 'bin', 'python'))


def install_package(pkg: str,
                    python_executable: str = sys.executable,
                    use_pip: bool = True,
                    editable: bool = False) -> None:
    """Install a package from test-data/packages/pkg/"""
    working_dir = os.path.join(package_path, pkg)
    with tempfile.TemporaryDirectory() as dir:
        if use_pip:
            install_cmd = [python_executable, '-m', 'pip', 'install']
            if editable:
                install_cmd.append('-e')
            install_cmd.append('.')
        else:
            install_cmd = [python_executable, 'setup.py']
            if editable:
                install_cmd.append('develop')
            else:
                install_cmd.append('install')
        # Note that newer versions of pip (21.3+) don't
        # follow this env variable, but this is for compatibility
        env = {'PIP_BUILD': dir}
        # Inherit environment for Windows
        env.update(os.environ)
        proc = subprocess.run(install_cmd,
                              cwd=working_dir,
                              stdout=PIPE,
                              stderr=PIPE,
                              env=env)
    if proc.returncode != 0:
        raise Exception(proc.stdout.decode('utf-8') + proc.stderr.decode('utf-8'))


def test_pep561(testcase: DataDrivenTestCase) -> None:
    """Test running mypy on files that depend on PEP 561 packages."""
    assert testcase.old_cwd is not None, "test was not properly set up"
    if 'python2' in testcase.name.lower():
        python = try_find_python2_interpreter()
        if python is None:
            pytest.skip()
    else:
        python = sys.executable

    assert python is not None, "Should be impossible"
    pkgs, pip_args = parse_pkgs(testcase.input[0])
    mypy_args = parse_mypy_args(testcase.input[1])
    use_pip = True
    editable = False
    for arg in pip_args:
        if arg == 'no-pip':
            use_pip = False
        elif arg == 'editable':
            editable = True
    assert pkgs != [], "No packages to install for PEP 561 test?"
    with virtualenv(python) as venv:
        venv_dir, python_executable = venv
        for pkg in pkgs:
            install_package(pkg, python_executable, use_pip, editable)

        cmd_line = list(mypy_args)
        has_program = not ('-p' in cmd_line or '--package' in cmd_line)
        if has_program:
            program = testcase.name + '.py'
            with open(program, 'w', encoding='utf-8') as f:
                for s in testcase.input:
                    f.write('{}\n'.format(s))
            cmd_line.append(program)

        cmd_line.extend(['--no-error-summary'])
        if python_executable != sys.executable:
            cmd_line.append('--python-executable={}'.format(python_executable))

        steps = testcase.find_steps()
        if steps != [[]]:
            steps = [[]] + steps  # type: ignore[operator,assignment]

        for i, operations in enumerate(steps):
            perform_file_operations(operations)

            output = []
            # Type check the module
            out, err, returncode = mypy.api.run(cmd_line)

            # split lines, remove newlines, and remove directory of test case
            for line in (out + err).splitlines():
                if line.startswith(test_temp_dir + os.sep):
                    output.append(line[len(test_temp_dir + os.sep):].rstrip("\r\n"))
                else:
                    # Normalize paths so that the output is the same on Windows and Linux/macOS.
                    line = line.replace(test_temp_dir + os.sep, test_temp_dir + '/')
                    output.append(line.rstrip("\r\n"))
            iter_count = '' if i == 0 else ' on iteration {}'.format(i + 1)
            expected = testcase.output if i == 0 else testcase.output2.get(i + 1, [])

            assert_string_arrays_equal(expected, output,
                               'Invalid output ({}, line {}){}'.format(
                                   testcase.file, testcase.line, iter_count))

        if has_program:
            os.remove(program)


def parse_pkgs(comment: str) -> Tuple[List[str], List[str]]:
    if not comment.startswith('# pkgs:'):
        return ([], [])
    else:
        pkgs_str, *args = comment[7:].split(';')
        return ([pkg.strip() for pkg in pkgs_str.split(',')], [arg.strip() for arg in args])


def parse_mypy_args(line: str) -> List[str]:
    m = re.match('# flags: (.*)$', line)
    if not m:
        return []  # No args; mypy will spit out an error.
    return m.group(1).split()


def test_mypy_path_is_respected() -> None:
    assert False
    packages = 'packages'
    pkg_name = 'a'
    with tempfile.TemporaryDirectory() as temp_dir:
        old_dir = os.getcwd()
        os.chdir(temp_dir)
        try:
            # Create the pkg for files to go into
            full_pkg_name = os.path.join(temp_dir, packages, pkg_name)
            os.makedirs(full_pkg_name)

            # Create the empty __init__ file to declare a package
            pkg_init_name = os.path.join(temp_dir, packages, pkg_name, '__init__.py')
            open(pkg_init_name, 'w', encoding='utf8').close()

            mypy_config_path = os.path.join(temp_dir, 'mypy.ini')
            with open(mypy_config_path, 'w') as mypy_file:
                mypy_file.write('[mypy]\n')
                mypy_file.write('mypy_path = ./{}\n'.format(packages))

            with virtualenv() as venv:
                venv_dir, python_executable = venv

                cmd_line_args = []
                if python_executable != sys.executable:
                    cmd_line_args.append('--python-executable={}'.format(python_executable))
                cmd_line_args.extend(['--config-file', mypy_config_path,
                                      '--package', pkg_name])

                out, err, returncode = mypy.api.run(cmd_line_args)
                assert returncode == 0
        finally:
            os.chdir(old_dir)
