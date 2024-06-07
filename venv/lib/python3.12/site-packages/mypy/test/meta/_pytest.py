import shlex
import subprocess
import sys
import textwrap
import uuid
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from mypy.test.config import test_data_prefix


@dataclass
class PytestResult:
    input: str
    input_updated: str  # any updates made by --update-data
    stdout: str
    stderr: str


def dedent_docstring(s: str) -> str:
    return textwrap.dedent(s).lstrip()


def run_pytest_data_suite(
    data_suite: str,
    *,
    data_file_prefix: str = "check",
    pytest_node_prefix: str = "mypy/test/testcheck.py::TypeCheckSuite",
    extra_args: Iterable[str],
    max_attempts: int,
) -> PytestResult:
    """
    Runs a suite of data test cases through pytest until either tests pass
    or until a maximum number of attempts (needed for incremental tests).

    :param data_suite: the actual "suite" i.e. the contents of a .test file
    """
    p_test_data = Path(test_data_prefix)
    p_root = p_test_data.parent.parent
    p = p_test_data / f"{data_file_prefix}-meta-{uuid.uuid4()}.test"
    assert not p.exists()
    data_suite = dedent_docstring(data_suite)
    try:
        p.write_text(data_suite)

        test_nodeid = f"{pytest_node_prefix}::{p.name}"
        extra_args = [
            sys.executable,
            "-m",
            "pytest",
            "-n",
            "0",
            "-s",
            *extra_args,
            test_nodeid,
        ]
        cmd = shlex.join(extra_args)
        for i in range(max_attempts - 1, -1, -1):
            print(f">> {cmd}")
            proc = subprocess.run(
                extra_args, capture_output=True, check=False, cwd=p_root
            )
            if proc.returncode == 0:
                break
            prefix = "NESTED PYTEST STDOUT"
            for line in proc.stdout.decode().splitlines():
                print(f"{prefix}: {line}")
                prefix = " " * len(prefix)
            prefix = "NESTED PYTEST STDERR"
            for line in proc.stderr.decode().splitlines():
                print(f"{prefix}: {line}")
                prefix = " " * len(prefix)
            print(f"Exit code {proc.returncode} ({i} attempts remaining)")

        return PytestResult(
            input=data_suite,
            input_updated=p.read_text(),
            stdout=proc.stdout.decode(),
            stderr=proc.stderr.decode(),
        )
    finally:
        p.unlink()
