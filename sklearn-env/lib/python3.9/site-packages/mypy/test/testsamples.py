import os.path
from typing import List, Set

from mypy.test.helpers import Suite, run_mypy


class SamplesSuite(Suite):
    """Test that we can type check some sample code."""

    def test_samples(self) -> None:
        for f in find_files(os.path.join('test-data', 'samples'), suffix='.py'):
            mypy_args = ['--no-strict-optional']
            if f == os.path.join('test-data', 'samples', 'crawl2.py'):
                # This test requires 3.5 for async functions
                mypy_args.append('--python-version=3.5')
            run_mypy(mypy_args + [f])

    def test_stdlibsamples(self) -> None:
        seen: Set[str] = set()
        stdlibsamples_dir = os.path.join('test-data', 'stdlib-samples', '3.2', 'test')
        modules: List[str] = []
        for f in find_files(stdlibsamples_dir, prefix='test_', suffix='.py'):
            if f not in seen:
                seen.add(f)
                modules.append(f)
        if modules:
            # TODO: Remove need for --no-strict-optional
            run_mypy(['--no-strict-optional', '--platform=linux'] + modules)


def find_files(base: str, prefix: str = '', suffix: str = '') -> List[str]:
    return [os.path.join(root, f)
            for root, dirs, files in os.walk(base)
            for f in files
            if f.startswith(prefix) and f.endswith(suffix)]
