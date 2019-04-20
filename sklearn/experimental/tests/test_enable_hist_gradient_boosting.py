"""Tests for making sure experimental imports work as expected."""

import sys
import os
import os.path as op
import textwrap
import tempfile
from subprocess import check_output, STDOUT, CalledProcessError
from subprocess import TimeoutExpired


def assert_run_python_script(source_code, timeout=60):
    """Utility to check assertions in an independent Python subprocess.

    The script provided in the source code should return 0 and not print
    anything on stderr or stdout.

    This is a port from cloudpickle https://github.com/cloudpipe/cloudpickle

    Parameters
    ----------
    source_code : str
        The Python source code to execute.
    timeout : int
        Time in seconds before timeout.
    """
    def _make_cwd_env():
        """Helper to prepare environment for the child processes"""
        sklearn_repo_folder = op.normpath(
            op.join(op.dirname(__file__), '../..'))
        env = os.environ.copy()
        pythonpath = "{src}{sep}tests{pathsep}{src}".format(
            src=sklearn_repo_folder, sep=os.sep, pathsep=os.pathsep)
        env['PYTHONPATH'] = pythonpath
        return sklearn_repo_folder, env

    fd, source_file = tempfile.mkstemp(suffix='_src_test_sklearn.py')
    os.close(fd)
    try:
        with open(source_file, 'wb') as f:
            f.write(source_code.encode('utf-8'))
        cmd = [sys.executable, source_file]
        cwd, env = _make_cwd_env()
        kwargs = {
            'cwd': cwd,
            'stderr': STDOUT,
            'env': env,
        }
        # If coverage is running, pass the config file to the subprocess
        coverage_rc = os.environ.get("COVERAGE_PROCESS_START")
        if coverage_rc:
            kwargs['env']['COVERAGE_PROCESS_START'] = coverage_rc

        kwargs['timeout'] = timeout
        try:
            try:
                out = check_output(cmd, **kwargs)
            except CalledProcessError as e:
                raise RuntimeError(u"script errored with output:\n%s"
                                   % e.output.decode('utf-8'))
            if out != b"":
                raise AssertionError(out.decode('utf-8'))
        except TimeoutExpired as e:
            raise RuntimeError(u"script timeout, output so far:\n%s"
                               % e.output.decode('utf-8'))
    finally:
        os.unlink(source_file)


def test_imports_strategies():
    # Make sure different import strategies work or fail as expected.

    # Since Python caches the imported modules, we need to run a child process
    # for every test case. Else, the tests would not be independent
    # (manually removing the imports from the cache (sys.modules) is not
    # recommended and can lead to many complications).

    good_import = """
    from sklearn.experimental import enable_hist_gradient_boosting
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.ensemble import GradientBoostingRegressor
    """
    assert_run_python_script(textwrap.dedent(good_import))

    good_import_with_ensemble_first = """
    import sklearn.ensemble
    from sklearn.experimental import enable_hist_gradient_boosting
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.ensemble import GradientBoostingRegressor
    """
    assert_run_python_script(textwrap.dedent(good_import_with_ensemble_first))

    bad_imports = """
    import pytest

    with pytest.raises(ImportError):
        from sklearn.ensemble import HistGradientBoostingClassifier

    with pytest.raises(ImportError):
        from sklearn.ensemble._hist_gradient_boosting import (
            HistGradientBoostingClassifier)

    import sklearn.experimental
    with pytest.raises(ImportError):
        from sklearn.ensemble import HistGradientBoostingClassifier
    """
    assert_run_python_script(textwrap.dedent(bad_imports))
