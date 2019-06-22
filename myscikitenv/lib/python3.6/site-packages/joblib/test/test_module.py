import sys
import joblib
import pytest
from joblib.testing import check_subprocess_call


def test_version():
    assert hasattr(joblib, '__version__'), (
        "There are no __version__ argument on the joblib module")


@pytest.mark.skipif(sys.version_info < (3, 3), reason="Need python3.3+")
def test_no_start_method_side_effect_on_import():
    # check that importing joblib does not implicitly set the global
    # start_method for multiprocessing.
    code = """if True:
        import joblib
        import multiprocessing as mp
        # The following line would raise RuntimeError if the
        # start_method is already set.
        mp.set_start_method("loky")
    """
    check_subprocess_call([sys.executable, '-c', code])
