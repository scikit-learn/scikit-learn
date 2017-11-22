from sklearn.utils.testing import SkipTest


def setup_feature_extraction():
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("Skipping feature_extraction.rst, pandas not installed")


def pytest_runtest_setup(item):
    fname = item.fspath.strpath
    if fname.endswith('modules/feature_extraction.rst'):
        setup_feature_extraction()
