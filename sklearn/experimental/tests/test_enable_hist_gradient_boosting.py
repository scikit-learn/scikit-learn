import pytest
import sys


@pytest.fixture
def clean_imports():
    # Removes the relevant scikit-learn related imports (also removes from the
    # cache). This is needed to keep the individual tests functions
    # independent.
    modules_to_delete = (
        'experimental',
        'enable_hist_gradient_boosting',
        'ensemble',
    )
    modules = list(sys.modules.keys())
    for module in modules:
        if any(mod_to_delete in module for mod_to_delete in modules_to_delete):
            del sys.modules[module]


def test_valid_import(clean_imports):
    # recommended way
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    from sklearn.ensemble import HistGradientBoostingClassifier


def test_valid_import_2(clean_imports):
    # recommended way, making sure ensemble can be imported before
    import sklearn.ensemble
    from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    from sklearn.ensemble import HistGradientBoostingClassifier


def test_import_failure(clean_imports):
    # missing enable_hist_gradient_boosting

    with pytest.raises(ImportError):
        from sklearn.ensemble import HistGradientBoostingClassifier

    with pytest.raises(ImportError):
        from sklearn.ensemble._hist_gradient_boosting import (
            HistGradientBoostingClassifier)

    import sklearn.experimental
    with pytest.raises(ImportError):
        from sklearn.ensemble import HistGradientBoostingClassifier
