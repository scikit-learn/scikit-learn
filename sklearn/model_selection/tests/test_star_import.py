"""Regression test for GH #35038: star-import of sklearn.model_selection.

The experimental Halving estimators must stay out of ``__all__``: star-import
resolves every name in ``__all__``, and the module's ``__getattr__`` raises
ImportError for them until ``enable_halving_search_cv`` is imported. With the
names listed, ``from sklearn.model_selection import *`` failed on a fresh
interpreter even though the module itself imported fine.
"""

import subprocess

import pytest


def test_star_import_does_not_trigger_experimental_guard():
    code = "from sklearn.model_selection import *"
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "HalvingGridSearchCV" not in result.stderr
    assert "is experimental" not in result.stderr


def test_experimental_names_absent_from_all():
    import sklearn.model_selection as module

    assert "HalvingGridSearchCV" not in module.__all__
    assert "HalvingRandomSearchCV" not in module.__all__


def test_experimental_names_still_guarded():
    # Explicit attribute access keeps the guidance error for users who want
    # the experimental estimators.
    import sklearn.model_selection as module

    with pytest.raises(ImportError, match="enable_halving_search_cv"):
        module.HalvingGridSearchCV
