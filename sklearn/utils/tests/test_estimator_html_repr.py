# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import importlib
import sys

import pytest


# TODO(1.8): Remove the entire file
def test_estimator_html_repr_warning():
    with pytest.warns(FutureWarning):
        # Make sure that we check for the warning when loading the module (reloading it
        # if needed).
        module_name = "sklearn.utils._estimator_html_repr"
        if module_name in sys.modules:
            importlib.reload(sys.modules[module_name])
        else:
            importlib.import_module(module_name)

    assert sys.modules[module_name] is not None
