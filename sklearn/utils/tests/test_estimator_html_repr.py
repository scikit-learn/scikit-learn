# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest


def test_estimator_html_repr_warning():
    with pytest.warns(FutureWarning):
        from sklearn.utils._estimator_html_repr import estimator_html_repr

    assert estimator_html_repr is not None
