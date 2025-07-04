# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest


# TODO(1.8): Remove the entire file
def test_estimator_html_repr_warning():
    with pytest.warns(FutureWarning):
        from sklearn.utils._estimator_html_repr import estimator_html_repr

    assert estimator_html_repr is not None
