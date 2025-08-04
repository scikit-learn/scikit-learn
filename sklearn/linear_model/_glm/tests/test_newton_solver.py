# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# Most of the solvers are tested elsewhere. Here, we only test certrain low level
# aspects.
import numpy as np

from sklearn._loss import HalfMultinomialLoss, HalfSquaredError
from sklearn.linear_model._glm._newton_solver import (
    min_norm_subgradient,
)
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.utils._testing import assert_allclose


def test_min_norm_subgradient():
    coef = np.array([-2.0, 1, 0, 1])  # coef[-1] might be the intercept

    # The test is mostly independent on the choice of base_loss, with the exception of
    # shapes.
    linear_loss = LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=True)
    alpha = 8.0  # l1_reg_strength
    grad = np.full_like(coef, 2.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [2 - 8, 2 + 8, 0, 2])

    grad = np.full_like(coef, -9.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [-9 - 8, -9 + 8, -1 * (9 - 8), -9])

    linear_loss = LinearModelLoss(base_loss=HalfSquaredError(), fit_intercept=False)
    alpha = 1.0
    grad = np.full_like(coef, 2.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [2 - 1, 2 + 1, +1 * (2 - 1), 2 + 1])

    mnsg = min_norm_subgradient(alpha, -grad, coef, linear_loss)
    assert_allclose(mnsg, [-2 - 1, -2 + 1, -1 * (2 - 1), -2 + 1])

    grad = np.full_like(coef, -9.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, [-9 - 1, -9 + 1, -1 * (9 - 1), -9 + 1])

    n_classes = 3
    coef = np.repeat(coef, n_classes)
    linear_loss = LinearModelLoss(base_loss=HalfMultinomialLoss(), fit_intercept=True)
    alpha = 1.0
    grad = np.full_like(coef, 2.0)
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    expected = np.repeat([2 - 1, 2 + 1, +1 * (2 - 1), 2], 3)
    assert_allclose(mnsg, expected)

    # Same as above but with 2-dim coef, coef.shape = (n_classes, n_features)
    coef = coef.reshape(n_classes, -1, order="F")
    grad = grad.reshape(n_classes, -1, order="F")
    expected = expected.reshape(n_classes, -1, order="F")
    mnsg = min_norm_subgradient(alpha, grad, coef, linear_loss)
    assert_allclose(mnsg, expected)

    # alpha = 0
    mnsg = min_norm_subgradient(0, grad, coef, linear_loss)
    assert mnsg is grad
