"""
Tests for LinearModelLoss

Note that correctness of losses (which compose LinearModelLoss) is already well
covered in the _loss module.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy import linalg, optimize

from sklearn._loss.loss import (
    HalfBinomialLoss,
    HalfMultinomialLoss,
    HalfPoissonLoss,
)
from sklearn.datasets import make_low_rank_matrix
from sklearn.linear_model._linear_loss import LinearModelLoss
from sklearn.utils.extmath import squared_norm
from sklearn.utils.fixes import CSR_CONTAINERS

# We do not need to test all losses, just what LinearModelLoss does on top of the
# base losses.
LOSSES = [HalfBinomialLoss, HalfMultinomialLoss, HalfPoissonLoss]


def random_X_y_coef(
    linear_model_loss, n_samples, n_features, coef_bound=(-2, 2), seed=42
):
    """Random generate y, X and coef in valid range."""
    rng = np.random.RandomState(seed)
    n_dof = n_features + linear_model_loss.fit_intercept
    X = make_low_rank_matrix(
        n_samples=n_samples,
        n_features=n_features,
        random_state=rng,
    )
    coef = linear_model_loss.init_zero_coef(X)

    if linear_model_loss.base_loss.is_multiclass:
        n_classes = linear_model_loss.base_loss.n_classes
        coef.flat[:] = rng.uniform(
            low=coef_bound[0],
            high=coef_bound[1],
            size=n_classes * n_dof,
        )
        if linear_model_loss.fit_intercept:
            raw_prediction = X @ coef[:, :-1].T + coef[:, -1]
        else:
            raw_prediction = X @ coef.T
        proba = linear_model_loss.base_loss.link.inverse(raw_prediction)

        # y = rng.choice(np.arange(n_classes), p=proba) does not work.
        # See https://stackoverflow.com/a/34190035/16761084
        def choice_vectorized(items, p):
            s = p.cumsum(axis=1)
            r = rng.rand(p.shape[0])[:, None]
            k = (s < r).sum(axis=1)
            return items[k]

        y = choice_vectorized(np.arange(n_classes), p=proba).astype(np.float64)
    else:
        coef.flat[:] = rng.uniform(
            low=coef_bound[0],
            high=coef_bound[1],
            size=n_dof,
        )
        if linear_model_loss.fit_intercept:
            raw_prediction = X @ coef[:-1] + coef[-1]
        else:
            raw_prediction = X @ coef
        y = linear_model_loss.base_loss.link.inverse(
            raw_prediction + rng.uniform(low=-1, high=1, size=n_samples)
        )

    return X, y, coef


@pytest.mark.parametrize("base_loss", LOSSES)
@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("n_features", [0, 1, 10])
@pytest.mark.parametrize("dtype", [None, np.float32, np.float64, np.int64])
def test_init_zero_coef(
    base_loss, fit_intercept, n_features, dtype, global_random_seed
):
    """Test that init_zero_coef initializes coef correctly."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    rng = np.random.RandomState(global_random_seed)
    X = rng.normal(size=(5, n_features))
    coef = loss.init_zero_coef(X, dtype=dtype)
    if loss.base_loss.is_multiclass:
        n_classes = loss.base_loss.n_classes
        assert coef.shape == (n_classes, n_features + fit_intercept)
        assert coef.flags["F_CONTIGUOUS"]
    else:
        assert coef.shape == (n_features + fit_intercept,)

    if dtype is None:
        assert coef.dtype == X.dtype
    else:
        assert coef.dtype == dtype

    assert np.count_nonzero(coef) == 0


@pytest.mark.parametrize("base_loss", LOSSES)
@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("sample_weight", [None, "range"])
@pytest.mark.parametrize("l2_reg_strength", [0, 1])
@pytest.mark.parametrize("csr_container", CSR_CONTAINERS)
def test_loss_grad_hess_are_the_same(
    base_loss,
    fit_intercept,
    sample_weight,
    l2_reg_strength,
    csr_container,
    global_random_seed,
):
    """Test that loss and gradient are the same across different functions."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss, n_samples=10, n_features=5, seed=global_random_seed
    )
    X_old, y_old, coef_old = X.copy(), y.copy(), coef.copy()

    if sample_weight == "range":
        sample_weight = np.linspace(1, y.shape[0], num=y.shape[0])

    l1 = loss.loss(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g1 = loss.gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    l2, g2 = loss.loss_gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g3, h3 = loss.gradient_hessian_product(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g4, h4, _ = loss.gradient_hessian(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    assert_allclose(l1, l2)
    assert_allclose(g1, g2)
    assert_allclose(g1, g3)
    assert_allclose(g1, g4)
    # The ravelling only takes effect for multiclass.
    assert_allclose(h4 @ g4.ravel(order="F"), h3(g3).ravel(order="F"))
    # Test that gradient_out and hessian_out are considered properly.
    g_out = np.empty_like(coef)
    h_out = np.empty_like(coef, shape=(coef.size, coef.size))
    g5, h5, _ = loss.gradient_hessian(
        coef,
        X,
        y,
        sample_weight=sample_weight,
        l2_reg_strength=l2_reg_strength,
        gradient_out=g_out,
        hessian_out=h_out,
    )
    assert np.shares_memory(g5, g_out)
    assert np.shares_memory(h5, h_out)
    assert_allclose(g5, g_out)
    assert_allclose(h5, h_out)
    assert_allclose(g1, g5)
    assert_allclose(h5, h4)

    # same for sparse X
    Xs = csr_container(X)
    l1_sp = loss.loss(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g1_sp = loss.gradient(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    l2_sp, g2_sp = loss.loss_gradient(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g3_sp, h3_sp = loss.gradient_hessian_product(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g4_sp, h4_sp, _ = loss.gradient_hessian(
        coef, Xs, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    assert_allclose(l1, l1_sp)
    assert_allclose(l1, l2_sp)
    assert_allclose(g1, g1_sp)
    assert_allclose(g1, g2_sp)
    assert_allclose(g1, g3_sp)
    assert_allclose(h3(g1), h3_sp(g1_sp))
    assert_allclose(g1, g4_sp)
    assert_allclose(h4, h4_sp)

    # X, y and coef should not have changed
    assert_allclose(X, X_old)
    assert_allclose(Xs.toarray(), X_old)
    assert_allclose(y, y_old)
    assert_allclose(coef, coef_old)


@pytest.mark.parametrize("base_loss", LOSSES)
@pytest.mark.parametrize("sample_weight", [None, "range"])
@pytest.mark.parametrize("l2_reg_strength", [0, 1])
@pytest.mark.parametrize("X_container", CSR_CONTAINERS + [None])
def test_loss_gradients_hessp_intercept(
    base_loss, sample_weight, l2_reg_strength, X_container, global_random_seed
):
    """Test that loss and gradient handle intercept correctly."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=False)
    loss_inter = LinearModelLoss(base_loss=base_loss(), fit_intercept=True)
    n_samples, n_features = 10, 5
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss,
        n_samples=n_samples,
        n_features=n_features,
        seed=global_random_seed,
    )

    X[:, -1] = 1  # make last column of 1 to mimic intercept term
    X_inter = X[
        :, :-1
    ]  # exclude intercept column as it is added automatically by loss_inter

    if X_container is not None:
        X = X_container(X)

    if sample_weight == "range":
        sample_weight = np.linspace(1, y.shape[0], num=y.shape[0])

    l, g = loss.loss_gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    _, hessp = loss.gradient_hessian_product(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    l_inter, g_inter = loss_inter.loss_gradient(
        coef, X_inter, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    _, hessp_inter = loss_inter.gradient_hessian_product(
        coef, X_inter, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )

    # Note, that intercept gets no L2 penalty.
    assert l == pytest.approx(
        l_inter + 0.5 * l2_reg_strength * squared_norm(coef.T[-1])
    )

    g_inter_corrected = g_inter
    g_inter_corrected.T[-1] += l2_reg_strength * coef.T[-1]
    assert_allclose(g, g_inter_corrected)

    s = np.random.RandomState(global_random_seed).randn(*coef.shape)
    h = hessp(s)
    h_inter = hessp_inter(s)
    h_inter_corrected = h_inter
    h_inter_corrected.T[-1] += l2_reg_strength * s.T[-1]
    assert_allclose(h, h_inter_corrected)


@pytest.mark.parametrize("base_loss", LOSSES)
@pytest.mark.parametrize("fit_intercept", [False, True])
@pytest.mark.parametrize("sample_weight", [None, "range"])
@pytest.mark.parametrize("l2_reg_strength", [0, 1])
def test_gradients_hessians_numerically(
    base_loss, fit_intercept, sample_weight, l2_reg_strength, global_random_seed
):
    """Test gradients and hessians with numerical derivatives.

    Gradient should equal the numerical derivatives of the loss function.
    Hessians should equal the numerical derivatives of gradients.
    """
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    n_samples, n_features = 10, 5
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss,
        n_samples=n_samples,
        n_features=n_features,
        seed=global_random_seed,
    )
    coef = coef.ravel(order="F")  # this is important only for multinomial loss

    if sample_weight == "range":
        sample_weight = np.linspace(1, y.shape[0], num=y.shape[0])

    # 1. Check gradients numerically
    eps = 1e-6
    g, hessp = loss.gradient_hessian_product(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    # Use a trick to get central finite difference of accuracy 4 (five-point stencil)
    # https://en.wikipedia.org/wiki/Numerical_differentiation
    # https://en.wikipedia.org/wiki/Finite_difference_coefficient
    # approx_g1 = (f(x + eps) - f(x - eps)) / (2*eps)
    approx_g1 = optimize.approx_fprime(
        coef,
        lambda coef: loss.loss(
            coef - eps,
            X,
            y,
            sample_weight=sample_weight,
            l2_reg_strength=l2_reg_strength,
        ),
        2 * eps,
    )
    # approx_g2 = (f(x + 2*eps) - f(x - 2*eps)) / (4*eps)
    approx_g2 = optimize.approx_fprime(
        coef,
        lambda coef: loss.loss(
            coef - 2 * eps,
            X,
            y,
            sample_weight=sample_weight,
            l2_reg_strength=l2_reg_strength,
        ),
        4 * eps,
    )
    # Five-point stencil approximation
    # See: https://en.wikipedia.org/wiki/Five-point_stencil#1D_first_derivative
    approx_g = (4 * approx_g1 - approx_g2) / 3
    assert_allclose(g, approx_g, rtol=1e-2, atol=1e-8)

    # 2. Check hessp numerically along the second direction of the gradient
    vector = np.zeros_like(g)
    vector[1] = 1
    hess_col = hessp(vector)
    # Computation of the Hessian is particularly fragile to numerical errors when doing
    # simple finite differences. Here we compute the grad along a path in the direction
    # of the vector and then use a least-square regression to estimate the slope
    eps = 1e-3
    d_x = np.linspace(-eps, eps, 30)
    d_grad = np.array(
        [
            loss.gradient(
                coef + t * vector,
                X,
                y,
                sample_weight=sample_weight,
                l2_reg_strength=l2_reg_strength,
            )
            for t in d_x
        ]
    )
    d_grad -= d_grad.mean(axis=0)
    approx_hess_col = linalg.lstsq(d_x[:, np.newaxis], d_grad)[0].ravel()
    assert_allclose(approx_hess_col, hess_col, rtol=1e-3)


@pytest.mark.parametrize("fit_intercept", [False, True])
def test_multinomial_coef_shape(fit_intercept, global_random_seed):
    """Test that multinomial LinearModelLoss respects shape of coef."""
    loss = LinearModelLoss(base_loss=HalfMultinomialLoss(), fit_intercept=fit_intercept)
    n_samples, n_features = 10, 5
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss,
        n_samples=n_samples,
        n_features=n_features,
        seed=global_random_seed,
    )
    s = np.random.RandomState(global_random_seed).randn(*coef.shape)

    l, g = loss.loss_gradient(coef, X, y)
    g1 = loss.gradient(coef, X, y)
    g2, hessp = loss.gradient_hessian_product(coef, X, y)
    h = hessp(s)
    assert g.shape == coef.shape
    assert h.shape == coef.shape
    assert_allclose(g, g1)
    assert_allclose(g, g2)
    g3, hess, _ = loss.gradient_hessian(coef, X, y)
    assert g3.shape == coef.shape
    # But full hessian is always 2d.
    assert hess.shape == (coef.size, coef.size)

    coef_r = coef.ravel(order="F")
    s_r = s.ravel(order="F")
    l_r, g_r = loss.loss_gradient(coef_r, X, y)
    g1_r = loss.gradient(coef_r, X, y)
    g2_r, hessp_r = loss.gradient_hessian_product(coef_r, X, y)
    h_r = hessp_r(s_r)
    assert g_r.shape == coef_r.shape
    assert h_r.shape == coef_r.shape
    assert_allclose(g_r, g1_r)
    assert_allclose(g_r, g2_r)

    assert_allclose(g, g_r.reshape(loss.base_loss.n_classes, -1, order="F"))
    assert_allclose(h, h_r.reshape(loss.base_loss.n_classes, -1, order="F"))


@pytest.mark.parametrize("sample_weight", [None, "range"])
def test_multinomial_hessian_3_classes(sample_weight, global_random_seed):
    """Test multinomial hessian for 3 classes and 2 points.

    For n_classes = 3 and n_samples = 2, we have
      p0 = [p0_0, p0_1]
      p1 = [p1_0, p1_1]
      p2 = [p2_0, p2_1]
    and with 2 x 2 diagonal subblocks
      H = [p0 * (1-p0),    -p0 * p1,    -p0 * p2]
          [   -p0 * p1, p1 * (1-p1),    -p1 * p2]
          [   -p0 * p2,    -p1 * p2, p2 * (1-p2)]
      hess = X' H X
    """
    n_samples, n_features, n_classes = 2, 5, 3
    loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=n_classes), fit_intercept=False
    )
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss,
        n_samples=n_samples,
        n_features=n_features,
        seed=global_random_seed,
    )
    coef = coef.ravel(order="F")  # this is important only for multinomial loss

    if sample_weight == "range":
        sample_weight = np.linspace(1, y.shape[0], num=y.shape[0])

    grad, hess, _ = loss.gradient_hessian(
        coef,
        X,
        y,
        sample_weight=sample_weight,
        l2_reg_strength=0,
    )
    # Hessian must be a symmetrix matrix.
    assert_allclose(hess, hess.T)

    weights, intercept, raw_prediction = loss.weight_intercept_raw(coef, X)
    grad_pointwise, proba = loss.base_loss.gradient_proba(
        y_true=y,
        raw_prediction=raw_prediction,
        sample_weight=sample_weight,
    )
    p0d, p1d, p2d, oned = (
        np.diag(proba[:, 0]),
        np.diag(proba[:, 1]),
        np.diag(proba[:, 2]),
        np.diag(np.ones(2)),
    )
    h = np.block(
        [
            [p0d * (oned - p0d), -p0d * p1d, -p0d * p2d],
            [-p0d * p1d, p1d * (oned - p1d), -p1d * p2d],
            [-p0d * p2d, -p1d * p2d, p2d * (oned - p2d)],
        ]
    )
    h = h.reshape((n_classes, n_samples, n_classes, n_samples))
    if sample_weight is None:
        h /= n_samples
    else:
        h *= sample_weight / np.sum(sample_weight)
    # hess_expected.shape = (n_features, n_classes, n_classes, n_features)
    hess_expected = np.einsum("ij, mini, ik->jmnk", X, h, X)
    hess_expected = np.moveaxis(hess_expected, 2, 3)
    hess_expected = hess_expected.reshape(
        n_classes * n_features, n_classes * n_features, order="C"
    )
    assert_allclose(hess_expected, hess_expected.T)
    assert_allclose(hess, hess_expected)


def test_linear_loss_gradient_hessian_raises_wrong_out_parameters():
    """Test that wrong gradient_out and hessian_out raises errors."""
    n_samples, n_features, n_classes = 5, 2, 3
    loss = LinearModelLoss(base_loss=HalfBinomialLoss(), fit_intercept=False)
    X = np.ones((n_samples, n_features))
    y = np.ones(n_samples)
    coef = loss.init_zero_coef(X)
    gradient_out = np.zeros(1)
    with pytest.raises(
        ValueError, match="gradient_out is required to have shape coef.shape"
    ):
        loss.gradient_hessian(
            coef=coef,
            X=X,
            y=y,
            gradient_out=gradient_out,
            hessian_out=None,
        )
    hessian_out = np.zeros(1)
    with pytest.raises(ValueError, match="hessian_out is required to have shape"):
        loss.gradient_hessian(
            coef=coef,
            X=X,
            y=y,
            gradient_out=None,
            hessian_out=hessian_out,
        )

    loss = LinearModelLoss(base_loss=HalfMultinomialLoss(), fit_intercept=False)
    coef = loss.init_zero_coef(X)
    gradient_out = np.zeros((2 * n_classes, n_features))[::2]
    with pytest.raises(ValueError, match="gradient_out must be F-contiguous"):
        loss.gradient_hessian(
            coef=coef,
            X=X,
            y=y,
            gradient_out=gradient_out,
        )
    hessian_out = np.zeros((2 * n_classes * n_features, n_classes * n_features))[::2]
    with pytest.raises(ValueError, match="hessian_out must be contiguous"):
        loss.gradient_hessian(
            coef=coef,
            X=X,
            y=y,
            gradient_out=None,
            hessian_out=hessian_out,
        )
