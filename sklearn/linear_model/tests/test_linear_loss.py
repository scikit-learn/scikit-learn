"""
Tests for LinearModelLoss

Note that correctness of losses (which compose LinearModelLoss) is already well
covered in the _loss module.
"""
import pytest
import numpy as np
from numpy.testing import assert_allclose
from scipy import linalg, optimize, sparse
from scipy.special import logit

from sklearn._loss.loss import (
    HalfBinomialLoss,
    HalfMultinomialLoss,
    HalfPoissonLoss,
)
from sklearn.datasets import make_low_rank_matrix
from sklearn.linear_model._linear_loss import (
    LinearModelLoss,
    Multinomial_LDL_Decomposition,
)
from sklearn.utils.extmath import squared_norm


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
def test_init_zero_coef(base_loss, fit_intercept, n_features, dtype):
    """Test that init_zero_coef initializes coef correctly."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    rng = np.random.RandomState(42)
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
def test_loss_grad_hess_are_the_same(
    base_loss, fit_intercept, sample_weight, l2_reg_strength
):
    """Test that loss and gradient are the same across different functions."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss, n_samples=10, n_features=5, seed=42
    )

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
    if not base_loss.is_multiclass:
        g4, h4, _ = loss.gradient_hessian(
            coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
        )
    else:
        with pytest.raises(NotImplementedError):
            loss.gradient_hessian(
                coef,
                X,
                y,
                sample_weight=sample_weight,
                l2_reg_strength=l2_reg_strength,
            )

    assert_allclose(l1, l2)
    assert_allclose(g1, g2)
    assert_allclose(g1, g3)
    if not base_loss.is_multiclass:
        assert_allclose(g1, g4)
        assert_allclose(h4 @ g4, h3(g3))

    # same for sparse X
    X = sparse.csr_matrix(X)
    l1_sp = loss.loss(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g1_sp = loss.gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    l2_sp, g2_sp = loss.loss_gradient(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    g3_sp, h3_sp = loss.gradient_hessian_product(
        coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
    )
    if not base_loss.is_multiclass:
        g4_sp, h4_sp, _ = loss.gradient_hessian(
            coef, X, y, sample_weight=sample_weight, l2_reg_strength=l2_reg_strength
        )

    assert_allclose(l1, l1_sp)
    assert_allclose(l1, l2_sp)
    assert_allclose(g1, g1_sp)
    assert_allclose(g1, g2_sp)
    assert_allclose(g1, g3_sp)
    assert_allclose(h3(g1), h3_sp(g1_sp))
    if not base_loss.is_multiclass:
        assert_allclose(g1, g4_sp)
        assert_allclose(h4 @ g4, h4_sp @ g1_sp)


@pytest.mark.parametrize("base_loss", LOSSES)
@pytest.mark.parametrize("sample_weight", [None, "range"])
@pytest.mark.parametrize("l2_reg_strength", [0, 1])
@pytest.mark.parametrize("X_sparse", [False, True])
def test_loss_gradients_hessp_intercept(
    base_loss, sample_weight, l2_reg_strength, X_sparse
):
    """Test that loss and gradient handle intercept correctly."""
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=False)
    loss_inter = LinearModelLoss(base_loss=base_loss(), fit_intercept=True)
    n_samples, n_features = 10, 5
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss, n_samples=n_samples, n_features=n_features, seed=42
    )

    X[:, -1] = 1  # make last column of 1 to mimic intercept term
    X_inter = X[
        :, :-1
    ]  # exclude intercept column as it is added automatically by loss_inter

    if X_sparse:
        X = sparse.csr_matrix(X)

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

    s = np.random.RandomState(42).randn(*coef.shape)
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
    base_loss, fit_intercept, sample_weight, l2_reg_strength
):
    """Test gradients and hessians with numerical derivatives.

    Gradient should equal the numerical derivatives of the loss function.
    Hessians should equal the numerical derivatives of gradients.
    """
    loss = LinearModelLoss(base_loss=base_loss(), fit_intercept=fit_intercept)
    n_samples, n_features = 10, 5
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss, n_samples=n_samples, n_features=n_features, seed=42
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
def test_multinomial_coef_shape(fit_intercept):
    """Test that multinomial LinearModelLoss respects shape of coef."""
    loss = LinearModelLoss(base_loss=HalfMultinomialLoss(), fit_intercept=fit_intercept)
    n_samples, n_features = 10, 5
    X, y, coef = random_X_y_coef(
        linear_model_loss=loss, n_samples=n_samples, n_features=n_features, seed=42
    )
    s = np.random.RandomState(42).randn(*coef.shape)

    l, g = loss.loss_gradient(coef, X, y)
    g1 = loss.gradient(coef, X, y)
    g2, hessp = loss.gradient_hessian_product(coef, X, y)
    h = hessp(s)
    assert g.shape == coef.shape
    assert h.shape == coef.shape
    assert_allclose(g, g1)
    assert_allclose(g, g2)

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

    assert_allclose(g, g_r.reshape((loss.base_loss.n_classes, -1), order="F"))
    assert_allclose(h, h_r.reshape((loss.base_loss.n_classes, -1), order="F"))


def test_multinomial_identifiability_properties(global_random_seed):
    """Test that multinomial LinearModelLoss fulfills identifiability properties.

    In our symmetrical overdetermined parametrization, we always expect
    np.sum(coef, axis=0) = 0. But if we add the same value cj to all classes, i.e.
    coef[:, j] += cj, certain properties hold:
      1. The unpenalized loss in invariant.
      2. The unpenalized Hessian @ c gives zero.
      3. The unpenalized Gradient @ c gives zero.
    """
    n_samples, n_features, n_classes = 30, 4, 6
    l2_reg_strength = 0
    rng = np.random.RandomState(global_random_seed)
    X = rng.standard_normal((n_samples, n_features))
    y = rng.randint(low=0, high=n_classes, size=(n_samples)).astype(float)
    coef = rng.standard_normal((n_classes, n_features))
    coef -= np.mean(coef, axis=0)
    assert_allclose(np.mean(coef, axis=0), 0, atol=1e-15)

    multinomial_loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=n_classes),
        fit_intercept=False,
    )
    gradient, hessp = multinomial_loss.gradient_hessian_product(
        coef=coef,
        X=X,
        y=y,
        l2_reg_strength=l2_reg_strength,
    )

    # Construct a coef-like array with same values for all classes, e.g.
    # c = [[1, 2, 3],
    #      [1, 2, 3]]
    c = np.tile(rng.standard_normal(n_features), (n_classes, 1))
    assert_allclose(hessp(c), 0, atol=1e-14)
    assert_allclose(gradient.ravel() @ c.ravel(), 0, atol=1e-13)

    loss1 = multinomial_loss.loss(coef=coef, X=X, y=y, l2_reg_strength=l2_reg_strength)
    loss2 = multinomial_loss.loss(
        coef=coef + c, X=X, y=y, l2_reg_strength=l2_reg_strength
    )
    assert_allclose(loss1, loss2)


def test_multinomial_LDL_decomposition_operates_inplace(global_random_seed):
    n_samples, n_classes = 3, 5
    rng = np.random.RandomState(global_random_seed)
    p = rng.uniform(low=0, high=1, size=(n_samples, n_classes))
    p /= np.sum(p, axis=1)[:, None]
    assert_allclose(np.sum(p, axis=1), np.ones(shape=n_samples))

    LDL = Multinomial_LDL_Decomposition(proba=p)
    v = rng.standard_normal(size=(n_samples, n_classes))
    res = LDL.sqrt_D_Lt_matmul(v)
    assert res is v

    v = rng.standard_normal(size=(n_samples, n_classes))
    res = LDL.L_sqrt_D_matmul(v)
    assert res is v

    v = rng.standard_normal(size=(n_samples, n_classes))
    res = LDL.inverse_L_sqrt_D_matmul(v)
    assert res is v


def test_multinomial_LDL_decomposition_binomial_single_point():
    """Test LDL' decomposition of multinomial hessian for simple cases.

    For the binomial case, we have p0 = 1 - p1
    LDL = [p0 * (1 - p0),      -p0 * p1] = p0 * (1 - p0) * [1, -1]
          [     -p0 * p1, p1 * (1 - p1)]                   [-1, 1]

    L = [ 1, 0]    D = [p0 * (1 - p0), 0]
        [-1, 1]        [            0, 0]
    """
    p0, p1 = 0.2, 0.8
    p = np.array([[p0, p1]])
    # Note that LDL sets q = 1 whenever q = 0 for easier handling of divisions by zero.
    # We compare those values to "zero", to make that clear.
    zero = 1

    LDL = Multinomial_LDL_Decomposition(proba=p)
    assert_allclose(1 / LDL.q_inv[0, :], [1 - p0, zero])
    assert_allclose(LDL.sqrt_d[0, :] ** 2, [p0 * (1 - p0), 0])

    # C = diag(D) L' x with x = [[1, 0]] (n_samples=1, n_classes=2)
    C = LDL.sqrt_D_Lt_matmul(np.array([[1.0, 0.0]]))
    assert_allclose(C[0, :], [np.sqrt(p0 * (1 - p0)), 0])

    # C = diag(D) L' x with x = [[0, 1]] (n_samples=1, n_classes=2)
    C = LDL.sqrt_D_Lt_matmul(np.array([[0.0, 1.0]]))
    assert_allclose(C[0, :], [-np.sqrt(p1 * (1 - p1)), 0])

    # Test with hessian product (hessp) of LinearModelLoss with X = [[1]] and coef such
    # that probabilities are again (p0, p1).
    # Hessian = X' LDL X
    loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=2),
        fit_intercept=False,
    )
    # Note that y has no effect on the hessian.
    coef = 0.5 * np.array([[logit(p0)], [logit(p1)]])  # tested below
    X = np.array([[1.0]])
    raw = X @ coef.T  # raw.shape = (n_samples, n_classes) = (1, 2)
    assert_allclose(loss.base_loss.predict_proba(raw), [[p0, p1]])
    C = LDL.sqrt_D_Lt_matmul(raw.copy())
    grad, hessp = loss.gradient_hessian_product(
        coef=coef,
        X=X,
        y=np.array([0.0]),
        l2_reg_strength=0.0,
    )
    # Note: hessp(coef).shape = (n_classes, n_features) = (2, 1)
    assert_allclose(
        hessp(coef),
        p0
        * (1 - p0)
        * np.array([raw[0, 0] - raw[0, 1], -raw[0, 0] + raw[0, 1]])[:, None],
    )
    assert_allclose(C @ C.T, coef.T @ hessp(coef))


def test_multinomial_LDL_decomposition_3_classes():
    """Test LDL' decomposition of multinomial hessian for 3 classes and 2 points.

    For n_classes = 3 and n_samples = 2, we have
      p0 = [p0_0, p0_1]
      p1 = [p1_0, p1_1]
      p2 = [p2_0, p2_1]
    and with 2 x 2 diagonal subblocks
      LDL = [p0 * (1-p0),    -p0 * p1,    -p0 * p2]
            [   -p0 * p1, p1 * (1-p1),    -p1 * p2]
            [   -p0 * p2,    -p1 * p2, p2 * (1-p2)]

      L = [         1,                 0, 0]
          [-p1 / (1-p0),               1, 0]
          [-p2 / (1-p0), -p2 / (1-p0-p1), 1]

      D = [p0 * (1-p0),                       0, 0]
          [          0, p1 * (1-p0-p1) / (1-p0), 0]
          [          0,                       0, 0]
    """
    n = 2 * 3  # n_samples * n_classes
    p0 = np.array([0.7, 0.6])
    p1 = np.array([0.2, 0.25])
    p2 = np.array([0.1, 0.15])
    one = np.ones(2)
    zero = np.zeros(2)
    p0d, p1d, p2d, oned, zerod = (
        np.diag(p0),
        np.diag(p1),
        np.diag(p2),
        np.diag(one),
        np.diag(zero),
    )
    H = np.block(
        [
            [p0d * (oned - p0d), -p0d * p1d, -p0d * p2d],
            [-p0d * p1d, p1d * (oned - p1d), -p1d * p2d],
            [-p0d * p2d, -p1d * p2d, p2d * (oned - p2d)],
        ]
    )
    L = np.block(
        [
            [oned, zerod, zerod],
            [np.diag(-p1 / (one - p0)), oned, zerod],
            [np.diag(-p2 / (one - p0)), np.diag(-p2 / (one - p0 - p1)), oned],
        ]
    )
    D = np.diag(np.r_[p0 * (1 - p0), p1 * (1 - p0 - p1) / (1 - p0), zero])
    L_inv = np.block(
        [
            [oned, zerod, zerod],
            [np.diag(p1 / (one - p0)), oned, zerod],
            [np.diag(p2 / (one - p0 - p1)), np.diag(p2 / (one - p0 - p1)), oned],
        ]
    )
    D_inv = np.zeros_like(D)
    D_inv[D > 0] = 1 / np.sqrt(D[D > 0])

    lu, d, perm = linalg.ldl(H)
    assert_allclose(lu @ d @ lu.T, H)
    assert_allclose(L @ D @ L.T, H)
    assert_allclose(L, lu)
    assert_allclose(D, d, atol=1e-14)
    assert_allclose(L_inv @ L, np.eye(n), atol=1e-14)
    assert_allclose(L @ L_inv, np.eye(n), atol=1e-14)

    LDL = Multinomial_LDL_Decomposition(proba=np.c_[p0, p1, p2])
    v = 1.0 + np.arange(2 * 3)
    DLt_v = LDL.sqrt_D_Lt_matmul(v.copy().reshape((2, 3), order="F"))
    assert_allclose(DLt_v.ravel(order="F"), np.sqrt(D) @ L.T @ v)
    LDL_v = LDL.L_sqrt_D_matmul(DLt_v.copy())
    assert_allclose(LDL_v.ravel(order="F"), H @ v)

    LD_v = LDL.L_sqrt_D_matmul(v.copy().reshape((2, 3), order="F"))
    assert_allclose(LD_v.ravel(order="F"), L @ np.sqrt(D) @ v)

    x1 = LDL.inverse_L_sqrt_D_matmul(v.copy().reshape((2, 3), order="F"))
    # This does not work as L @ D is singular.
    #   x2 = linalg.solve(L @ np.sqrt(D), v)
    # We can just neglect the last class as it is redundant.
    x2 = linalg.solve(L[:-2, :-2] @ np.sqrt(D[:-2, :-2]), v[:-2])
    assert_allclose(x1.ravel(order="F")[:-2], x2)
    assert_allclose(x1.ravel(order="F"), D_inv @ L_inv @ v)

    # Test consistency of L_sqrt_D_matmul, sqrt_D_Lt_matmul and inverse_L_sqrt_D_matmul,
    # i.e. reconstructing the matrices based on X @ unit_vector and compare with
    # explicit matrices D and L.
    unit_vec = np.zeros(n)
    for op, result in (
        (LDL.sqrt_D_Lt_matmul, np.sqrt(D) @ L.T),
        (LDL.L_sqrt_D_matmul, L @ np.sqrt(D)),
        (LDL.inverse_L_sqrt_D_matmul, D_inv @ L_inv),
    ):
        X = np.zeros((n, n))
        for i in range(n):
            unit_vec[i] = 1
            X[:, i] = op(unit_vec.copy().reshape((2, 3), order="F")).ravel(order="F")
            unit_vec[i] = 0
        assert_allclose(X, result)


def test_multinomial_LDL_decomposition(global_random_seed):
    """Test LDL' decomposition of multinomial hessian."""
    n_samples, n_features, n_classes = 3, 4, 5
    rng = np.random.RandomState(global_random_seed)
    X = rng.standard_normal(size=(n_samples, n_features))
    coef = rng.standard_normal(size=(n_classes, n_features))
    raw_prediction = X @ coef.T  # shape = (n_samples, n_classes)
    loss = LinearModelLoss(
        base_loss=HalfMultinomialLoss(n_classes=n_classes),
        fit_intercept=False,
    )
    p = loss.base_loss.predict_proba(raw_prediction=raw_prediction)
    assert_allclose(np.sum(p, axis=1), np.ones(shape=n_samples))

    LDL = Multinomial_LDL_Decomposition(proba=p)
    # Note that LDL sets q = 1 whenever q = 0 for easier handling of divisions by zero.
    # As q[:, -1] = 0 if p sums to 1, we do not compare the last values corresponding
    # to the last class.
    assert_allclose(1 / LDL.q_inv[:, :-1], 1 - np.cumsum(LDL.p, axis=1)[:, :-1])

    # C = diag(D) L' x with x = X @ coef = raw_prediction
    C = LDL.sqrt_D_Lt_matmul(raw_prediction.copy())

    # Note that y has no effect on the hessian.
    grad, hessp = loss.gradient_hessian_product(
        coef=coef,
        X=X,
        y=rng.randint(low=0, high=n_classes, size=n_samples).astype(X.dtype),
        l2_reg_strength=0.0,
    )
    # C.shape = (n_samples, n_classes), hessp(coef).shape = (n_classes, n_features)
    assert_allclose(np.tensordot(C, C, axes=2), np.tensordot(coef, hessp(coef), axes=2))
    CtC = LDL.L_sqrt_D_matmul(C.copy())
    CtC = np.tensordot(raw_prediction, CtC, axes=2)
    assert_allclose(CtC, np.tensordot(C, C, axes=2))


def test_multinomial_LDL_inverse_sqrt_D_Lt_matmul(global_random_seed):
    """Test inverse of LDL' decomposition of multinomial hessian."""
    n_samples, n_classes = 4, 5
    rng = np.random.RandomState(global_random_seed)
    p = rng.uniform(low=0, high=1, size=(n_samples, n_classes))
    p /= np.sum(p, axis=1)[:, None]
    LDL = Multinomial_LDL_Decomposition(proba=p)
    x = rng.standard_normal(size=(n_samples, n_classes))
    # Without centering, the last class, x[:, -1] would fail in the assert below.
    x -= np.mean(x, axis=1)[:, None]

    invDL_x = LDL.inverse_L_sqrt_D_matmul(x.copy())
    LD_invDL_x = LDL.L_sqrt_D_matmul(invDL_x.copy())
    assert_allclose(LD_invDL_x, x)


def test_multinomial_LDL_warning():
    """Test that LDL warns if probabilities do not sum to 1."""
    p = np.arange(4, dtype=float).reshape((2, 2))
    msg = "Probabilities proba are assumed to sum to 1, but they don't."
    with pytest.warns(UserWarning, match=msg):
        Multinomial_LDL_Decomposition(proba=p, proba_sum_to_1=True)
