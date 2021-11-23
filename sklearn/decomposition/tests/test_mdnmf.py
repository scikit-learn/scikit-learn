import copy
import numpy as np
import pytest

from sklearn.decomposition import MDNMF
from sklearn.decomposition._mdnmf import (
    calc_nn_graph_laplacian,
    fit_coordinate_descent,
    fit_multiplicative_update,
    mdnmf_loss,
)
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._testing import assert_array_equal, assert_array_almost_equal


@pytest.mark.parametrize("solver", ("cd", "mu"))
def test_convergence_warning(solver):
    convergence_warning = (
        "Maximum number of iterations 1 reached. Increase it to improve convergence."
    )
    A = np.ones((2, 2))
    with pytest.warns(ConvergenceWarning, match=convergence_warning):
        MDNMF(solver=solver, max_iter=1).fit(A)


def test_parameter_checking():
    A = np.ones((2, 2))
    name = "spam"
    msg = "Invalid solver parameter: got 'spam' instead of one of"
    with pytest.raises(ValueError, match=msg):
        MDNMF(solver=name).fit(A)

    msg = "Negative values in data passed to"
    with pytest.raises(ValueError, match=msg):
        MDNMF().fit(-A)
    clf = MDNMF(2, tol=0.1).fit(A)
    with pytest.raises(ValueError, match=msg):
        clf.transform(-A)


def test_feature_names_out():
    random_state = np.random.RandomState(0)
    X = np.abs(random_state.randn(10, 4))
    mdnmf = MDNMF(n_components=3).fit(X)

    names = mdnmf.get_feature_names_out()
    assert_array_equal([f"mdnmf{i}" for i in range(3)], names)


@pytest.mark.parametrize("solver", ("cd", "mu"))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_mdnmf_fit_nn_output(solver, alpha, beta, gamma):
    # Test that the decomposition does not contain negative values
    A = np.c_[5.0 - np.arange(1, 6), 5.0 + np.arange(1, 6)]
    model = MDNMF(
        n_components=2,
        solver=solver,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        random_state=0,
    )
    transf = model.fit_transform(A)
    assert not ((model.components_ < 0).any() or (transf < 0).any())


def test_n_components_greater_n_features():
    # Smoke test for the case of more components than features.
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(30, 10))
    y = rng.randint(0, 2, 30)
    X_y = np.c_[X, y]
    MDNMF(n_components=15, random_state=0, tol=1e-2).fit(X_y)


@pytest.mark.parametrize("dist_unit", (0.1, 1, 10))
def test_nn_graph_calc_one_class(dist_unit):
    # Given
    X = (
        np.array(
            [
                [1, 0, 0],
                [2, 0, 0],
                [3, 0, 0],
                [0, 1, 0],
                [0, 2, 0],
            ]
        )
        * dist_unit
    )
    y = [0, 0, 0, 0, 0]

    # When
    Sg, Dg, Lg, Sc, Dc, Lc = calc_nn_graph_laplacian(X, y, 2, 2, return_all=True)

    # Then
    assert_array_equal(
        Sg,
        np.array(
            [
                [0, 1, 0, 1, 0],
                [1, 0, 1, 0, 0],
                [1, 1, 0, 0, 0],
                [1, 0, 0, 0, 1],
                [1, 0, 0, 1, 0],
            ]
        ),
    )
    assert_array_equal(Dg, np.diag(np.ones(5) * 2))
    assert_array_almost_equal(
        Lg,
        np.array(
            [
                [2, -1, 0, -1, 0],
                [-1, 2, -1, 0, 0],
                [-1, -1, 2, 0, 0],
                [-1, 0, 0, 2, -1],
                [-1, 0, 0, -1, 2],
            ]
        ),
        decimal=2,
    )
    assert_array_equal(Sc, np.zeros((5, 5)))
    assert_array_equal(Dc, np.zeros((5, 5)))
    assert_array_almost_equal(Lc, np.zeros((5, 5)), decimal=2)


@pytest.mark.parametrize("dist_unit", (0.1, 1, 10))
def test_nn_graph_calc_two_classes(dist_unit):
    # Given
    X = (
        np.array(
            [
                [10, 0, 0],
                [0, 10, 0],
                [0, 0, 10],
                [5, 0, 0],
                [0, 5, 0],
                [0, 0, 4],
                [7, 0, 0],
                [0, 0, 3],
            ]
        )
        * dist_unit
    )
    y = [0, 0, 0, 0, 1, 0, 0, 1]

    # When
    Sg, Dg, Lg, Sc, Dc, Lc = calc_nn_graph_laplacian(X, y, 2, 2, return_all=True)

    # Then
    assert_array_equal(
        Sg,
        np.array(
            [
                [0, 0, 0, 1, 0, 0, 1, 0],
                [0, 0, 0, 1, 0, 1, 0, 0],
                [0, 0, 0, 1, 0, 1, 0, 0],
                [1, 0, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 0, 0, 0, 1],
                [0, 0, 1, 1, 0, 0, 0, 0],
                [1, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
            ]
        ),
    )
    assert_array_equal(
        Dg,
        np.array(
            [
                [2, 0, 0, 0, 0, 0, 0, 0],
                [0, 2, 0, 0, 0, 0, 0, 0],
                [0, 0, 2, 0, 0, 0, 0, 0],
                [0, 0, 0, 2, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 0, 0, 2, 0, 0],
                [0, 0, 0, 0, 0, 0, 2, 0],
                [0, 0, 0, 0, 0, 0, 0, 1],
            ]
        ),
    )
    assert_array_almost_equal(
        Lg,
        np.array(
            [
                [2, 0, 0, -1, 0, 0, -1, 0],
                [0, 2, 0, -1, 0, -1, 0, 0],
                [0, 0, 2, -1, 0, -1, 0, 0],
                [-1, 0, 0, 2, 0, 0, -1, 0],
                [0, 0, 0, 0, 1, 0, 0, -1],
                [0, 0, -1, -1, 0, 2, 0, 0],
                [-1, 0, 0, -1, 0, 0, 2, 0],
                [0, 0, 0, 0, -1, 0, 0, 1],
            ]
        ),
        decimal=2,
    )
    assert_array_equal(
        Sc,
        np.array(
            [
                [0, 0, 0, 0, 1, 0, 0, 1],
                [0, 0, 0, 0, 1, 0, 0, 1],
                [0, 0, 0, 0, 1, 0, 0, 1],
                [0, 0, 0, 0, 1, 0, 0, 1],
                [0, 1, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 1],
                [0, 0, 0, 0, 1, 0, 0, 1],
                [0, 0, 0, 1, 0, 1, 0, 0],
            ]
        ),
    )
    assert_array_equal(Dc, np.diag(np.ones(8) * 2))
    assert_array_almost_equal(
        Lc,
        np.array(
            [
                [2, 0, 0, 0, -1, 0, 0, -1],
                [0, 2, 0, 0, -1, 0, 0, -1],
                [0, 0, 2, 0, -1, 0, 0, -1],
                [0, 0, 0, 2, -1, 0, 0, -1],
                [0, -1, 0, 0, 2, -1, 0, 0],
                [0, 0, 0, 0, -1, 2, 0, -1],
                [0, 0, 0, 0, -1, 0, 2, -1],
                [0, 0, 0, -1, 0, -1, 0, 2],
            ]
        ),
        decimal=2,
    )


def test_nn_graph_calc_no_labels():
    # Given
    X = np.random.randn(3, 5)

    # When
    Sg1, Dg1, Lg1, Lc1 = calc_nn_graph_laplacian(X, None, 2, 2)
    Sg2, Dg2, Lg2, Sc2, Dc2, Lc2 = calc_nn_graph_laplacian(
        X, None, 2, 2, return_all=True
    )

    # Then
    assert_array_equal(Sg1, np.zeros((3, 3)))
    assert_array_equal(Dg1, np.zeros((3, 3)))
    assert_array_equal(Lg1, np.zeros((3, 3)))
    assert_array_equal(Lc1, np.zeros((3, 3)))
    assert_array_equal(Sg2, np.zeros((3, 3)))
    assert_array_equal(Dg2, np.zeros((3, 3)))
    assert_array_equal(Lg2, np.zeros((3, 3)))
    assert_array_equal(Sc2, np.zeros((3, 3)))
    assert_array_equal(Dc2, np.zeros((3, 3)))
    assert_array_equal(Lc2, np.zeros((3, 3)))


@pytest.mark.parametrize("alpha", (0.001, 0.01, 0.1))
@pytest.mark.parametrize("beta", (0.01, 0.1, 1.0))
@pytest.mark.parametrize("gamma", (10, 100, 1000))
def test_mdnmf_loss_smaller_for_better_reconstruction(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 10))
    H_rand = np.abs(rng.randn(10, 20))
    H_calc = np.linalg.inv(W) @ X  # H with low reconstruction error
    _, _, Lg, Lc = calc_nn_graph_laplacian(X, y, 2, 2)
    e = np.ones((10, 10)) - np.eye(10)

    # When
    loss_rand = mdnmf_loss(X.T, H_rand.T, W.T, Lg, Lc, e, alpha, beta, gamma)
    loss_calc = mdnmf_loss(X.T, H_calc.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # Then
    assert loss_rand > loss_calc


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_mdnmf_loss_no_neighbor(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    W = np.abs(rng.randn(10, 10))
    H_rand = np.abs(rng.randn(10, 20))
    H_calc = np.linalg.inv(W) @ X  # H with low reconstruction error
    Lg = Lc = np.zeros((10, 10))
    e = np.ones((10, 10)) - np.eye(10)

    # When
    loss_rand = mdnmf_loss(X.T, H_rand.T, W.T, Lg, Lc, e, alpha, beta, gamma)
    loss_calc = mdnmf_loss(X.T, H_calc.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # Then
    assert loss_rand > loss_calc


@pytest.mark.parametrize("alpha", (0, 0.001, 0.01, 0.1))
@pytest.mark.parametrize("beta", (0, 0.01, 0.1, 1.0))
@pytest.mark.parametrize("gamma", (0, 10, 100, 1000))
def test_fit_mu_update_w_h(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_multiplicative_update(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=20,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        verbose=1,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_mu_update_w_h_no_neighbor(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg = Dg = Lg = Lc = np.zeros((10, 10))
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_multiplicative_update(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=20,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_mu_fix_w(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)

    X_fit = X.T
    W_original = np.abs(rng.randn(20, 5))
    W_fit = copy.deepcopy(W_original)
    H_fit = np.abs(rng.randn(5, 10))
    E = np.ones((20, 10))
    e = np.ones((5, 5)) - np.eye(5)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_multiplicative_update(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=20,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_W=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(W_fit, W_original)


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_mu_fix_h(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(12, 16))
    y = rng.randint(0, 3, 12)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 4, 2)

    X_fit = X.T
    W_fit = np.abs(rng.randn(16, 7))
    H_original = np.abs(rng.randn(7, 12))
    H_fit = copy.deepcopy(H_original)
    E = np.ones((16, 12))
    e = np.ones((7, 7)) - np.eye(7)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_multiplicative_update(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=20,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_H=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(H_fit, H_original)


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_mu_fix_h_no_neighbor(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(12, 16))
    Sg = Dg = Lg = Lc = np.zeros((12, 12))

    X_fit = X.T
    W_fit = np.abs(rng.randn(16, 7))
    H_original = np.abs(rng.randn(7, 12))
    H_fit = copy.deepcopy(H_original)
    E = np.ones((16, 12))
    e = np.ones((7, 7)) - np.eye(7)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_multiplicative_update(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=20,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_H=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(H_fit, H_original)


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_mu_fix_w_h(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(14, 9))
    y = rng.randint(0, 5, 14)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 1, 5)

    X_fit = X.T
    W_original = np.abs(rng.randn(9, 5))
    W_fit = copy.deepcopy(W_original)
    H_original = np.abs(rng.randn(5, 14))
    H_fit = copy.deepcopy(H_original)
    E = np.ones((9, 14))
    e = np.ones((5, 5)) - np.eye(5)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_multiplicative_update(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=20,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_W=False,
        update_H=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 == loss_0
    assert_array_equal(W_fit, W_original)
    assert_array_equal(H_fit, H_original)


@pytest.mark.parametrize("omega", (0.1, 0.3))
@pytest.mark.parametrize("alpha", (0, 0.01, 0.1, 1))
@pytest.mark.parametrize("beta", (0, 0.1, 1))
@pytest.mark.parametrize("gamma", (0, 10, 100))
def test_fit_cd_update_w_h_sor_newton(omega, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=5,
        optimizer_config=dict(omega=omega),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("omega", (0.1,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_update_w_h_sor_newton_no_neighbor(omega, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg = Dg = Lg = Lc = np.zeros((10, 10))
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=5,
        optimizer_config=dict(omega=omega),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("omega", (0.1,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_fix_w_sor_newton(omega, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)

    X_fit = X.T
    W_original = np.abs(rng.randn(20, 5))
    W_fit = copy.deepcopy(W_original)
    H_fit = np.abs(rng.randn(5, 10))
    E = np.ones((20, 10))
    e = np.ones((5, 5)) - np.eye(5)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_coordinate_descent(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=5,
        optimizer_config=dict(omega=omega),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_W=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(W_fit, W_original)


@pytest.mark.parametrize("omega", (0.1,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_fix_w_sor_newton_no_neighbor(omega, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    Sg = Dg = Lg = Lc = np.zeros((10, 10))

    X_fit = X.T
    W_original = np.abs(rng.randn(20, 5))
    W_fit = copy.deepcopy(W_original)
    H_fit = np.abs(rng.randn(5, 10))
    E = np.ones((20, 10))
    e = np.ones((5, 5)) - np.eye(5)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_coordinate_descent(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=5,
        optimizer_config=dict(omega=omega),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_W=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(W_fit, W_original)


@pytest.mark.parametrize("omega", (0.1,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_fix_h_sor_newton(omega, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(12, 16))
    y = rng.randint(0, 3, 12)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 4, 2)

    X_fit = X.T
    W_fit = np.abs(rng.randn(16, 7))
    H_original = np.abs(rng.randn(7, 12))
    H_fit = copy.deepcopy(H_original)
    E = np.ones((16, 12))
    e = np.ones((7, 7)) - np.eye(7)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_coordinate_descent(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=5,
        optimizer_config=dict(omega=omega),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_H=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(H_fit, H_original)


@pytest.mark.parametrize("omega", (0.1,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_fix_h_sor_newton_no_neighbor(omega, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(12, 16))
    Sg = Dg = Lg = Lc = np.zeros((12, 12))

    X_fit = X.T
    W_fit = np.abs(rng.randn(16, 7))
    H_original = np.abs(rng.randn(7, 12))
    H_fit = copy.deepcopy(H_original)
    E = np.ones((16, 12))
    e = np.ones((7, 7)) - np.eye(7)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_coordinate_descent(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=5,
        optimizer_config=dict(omega=omega),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_H=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert_array_equal(H_fit, H_original)


@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_fix_w_h(alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(14, 9))
    y = rng.randint(0, 5, 14)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 1, 5)

    X_fit = X.T
    W_original = np.abs(rng.randn(9, 5))
    W_fit = copy.deepcopy(W_original)
    H_original = np.abs(rng.randn(5, 14))
    H_fit = copy.deepcopy(H_original)
    E = np.ones((9, 14))
    e = np.ones((5, 5)) - np.eye(5)
    loss_0 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)

    # When
    W_fit, H_fit, _ = fit_coordinate_descent(
        X=X_fit,
        W=W_fit,
        H=H_fit,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E,
        e=e,
        max_iter=5,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        update_W=False,
        update_H=False,
    )

    # Then
    loss_1 = mdnmf_loss(X_fit, W_fit, H_fit, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 == loss_0
    assert_array_equal(W_fit, W_original)
    assert_array_equal(H_fit, H_original)


@pytest.mark.parametrize("opt_gamma", (0.9,))
@pytest.mark.parametrize("eta", (0.1,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (100,))
def test_fit_cd_momentum(opt_gamma, eta, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=20,
        optimizer_config=dict(name="momentum", gamma=opt_gamma, eta=eta),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("opt_gamma", (0.9,))
@pytest.mark.parametrize("eta", (0.01,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (100,))
def test_fit_cd_nesterov(opt_gamma, eta, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=20,
        optimizer_config=dict(name="nesterov", gamma=opt_gamma, eta=eta),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("beta_1", (0.9,))
@pytest.mark.parametrize("beta_2", (0.999,))
@pytest.mark.parametrize("epsilon", (1e-8,))
@pytest.mark.parametrize("eta", (0.01,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_adam(beta_1, beta_2, epsilon, eta, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=5,
        optimizer_config=dict(
            name="adam", beta_1=beta_1, beta_2=beta_2, epsilon=epsilon, eta=eta
        ),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("beta_1", (0.9,))
@pytest.mark.parametrize("beta_2", (0.999,))
@pytest.mark.parametrize("epsilon", (1e-8,))
@pytest.mark.parametrize("eta", (0.01,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_adamax(beta_1, beta_2, epsilon, eta, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=5,
        optimizer_config=dict(
            name="adamax", beta_1=beta_1, beta_2=beta_2, epsilon=epsilon, eta=eta
        ),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("beta_1", (0.9,))
@pytest.mark.parametrize("beta_2", (0.999,))
@pytest.mark.parametrize("epsilon", (1e-8,))
@pytest.mark.parametrize("eta", (0.01,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_nadam(beta_1, beta_2, epsilon, eta, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=5,
        optimizer_config=dict(
            name="nadam", beta_1=beta_1, beta_2=beta_2, epsilon=epsilon, eta=eta
        ),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("beta_1", (0.9,))
@pytest.mark.parametrize("beta_2", (0.999,))
@pytest.mark.parametrize("epsilon", (1e-8,))
@pytest.mark.parametrize("eta", (0.01,))
@pytest.mark.parametrize("alpha", (0, 0.01))
@pytest.mark.parametrize("beta", (0, 0.1))
@pytest.mark.parametrize("gamma", (0, 100))
def test_fit_cd_amsgrad(beta_1, beta_2, epsilon, eta, alpha, beta, gamma):
    # Given
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(10, 20))
    y = rng.randint(0, 2, 10)
    W = np.abs(rng.randn(10, 5))
    H = np.abs(rng.randn(5, 20))
    E = np.ones((10, 20))
    e = np.ones((5, 5)) - np.eye(5)
    Sg, Dg, Lg, Lc = calc_nn_graph_laplacian(X, y, 3, 3)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    W, H, _ = fit_coordinate_descent(
        X=X.T,
        W=H.T,
        H=W.T,
        Sg=Sg,
        Dg=Dg,
        Lg=Lg,
        Lc=Lc,
        E=E.T,
        e=e,
        max_iter=5,
        optimizer_config=dict(
            name="amsgrad", beta_1=beta_1, beta_2=beta_2, epsilon=epsilon, eta=eta
        ),
        alpha=alpha,
        beta=beta,
        gamma=gamma,
    )

    # Then
    loss_1 = mdnmf_loss(X.T, W, H, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0


@pytest.mark.parametrize("alpha", (0, 0.001, 0.01, 0.1))
@pytest.mark.parametrize("beta", (0, 0.01, 0.1, 1.0))
@pytest.mark.parametrize("gamma", (0, 10, 100, 1000))
def test_mdnmf_fit_transform_mu(alpha, beta, gamma):
    # Given
    n_samples = 10
    n_features = 20
    n_components = 5
    n_classes = 3
    random_state = 42
    k1 = 3
    k2 = 3

    rng = np.random.mtrand.RandomState(random_state)
    X = np.abs(rng.randn(n_samples, n_features))
    y = rng.randint(0, n_classes - 1, n_samples)
    X_y = np.c_[X, y]
    W = np.abs(rng.randn(n_samples, n_components))
    H = np.abs(rng.randn(n_components, n_features))
    e = np.ones((n_components, n_components)) - np.eye(n_components)
    _, _, Lg, Lc = calc_nn_graph_laplacian(X, y, k1, k2)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    mdnmf = MDNMF(
        n_components=n_components,
        solver="mu",
        max_iter=5,
        random_state=random_state,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        k1=k1,
        k2=k2,
    )
    W_ft = mdnmf.fit_transform(X_y, None, W, H)
    H_ft = mdnmf.components_
    mdnmf = MDNMF(
        n_components=n_components,
        solver="mu",
        max_iter=5,
        random_state=random_state,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        k1=k1,
        k2=k2,
    ).fit(X_y, None, W, H)
    W_f_and_t = mdnmf.transform(X_y)
    H_f_and_t = mdnmf.components_

    # Then
    loss_1 = mdnmf_loss(X.T, H_ft.T, W_ft.T, Lg, Lc, e, alpha, beta, gamma)
    loss_2 = mdnmf_loss(X.T, H_f_and_t.T, W_f_and_t.T, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert loss_2 < loss_0
    assert_array_equal(W_ft, W_f_and_t)


@pytest.mark.parametrize("alpha", (0, 0.001, 0.01, 0.1))
@pytest.mark.parametrize("beta", (0, 0.01, 0.1, 1.0))
@pytest.mark.parametrize("gamma", (0, 10, 100, 1000))
def test_mdnmf_fit_transform_cd(alpha, beta, gamma):
    # Given
    n_samples = 15
    n_features = 25
    n_components = 6
    n_classes = 4
    random_state = 42
    k1 = 4
    k2 = 2

    rng = np.random.mtrand.RandomState(random_state)
    X = np.abs(rng.randn(n_samples, n_features))
    y = rng.randint(0, n_classes - 1, n_samples)
    X_y = np.c_[X, y]
    W = np.abs(rng.randn(n_samples, n_components))
    H = np.abs(rng.randn(n_components, n_features))
    e = np.ones((n_components, n_components)) - np.eye(n_components)
    _, _, Lg, Lc = calc_nn_graph_laplacian(X, y, k1, k2)
    loss_0 = mdnmf_loss(X.T, H.T, W.T, Lg, Lc, e, alpha, beta, gamma)

    # When
    mdnmf = MDNMF(
        n_components=n_components,
        solver="cd",
        max_iter=5,
        random_state=random_state,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        k1=k1,
        k2=k2,
    )
    W_ft = mdnmf.fit_transform(X_y, None, W, H)
    H_ft = mdnmf.components_
    mdnmf = MDNMF(
        n_components=n_components,
        solver="cd",
        max_iter=5,
        random_state=random_state,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        k1=k1,
        k2=k2,
    ).fit(X_y, None, W, H)
    W_f_and_t = mdnmf.transform(X_y)
    H_f_and_t = mdnmf.components_

    # Then
    loss_1 = mdnmf_loss(X.T, H_ft.T, W_ft.T, Lg, Lc, e, alpha, beta, gamma)
    loss_2 = mdnmf_loss(X.T, H_f_and_t.T, W_f_and_t.T, Lg, Lc, e, alpha, beta, gamma)
    assert loss_1 < loss_0
    assert loss_2 < loss_0
    assert_array_equal(W_ft, W_f_and_t)
