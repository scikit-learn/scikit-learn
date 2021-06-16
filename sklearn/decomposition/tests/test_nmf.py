import re

import numpy as np
import scipy.sparse as sp

from scipy import linalg
from sklearn.decomposition import NMF, MiniBatchNMF
from sklearn.decomposition import non_negative_factorization
from sklearn.decomposition import _nmf as nmf  # For testing internals
from scipy.sparse import csc_matrix

import pytest

from sklearn.utils._testing import assert_array_equal
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils._testing import assert_almost_equal
from sklearn.utils._testing import assert_allclose
from sklearn.utils._testing import ignore_warnings
from sklearn.utils.extmath import squared_norm
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize('regularization',
                         [None, 'both', 'components', 'transformation'])
def test_convergence_warning(Estimator, solver, regularization):
    convergence_warning = ("Maximum number of iterations 1 reached. "
                           "Increase it to improve convergence.")
    A = np.ones((2, 2))
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    with pytest.warns(ConvergenceWarning, match=convergence_warning):
        Estimator(
            solver=solver, regularization=regularization,
            max_iter=1, init=init
        ).fit(A)


def test_initialize_nn_output():
    # Test that initialization does not return negative values
    rng = np.random.mtrand.RandomState(42)
    data = np.abs(rng.randn(10, 10))
    for init in ('random', 'nndsvd', 'nndsvda', 'nndsvdar'):
        W, H = nmf._initialize_nmf(data, 10, init=init, random_state=0)
        assert not ((W < 0).any() or (H < 0).any())


def test_parameter_checking():
    A = np.ones((2, 2))
    name = 'spam'
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    msg = "Invalid solver parameter: got 'spam' instead of one of"
    with pytest.raises(ValueError, match=msg):
        NMF(solver=name, init=init).fit(A)
    with pytest.raises(ValueError, match=msg):
        MiniBatchNMF(solver=name).fit(A)
    msg = "Invalid init parameter: got 'spam' instead of one of"
    with pytest.raises(ValueError, match=msg):
        NMF(init=name).fit(A)
    msg = "Invalid regularization parameter: got 'spam' instead of one of"
    with pytest.raises(ValueError, match=msg):
        NMF(regularization=name, init=init).fit(A)
    msg = "Invalid beta_loss parameter: got 'spam' instead of one"
    with pytest.raises(ValueError, match=msg):
        NMF(solver='mu', init=init, beta_loss=name).fit(A)
    with pytest.raises(ValueError, match=msg):
        MiniBatchNMF(solver='mu', beta_loss=name).fit(A)
    msg = ("Invalid beta_loss parameter: solver 'cd' does not handle "
           "beta_loss = 1.0")
    with pytest.raises(ValueError, match=msg):
        NMF(solver='cd', init=init, beta_loss=1.0).fit(A)
    msg = ("Invalid solver 'cd' not supported "
           "when batch_size is not None.")
    with pytest.raises(ValueError, match=msg):
        MiniBatchNMF(solver='cd', beta_loss='frobenius').fit(A)
    msg = "Negative values in data passed to"
    with pytest.raises(ValueError, match=msg):
        NMF(init=init).fit(-A)
    with pytest.raises(ValueError, match=msg):
        MiniBatchNMF().fit(-A)
    with pytest.raises(ValueError, match=msg):
        nmf._initialize_nmf(-A, 2, 'nndsvd')
    clf = NMF(2, tol=0.1, init=init).fit(A)
    with pytest.raises(ValueError, match=msg):
        clf.transform(-A)

    for init in ['nndsvd', 'nndsvda', 'nndsvdar']:
        msg = re.escape(
            "init = '{}' can only be used when "
            "n_components <= min(n_samples, n_features)"
            .format(init)
        )
        with pytest.raises(ValueError, match=msg):
            NMF(3, init=init).fit(A)
        with pytest.raises(ValueError, match=msg):
            MiniBatchNMF(3, init=init).fit(A)
        with pytest.raises(ValueError, match=msg):
            nmf._initialize_nmf(A, 3, init)


def test_initialize_close():
    # Test NNDSVD error
    # Test that _initialize_nmf error is less than the standard deviation of
    # the entries in the matrix.
    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(10, 10))
    W, H = nmf._initialize_nmf(A, 10, init='nndsvd')
    error = linalg.norm(np.dot(W, H) - A)
    sdev = linalg.norm(A - A.mean())
    assert error <= sdev


def test_initialize_variants():
    # Test NNDSVD variants correctness
    # Test that the variants 'nndsvda' and 'nndsvdar' differ from basic
    # 'nndsvd' only where the basic version has zeros.
    rng = np.random.mtrand.RandomState(42)
    data = np.abs(rng.randn(10, 10))
    W0, H0 = nmf._initialize_nmf(data, 10, init='nndsvd')
    Wa, Ha = nmf._initialize_nmf(data, 10, init='nndsvda')
    War, Har = nmf._initialize_nmf(data, 10, init='nndsvdar',
                                   random_state=0)

    for ref, evl in ((W0, Wa), (W0, War), (H0, Ha), (H0, Har)):
        assert_almost_equal(evl[ref != 0], ref[ref != 0])


# ignore UserWarning raised when both solver='mu' and init='nndsvd'
@ignore_warnings(category=UserWarning)
@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize('init',
                         (None, 'nndsvd', 'nndsvda', 'nndsvdar', 'random'))
@pytest.mark.parametrize('regularization',
                         (None, 'both', 'components', 'transformation'))
def test_nmf_fit_nn_output(Estimator, solver, init, regularization):
    # Test that the decomposition does not contain negative values
    A = np.c_[5. - np.arange(1, 6),
              5. + np.arange(1, 6)]
    model = Estimator(n_components=2, solver=solver, init=init,
                      regularization=regularization, random_state=0)
    transf = model.fit_transform(A)
    assert not((model.components_ < 0).any() or
               (transf < 0).any())


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize('regularization',
                         (None, 'both', 'components', 'transformation'))
def test_nmf_fit_close(Estimator, solver, regularization):
    rng = np.random.mtrand.RandomState(42)
    # Test that the fit is not too far away
    pnmf = Estimator(5, solver=solver, init='nndsvdar', random_state=0,
                     regularization=regularization, max_iter=600)
    X = np.abs(rng.randn(6, 5))
    assert pnmf.fit(X).reconstruction_err_ < 0.1


@pytest.mark.parametrize('regularization',
                         (None, 'both', 'components', 'transformation'))
def test_nmf_true_reconstruction(regularization):
    # Test that the fit is not too far away from an exact solution
    # (by construction)
    n_samples = 15
    n_components = 5
    n_features = 10
    beta_loss = 1
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    batch_size = 3
    max_iter = 1000

    rng = np.random.mtrand.RandomState(42)
    W_true = np.zeros([n_samples, n_components])
    W_array = np.abs(rng.randn(n_samples))
    for j in range(n_components):
        W_true[j % n_samples, j] = W_array[j % n_samples]
    H_true = np.zeros([n_components, n_features])
    H_array = np.abs(rng.randn(n_components))
    for j in range(n_features):
        H_true[j % n_components, j] = H_array[j % n_components]
    X = np.dot(W_true, H_true)

    model = NMF(n_components=n_components, solver='mu',
                init=init, beta_loss=beta_loss, max_iter=max_iter,
                regularization=regularization, random_state=0)
    transf = model.fit_transform(X)
    X_calc = np.dot(transf, model.components_)

    assert model.reconstruction_err_ < 0.1
    assert_array_almost_equal(X, X_calc)

    mbmodel = MiniBatchNMF(n_components=n_components, solver='mu',
                           init=init, beta_loss=beta_loss,
                           batch_size=batch_size, forget_factor=0.3,
                           regularization=regularization, random_state=0,
                           max_iter=max_iter)
    transf = mbmodel.fit_transform(X)
    X_calc = np.dot(transf, mbmodel.components_)

    assert mbmodel.reconstruction_err_ < 0.1
    assert_array_almost_equal(X, X_calc, decimal=1)


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize('regularization',
                         (None, 'both', 'components', 'transformation'))
def test_nmf_transform(Estimator, solver, regularization):
    # Test that NMF.transform returns close values
    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(6, 5))
    m = Estimator(solver=solver, n_components=3, init='random',
                  regularization=regularization, random_state=0, tol=1e-6)
    ft = m.fit_transform(A)
    t = m.transform(A)
    assert_array_almost_equal(ft, t, decimal=2)


@pytest.mark.parametrize('Estimator', [NMF, MiniBatchNMF])
def test_nmf_transform_custom_init(Estimator):
    # Smoke test that checks if NMF.transform works with custom initialization
    random_state = np.random.RandomState(0)
    A = np.abs(random_state.randn(6, 5))
    n_components = 4
    avg = np.sqrt(A.mean() / n_components)
    H_init = np.abs(avg * random_state.randn(n_components, 5))
    W_init = np.abs(avg * random_state.randn(6, n_components))

    m = Estimator(solver='mu', n_components=n_components, init='custom',
                  random_state=0)
    m.fit_transform(A, W=W_init, H=H_init)
    m.transform(A)


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize('regularization',
                         (None, 'both', 'components', 'transformation'))
def test_nmf_inverse_transform(Estimator, solver, regularization):
    # Test that NMF.inverse_transform returns close values
    random_state = np.random.RandomState(0)
    A = np.abs(random_state.randn(6, 4))
    m = Estimator(solver=solver, n_components=4, init='random', random_state=0,
                  regularization=regularization, max_iter=5000, tol=1e-6)
    ft = m.fit_transform(A)
    A_new = m.inverse_transform(ft)
    assert_array_almost_equal(A, A_new, decimal=2)


@pytest.mark.parametrize('Estimator', [NMF, MiniBatchNMF])
def test_n_components_greater_n_features(Estimator):
    # Smoke test for the case of more components than features.
    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(30, 10))
    init = 'random'  # FIXME : should be removed in 1.1
    Estimator(n_components=15, random_state=0, tol=1e-2, init=init).fit(A)


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize('regularization',
                         [None, 'both', 'components', 'transformation'])
def test_nmf_sparse_input(Estimator, solver, regularization):
    # Test that sparse matrices are accepted as input
    from scipy.sparse import csc_matrix

    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0
    A_sparse = csc_matrix(A)

    est1 = Estimator(solver=solver, n_components=5, init='random',
                     regularization=regularization, random_state=0,
                     tol=1e-2)
    est2 = clone(est1)

    W1 = est1.fit_transform(A)
    W2 = est2.fit_transform(A_sparse)
    H1 = est1.components_
    H2 = est2.components_

    assert_array_almost_equal(W1, W2)
    assert_array_almost_equal(H1, H2)


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
def test_nmf_sparse_transform(Estimator, solver):
    # Test that transform works on sparse data.  Issue #2124
    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(3, 2))
    A[1, 1] = 0
    A = csc_matrix(A)

    init = 'nndsvd'  # FIXME : should be removed in 1.1

    model = Estimator(solver=solver, random_state=0, n_components=2,
                      max_iter=400, init=init)
    A_fit_tr = model.fit_transform(A)
    A_tr = model.transform(A)
    assert_array_almost_equal(A_fit_tr, A_tr, decimal=1)


@pytest.mark.parametrize('init', ['random', 'nndsvd'])
@pytest.mark.parametrize(['Estimator', 'solver', 'batch_size',
                          'forget_factor'],
                         [[NMF, 'cd', None, None],
                          [NMF, 'mu', None, None],
                          [MiniBatchNMF, 'mu', 10, 0.7]])
@pytest.mark.parametrize('regularization',
                         (None, 'both', 'components', 'transformation'))
def test_non_negative_factorization_consistency(Estimator, init,
                                                solver, regularization,
                                                batch_size, forget_factor):
    # Test that the function is called in the same way, either directly
    # or through the NMF class
    max_iter = 500
    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0

    W_nmf, H, *_ = non_negative_factorization(
        A, init=init, solver=solver, max_iter=max_iter,
        regularization=regularization, random_state=1, tol=1e-2,
        batch_size=batch_size, forget_factor=forget_factor)
    W_nmf_2, *_ = non_negative_factorization(
        A, H=H, update_H=False, init=init, solver=solver,
        max_iter=max_iter, batch_size=batch_size, forget_factor=forget_factor,
        regularization=regularization, random_state=1, tol=1e-2)

    model_class = Estimator(init=init, solver=solver,
                            regularization=regularization, max_iter=max_iter,
                            random_state=1, tol=1e-2)
    W_cls = model_class.fit_transform(A)
    W_cls_2 = model_class.transform(A)

    assert_array_almost_equal(W_nmf, W_cls, decimal=10)
    assert_array_almost_equal(W_nmf_2, W_cls_2, decimal=10)


def test_non_negative_factorization_checking():
    A = np.ones((2, 2))
    # Test parameters checking is public function
    nnmf = non_negative_factorization
    msg = re.escape(
        "Number of components must be a positive integer; "
        "got (n_components=1.5)"
    )
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, A, 1.5, init='random')
    msg = re.escape(
        "Number of components must be a positive integer; "
        "got (n_components='2')"
    )
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, A, '2', init='random')
    msg = re.escape("Negative values in data passed to NMF (input H)")
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, -A, 2, init='custom')
    msg = re.escape("Negative values in data passed to NMF (input W)")
    with pytest.raises(ValueError, match=msg):
        nnmf(A, -A, A, 2, init='custom')
    msg = re.escape("Array passed to NMF (input H) is full of zeros")
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, 0 * A, 2, init='custom')
    msg = "Invalid regularization parameter: got 'spam' instead of one of"
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, 0 * A, 2, init='custom', regularization='spam')
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    msg = ("Number of samples per batch must be a positive integer; "
           "got (batch_size=0.5)")
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, A, 2, batch_size=0.5, init=init, solver='mu', beta_loss=1)
    msg = ("Number of samples per batch must be a positive integer; "
           "got (batch_size='3')")
    with pytest.raises(ValueError, match=msg):
        nnmf(A, A, A, 2, batch_size='3', init=init, solver='mu', beta_loss=1)


def _beta_divergence_dense(X, W, H, beta):
    """Compute the beta-divergence of X and W.H for dense array only.

    Used as a reference for testing nmf._beta_divergence.
    """
    WH = np.dot(W, H)

    if beta == 2:
        return squared_norm(X - WH) / 2

    WH_Xnonzero = WH[X != 0]
    X_nonzero = X[X != 0]
    np.maximum(WH_Xnonzero, 1e-9, out=WH_Xnonzero)

    if beta == 1:
        res = np.sum(X_nonzero * np.log(X_nonzero / WH_Xnonzero))
        res += WH.sum() - X.sum()

    elif beta == 0:
        div = X_nonzero / WH_Xnonzero
        res = np.sum(div) - X.size - np.sum(np.log(div))
    else:
        res = (X_nonzero ** beta).sum()
        res += (beta - 1) * (WH ** beta).sum()
        res -= beta * (X_nonzero * (WH_Xnonzero ** (beta - 1))).sum()
        res /= beta * (beta - 1)

    return res


def test_beta_divergence():
    # Compare _beta_divergence with the reference _beta_divergence_dense
    n_samples = 20
    n_features = 10
    n_components = 5
    beta_losses = [0., 0.5, 1., 1.5, 2.]

    # initialization
    rng = np.random.mtrand.RandomState(42)
    X = rng.randn(n_samples, n_features)
    np.clip(X, 0, None, out=X)
    X_csr = sp.csr_matrix(X)
    W, H = nmf._initialize_nmf(X, n_components, init='random',
                               random_state=42)

    for beta in beta_losses:
        ref = _beta_divergence_dense(X, W, H, beta)
        loss = nmf._beta_divergence(X, W, H, beta)
        loss_csr = nmf._beta_divergence(X_csr, W, H, beta)

        assert_almost_equal(ref, loss, decimal=7)
        assert_almost_equal(ref, loss_csr, decimal=7)


def test_special_sparse_dot():
    # Test the function that computes np.dot(W, H), only where X is non zero.
    n_samples = 10
    n_features = 5
    n_components = 3
    rng = np.random.mtrand.RandomState(42)
    X = rng.randn(n_samples, n_features)
    np.clip(X, 0, None, out=X)
    X_csr = sp.csr_matrix(X)

    W = np.abs(rng.randn(n_samples, n_components))
    H = np.abs(rng.randn(n_components, n_features))

    WH_safe = nmf._special_sparse_dot(W, H, X_csr)
    WH = nmf._special_sparse_dot(W, H, X)

    # test that both results have same values, in X_csr nonzero elements
    ii, jj = X_csr.nonzero()
    WH_safe_data = np.asarray(WH_safe[ii, jj]).ravel()
    assert_array_almost_equal(WH_safe_data, WH[ii, jj], decimal=10)

    # test that WH_safe and X_csr have the same sparse structure
    assert_array_equal(WH_safe.indices, X_csr.indices)
    assert_array_equal(WH_safe.indptr, X_csr.indptr)
    assert_array_equal(WH_safe.shape, X_csr.shape)


@ignore_warnings(category=ConvergenceWarning)
@pytest.mark.parametrize('forget_factor', [None, 0.7])
def test_nmf_multiplicative_update_sparse(forget_factor):
    # Compare sparse and dense input in multiplicative update NMF
    # Also test continuity of the results with respect to beta_loss parameter
    n_samples = 20
    n_features = 10
    n_components = 5
    alpha = 0.1
    l1_ratio = 0.5
    n_iter = 20

    # initialization
    rng = np.random.mtrand.RandomState(1337)
    X = rng.randn(n_samples, n_features)
    X = np.abs(X)
    X_csr = sp.csr_matrix(X)
    W0, H0 = nmf._initialize_nmf(X, n_components, init='random',
                                 random_state=42)

    for beta_loss in (-1.2, 0, 0.2, 1., 2., 2.5):
        # Reference with dense array X
        W, H = W0.copy(), H0.copy()
        W1, H1, *_ = non_negative_factorization(
            X, W, H, n_components, init='custom', update_H=True,
            solver='mu', beta_loss=beta_loss, max_iter=n_iter, alpha=alpha,
            l1_ratio=l1_ratio, regularization='both', random_state=42,
            forget_factor=forget_factor)

        # Compare with sparse X
        W, H = W0.copy(), H0.copy()
        W2, H2, *_ = non_negative_factorization(
            X_csr, W, H, n_components, init='custom', update_H=True,
            solver='mu', beta_loss=beta_loss, max_iter=n_iter, alpha=alpha,
            l1_ratio=l1_ratio, regularization='both', random_state=42,
            forget_factor=forget_factor)

        assert_array_almost_equal(W1, W2, decimal=7)
        assert_array_almost_equal(H1, H2, decimal=7)

        # Compare with almost same beta_loss, since some values have a specific
        # behavior, but the results should be continuous w.r.t beta_loss
        beta_loss -= 1.e-5
        W, H = W0.copy(), H0.copy()
        W3, H3, *_ = non_negative_factorization(
            X_csr, W, H, n_components, init='custom', update_H=True,
            solver='mu', beta_loss=beta_loss, max_iter=n_iter, alpha=alpha,
            l1_ratio=l1_ratio, regularization='both', random_state=42,
            forget_factor=forget_factor)

        assert_array_almost_equal(W1, W3, decimal=4)
        assert_array_almost_equal(H1, H3, decimal=4)


@pytest.mark.parametrize('forget_factor', [None, 0.7])
def test_nmf_negative_beta_loss(forget_factor):
    # Test that an error is raised if beta_loss < 0 and X contains zeros.
    # Test that the output has not NaN values when the input contains zeros.
    n_samples = 6
    n_features = 5
    n_components = 3

    rng = np.random.mtrand.RandomState(42)
    X = rng.randn(n_samples, n_features)
    np.clip(X, 0, None, out=X)
    X_csr = sp.csr_matrix(X)

    def _assert_nmf_no_nan(X, beta_loss):
        W, H, *_ = non_negative_factorization(
            X, init='random', n_components=n_components, solver='mu',
            beta_loss=beta_loss, random_state=0, max_iter=1000,
            forget_factor=forget_factor)
        assert not np.any(np.isnan(W))
        assert not np.any(np.isnan(H))

    msg = "When beta_loss <= 0 and X contains zeros, the solver may diverge."
    for beta_loss in (-0.6, 0.):
        with pytest.raises(ValueError, match=msg):
            _assert_nmf_no_nan(X, beta_loss)
        _assert_nmf_no_nan(X + 1e-9, beta_loss)

    for beta_loss in (0.2, 1., 1.2, 2., 2.5):
        _assert_nmf_no_nan(X, beta_loss)
        _assert_nmf_no_nan(X_csr, beta_loss)


@pytest.mark.parametrize(['Estimator', 'solver', 'beta_loss'],
                         [[NMF, 'cd', 2], [NMF, 'mu', 2],
                          [MiniBatchNMF, 'mu', 1]])
def test_nmf_regularization(Estimator, solver, beta_loss):
    # Test the effect of L1 and L2 regularizations
    n_samples = 6
    n_features = 5
    n_components = 3
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(n_samples, n_features))

    init = 'nndsvdar'
    # L1 regularization should increase the number of zeros
    l1_ratio = 1.
    max_iter = 500
    regul = Estimator(n_components=n_components, solver=solver,
                      alpha=0.5, l1_ratio=l1_ratio, random_state=42,
                      init=init, max_iter=max_iter, beta_loss=beta_loss)
    model = Estimator(n_components=n_components, solver=solver,
                      alpha=0., l1_ratio=l1_ratio, random_state=42,
                      init=init, max_iter=max_iter, beta_loss=beta_loss)

    W_regul = regul.fit_transform(X)
    W_model = model.fit_transform(X)

    H_regul = regul.components_
    H_model = model.components_

    W_regul_n_zeros = W_regul[W_regul == 0].size
    W_model_n_zeros = W_model[W_model == 0].size
    H_regul_n_zeros = H_regul[H_regul == 0].size
    H_model_n_zeros = H_model[H_model == 0].size

    assert W_regul_n_zeros > W_model_n_zeros
    assert H_regul_n_zeros > H_model_n_zeros

    # L2 regularization should decrease the sum of the squared norm
    # of the matrices
    l1_ratio = 0.
    regul = Estimator(n_components=n_components, solver=solver,
                      alpha=0.5, l1_ratio=l1_ratio, random_state=42,
                      init=init, max_iter=max_iter)
    model = Estimator(n_components=n_components, solver=solver,
                      alpha=0., l1_ratio=l1_ratio, random_state=42,
                      init=init, max_iter=max_iter)

    W_regul = regul.fit_transform(X)
    W_model = model.fit_transform(X)

    H_regul = regul.components_
    H_model = model.components_

    assert (linalg.norm(W_model))**2. + (linalg.norm(H_model))**2. > \
           (linalg.norm(W_regul))**2. + (linalg.norm(H_regul))**2.


@ignore_warnings(category=ConvergenceWarning)
@pytest.mark.parametrize('forget_factor', [None, 0.7])
def test_nmf_decreasing(forget_factor):
    # test that the objective function is decreasing at each iteration
    n_samples = 20
    n_features = 15
    n_components = 10
    alpha = 0.1
    l1_ratio = 0.5
    tol = 0.

    # initialization
    rng = np.random.mtrand.RandomState(42)
    X = rng.randn(n_samples, n_features)
    np.abs(X, X)
    W0, H0 = nmf._initialize_nmf(X, n_components, init='random',
                                 random_state=42)

    for beta_loss in (-1.2, 0, 0.2, 1., 2., 2.5):
        for solver in ('cd', 'mu'):
            if solver != 'mu' and beta_loss != 2:
                # not implemented
                continue
            if solver == 'cd' and forget_factor is not None:
                # not allowed
                continue
            W, H = W0.copy(), H0.copy()
            previous_loss = None
            for _ in range(30):
                # one more iteration starting from the previous results
                W, H, *_ = non_negative_factorization(
                    X, W, H, beta_loss=beta_loss, init='custom',
                    forget_factor=forget_factor,
                    n_components=n_components, max_iter=1, alpha=alpha,
                    solver=solver, tol=tol, l1_ratio=l1_ratio, verbose=0,
                    regularization='both', random_state=0, update_H=True)

                loss = nmf._beta_divergence(X, W, H, beta_loss)
                if previous_loss is not None:
                    assert previous_loss > loss
                previous_loss = loss


def test_nmf_underflow():
    # Regression test for an underflow issue in _beta_divergence
    rng = np.random.RandomState(0)
    n_samples, n_features, n_components = 10, 2, 2
    X = np.abs(rng.randn(n_samples, n_features)) * 10
    W = np.abs(rng.randn(n_samples, n_components)) * 10
    H = np.abs(rng.randn(n_components, n_features))

    X[0, 0] = 0
    ref = nmf._beta_divergence(X, W, H, beta=1.0)
    X[0, 0] = 1e-323
    res = nmf._beta_divergence(X, W, H, beta=1.0)
    assert_almost_equal(res, ref)


@pytest.mark.parametrize("dtype_in, dtype_out", [
    (np.float32, np.float32),
    (np.float64, np.float64),
    (np.int32, np.float64),
    (np.int64, np.float64)])
@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize("regularization",
                         (None, "both", "components", "transformation"))
def test_nmf_dtype_match(Estimator, dtype_in, dtype_out,
                         solver, regularization):
    # Check that NMF preserves dtype (float32 and float64)
    X = np.random.RandomState(0).randn(20, 15).astype(dtype_in, copy=False)
    np.abs(X, out=X)
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    nmf = Estimator(solver=solver, regularization=regularization, init=init)

    assert nmf.fit(X).transform(X).dtype == dtype_out
    assert nmf.fit_transform(X).dtype == dtype_out
    assert nmf.components_.dtype == dtype_out


@pytest.mark.parametrize(['Estimator', 'solver'],
                         [[NMF, 'cd'], [NMF, 'mu'],
                          [MiniBatchNMF, 'mu']])
@pytest.mark.parametrize("regularization",
                         (None, "both", "components", "transformation"))
def test_nmf_float32_float64_consistency(Estimator, solver, regularization):
    # Check that the result of NMF is the same between float32 and float64
    X = np.random.RandomState(0).randn(50, 7)
    np.abs(X, out=X)
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    tol = 1e-6
    nmf32 = Estimator(solver=solver, regularization=regularization,
                      random_state=0, init=init, tol=tol)
    W32 = nmf32.fit_transform(X.astype(np.float32))
    nmf64 = Estimator(solver=solver, regularization=regularization,
                      random_state=0, init=init, tol=tol)
    W64 = nmf64.fit_transform(X)

    assert_allclose(W32, W64, rtol=1e-6, atol=1e-5)


@pytest.mark.parametrize('Estimator', [NMF, MiniBatchNMF])
def test_nmf_custom_init_dtype_error(Estimator):
    # Check that an error is raise if custom H and/or W don't have the same
    # dtype as X.
    rng = np.random.RandomState(0)
    X = rng.random_sample((20, 15))
    H = rng.random_sample((15, 15)).astype(np.float32)
    W = rng.random_sample((20, 15))

    with pytest.raises(TypeError, match="should have the same dtype as X"):
        Estimator(init='custom').fit(X, H=H, W=W)

    with pytest.raises(TypeError, match="should have the same dtype as X"):
        non_negative_factorization(X, H=H, update_H=False)


def test_nmf_is_minibatch_nmf():
    # Test that the standard nmf is the minibatch nmf after 1 iteration
    # with batch_size = n_samples and forget_factor 0.0
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(48, 5))
    max_iter = 1
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    nmf = NMF(5, solver='mu', init=init, random_state=0,
              max_iter=max_iter,)
    mbnmf = MiniBatchNMF(5, solver='mu', init=init, random_state=0,
                         max_iter=max_iter,
                         batch_size=X.shape[0], forget_factor=0.0)
    W = nmf.fit_transform(X)
    mbW = mbnmf.fit_transform(X)
    assert_array_almost_equal(W, mbW)


@pytest.mark.parametrize('batch_size', [24, 32, 48])
def test_nmf_close_minibatch_nmf(batch_size):
    # Test that the decomposition with standard and minibatch nmf
    # gives close results
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(48, 5))
    max_iter = 5000
    solver = 'mu'
    beta_loss = 'kullback-leibler'
    init = 'nndsvda'  # FIXME : should be removed in 1.1
    nmf = NMF(5, solver=solver, init=init, random_state=0,
              max_iter=max_iter, beta_loss=beta_loss)
    mbnmf = MiniBatchNMF(5, solver=solver, init=init, random_state=0,
                         max_iter=max_iter, batch_size=batch_size,
                         beta_loss=beta_loss)
    W = nmf.fit_transform(X)
    mbW = mbnmf.fit_transform(X)
    assert_array_almost_equal(W, mbW, decimal=1)


def test_minibatch_nmf_partial_fit():
    rng = np.random.mtrand.RandomState(42)
    X = np.abs(rng.randn(48, 5))
    mbnmf1 = MiniBatchNMF(5, solver='mu', init='nndsvdar', random_state=0,
                          max_iter=200, batch_size=24)
    mbnmf2 = MiniBatchNMF(5, solver='mu', init='nndsvdar', random_state=0,
                          max_iter=1, batch_size=24)

    mbnmf1.fit(X)
    for i in range(mbnmf1.n_iter_):
        mbnmf2.partial_fit(X)

    assert mbnmf1.n_iter_ == mbnmf2.n_iter_
    assert_array_almost_equal(mbnmf1.components_, mbnmf2.components_,
                              decimal=0)


# FIXME : should be removed in 1.1
def test_init_default_deprecation():
    # Test FutureWarning on init default
    msg = (r"The 'init' value, when 'init=None' and "
           r"n_components is less than n_samples and "
           r"n_features, will be changed from 'nndsvd' to "
           r"'nndsvda' in 1.1 \(renaming of 0.26\).")
    rng = np.random.mtrand.RandomState(42)
    A = np.abs(rng.randn(6, 5))
    with pytest.warns(FutureWarning, match=msg):
        nmf._initialize_nmf(A, 3)
    with pytest.warns(FutureWarning, match=msg):
        NMF().fit(A)
    with pytest.warns(FutureWarning, match=msg):
        non_negative_factorization(A)
