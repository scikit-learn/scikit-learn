import numpy as np
from scipy import linalg
from sklearn.decomposition import (NMF, ProjectedGradientNMF,
                                   non_negative_factorization)
from sklearn.decomposition import nmf   # For testing internals
from scipy.sparse import csc_matrix

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import ignore_warnings
from sklearn.base import clone


random_state = np.random.mtrand.RandomState(0)


def test_initialize_nn_output():
    # Test that initialization does not return negative values
    data = np.abs(random_state.randn(10, 10))
    for init in ('random', 'nndsvd', 'nndsvda', 'nndsvdar'):
        W, H = nmf._initialize_nmf(data, 10, init=init, random_state=0)
        assert_false((W < 0).any() or (H < 0).any())


@ignore_warnings
def test_parameter_checking():
    A = np.ones((2, 2))
    name = 'spam'
    msg = "Invalid solver parameter: got 'spam' instead of one of"
    assert_raise_message(ValueError, msg, NMF(solver=name).fit, A)
    msg = "Invalid init parameter: got 'spam' instead of one of"
    assert_raise_message(ValueError, msg, NMF(init=name).fit, A)
    msg = "Invalid sparseness parameter: got 'spam' instead of one of"
    assert_raise_message(ValueError, msg, NMF(sparseness=name).fit, A)

    msg = "Negative values in data passed to"
    assert_raise_message(ValueError, msg, NMF().fit, -A)
    assert_raise_message(ValueError, msg, nmf._initialize_nmf, -A,
                         2, 'nndsvd')
    clf = NMF(2, tol=0.1).fit(A)
    assert_raise_message(ValueError, msg, clf.transform, -A)


def test_initialize_close():
    # Test NNDSVD error
    # Test that _initialize_nmf error is less than the standard deviation of
    # the entries in the matrix.
    A = np.abs(random_state.randn(10, 10))
    W, H = nmf._initialize_nmf(A, 10, init='nndsvd')
    error = linalg.norm(np.dot(W, H) - A)
    sdev = linalg.norm(A - A.mean())
    assert_true(error <= sdev)


def test_initialize_variants():
    # Test NNDSVD variants correctness
    # Test that the variants 'nndsvda' and 'nndsvdar' differ from basic
    # 'nndsvd' only where the basic version has zeros.
    data = np.abs(random_state.randn(10, 10))
    W0, H0 = nmf._initialize_nmf(data, 10, init='nndsvd')
    Wa, Ha = nmf._initialize_nmf(data, 10, init='nndsvda')
    War, Har = nmf._initialize_nmf(data, 10, init='nndsvdar',
                                   random_state=0)

    for ref, evl in ((W0, Wa), (W0, War), (H0, Ha), (H0, Har)):
        assert_true(np.allclose(evl[ref != 0], ref[ref != 0]))


@ignore_warnings
def test_nmf_fit_nn_output():
    # Test that the decomposition does not contain negative values
    A = np.c_[5 * np.ones(5) - np.arange(1, 6),
              5 * np.ones(5) + np.arange(1, 6)]
    for solver in ('pg', 'cd'):
        for init in (None, 'nndsvd', 'nndsvda', 'nndsvdar'):
            model = NMF(n_components=2, solver=solver, init=init,
                        random_state=0)
            transf = model.fit_transform(A)
            assert_false((model.components_ < 0).any() or
                         (transf < 0).any())


@ignore_warnings
def test_nmf_fit_close():
    # Test that the fit is not too far away
    for solver in ('pg', 'cd'):
        pnmf = NMF(5, solver=solver, init='nndsvd', random_state=0)
        X = np.abs(random_state.randn(6, 5))
        assert_less(pnmf.fit(X).reconstruction_err_, 0.05)


def test_nls_nn_output():
    # Test that NLS solver doesn't return negative values
    A = np.arange(1, 5).reshape(1, -1)
    Ap, _, _ = nmf._nls_subproblem(np.dot(A.T, -A), A.T, A, 0.001, 100)
    assert_false((Ap < 0).any())


def test_nls_close():
    # Test that the NLS results should be close
    A = np.arange(1, 5).reshape(1, -1)
    Ap, _, _ = nmf._nls_subproblem(np.dot(A.T, A), A.T, np.zeros_like(A),
                                   0.001, 100)
    assert_true((np.abs(Ap - A) < 0.01).all())


@ignore_warnings
def test_nmf_transform():
    # Test that NMF.transform returns close values
    A = np.abs(random_state.randn(6, 5))
    for solver in ('pg', 'cd'):
        m = NMF(solver=solver, n_components=4, init='nndsvd', random_state=0)
        ft = m.fit_transform(A)
        t = m.transform(A)
        assert_array_almost_equal(ft, t, decimal=2)


def test_nmf_transform_custom_init():
    # Smoke test that checks if NMF.transform works with custom initialization
    A = np.abs(random_state.randn(6, 5))
    n_components = 4
    avg = np.sqrt(A.mean() / n_components)
    H_init = np.abs(avg * random_state.randn(n_components, 5))
    W_init = np.abs(avg * random_state.randn(6, n_components))

    m = NMF(solver='cd', n_components=n_components, init='custom', random_state=0)
    ft = m.fit_transform(A, W=W_init, H=H_init)
    t = m.transform(A)



@ignore_warnings
def test_nmf_inverse_transform():
    # Test that NMF.inverse_transform returns close values
    random_state = np.random.RandomState(0)
    A = np.abs(random_state.randn(6, 4))
    for solver in ('pg', 'cd'):
        m = NMF(solver=solver, n_components=4, init='random', random_state=0)
        ft = m.fit_transform(A)
        t = m.transform(A)
        A_new = m.inverse_transform(t)
        assert_array_almost_equal(A, A_new, decimal=2)


@ignore_warnings
def test_n_components_greater_n_features():
    # Smoke test for the case of more components than features.
    A = np.abs(random_state.randn(30, 10))
    NMF(n_components=15, random_state=0, tol=1e-2).fit(A)


@ignore_warnings
def test_projgrad_nmf_sparseness():
    # Test sparseness
    # Test that sparsity constraints actually increase sparseness in the
    # part where they are applied.
    tol = 1e-2
    A = np.abs(random_state.randn(10, 10))
    m = ProjectedGradientNMF(n_components=5, random_state=0, tol=tol).fit(A)
    data_sp = ProjectedGradientNMF(n_components=5, sparseness='data',
                                   random_state=0,
                                   tol=tol).fit(A).data_sparseness_
    comp_sp = ProjectedGradientNMF(n_components=5, sparseness='components',
                                   random_state=0,
                                   tol=tol).fit(A).comp_sparseness_
    assert_greater(data_sp, m.data_sparseness_)
    assert_greater(comp_sp, m.comp_sparseness_)


@ignore_warnings
def test_sparse_input():
    # Test that sparse matrices are accepted as input
    from scipy.sparse import csc_matrix

    A = np.abs(random_state.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0
    A_sparse = csc_matrix(A)

    for solver in ('pg', 'cd'):
        est1 = NMF(solver=solver, n_components=5, init='random',
                   random_state=0, tol=1e-2)
        est2 = clone(est1)

        W1 = est1.fit_transform(A)
        W2 = est2.fit_transform(A_sparse)
        H1 = est1.components_
        H2 = est2.components_

        assert_array_almost_equal(W1, W2)
        assert_array_almost_equal(H1, H2)


@ignore_warnings
def test_sparse_transform():
    # Test that transform works on sparse data.  Issue #2124

    A = np.abs(random_state.randn(3, 2))
    A[A > 1.0] = 0
    A = csc_matrix(A)

    for solver in ('pg', 'cd'):
        model = NMF(solver=solver, random_state=0, tol=1e-4, n_components=2)
        A_fit_tr = model.fit_transform(A)
        A_tr = model.transform(A)
        assert_array_almost_equal(A_fit_tr, A_tr, decimal=1)


@ignore_warnings
def test_non_negative_factorization_consistency():
    # Test that the function is called in the same way, either directly
    # or through the NMF class
    A = np.abs(random_state.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0

    for solver in ('pg', 'cd'):
        W_nmf, H, _ = non_negative_factorization(
            A, solver=solver, random_state=1, tol=1e-2)
        W_nmf_2, _, _ = non_negative_factorization(
            A, H=H, update_H=False, solver=solver, random_state=1, tol=1e-2)

        model_class = NMF(solver=solver, random_state=1, tol=1e-2)
        W_cls = model_class.fit_transform(A)
        W_cls_2 = model_class.transform(A)
        assert_array_almost_equal(W_nmf, W_cls, decimal=10)
        assert_array_almost_equal(W_nmf_2, W_cls_2, decimal=10)


@ignore_warnings
def test_non_negative_factorization_checking():
    A = np.ones((2, 2))
    # Test parameters checking is public function
    nnmf = non_negative_factorization
    msg = "Number of components must be positive; got (n_components='2')"
    assert_raise_message(ValueError, msg, nnmf, A, A, A, '2')
    msg = "Negative values in data passed to NMF (input H)"
    assert_raise_message(ValueError, msg, nnmf, A, A, -A, 2, 'custom')
    msg = "Negative values in data passed to NMF (input W)"
    assert_raise_message(ValueError, msg, nnmf, A, -A, A, 2, 'custom')
    msg = "Array passed to NMF (input H) is full of zeros"
    assert_raise_message(ValueError, msg, nnmf, A, A, 0 * A, 2, 'custom')


def test_safe_compute_error():
    A = np.abs(random_state.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0
    A_sparse = csc_matrix(A)

    W, H = nmf._initialize_nmf(A, 5, init='random', random_state=0)

    error = nmf._safe_compute_error(A, W, H)
    error_sparse = nmf._safe_compute_error(A_sparse, W, H)

    assert_almost_equal(error, error_sparse)
