import numpy as np
from scipy import linalg
from sklearn.decomposition import nmf

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import raises
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less


random_state = np.random.mtrand.RandomState(0)


@raises(ValueError)
def test_initialize_nn_input():
    """Test NNDSVD behaviour on negative input"""
    nmf._initialize_nmf(-np.ones((2, 2)), 2)


def test_initialize_nn_output():
    """Test that NNDSVD does not return negative values"""
    data = np.abs(random_state.randn(10, 10))
    for var in (None, 'a', 'ar'):
        W, H = nmf._initialize_nmf(data, 10, random_state=0)
        assert_false((W < 0).any() or (H < 0).any())


def test_initialize_close():
    """Test NNDSVD error

    Test that _initialize_nmf error is less than the standard deviation of the
    entries in the matrix.
    """
    A = np.abs(random_state.randn(10, 10))
    W, H = nmf._initialize_nmf(A, 10)
    error = linalg.norm(np.dot(W, H) - A)
    sdev = linalg.norm(A - A.mean())
    assert_true(error <= sdev)


def test_initialize_variants():
    """Test NNDSVD variants correctness

    Test that the variants 'a' and 'ar' differ from basic NNDSVD only where
    the basic version has zeros.
    """
    data = np.abs(random_state.randn(10, 10))
    W0, H0 = nmf._initialize_nmf(data, 10, variant=None)
    Wa, Ha = nmf._initialize_nmf(data, 10, variant='a')
    War, Har = nmf._initialize_nmf(data, 10, variant='ar', random_state=0)

    for ref, evl in ((W0, Wa), (W0, War), (H0, Ha), (H0, Har)):
        assert_true(np.allclose(evl[ref != 0], ref[ref != 0]))


@raises(ValueError)
def test_projgrad_nmf_fit_nn_input():
    """Test model fit behaviour on negative input"""
    A = -np.ones((2, 2))
    m = nmf.ProjectedGradientNMF(n_components=2, init=None, random_state=0)
    m.fit(A)


def test_projgrad_nmf_fit_nn_output():
    """Test that the decomposition does not contain negative values"""
    A = np.c_[5 * np.ones(5) - np.arange(1, 6),
              5 * np.ones(5) + np.arange(1, 6)]
    for init in (None, 'nndsvd', 'nndsvda', 'nndsvdar'):
        model = nmf.ProjectedGradientNMF(n_components=2, init=init,
                                         random_state=0)
        transf = model.fit_transform(A)
        assert_false((model.components_ < 0).any() or
                     (transf < 0).any())


def test_projgrad_nmf_fit_close():
    """Test that the fit is not too far away"""
    X = np.abs(random_state.randn(6, 5))
    for solver in ('pg', 'lbfgs'):
        pnmf = nmf.NMF(n_components=5, init='nndsvda', nnls_solver=solver,
                       random_state=0)
        print pnmf
        assert_less(pnmf.fit(X).reconstruction_err_, 0.05)


@raises(ValueError)
def test_nls_nn_input():
    """Test NLS solver's behaviour on negative input"""
    A = np.ones((2, 2))
    for solver in ('pg', 'lbfgs'):
        model = nmf.NMF(nnls_solver=solver)
        model._nnls(A, A, -A, 0.001, 20)


def test_nls_nn_output():
    """Test that NLS solver doesn't return negative values"""
    A = np.arange(1, 5).reshape(1, -1)
    for solver in ('pg', 'lbfgs'):
        model = nmf.NMF(nnls_solver=solver)
        Ap, _, _ = model._nnls(A.T, np.dot(A.T, -A),  A, 0.001, 100)
    assert_false((Ap < 0).any())


def test_nls_close():
    """Test that the NLS results should be close"""
    A = np.arange(1, 5).reshape(1, -1)
    for solver in ('pg', 'lbfgs'):
        model = nmf.NMF(nnls_solver=solver)
        Ap, _, _ = model._nnls(A.T, np.dot(A.T, A), np.zeros_like(A),
                               0.001, 100)
        print model
        print Ap
        assert_array_almost_equal(Ap, A, decimal=3)


def test_projgrad_nmf_transform():
    """Test that NMF.transform returns close values

    (transform uses scipy.optimize.nnls for now)
    """
    A = np.abs(random_state.randn(6, 5))
    for solver in ('pg', 'lbfgs'):
        m = nmf.NMF(n_components=5, init='nndsvd', random_state=0,
                    nnls_solver=solver)
        transf = m.fit_transform(A)
        assert_array_almost_equal(transf, m.transform(A), decimal=2)


def test_n_components_greater_n_features():
    """Smoke test for the case of more components than features."""
    A = np.abs(random_state.randn(30, 10))
    for solver in ('pg', 'lbfgs'):
        nmf.NMF(n_components=15, sparseness='data', random_state=0,
                nnls_solver=solver).fit(A)


def test_projgrad_nmf_sparseness():
    """Test sparseness

    Test that sparsity constraints actually increase sparseness in the
    part where they are applied.
    """
    A = np.abs(random_state.randn(10, 10))
    m = nmf.ProjectedGradientNMF(n_components=5, random_state=0).fit(A)
    data_sp = nmf.ProjectedGradientNMF(n_components=5, sparseness='data',
                                       random_state=0).fit(A).data_sparseness_
    comp_sp = nmf.ProjectedGradientNMF(n_components=5, sparseness='components',
                                       random_state=0).fit(A).comp_sparseness_
    assert_greater(data_sp, m.data_sparseness_)
    assert_greater(comp_sp, m.comp_sparseness_)


def test_sparse_input():
    """Test that sparse matrices are accepted as input"""
    from scipy.sparse import csr_matrix

    A = np.abs(random_state.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0
    T1 = nmf.ProjectedGradientNMF(n_components=5, init='random',
                                  random_state=999).fit_transform(A)

    A_sparse = csr_matrix(A)
    pg_nmf = nmf.ProjectedGradientNMF(n_components=5, init='random',
                                      random_state=999)
    T2 = pg_nmf.fit_transform(A_sparse)
    assert_array_almost_equal(pg_nmf.reconstruction_err_,
                              linalg.norm(A - np.dot(T2, pg_nmf.components_),
                                          'fro'))
    assert_array_almost_equal(T1, T2)

    # same with sparseness

    T2 = nmf.ProjectedGradientNMF(
        n_components=5, init='random', sparseness='data',
        random_state=999).fit_transform(A_sparse)
    T1 = nmf.ProjectedGradientNMF(
        n_components=5, init='random', sparseness='data',
        random_state=999).fit_transform(A)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
