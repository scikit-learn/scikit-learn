import numpy as np
from .. import nmf
from nose.tools import assert_true, assert_false, raises
from numpy.testing import assert_array_almost_equal

from sklearn.utils.testing import assert_greater

random_state = np.random.mtrand.RandomState(0)


@raises(ValueError)
def test_initialize_nn_input():
    """Test NNDSVD behaviour on negative input"""
    nmf._initialize_nmf(-np.ones((2, 2)), 2)


def test_initialize_nn_output():
    """Test that NNDSVD does not return negative values"""
    data = np.abs(random_state.randn(10, 10))
    for var in (None, 'a', 'ar'):
        W, H = nmf._initialize_nmf(data, 10)
        assert_false((W < 0).any() or (H < 0).any())


def test_initialize_close():
    """Test NNDSVD error

    Test that _initialize_nmf error is less than the standard deviation of the
    entries in the matrix.
    """
    A = np.abs(random_state.randn(10, 10))
    W, H = nmf._initialize_nmf(A, 10)
    error = np.linalg.norm(np.dot(W, H) - A)
    sdev = np.linalg.norm(A - A.mean())
    assert_true(error <= sdev)


def test_initialize_variants():
    """Test NNDSVD variants correctness

    Test that the variants 'a' and 'ar' differ from basic NNDSVD only where
    the basic version has zeros.
    """
    data = np.abs(random_state.randn(10, 10))
    W0, H0 = nmf._initialize_nmf(data, 10, variant=None)
    Wa, Ha = nmf._initialize_nmf(data, 10, variant='a')
    War, Har = nmf._initialize_nmf(data, 10, variant='ar')

    for ref, evl in ((W0, Wa), (W0, War), (H0, Ha), (H0, Har)):
        assert_true(np.allclose(evl[ref != 0], ref[ref != 0]))


@raises(ValueError)
def test_projgrad_nmf_fit_nn_input():
    """Test model fit behaviour on negative input"""
    A = -np.ones((2, 2))
    m = nmf.ProjectedGradientNMF(n_components=2, init=None)
    m.fit(A)


def test_projgrad_nmf_fit_nn_output():
    """Test that the decomposition does not contain negative values"""
    A = np.c_[5 * np.ones(5) - xrange(1, 6),
              5 * np.ones(5) + xrange(1, 6)]
    for init in (None, 'nndsvd', 'nndsvda', 'nndsvdar'):
        model = nmf.ProjectedGradientNMF(n_components=2, init=init)
        transf = model.fit_transform(A)
        assert_false((model.components_ < 0).any() or
                     (transf < 0).any())


def test_projgrad_nmf_fit_close():
    """Test that the fit is not too far away"""
    assert_true(nmf.ProjectedGradientNMF(5, init='nndsvda').fit(np.abs(
      random_state.randn(6, 5))).reconstruction_err_ < 0.05)


@raises(ValueError)
def test_nls_nn_input():
    """Test NLS solver's behaviour on negative input"""
    A = np.ones((2, 2))
    nmf._nls_subproblem(A, A, -A, 0.001, 20)


def test_nls_nn_output():
    """Test that NLS solver doesn't return negative values"""
    A = np.atleast_2d(range(1, 5))
    Ap, _, _ = nmf._nls_subproblem(np.dot(A.T, -A), A.T, A, 0.001, 100)
    assert_false((Ap < 0).any())


def test_nls_close():
    """Test that the NLS results should be close"""
    A = np.atleast_2d(range(1, 5))
    Ap, _, _ = nmf._nls_subproblem(np.dot(A.T, A), A.T, np.zeros_like(A),
                                   0.001, 100)
    assert_true((np.abs(Ap - A) < 0.01).all())


def test_projgrad_nmf_transform():
    """Test that NMF.transform returns close values

    (transform uses scipy.optimize.nnls for now)
    """
    A = np.abs(random_state.randn(6, 5))
    m = nmf.ProjectedGradientNMF(n_components=5, init='nndsvd')
    transf = m.fit_transform(A)
    assert_true(np.allclose(transf, m.transform(A), atol=1e-2, rtol=0))


def test_projgrad_nmf_sparseness():
    """Test sparseness

    Test that sparsity contraints actually increase sparseness in the
    part where they are applied.
    """

    A = np.abs(random_state.randn(10, 10))
    m = nmf.ProjectedGradientNMF(n_components=5).fit(A)
    data_sp = nmf.ProjectedGradientNMF(n_components=5,
                  sparseness='data').fit(A).data_sparseness_
    comp_sp = nmf.ProjectedGradientNMF(n_components=5,
                  sparseness='components').fit(A).comp_sparseness_
    assert_greater(data_sp, m.data_sparseness_)
    assert_greater(comp_sp, m.comp_sparseness_)


def test_sparse_input():
    """Test that sparse matrices are accepted as input"""
    from scipy.sparse import csr_matrix

    A = np.abs(random_state.randn(10, 10))
    A[:, 2 * np.arange(5)] = 0
    T1 = nmf.ProjectedGradientNMF(n_components=5, init=999).fit_transform(A)

    A = csr_matrix(A)
    T2 = nmf.ProjectedGradientNMF(n_components=5, init=999).fit_transform(A)
    assert_array_almost_equal(T1, T2)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
