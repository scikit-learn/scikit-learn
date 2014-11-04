import numpy as np

from sklearn.utils.testing import (assert_array_almost_equal,
                                   assert_equal,
                                   assert_less_equal,
                                   assert_true)

from sklearn.decomposition import DictionaryLearningKSVD, learn_dictionary_ksvd
from ...utils.extmath import norm


random_state = np.random.RandomState(0)
data5x3 = random_state.randn(5, 3)
data10x5 = random_state.randn(10, 5)


def check_learned_shape(data, n_features, n_components):
    coder = DictionaryLearningKSVD(n_components=n_components,
                                   random_state=random_state).fit(data)
    assert_true(coder.components_.shape == (n_components, n_features))


def test_learned_dictionary_shape():
    check_learned_shape(data5x3, 3, 4)
    check_learned_shape(data10x5, 5, 4)
    check_learned_shape(data10x5, 5, 7)


def check_reconstruction(data, n_components, decimal,
                         **ksvd_parameters):
    """Check that learned dictionary gives good reconstruction."""

    coder = DictionaryLearningKSVD(n_components=n_components,
                                   random_state=random_state,
                                   **ksvd_parameters)
    code = coder.fit(data).transform(data)
    reconstructed = np.dot(code, coder.components_)
    assert_array_almost_equal(data, reconstructed, decimal)


def test_reconstruction():
    check_reconstruction(data10x5, 7, 0)
    check_reconstruction(data10x5, 7, 0, approximate_svd=False)
    check_reconstruction(data10x5, 10, 6)
    check_reconstruction(data10x5, 10, 6, approximate_svd=False)
    check_reconstruction(data10x5, 12, 6)
    check_reconstruction(data10x5, 12, 6, approximate_svd=False)


def test_n_nonzero_coefs():
    data, n_samples, n_features, n_components = data10x5, 10, 5, 7
    n_nonzero_coefs = 2
    sparse_codes, _, _ = learn_dictionary_ksvd(data,
                                               n_nonzero_coefs=n_nonzero_coefs,
                                               iteration_count=10,
                                               n_components=n_components,
                                               random_state=random_state)
    for sample_index in range(n_samples):
        assert_less_equal(len(np.flatnonzero(sparse_codes[sample_index, :])),
                          n_nonzero_coefs)


def test_error_path():
    """Check that errors decrease and the last error is correct."""

    data, n_samples, n_features, n_components = data10x5, 10, 5, 7
    n_nonzero_coefs = 2
    codes, dictionary, errors = learn_dictionary_ksvd(
        data, n_nonzero_coefs=n_nonzero_coefs, iteration_count=10,
        n_components=n_components, random_state=random_state)

    for prev_error, error in zip(errors[:-1], errors[1:]):
        assert_less_equal(error, prev_error)
    last_residue = data - np.dot(codes, dictionary)
    assert_equal(errors[-1], norm(last_residue))


def test_params():
    """Test passing different params."""

    init_dict = random_state.randn(4, 3)
    _, _, _ = learn_dictionary_ksvd(data5x3, n_nonzero_coefs=2,
                                    iteration_count=10,
                                    init_dictionary=init_dict,
                                    n_components=4,
                                    random_state=random_state)

    _, _, _ = learn_dictionary_ksvd(data5x3, n_nonzero_coefs=2,
                                    iteration_count=10,
                                    n_components=4,
                                    random_state=random_state,
                                    callback=lambda l: None)

    coder = DictionaryLearningKSVD(n_components=4,
                                   fit_n_nonzero_coefs=2,
                                   random_state=random_state)
    coder.fit(data5x3)
