import numpy as np
from scipy.linalg import block_diag
from scipy.sparse import csr_matrix
from scipy.special import psi

from sklearn.decomposition import LatentDirichletAllocation
from sklearn.decomposition._online_lda import (_dirichlet_expectation_1d,
                                               _dirichlet_expectation_2d)

from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import if_safe_multiprocessing_with_blas

from sklearn.utils.validation import NotFittedError
from sklearn.externals.six.moves import xrange


def _build_sparse_mtx():
    # Create 3 topics and each topic has 3 disticnt words.
    # (Each word only belongs to a single topic.)
    n_topics = 3
    block = n_topics * np.ones((3, 3))
    blocks = [block] * n_topics
    X = block_diag(*blocks)
    X = csr_matrix(X)
    return (n_topics, X)


def test_lda_default_prior_params():
    # default prior parameter should be `1 / topics`
    # and verbose params should not affect result
    n_topics, X = _build_sparse_mtx()
    prior = 1. / n_topics
    lda_1 = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=prior,
                                      topic_word_prior=prior, random_state=0)
    lda_2 = LatentDirichletAllocation(n_topics=n_topics, random_state=0)

    topic_distr_1 = lda_1.fit_transform(X)
    topic_distr_2 = lda_2.fit_transform(X)
    assert_almost_equal(topic_distr_1, topic_distr_2)


def test_lda_fit_batch():
    # Test LDA batch learning_offset (`fit` method with 'batch' learning)
    rng = np.random.RandomState(0)
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, evaluate_every=1,
                                    learning_method='batch', random_state=rng)
    lda.fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for component in lda.components_:
        # Find top 3 words in each LDA component
        top_idx = set(component.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_fit_online():
    # Test LDA online learning (`fit` method with 'online' learning)
    rng = np.random.RandomState(0)
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, learning_offset=10.,
                                    evaluate_every=1, learning_method='online',
                                    random_state=rng)
    lda.fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for component in lda.components_:
        # Find top 3 words in each LDA component
        top_idx = set(component.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_partial_fit():
    # Test LDA online learning (`partial_fit` method)
    # (same as test_lda_batch)
    rng = np.random.RandomState(0)
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, learning_offset=10.,
                                    total_samples=100, random_state=rng)
    for i in xrange(3):
        lda.partial_fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for c in lda.components_:
        top_idx = set(c.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_dense_input():
    # Test LDA with dense input.
    rng = np.random.RandomState(0)
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, learning_method='batch',
                                    random_state=rng)
    lda.fit(X.toarray())

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for component in lda.components_:
        # Find top 3 words in each LDA component
        top_idx = set(component.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_transform():
    # Test LDA transform.
    # Transform result cannot be negative
    rng = np.random.RandomState(0)
    X = rng.randint(5, size=(20, 10))
    n_topics = 3
    lda = LatentDirichletAllocation(n_topics=n_topics, random_state=rng)
    X_trans = lda.fit_transform(X)
    assert_true((X_trans > 0.0).any())


def test_lda_fit_transform():
    # Test LDA fit_transform & transform
    # fit_transform and transform result should be the same
    for method in ('online', 'batch'):
        rng = np.random.RandomState(0)
        X = rng.randint(10, size=(50, 20))
        lda = LatentDirichletAllocation(n_topics=5, learning_method=method,
                                        random_state=rng)
        X_fit = lda.fit_transform(X)
        X_trans = lda.transform(X)
        assert_array_almost_equal(X_fit, X_trans, 4)


def test_lda_partial_fit_dim_mismatch():
    # test `n_features` mismatch in `partial_fit`
    rng = np.random.RandomState(0)
    n_topics = rng.randint(3, 6)
    n_col = rng.randint(6, 10)
    X_1 = np.random.randint(4, size=(10, n_col))
    X_2 = np.random.randint(4, size=(10, n_col + 1))
    lda = LatentDirichletAllocation(n_topics=n_topics, learning_offset=5.,
                                    total_samples=20, random_state=rng)
    lda.partial_fit(X_1)
    assert_raises_regexp(ValueError, r"^The provided data has",
                         lda.partial_fit, X_2)


def test_invalid_params():
    # test `_check_params` method
    X = np.ones((5, 10))

    invalid_models = (
        ('n_topics', LatentDirichletAllocation(n_topics=0)),
        ('learning_method',
         LatentDirichletAllocation(learning_method='unknown')),
        ('total_samples', LatentDirichletAllocation(total_samples=0)),
        ('learning_offset', LatentDirichletAllocation(learning_offset=-1)),
    )
    for param, model in invalid_models:
        regex = r"^Invalid %r parameter" % param
        assert_raises_regexp(ValueError, regex, model.fit, X)


def test_lda_negative_input():
    # test pass dense matrix with sparse negative input.
    X = -np.ones((5, 10))
    lda = LatentDirichletAllocation()
    regex = r"^Negative values in data passed"
    assert_raises_regexp(ValueError, regex, lda.fit, X)


def test_lda_no_component_error():
    # test `transform` and `perplexity` before `fit`
    rng = np.random.RandomState(0)
    X = rng.randint(4, size=(20, 10))
    lda = LatentDirichletAllocation()
    regex = r"^no 'components_' attribute"
    assert_raises_regexp(NotFittedError, regex, lda.transform, X)
    assert_raises_regexp(NotFittedError, regex, lda.perplexity, X)


def test_lda_transform_mismatch():
    # test `n_features` mismatch in partial_fit and transform
    rng = np.random.RandomState(0)
    X = rng.randint(4, size=(20, 10))
    X_2 = rng.randint(4, size=(10, 8))

    n_topics = rng.randint(3, 6)
    lda = LatentDirichletAllocation(n_topics=n_topics, random_state=rng)
    lda.partial_fit(X)
    assert_raises_regexp(ValueError, r"^The provided data has",
                         lda.partial_fit, X_2)


@if_safe_multiprocessing_with_blas
def test_lda_multi_jobs():
    n_topics, X = _build_sparse_mtx()
    # Test LDA batch training with multi CPU
    for method in ('online', 'batch'):
        rng = np.random.RandomState(0)
        lda = LatentDirichletAllocation(n_topics=n_topics, n_jobs=2,
                                        learning_method=method,
                                        evaluate_every=1,
                                        random_state=rng)
        lda.fit(X)

        correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
        for c in lda.components_:
            top_idx = set(c.argsort()[-3:][::-1])
            assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


@if_safe_multiprocessing_with_blas
def test_lda_partial_fit_multi_jobs():
    # Test LDA online training with multi CPU
    rng = np.random.RandomState(0)
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, n_jobs=2,
                                    learning_offset=5., total_samples=30,
                                    random_state=rng)
    for i in range(2):
        lda.partial_fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for c in lda.components_:
        top_idx = set(c.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_preplexity_mismatch():
    # test dimension mismatch in `perplexity` method
    rng = np.random.RandomState(0)
    n_topics = rng.randint(3, 6)
    n_samples = rng.randint(6, 10)
    X = np.random.randint(4, size=(n_samples, 10))
    lda = LatentDirichletAllocation(n_topics=n_topics, learning_offset=5.,
                                    total_samples=20, random_state=rng)
    lda.fit(X)
    # invalid samples
    invalid_n_samples = rng.randint(4, size=(n_samples + 1, n_topics))
    assert_raises_regexp(ValueError, r'Number of samples', lda.perplexity, X,
                         invalid_n_samples)
    # invalid topic number
    invalid_n_topics = rng.randint(4, size=(n_samples, n_topics + 1))
    assert_raises_regexp(ValueError, r'Number of topics', lda.perplexity, X,
                         invalid_n_topics)


def test_lda_perplexity():
    # Test LDA perplexity for batch training
    # perplexity should be lower after each iteration
    n_topics, X = _build_sparse_mtx()
    for method in ('online', 'batch'):
        lda_1 = LatentDirichletAllocation(n_topics=n_topics, max_iter=1,
                                          learning_method=method,
                                          total_samples=100, random_state=0)
        lda_2 = LatentDirichletAllocation(n_topics=n_topics, max_iter=10,
                                          learning_method=method,
                                          total_samples=100, random_state=0)
        distr_1 = lda_1.fit_transform(X)
        perp_1 = lda_1.perplexity(X, distr_1, sub_sampling=False)

        distr_2 = lda_2.fit_transform(X)
        perp_2 = lda_2.perplexity(X, distr_2, sub_sampling=False)
        assert_greater_equal(perp_1, perp_2)

        perp_1_subsampling = lda_1.perplexity(X, distr_1, sub_sampling=True)
        perp_2_subsampling = lda_2.perplexity(X, distr_2, sub_sampling=True)
        assert_greater_equal(perp_1_subsampling, perp_2_subsampling)


def test_lda_score():
    # Test LDA score for batch training
    # score should be higher after each iteration
    n_topics, X = _build_sparse_mtx()
    for method in ('online', 'batch'):
        lda_1 = LatentDirichletAllocation(n_topics=n_topics, max_iter=1,
                                          learning_method=method,
                                          total_samples=100, random_state=0)
        lda_2 = LatentDirichletAllocation(n_topics=n_topics, max_iter=10,
                                          learning_method=method,
                                          total_samples=100, random_state=0)
        lda_1.fit_transform(X)
        score_1 = lda_1.score(X)

        lda_2.fit_transform(X)
        score_2 = lda_2.score(X)
        assert_greater_equal(score_2, score_1)


def test_perplexity_input_format():
    # Test LDA perplexity for sparse and dense input
    # score should be the same for both dense and sparse input
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, max_iter=1,
                                    learning_method='batch',
                                    total_samples=100, random_state=0)
    distr = lda.fit_transform(X)
    perp_1 = lda.perplexity(X)
    perp_2 = lda.perplexity(X, distr)
    perp_3 = lda.perplexity(X.toarray(), distr)
    assert_almost_equal(perp_1, perp_2)
    assert_almost_equal(perp_1, perp_3)


def test_lda_score_perplexity():
    # Test the relationship between LDA score and perplexity
    n_topics, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, max_iter=10,
                                    random_state=0)
    distr = lda.fit_transform(X)
    perplexity_1 = lda.perplexity(X, distr, sub_sampling=False)

    score = lda.score(X)
    perplexity_2 = np.exp(-1. * (score / np.sum(X.data)))
    assert_almost_equal(perplexity_1, perplexity_2)


def test_lda_empty_docs():
    """Test LDA on empty document (all-zero rows)."""
    Z = np.zeros((5, 4))
    for X in [Z, csr_matrix(Z)]:
        lda = LatentDirichletAllocation(max_iter=750).fit(X)
        assert_almost_equal(lda.components_.sum(axis=0),
                            np.ones(lda.components_.shape[1]))


def test_dirichlet_expectation():
    """Test Cython version of Dirichlet expectation calculation."""
    x = np.logspace(-100, 10, 10000)
    expectation = np.empty_like(x)
    _dirichlet_expectation_1d(x, 0, expectation)
    assert_allclose(expectation, np.exp(psi(x) - psi(np.sum(x))),
                    atol=1e-19)

    x = x.reshape(100, 100)
    assert_allclose(_dirichlet_expectation_2d(x),
                    psi(x) - psi(np.sum(x, axis=1)[:, np.newaxis]),
                    rtol=1e-11, atol=3e-9)
