import numpy as np
from scipy.sparse import csr_matrix
from scipy.linalg import block_diag
from sklearn.decomposition import LatentDirichletAllocation

from sklearn.utils.testing import assert_true
from sklearn.utils.testing import raises
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import if_not_mac_os

from sklearn.utils.validation import NotFittedError
from sklearn.externals.six.moves import xrange


def _build_sparse_mtx():
    # Create 3 topics and each have 3 words
    n_topics = 3
    doc_topic_prior = 1. / n_topics
    topic_word_prior = 1. / n_topics
    block = n_topics * np.ones((3, 3))
    blocks = [block] * n_topics
    X = block_diag(*blocks)
    X = csr_matrix(X)
    return (n_topics, doc_topic_prior, topic_word_prior, X)


def test_lda_fit_batch():
    """
    Test LDA batch learning_offset (`fit` method with 'batch' learning)
    """
    rng = np.random.RandomState(0)
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior,
                                    learning_method='batch', random_state=rng)
    lda.fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for c in lda.components_:
        top_idx = set(c.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_fit_online():
    """
    Test LDA online learning (`fit` method with 'online' learning)
    """
    rng = np.random.RandomState(0)
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, learning_offset=30.,
                                    learning_method='online', random_state=rng)
    lda.fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for c in lda.components_:
        top_idx = set(c.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_partial_fit():
    """
    Test LDA online learning (`partial_fit` method)
    (same as test_lda_batch)
    """
    rng = np.random.RandomState(0)
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, learning_offset=30.,
                                    random_state=rng)
    for i in xrange(3):
        lda.partial_fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for c in lda.components_:
        top_idx = set(c.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_dense_input():
    """
    Test LDA with dense input.
    """
    rng = np.random.RandomState(0)
    X = rng.randint(5, size=(20, 10))
    n_topics = 3
    doc_topic_prior = 1. / n_topics
    topic_word_prior = 1. / n_topics
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, random_state=rng)
    X_trans = lda.fit_transform(X)
    assert_true((X_trans > 0.0).any())


def test_lda_fit_transform():
    """
    Test LDA fit_transform & transform
    fit_transform and transform result should be the same
    """
    for method in ('online', 'batch'):
        rng = np.random.RandomState(0)
        n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
        lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                        topic_word_prior=topic_word_prior,
                                        learning_method=method, random_state=rng)
        X_fit = lda.fit_transform(X)
        X_trans = lda.transform(X)
        assert_array_almost_equal(X_fit, X_trans, 4)


def test_lda_normalize_docs():
    """
    test sum of topic distribution equals to 1 for each doc
    """
    rng = np.random.RandomState(0)
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, normalize_doc=True,
                                    random_state=rng)
    X_fit = lda.fit_transform(X)
    assert_array_almost_equal(X_fit.sum(axis=1), np.ones(X.shape[0]))


@raises(ValueError)
def test_lda_partial_fit_dim_mismatch():
    """
    test `n_features` mismatch in `partial_fit`
    """
    rng = np.random.RandomState(0)
    n_topics = rng.randint(3, 6)
    doc_topic_prior = 1. / n_topics
    topic_word_prior = 1. / n_topics

    n_col = rng.randint(6, 10)
    X_1 = np.random.randint(4, size=(10, n_col))
    X_2 = np.random.randint(4, size=(10, n_col + 1))
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, learning_offset=5.,
                                    total_samples=20, random_state=rng)
    for X in [X_1, X_2]:
        lda.partial_fit(X)


@raises(ValueError)
def test_invalid_learning_method():
    """
    test invalid learing method
    """
    lda = LatentDirichletAllocation(learning_method='unknown')


@raises(ValueError)
def test_lda_negative_input():
    """
    test pass dense matrix with sparse negative input.
    """
    X = -np.ones((5, 10))
    lda = LatentDirichletAllocation()
    lda.fit(X)


@raises(NotFittedError)
def test_lda_transform_before_fit():
    """
    test `transform` before `fit`
    """
    rng = np.random.RandomState(0)
    X = rng.randint(4, size=(20, 10))
    lda = LatentDirichletAllocation()
    lda.transform(X)


@raises(ValueError)
def test_lda_transform_mismatch():
    """
    test `n_features` mismatch in fit and transform
    """
    rng = np.random.RandomState(0)
    X = rng.randint(4, size=(20, 10))
    X_2 = rng.randint(4, size=(10, 8))

    n_topics = rng.randint(3, 6)
    doc_topic_prior = 1. / n_topics
    topic_word_prior = 1. / n_topics

    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, random_state=rng)
    lda.partial_fit(X)
    lda.transform(X_2)


@if_not_mac_os()
def test_lda_multi_jobs():
    """
    Test LDA batch training with multi CPU
    """
    for method in ('online', 'batch'):
        rng = np.random.RandomState(0)
        n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
        lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                        topic_word_prior=topic_word_prior, n_jobs=3,
                                        learning_method=method, random_state=rng)
        lda.fit(X)

        correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
        for c in lda.components_:
            top_idx = set(c.argsort()[-3:][::-1])
            assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


@if_not_mac_os()
def test_lda_online_multi_jobs():
    """
    Test LDA online training with multi CPU
    """
    rng = np.random.RandomState(0)
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, n_jobs=2,
                                    learning_offset=5., total_samples=30, random_state=rng)
    for i in xrange(3):
        lda.partial_fit(X)

    correct_idx_grps = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    for c in lda.components_:
        top_idx = set(c.argsort()[-3:][::-1])
        assert_true(tuple(sorted(top_idx)) in correct_idx_grps)


def test_lda_perplexity():
    """
    Test LDA perplexity for batch training
    perplexity should be lower after each iteration
    """
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    for method in ('online', 'batch'):
        lda_1 = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                          topic_word_prior=topic_word_prior, max_iter=1,
                                          learning_method=method, total_samples=100, random_state=0)
        lda_2 = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                          topic_word_prior=topic_word_prior, max_iter=10,
                                          learning_method=method, total_samples=100, random_state=0)
        distr_1 = lda_1.fit_transform(X)
        perp_1 = lda_1.perplexity(X, distr_1, sub_sampling=False)

        distr_2 = lda_2.fit_transform(X)
        perp_2 = lda_2.perplexity(X, distr_2, sub_sampling=False)
        assert_greater_equal(perp_1, perp_2)


def test_lda_score():
    """
    Test LDA score for batch training
    score should be higher after each iteration
    """
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    for method in ('online', 'batch'):
        lda_1 = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                          topic_word_prior=topic_word_prior, max_iter=1,
                                          learning_method=method, total_samples=100, random_state=0)
        lda_2 = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                          topic_word_prior=topic_word_prior, max_iter=10,
                                          learning_method=method, total_samples=100, random_state=0)
        lda_1.fit_transform(X)
        score_1 = lda_1.score(X)

        lda_2.fit_transform(X)
        score_2 = lda_2.score(X)
        assert_greater_equal(score_2, score_1)


def test_lda_score_perplexity():
    """
    Test the relationship between LDA score and perplexity
    """
    n_topics, doc_topic_prior, topic_word_prior, X = _build_sparse_mtx()
    lda = LatentDirichletAllocation(n_topics=n_topics, doc_topic_prior=doc_topic_prior,
                                    topic_word_prior=topic_word_prior, max_iter=10,
                                    random_state=0)
    distr = lda.fit_transform(X)
    perplexity_1 = lda.perplexity(X, distr, sub_sampling=False)

    score = lda.score(X)
    perplexity_2 = np.exp(-1. * (score / np.sum(X.data)))
    assert_almost_equal(perplexity_1, perplexity_2)


if __name__ == '__main__':
    import nose
    nose.run(argv=['', __file__])
