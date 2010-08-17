from scikits.learn.features.text import strip_accents
from scikits.learn.features.text import SimpleAnalyzer
from scikits.learn.features.text import HashingVectorizer
from scikits.learn.features.text import SparseHashingVectorizer
from scikits.learn.logistic import LogisticRegression
from scikits.learn.sparse.svm import SVC
import numpy as np
from nose.tools import *
from numpy.testing import assert_array_almost_equal

JUNK_FOOD_DOCS = (
    "the pizza pizza beer",
    "the pizza pizza beer",
    "the the pizza beer beer",
    "the pizza beer beer",
    "the coke beer coke",
    "the coke pizza pizza",
)

NOTJUNK_FOOD_DOCS = (
    "the salad celeri",
    "the salad salad sparkling water",
    "the the celeri celeri",
    "the tomato tomato salad water",
    "the tomato salad water",
)

def test_strip_accents():
    # check some classical latin accentuated symbols
    a = u'\xe0\xe1\xe2\xe3\xe4\xe5\xe7\xe8\xe9\xea\xeb'
    expected = u'aaaaaaceeee'
    assert_equal(strip_accents(a), expected)

    a = u'\xec\xed\xee\xef\xf1\xf2\xf3\xf4\xf5\xf6\xf9\xfa\xfb\xfc\xfd'
    expected = u'iiiinooooouuuuy'
    assert_equal(strip_accents(a), expected)

    # check some arabic
    a = u'\u0625' # halef with a hamza below
    expected = u'\u0627' # simple halef
    assert_equal(strip_accents(a), expected)

    # mix letters accentuated and not
    a = u"this is \xe0 test"
    expected = u'this is a test'
    assert_equal(strip_accents(a), expected)


def test_simple_analyzer():
    sa = SimpleAnalyzer()

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u'ai', u'mange', u'du', u'kangourou', u'ce', u'midi',
                u'etait', u'pas', u'tres', u'bon']
    assert_equal(sa.analyze(text), expected)

    text = "This is a test, really.\n\n I met Harry yesterday."
    expected = [u'this', u'is', u'test', u'really', u'met', u'harry',
                u'yesterday']
    assert_equal(sa.analyze(text), expected)


def test_dense_tf_idf():
    hv = HashingVectorizer(dim=1000, probes=3)
    hv.vectorize(JUNK_FOOD_DOCS)
    hv.vectorize(NOTJUNK_FOOD_DOCS)

    # extract the TF-IDF data
    X = hv.get_tfidf()
    assert_equal(X.shape, (11, 1000))

    # label junk food as -1, the others as +1
    y = np.ones(X.shape[0])
    y[:6] = -1

    # train and test a classifier
    clf = LogisticRegression().fit(X[1:-1], y[1:-1])
    assert_equal(clf.predict([X[0]]), [-1])
    assert_equal(clf.predict([X[-1]]), [1])


def test_sparse_tf_idf():
    hv = SparseHashingVectorizer(dim=1000000, probes=3)
    hv.vectorize(JUNK_FOOD_DOCS)
    hv.vectorize(NOTJUNK_FOOD_DOCS)

    # extract the TF-IDF data
    X = hv.get_tfidf()
    assert_equal(X.shape, (11, 1000000))

    # label junk food as -1, the others as +1
    y = np.ones(X.shape[0])
    y[:6] = -1

    # train and test a classifier
    clf = SVC(kernel='linear', C=10).fit(X[1:-1], y[1:-1])
    assert_equal(clf.predict(X[0, :]), [-1])
    assert_equal(clf.predict(X[-1, :]), [1])

def test_dense_sparse_idf_sanity():
    hv = HashingVectorizer(dim=100, probes=3)
    shv = SparseHashingVectorizer(dim=100, probes=3)

    hv.vectorize(JUNK_FOOD_DOCS)
    shv.vectorize(JUNK_FOOD_DOCS)

    # check that running TF IDF estimates are the same
    dense_tf_idf = hv.get_tfidf()
    sparse_tfidf = shv.get_tfidf().todense()

    assert_array_almost_equal(dense_tf_idf, sparse_tfidf)

    # check that incremental behaviour stays the same
    hv.vectorize(NOTJUNK_FOOD_DOCS)
    shv.vectorize(NOTJUNK_FOOD_DOCS)

    dense_tf_idf = hv.get_tfidf()
    sparse_tfidf = shv.get_tfidf().todense()

    assert_array_almost_equal(dense_tf_idf, sparse_tfidf)


