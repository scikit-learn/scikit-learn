from scikits.learn.feature_extraction.text import CharNGramAnalyzer
from scikits.learn.feature_extraction.text import CountVectorizer
from scikits.learn.feature_extraction.text import HashingVectorizer
from scikits.learn.feature_extraction.text import TfidfTransformer
from scikits.learn.feature_extraction.text import Vectorizer
from scikits.learn.feature_extraction.sparse.text import SparseCountVectorizer
from scikits.learn.feature_extraction.sparse.text import SparseHashingVectorizer
from scikits.learn.feature_extraction.sparse.text import SparseTfidfTransformer
from scikits.learn.feature_extraction.sparse.text import SparseVectorizer
from scikits.learn.feature_extraction.text import WordNGramAnalyzer
from scikits.learn.feature_extraction.text import strip_accents
from scikits.learn.grid_search import GridSearchCV
from scikits.learn.pipeline import Pipeline
from scikits.learn.svm import LinearSVC as DenseLinearSVC
from scikits.learn.svm.sparse import LinearSVC as SparseLinearSVC

import numpy as np
import numpy.linalg as la
from nose.tools import *
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

import pickle

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


def test_word_analyzer_unigrams():
    wa = WordNGramAnalyzer(min_n=1, max_n=1)

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u'ai', u'mange', u'du', u'kangourou', u'ce', u'midi',
                u'etait', u'pas', u'tres', u'bon']
    assert_equal(wa.analyze(text), expected)

    text = "This is a test, really.\n\n I met Harry yesterday."
    expected = [u'this', u'is', u'test', u'really', u'met', u'harry',
                u'yesterday']
    assert_equal(wa.analyze(text), expected)


def test_word_analyzer_unigrams_and_bigrams():
    wa = WordNGramAnalyzer(min_n=1, max_n=2)

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u'ai', u'mange', u'du', u'kangourou', u'ce', u'midi', u'etait',
                u'pas', u'tres', u'bon', u'ai mange', u'mange du',
                u'du kangourou', u'kangourou ce', u'ce midi', u'midi etait',
                u'etait pas', u'pas tres', u'tres bon']
    assert_equal(wa.analyze(text), expected)


def test_char_ngram_analyzer():
    cnga = CharNGramAnalyzer(min_n=3, max_n=6)

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u"j'a", u"'ai", u'ai ', u'i m', u' ma']
    assert_equal(cnga.analyze(text)[:5], expected)
    expected = [u's tres', u' tres ', u'tres b', u'res bo', u'es bon']
    assert_equal(cnga.analyze(text)[-5:], expected)

    text = "This \n\tis a test, really.\n\n I met Harry yesterday."
    expected = [u'thi', u'his', u'is ', u's i', u' is']
    assert_equal(cnga.analyze(text)[:5], expected)
    expected = [u' yeste', u'yester', u'esterd', u'sterda', u'terday']
    assert_equal(cnga.analyze(text)[-5:], expected)


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
    clf = DenseLinearSVC(C=10).fit(X[1:-1], y[1:-1])
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
    clf = SparseLinearSVC(C=10).fit(X[1:-1], y[1:-1])
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

def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a

def _test_vectorizer(cv_class, tf_class, v_class):
    # results to be compared
    res = []

    # raw documents
    train_data = iter(JUNK_FOOD_DOCS[:-1])
    n_train = len(JUNK_FOOD_DOCS[:-1])
    test_data = [JUNK_FOOD_DOCS[-1]]

    # test without vocabulary
    v1 = cv_class()
    counts_train = v1.fit_transform(train_data)
    assert_equal(counts_train[0, v1.vocabulary["pizza"]], 2)

    v2 = cv_class(vocabulary=v1.vocabulary)

    # test with a pre-existing vocabulary
    for v in (v1, v2):
        counts_test = v.transform(test_data)
        assert_equal(counts_test[0, v.vocabulary["coke"]], 1)

    # test tf-idf
    t1 = tf_class()
    tfidf = toarray(t1.fit(counts_train).transform(counts_train))
    assert_equal(len(t1.idf), len(v1.vocabulary))
    assert_equal(tfidf.shape,
                 (n_train, len(v1.vocabulary)))

    res.append(tfidf)
    res.append(t1.idf)

    # test tf-idf with new data
    tfidf_test = toarray(t1.transform(counts_test))
    assert_equal(tfidf_test.shape,
                 (len(test_data), len(v1.vocabulary)))

    # test tf alone
    t2 = tf_class(use_idf=False)
    tf = toarray(t2.fit(counts_train).transform(counts_train))
    assert_equal(t2.idf, None)
    assert_array_almost_equal(np.sum(tf, axis=1),
                              [1.0] * n_train)

    # test normalization
    t3 = tf_class(normalize=True)
    tfidf_n = toarray(t3.fit(counts_train).transform(counts_train))
    assert_equal(la.norm(tfidf_n[0]), 1.0)

    # test the direct tfidf vectorizer
    # (equivalent to term count vectorizer + tfidf transformer)
    train_data = iter(JUNK_FOOD_DOCS[:-1])
    tv = v_class()
    tfidf2 = toarray(tv.fit_transform(train_data))
    assert_array_almost_equal(tfidf, tfidf2)

    # test the direct tfidf vectorizer with new data
    tfidf_test2 = toarray(tv.transform(test_data))
    assert_array_almost_equal(tfidf_test, tfidf_test2)

    return res

def test_vectorizer():
    res_dense =_test_vectorizer(CountVectorizer,
                                TfidfTransformer,
                                Vectorizer)
    res_sparse = _test_vectorizer(SparseCountVectorizer,
                                  SparseTfidfTransformer,
                                  SparseVectorizer)

    for i in xrange(len(res_sparse)):
        # check that the dense and sparse implementations
        # return the same results
        assert_array_equal(res_dense[i], res_sparse[i])

def test_dense_vectorizer_pipeline_grid_selection():
    # raw documents
    data = JUNK_FOOD_DOCS + NOTJUNK_FOOD_DOCS
    # simulate iterables
    train_data = iter(data[1:-1])
    test_data = iter([data[0], data[-1]])

    # label junk food as -1, the others as +1
    y = np.ones(len(data))
    y[:6] = -1
    y_train = y[1:-1]
    y_test = np.array([y[0],y[-1]])

    pipeline = Pipeline([('vect', CountVectorizer()),
                         ('svc', DenseLinearSVC())])

    parameters = {
        'vect__analyzer': (WordNGramAnalyzer(min_n=1, max_n=1),
                           WordNGramAnalyzer(min_n=1, max_n=2)),
        'svc__loss'  : ('l1', 'l2')
    }


    # find the best parameters for both the feature extraction and the
    # classifier
    clf = GridSearchCV(pipeline, parameters, n_jobs=1)

    # cross-validation doesn't work if the length of the data is not known,
    # hence use lists instead of iterators
    pred = clf.fit(list(train_data), y_train).predict(list(test_data))
    assert_array_equal(pred, y_test)

    # check that the bigram representation yields higher predictive accurracy
    # this test is unstable...
    #assert_equal(clf.best_estimator.steps[0][1].analyzer.max_n, 2)

def test_pickle():
    for obj in (CountVectorizer(), SparseCountVectorizer(),
                TfidfTransformer(), SparseTfidfTransformer(),
                Vectorizer(), SparseVectorizer()):

        s = pickle.dumps(obj)
        assert_equal(type(pickle.loads(s)), obj.__class__)
