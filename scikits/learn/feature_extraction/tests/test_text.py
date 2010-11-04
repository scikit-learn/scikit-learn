from scikits.learn.feature_extraction.text import CharNGramAnalyzer
from scikits.learn.feature_extraction.text import WordNGramAnalyzer
from scikits.learn.feature_extraction.text import strip_accents

from scikits.learn.feature_extraction.text import CountVectorizer
from scikits.learn.feature_extraction.text import TfidfTransformer
from scikits.learn.feature_extraction.text import Vectorizer

import scikits.learn.feature_extraction.text.sparse as st

SparseCountVectorizer = st.CountVectorizer
SparseTfidfTransformer = st.TfidfTransformer
SparseVectorizer = st.Vectorizer

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
    wa = WordNGramAnalyzer(min_n=1, max_n=1, stop_words=None)

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u'ai', u'mange', u'du', u'kangourou', u'ce', u'midi',
                u'etait', u'pas', u'tres', u'bon']
    assert_equal(wa.analyze(text), expected)

    text = "This is a test, really.\n\n I met Harry yesterday."
    expected = [u'this', u'is', u'test', u'really', u'met', u'harry',
                u'yesterday']
    assert_equal(wa.analyze(text), expected)


def test_word_analyzer_unigrams_and_bigrams():
    wa = WordNGramAnalyzer(min_n=1, max_n=2, stop_words=None)

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
    grid_search = GridSearchCV(pipeline, parameters, n_jobs=1)

    # cross-validation doesn't work if the length of the data is not known,
    # hence use lists instead of iterators
    pred = grid_search.fit(list(train_data), y_train).predict(list(test_data))
    assert_array_equal(pred, y_test)

    # on this toy dataset bigram representation which is used in the last of the
    # grid_search is considered the best estimator since they all converge to
    # 100% accurracy models
    assert_equal(grid_search.best_score, 1.0)
    best_vectorizer = grid_search.best_estimator.named_steps['vect']
    assert_equal(best_vectorizer.analyzer.max_n, 2)


def test_pickle():
    for obj in (CountVectorizer(), SparseCountVectorizer(),
                TfidfTransformer(), SparseTfidfTransformer(),
                Vectorizer(), SparseVectorizer()):

        s = pickle.dumps(obj)
        assert_equal(type(pickle.loads(s)), obj.__class__)

