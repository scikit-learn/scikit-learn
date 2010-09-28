from scikits.learn.features.text import strip_accents
from scikits.learn.features.text import WordNGramAnalyzer
from scikits.learn.features.text import CharNGramAnalyzer
from scikits.learn.features.text import TermCountVectorizer
from scikits.learn.features.text import TfidfTransformer
from scikits.learn.features.text import TfidfVectorizer
from scikits.learn.features.text import HashingVectorizer
from scikits.learn.features.text import SparseHashingVectorizer
from scikits.learn.svm import LinearSVC as DenseLinearSVC
from scikits.learn.svm.sparse import LinearSVC as SparseLinearSVC
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


def test_dense_vectorizer():
    wa = WordNGramAnalyzer()
    train_data = [wa.analyze(d) for d in JUNK_FOOD_DOCS[:-1]]
    test_data = [wa.analyze(JUNK_FOOD_DOCS[-1])]

    # test without vocabulary
    v1 = TermCountVectorizer()
    counts_train = v1.transform(train_data)
    assert_equal(counts_train[0, v1.vocabulary["pizza"]], 2)

    v2 = TermCountVectorizer(vocabulary=v1.vocabulary)

    # test with a pre-existing vocabulary
    for v in (v1, v2):
        counts_test = v.transform(test_data)
        assert_equal(counts_test[0, v.vocabulary["coke"]], 1)

    # test tf-idf
    t1 = TfidfTransformer()
    tfidf = t1.fit(counts_train).transform(counts_train)
    assert_equal(len(t1.idf), len(v1.vocabulary))
    assert_equal(tfidf.shape,
                 (len(train_data), len(v1.vocabulary)))

    # test tf-idf with new data
    tfidf_test = t1.transform(counts_test)
    assert_equal(tfidf_test.shape,
                 (len(test_data), len(v1.vocabulary)))

    # test tf alone
    t2 = TfidfTransformer(use_idf=False)
    tf = t2.fit(counts_train).transform(counts_train)
    assert_equal(t2.idf, None)
    assert_array_almost_equal(np.sum(tf, axis=1),
                              [1.0] * len(train_data))

    # test the direct tfidf vectorizer
    # (equivalent to term count vectorizer + tfidf transformer)
    tv = TfidfVectorizer()
    tfidf2 = tv.fit(train_data).transform(train_data)
    assert_array_almost_equal(tfidf, tfidf2)

    # test the direct tfidf vectorizer with new data
    tfidf_test2 = tv.transform(test_data)
    assert_array_almost_equal(tfidf_test, tfidf_test2)
