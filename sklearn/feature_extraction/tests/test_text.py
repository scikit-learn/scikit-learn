from sklearn.feature_extraction.text import strip_accents_unicode
from sklearn.feature_extraction.text import strip_accents_ascii

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.feature_extraction.text import Vectorizer

from sklearn.grid_search import GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.svm import LinearSVC

import numpy as np
from nose.tools import assert_equal, assert_equals, \
            assert_false, assert_not_equal, assert_true
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_raises

import pickle
from StringIO import StringIO

JUNK_FOOD_DOCS = (
    "the pizza pizza beer copyright",
    "the pizza burger beer copyright",
    "the the pizza beer beer copyright",
    "the burger beer beer copyright",
    "the coke burger coke copyright",
    "the coke burger burger",
)

NOTJUNK_FOOD_DOCS = (
    "the salad celeri copyright",
    "the salad salad sparkling water copyright",
    "the the celeri celeri copyright",
    "the tomato tomato salad water",
    "the tomato salad water copyright",
)

ALL_FOOD_DOCS = JUNK_FOOD_DOCS + NOTJUNK_FOOD_DOCS


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def test_strip_accents():
    # check some classical latin accentuated symbols
    a = u'\xe0\xe1\xe2\xe3\xe4\xe5\xe7\xe8\xe9\xea\xeb'
    expected = u'aaaaaaceeee'
    assert_equal(strip_accents_unicode(a), expected)

    a = u'\xec\xed\xee\xef\xf1\xf2\xf3\xf4\xf5\xf6\xf9\xfa\xfb\xfc\xfd'
    expected = u'iiiinooooouuuuy'
    assert_equal(strip_accents_unicode(a), expected)

    # check some arabic
    a = u'\u0625'  # halef with a hamza below
    expected = u'\u0627'  # simple halef
    assert_equal(strip_accents_unicode(a), expected)

    # mix letters accentuated and not
    a = u"this is \xe0 test"
    expected = u'this is a test'
    assert_equal(strip_accents_unicode(a), expected)


def test_to_ascii():
    # check some classical latin accentuated symbols
    a = u'\xe0\xe1\xe2\xe3\xe4\xe5\xe7\xe8\xe9\xea\xeb'
    expected = u'aaaaaaceeee'
    assert_equal(strip_accents_ascii(a), expected)

    a = u'\xec\xed\xee\xef\xf1\xf2\xf3\xf4\xf5\xf6\xf9\xfa\xfb\xfc\xfd'
    expected = u'iiiinooooouuuuy'
    assert_equal(strip_accents_ascii(a), expected)

    # check some arabic
    a = u'\u0625'  # halef with a hamza below
    expected = u''  # halef has no direct ascii match
    assert_equal(strip_accents_ascii(a), expected)

    # mix letters accentuated and not
    a = u"this is \xe0 test"
    expected = u'this is a test'
    assert_equal(strip_accents_ascii(a), expected)


def test_word_analyzer_unigrams():
    wa = CountVectorizer().build_analyzer()
    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u'ai', u'mange', u'du', u'kangourou', u'ce', u'midi',
                u'etait', u'pas', u'tres', u'bon']
    assert_equal(wa(text), expected)

    text = "This is a test, really.\n\n I met Harry yesterday."
    expected = [u'this', u'is', u'test', u'really', u'met', u'harry',
                u'yesterday']
    assert_equal(wa(text), expected)

    wa = CountVectorizer(input='file').build_analyzer()
    text = StringIO("This is a test with a file-like object!")
    expected = [u'this', u'is', u'test', u'with', u'file', u'like',
                u'object']
    assert_equal(wa(text), expected)


def test_word_analyzer_unigrams_and_bigrams():
    wa = CountVectorizer(tokenize="word", min_n=1, max_n=2).build_analyzer()

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    expected = [u'ai', u'mange', u'du', u'kangourou', u'ce', u'midi', u'etait',
                u'pas', u'tres', u'bon', u'ai mange', u'mange du',
                u'du kangourou', u'kangourou ce', u'ce midi', u'midi etait',
                u'etait pas', u'pas tres', u'tres bon']
    assert_equal(wa(text), expected)


def test_unicode_decode_error():
    # decode_error default to strict, so this should fail
    # First, encode (as bytes) a unicode string.
    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon."
    text_bytes = text.encode('utf-8')

    # Then let the Analyzer try to decode it as ascii. It should fail,
    # because we have given it an incorrect charset.
    wa = CountVectorizer(min_n=1, max_n=2, charset='ascii').build_analyzer()
    assert_raises(UnicodeDecodeError, wa, text_bytes)

    ca = CountVectorizer(tokenize='char', min_n=3, max_n=6,
                         charset='ascii').build_analyzer()
    assert_raises(UnicodeDecodeError, ca, text_bytes)


def test_char_ngram_analyzer():
    cnga = CountVectorizer(tokenize='char', min_n=3, max_n=6).build_analyzer()

    text = u"J'ai mang\xe9 du kangourou  ce midi, c'\xe9tait pas tr\xeas bon"
    expected = [u"j'a", u"'ai", u'ai ', u'i m', u' ma']
    assert_equal(cnga(text)[:5], expected)
    expected = [u's tres', u' tres ', u'tres b', u'res bo', u'es bon']
    assert_equal(cnga(text)[-5:], expected)

    text = "This \n\tis a test, really.\n\n I met Harry yesterday"
    expected = [u'thi', u'his', u'is ', u's i', u' is']
    assert_equal(cnga(text)[:5], expected)

    expected = [u' yeste', u'yester', u'esterd', u'sterda', u'terday']
    assert_equal(cnga(text)[-5:], expected)

    cnga = CountVectorizer(input='file', tokenize='char',
                           min_n=3, max_n=6).build_analyzer()
    text = StringIO("This is a test with a file-like object!")
    expected = [u'thi', u'his', u'is ', u's i', u' is']
    assert_equal(cnga(text)[:5], expected)


def test_countvectorizer_custom_vocabulary():
    what_we_like = ["pizza", "beer"]
    vect = CountVectorizer(fixed_vocabulary=what_we_like)
    vect.fit(JUNK_FOOD_DOCS)
    assert_equal(set(vect.fixed_vocabulary), set(what_we_like))
    X = vect.transform(JUNK_FOOD_DOCS)
    assert_equal(X.shape[1], len(what_we_like))

    # try again with a dict vocabulary
    vocab = {"pizza": 0, "beer": 1}
    vect = CountVectorizer(fixed_vocabulary=vocab)
    vect.fit(JUNK_FOOD_DOCS)
    assert_equal(vect.fixed_vocabulary, vocab)
    X = vect.transform(JUNK_FOOD_DOCS)
    assert_equal(X.shape[1], len(what_we_like))


def test_countvectorizer_custom_vocabulary_pipeline():
    what_we_like = ["pizza", "beer"]
    pipe = Pipeline([
        ('count', CountVectorizer(fixed_vocabulary=what_we_like)),
        ('tfidf', TfidfTransformer())])
    X = pipe.fit_transform(ALL_FOOD_DOCS)
    assert_equal(set(pipe.named_steps['count'].fixed_vocabulary),
                 set(what_we_like))
    assert_equal(X.shape[1], len(what_we_like))


def test_fit_countvectorizer_twice():
    cv = CountVectorizer()
    X1 = cv.fit_transform(ALL_FOOD_DOCS[:5])
    X2 = cv.fit_transform(ALL_FOOD_DOCS[5:])
    assert_not_equal(X1.shape[1], X2.shape[1])


def test_sublinear_tf():
    X = [[1], [2], [3]]
    tr = TfidfTransformer(sublinear_tf=True, use_idf=False, norm=None)
    tfidf = toarray(tr.fit_transform(X))
    assert_equal(tfidf[0], 1)
    assert_true(tfidf[1] > tfidf[0])
    assert_true(tfidf[2] > tfidf[1])
    assert_true(tfidf[1] < 2)
    assert_true(tfidf[2] < 3)


def test_vectorizer():
    # raw documents as an iterator
    train_data = iter(ALL_FOOD_DOCS[:-1])
    test_data = [ALL_FOOD_DOCS[-1]]
    n_train = len(ALL_FOOD_DOCS) - 1

    # test without vocabulary
    v1 = CountVectorizer(max_df=0.5)
    counts_train = v1.fit_transform(train_data)
    if hasattr(counts_train, 'tocsr'):
        counts_train = counts_train.tocsr()
    assert_equal(counts_train[0, v1.vocabulary_[u"pizza"]], 2)

    # build a vectorizer v1 with the same vocabulary as the one fitted by v1
    v2 = CountVectorizer(fixed_vocabulary=v1.vocabulary_)

    # compare that the two vectorizer give the same output on the test sample
    for v in (v1, v2):
        counts_test = v.transform(test_data)
        if hasattr(counts_test, 'tocsr'):
            counts_test = counts_test.tocsr()

        vocabulary = v.get_vocabulary()
        assert_equal(counts_test[0, vocabulary[u"salad"]], 1)
        assert_equal(counts_test[0, vocabulary[u"tomato"]], 1)
        assert_equal(counts_test[0, vocabulary[u"water"]], 1)

        # stop word from the fixed list
        assert_false(u"the" in vocabulary)

        # stop word found automatically by the vectorizer DF thresholding
        # words that are high frequent across the complete corpus are likely
        # to be not informative (either real stop words of extraction
        # artifacts)
        assert_false(u"copyright" in vocabulary)

        # not present in the sample
        assert_equal(counts_test[0, vocabulary[u"coke"]], 0)
        assert_equal(counts_test[0, vocabulary[u"burger"]], 0)
        assert_equal(counts_test[0, vocabulary[u"beer"]], 0)
        assert_equal(counts_test[0, vocabulary[u"pizza"]], 0)

    # test tf-idf
    t1 = TfidfTransformer(norm='l1')
    tfidf = toarray(t1.fit(counts_train).transform(counts_train))
    assert_equal(len(t1.idf_), len(v1.vocabulary_))
    assert_equal(tfidf.shape, (n_train, len(v1.vocabulary_)))

    # test tf-idf with new data
    tfidf_test = toarray(t1.transform(counts_test))
    assert_equal(tfidf_test.shape, (len(test_data), len(v1.vocabulary_)))

    # test tf alone
    t2 = TfidfTransformer(norm='l1', use_idf=False)
    tf = toarray(t2.fit(counts_train).transform(counts_train))
    assert_equal(t2.idf_, None)

    # L1-normalized term frequencies sum to one
    assert_array_almost_equal(np.sum(tf, axis=1), [1.0] * n_train)

    # test the direct tfidf vectorizer
    # (equivalent to term count vectorizer + tfidf transformer)
    train_data = iter(ALL_FOOD_DOCS[:-1])
    tv = Vectorizer(norm='l1')
    tv.max_df = v1.max_df
    tfidf2 = toarray(tv.fit_transform(train_data))
    assert_array_almost_equal(tfidf, tfidf2)

    # test the direct tfidf vectorizer with new data
    tfidf_test2 = toarray(tv.transform(test_data))
    assert_array_almost_equal(tfidf_test, tfidf_test2)

    # test transform on unfitted vectorizer with empty vocabulary
    v3 = CountVectorizer(fixed_vocabulary=None)
    assert_raises(ValueError, v3.transform, train_data)


def test_feature_names():
    cv = CountVectorizer(max_df=0.5)
    X = cv.fit_transform(ALL_FOOD_DOCS)

    n_samples, n_features = X.shape
    assert_equal(len(cv.vocabulary_), n_features)

    feature_names = cv.get_feature_names()
    assert_equal(feature_names.shape, (n_features,))
    assert_array_equal(['beer', 'burger', 'celeri', 'coke', 'pizza',
                        'salad', 'sparkling', 'tomato', 'water'],
                       feature_names)

    for idx, name in enumerate(feature_names):
        assert_equal(idx, cv.vocabulary_.get(name))


def test_vectorizer_max_features():
    vec_factories = (
        CountVectorizer,
        Vectorizer,
    )

    expected_vocabulary = set(['burger', 'beer', 'salad', 'pizza'])

    for vec_factory in vec_factories:
        # test bounded number of extracted features
        vectorizer = vec_factory(max_df=0.6, max_features=4)
        vectorizer.fit(ALL_FOOD_DOCS)
        assert_equals(set(vectorizer.vocabulary_), expected_vocabulary)


def test_vectorizer_max_df():
    test_data = [u'abc', u'dea']  # the letter a occurs in both strings
    vect = CountVectorizer(tokenize='char', min_n=1, max_n=1, max_df=1.0)
    vect.fit(test_data)
    assert_true(u'a' in vect.vocabulary_.keys())
    assert_equals(len(vect.vocabulary_.keys()), 5)

    vect.max_df = 0.5
    vect.fit(test_data)
    assert_true(u'a' not in vect.vocabulary_.keys())  # 'a' is ignored
    assert_equals(len(vect.vocabulary_.keys()), 4)  # the others remain


def test_binary_occurrences():
    # by default multiple occurrences are counted as longs
    test_data = [u'aaabc', u'abbde']
    vect = CountVectorizer(tokenize='char', min_n=1, max_n=1, max_df=1.0)
    X = vect.fit_transform(test_data).toarray()
    assert_array_equal(['a', 'b', 'c', 'd', 'e'], vect.get_feature_names())
    assert_array_equal([[3, 1, 1, 0, 0],
                        [1, 2, 0, 1, 1]], X)

    # using boolean features, we can fetch the binary occurrence info
    # instead.
    vect = CountVectorizer(tokenize='char', min_n=1, max_n=1, max_df=1.0,
                           binary=True)
    X = vect.fit_transform(test_data).toarray()
    assert_array_equal([[1, 1, 1, 0, 0],
                        [1, 1, 0, 1, 1]], X)

    # check the ability to change the dtype
    vect = CountVectorizer(tokenize='char', min_n=1, max_n=1, max_df=1.0,
                           binary=True, dtype=np.float32)
    X_sparse = vect.fit_transform(test_data)
    assert_equal(X_sparse.dtype, np.float32)


def test_vectorizer_inverse_transform():
    # raw documents
    data = ALL_FOOD_DOCS
    for vectorizer in (Vectorizer(), CountVectorizer()):
        transformed_data = vectorizer.fit_transform(data)
        inversed_data = vectorizer.inverse_transform(transformed_data)
        analyze = vectorizer.build_analyzer()
        for i, doc in enumerate(data):
            terms = np.sort(np.unique(analyze(doc)))
            inversed_terms = np.sort(np.unique(inversed_data[i]))
            assert_array_equal(terms, inversed_terms)

    # Test that inverse_transform also works with numpy arrays
    transformed_data = transformed_data.toarray()
    inversed_data2 = vectorizer.inverse_transform(transformed_data)
    for terms, terms2 in zip(inversed_data, inversed_data2):
        assert_array_equal(terms, terms2)


def test_vectorizer_pipeline_grid_selection():
    # raw documents
    data = JUNK_FOOD_DOCS + NOTJUNK_FOOD_DOCS
    # simulate iterables
    train_data = iter(data[1:-1])
    test_data = iter([data[0], data[-1]])

    # label junk food as -1, the others as +1
    y = np.ones(len(data))
    y[:6] = -1
    y_train = y[1:-1]
    y_test = np.array([y[0], y[-1]])

    pipeline = Pipeline([('vect', CountVectorizer()),
                         ('svc', LinearSVC())])

    parameters = {
        'vect__max_n': (1, 2),
        'svc__loss': ('l1', 'l2')
    }

    # find the best parameters for both the feature extraction and the
    # classifier
    grid_search = GridSearchCV(pipeline, parameters, n_jobs=1)

    # cross-validation doesn't work if the length of the data is not known,
    # hence use lists instead of iterators
    pred = grid_search.fit(list(train_data), y_train).predict(list(test_data))
    assert_array_equal(pred, y_test)

    # on this toy dataset bigram representation which is used in the last of
    # the grid_search is considered the best estimator since they all converge
    # to 100% accuracy models
    assert_equal(grid_search.best_score_, 1.0)
    best_vectorizer = grid_search.best_estimator_.named_steps['vect']
    assert_equal(best_vectorizer.max_n, 1)


def test_pickle():
    for obj in (CountVectorizer(), TfidfTransformer(), Vectorizer()):
        s = pickle.dumps(obj)
        assert_equal(type(pickle.loads(s)), obj.__class__)
