# -*- coding: utf-8 -*-
from __future__ import division

import re

import numpy as np
from scipy import sparse
import pytest

from sklearn.exceptions import NotFittedError
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_no_warnings

from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder


def toarray(a):
    if hasattr(a, "toarray"):
        a = a.toarray()
    return a


def test_one_hot_encoder_sparse():
    # Test OneHotEncoder's fit and transform.
    X = [[3, 2, 1], [0, 1, 1]]
    enc = OneHotEncoder()
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        # discover max values automatically
        X_trans = enc.fit_transform(X).toarray()
        assert_equal(X_trans.shape, (2, 5))
        assert_array_equal(enc.active_features_,
                           np.where([1, 0, 0, 1, 0, 1, 1, 0, 1])[0])
        assert_array_equal(enc.feature_indices_, [0, 4, 7, 9])

        # check outcome
        assert_array_equal(X_trans,
                           [[0., 1., 0., 1., 1.],
                            [1., 0., 1., 0., 1.]])

    # max value given as 3
    # enc = assert_warns(DeprecationWarning, OneHotEncoder, n_values=4)
    enc = OneHotEncoder(n_values=4)
    with ignore_warnings(category=DeprecationWarning):
        X_trans = enc.fit_transform(X)
        assert_equal(X_trans.shape, (2, 4 * 3))
        assert_array_equal(enc.feature_indices_, [0, 4, 8, 12])

    # max value given per feature
    # enc = assert_warns(DeprecationWarning, OneHotEncoder, n_values=[3, 2, 2])
    enc = OneHotEncoder(n_values=[3, 2, 2])
    with ignore_warnings(category=DeprecationWarning):
        X = [[1, 0, 1], [0, 1, 1]]
        X_trans = enc.fit_transform(X)
        assert_equal(X_trans.shape, (2, 3 + 2 + 2))
        assert_array_equal(enc.n_values_, [3, 2, 2])
    # check that testing with larger feature works:
    X = np.array([[2, 0, 1], [0, 1, 1]])
    enc.transform(X)

    # test that an error is raised when out of bounds:
    X_too_large = [[0, 2, 1], [0, 1, 1]]
    assert_raises(ValueError, enc.transform, X_too_large)
    error_msg = r"unknown categorical feature present \[2\] during transform"
    assert_raises_regex(ValueError, error_msg, enc.transform, X_too_large)
    with ignore_warnings(category=DeprecationWarning):
        assert_raises(
            ValueError,
            OneHotEncoder(n_values=2).fit_transform, X)

    # test that error is raised when wrong number of features
    assert_raises(ValueError, enc.transform, X[:, :-1])

    # test that error is raised when wrong number of features in fit
    # with prespecified n_values
    with ignore_warnings(category=DeprecationWarning):
        assert_raises(ValueError, enc.fit, X[:, :-1])
    # test exception on wrong init param
    with ignore_warnings(category=DeprecationWarning):
        assert_raises(
            TypeError, OneHotEncoder(n_values=np.int).fit, X)

    enc = OneHotEncoder()
    # test negative input to fit
    with ignore_warnings(category=FutureWarning):
        assert_raises(ValueError, enc.fit, [[0], [-1]])

    # test negative input to transform
    with ignore_warnings(category=FutureWarning):
        enc.fit([[0], [1]])
    assert_raises(ValueError, enc.transform, [[0], [-1]])


def test_one_hot_encoder_dense():
    # check for sparse=False
    X = [[3, 2, 1], [0, 1, 1]]
    enc = OneHotEncoder(sparse=False)
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        # discover max values automatically
        X_trans = enc.fit_transform(X)
        assert_equal(X_trans.shape, (2, 5))
        assert_array_equal(enc.active_features_,
                           np.where([1, 0, 0, 1, 0, 1, 1, 0, 1])[0])
        assert_array_equal(enc.feature_indices_, [0, 4, 7, 9])

    # check outcome
    assert_array_equal(X_trans,
                       np.array([[0., 1., 0., 1., 1.],
                                 [1., 0., 1., 0., 1.]]))


def test_one_hot_encoder_deprecationwarnings():
    for X in [[[3, 2, 1], [0, 1, 1]],
              [[3., 2., 1.], [0., 1., 1.]]]:
        enc = OneHotEncoder()
        assert_warns_message(FutureWarning, "handling of integer",
                             enc.fit, X)
        enc = OneHotEncoder()
        assert_warns_message(FutureWarning, "handling of integer",
                             enc.fit_transform, X)

        # check it still works correctly as well
        with ignore_warnings(category=FutureWarning):
            X_trans = enc.fit_transform(X).toarray()
        res = [[0., 1., 0., 1., 1.],
               [1., 0., 1., 0., 1.]]
        assert_array_equal(X_trans, res)

        # check deprecated attributes
        assert_warns(DeprecationWarning, lambda: enc.active_features_)
        assert_warns(DeprecationWarning, lambda: enc.feature_indices_)
        assert_warns(DeprecationWarning, lambda: enc.n_values_)

        # check no warning is raised if keyword is specified
        enc = OneHotEncoder(categories='auto')
        assert_no_warnings(enc.fit, X)
        enc = OneHotEncoder(categories='auto')
        assert_no_warnings(enc.fit_transform, X)
        X_trans = enc.fit_transform(X).toarray()
        assert_array_equal(X_trans, res)

        # check there is also a warning if the default is passed
        enc = OneHotEncoder(n_values='auto', handle_unknown='ignore')
        assert_warns(DeprecationWarning, enc.fit, X)

    X = np.array([['cat1', 'cat2']], dtype=object).T
    enc = OneHotEncoder(categorical_features='all')
    assert_warns(DeprecationWarning, enc.fit, X)


def test_one_hot_encoder_force_new_behaviour():
    # ambiguous integer case (non secutive range of categories)
    X = np.array([[1, 2]]).T
    X2 = np.array([[0, 1]]).T

    # without argument -> by default using legacy behaviour with warnings
    enc = OneHotEncoder()

    with ignore_warnings(category=FutureWarning):
        enc.fit(X)

    res = enc.transform(X2)
    exp = np.array([[0, 0], [1, 0]])
    assert_array_equal(res.toarray(), exp)

    # with explicit auto argument -> don't use legacy behaviour
    # (so will raise an error on unseen value within range)
    enc = OneHotEncoder(categories='auto')
    enc.fit(X)
    assert_raises(ValueError, enc.transform, X2)


def _run_one_hot(X, X2, cat):
    # enc = assert_warns(
    #     DeprecationWarning,
    #     OneHotEncoder, categorical_features=cat)
    enc = OneHotEncoder(categorical_features=cat)
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        Xtr = enc.fit_transform(X)
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        X2tr = enc.fit(X).transform(X2)
    return Xtr, X2tr


def _check_one_hot(X, X2, cat, n_features):
    ind = np.where(cat)[0]
    # With mask
    A, B = _run_one_hot(X, X2, cat)
    # With indices
    C, D = _run_one_hot(X, X2, ind)
    # Check shape
    assert_equal(A.shape, (2, n_features))
    assert_equal(B.shape, (1, n_features))
    assert_equal(C.shape, (2, n_features))
    assert_equal(D.shape, (1, n_features))
    # Check that mask and indices give the same results
    assert_array_equal(toarray(A), toarray(C))
    assert_array_equal(toarray(B), toarray(D))


def test_one_hot_encoder_categorical_features():
    X = np.array([[3, 2, 1], [0, 1, 1]])
    X2 = np.array([[1, 1, 1]])

    cat = [True, False, False]
    _check_one_hot(X, X2, cat, 4)

    # Edge case: all non-categorical
    cat = [False, False, False]
    _check_one_hot(X, X2, cat, 3)

    # Edge case: all categorical
    cat = [True, True, True]
    _check_one_hot(X, X2, cat, 5)

    # check error raised if also specifying categories
    oh = OneHotEncoder(categories=[range(3)],
                       categorical_features=[True, False, False])
    assert_raises(ValueError, oh.fit, X)


def test_one_hot_encoder_categorical_features_ignore_unknown():
    # GH12881 bug in combination of categorical_features with ignore
    X = np.array([[1, 2, 3], [4, 5, 6], [2, 3, 2]]).T
    oh = OneHotEncoder(categorical_features=[2], handle_unknown='ignore')

    with ignore_warnings(category=DeprecationWarning):
        res = oh.fit_transform(X)

    expected = np.array([[1, 0, 1], [0, 1, 0], [1, 2, 3], [4, 5, 6]]).T
    assert_array_equal(res.toarray(), expected)


def test_one_hot_encoder_handle_unknown():
    X = np.array([[0, 2, 1], [1, 0, 3], [1, 0, 2]])
    X2 = np.array([[4, 1, 1]])

    # Test that one hot encoder raises error for unknown features
    # present during transform.
    oh = OneHotEncoder(handle_unknown='error')
    assert_warns(FutureWarning, oh.fit, X)
    assert_raises(ValueError, oh.transform, X2)

    # Test the ignore option, ignores unknown features (giving all 0's)
    oh = OneHotEncoder(handle_unknown='ignore')
    oh.fit(X)
    X2_passed = X2.copy()
    assert_array_equal(
        oh.transform(X2_passed).toarray(),
        np.array([[0.,  0.,  0.,  0.,  1.,  0.,  0.]]))
    # ensure transformed data was not modified in place
    assert_allclose(X2, X2_passed)

    # Raise error if handle_unknown is neither ignore or error.
    oh = OneHotEncoder(handle_unknown='42')
    assert_raises(ValueError, oh.fit, X)


def test_one_hot_encoder_not_fitted():
    X = np.array([['a'], ['b']])
    enc = OneHotEncoder(categories=['a', 'b'])
    msg = ("This OneHotEncoder instance is not fitted yet. "
           "Call 'fit' with appropriate arguments before using this method.")
    with pytest.raises(NotFittedError, match=msg):
        enc.transform(X)


def test_one_hot_encoder_no_categorical_features():
    X = np.array([[3, 2, 1], [0, 1, 1]], dtype='float64')

    cat = [False, False, False]
    enc = OneHotEncoder(categorical_features=cat)
    with ignore_warnings(category=(DeprecationWarning, FutureWarning)):
        X_tr = enc.fit_transform(X)
    expected_features = np.array(list(), dtype='object')
    assert_array_equal(X, X_tr)
    assert_array_equal(enc.get_feature_names(), expected_features)
    assert enc.categories_ == []


def test_one_hot_encoder_handle_unknown_strings():
    X = np.array(['11111111', '22', '333', '4444']).reshape((-1, 1))
    X2 = np.array(['55555', '22']).reshape((-1, 1))
    # Non Regression test for the issue #12470
    # Test the ignore option, when categories are numpy string dtype
    # particularly when the known category strings are larger
    # than the unknown category strings
    oh = OneHotEncoder(handle_unknown='ignore')
    oh.fit(X)
    X2_passed = X2.copy()
    assert_array_equal(
        oh.transform(X2_passed).toarray(),
        np.array([[0.,  0.,  0.,  0.], [0.,  1.,  0.,  0.]]))
    # ensure transformed data was not modified in place
    assert_array_equal(X2, X2_passed)


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
@pytest.mark.parametrize("input_dtype", [np.int32, np.float32, np.float64])
def test_one_hot_encoder_dtype(input_dtype, output_dtype):
    X = np.asarray([[0, 1]], dtype=input_dtype).T
    X_expected = np.asarray([[1, 0], [0, 1]], dtype=output_dtype)

    oh = OneHotEncoder(categories='auto', dtype=output_dtype)
    assert_array_equal(oh.fit_transform(X).toarray(), X_expected)
    assert_array_equal(oh.fit(X).transform(X).toarray(), X_expected)

    oh = OneHotEncoder(categories='auto', dtype=output_dtype, sparse=False)
    assert_array_equal(oh.fit_transform(X), X_expected)
    assert_array_equal(oh.fit(X).transform(X), X_expected)


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
def test_one_hot_encoder_dtype_pandas(output_dtype):
    pd = pytest.importorskip('pandas')

    X_df = pd.DataFrame({'A': ['a', 'b'], 'B': [1, 2]})
    X_expected = np.array([[1, 0, 1, 0], [0, 1, 0, 1]], dtype=output_dtype)

    oh = OneHotEncoder(dtype=output_dtype)
    assert_array_equal(oh.fit_transform(X_df).toarray(), X_expected)
    assert_array_equal(oh.fit(X_df).transform(X_df).toarray(), X_expected)

    oh = OneHotEncoder(dtype=output_dtype, sparse=False)
    assert_array_equal(oh.fit_transform(X_df), X_expected)
    assert_array_equal(oh.fit(X_df).transform(X_df), X_expected)


def test_one_hot_encoder_set_params():
    X = np.array([[1, 2]]).T
    oh = OneHotEncoder()
    # set params on not yet fitted object
    oh.set_params(categories=[[0, 1, 2, 3]])
    assert oh.get_params()['categories'] == [[0, 1, 2, 3]]
    assert oh.fit_transform(X).toarray().shape == (2, 4)
    # set params on already fitted object
    oh.set_params(categories=[[0, 1, 2, 3, 4]])
    assert oh.fit_transform(X).toarray().shape == (2, 5)


def check_categorical_onehot(X):
    enc = OneHotEncoder(categories='auto')
    Xtr1 = enc.fit_transform(X)

    enc = OneHotEncoder(categories='auto', sparse=False)
    Xtr2 = enc.fit_transform(X)

    assert_allclose(Xtr1.toarray(), Xtr2)

    assert sparse.isspmatrix_csr(Xtr1)
    return Xtr1.toarray()


@pytest.mark.parametrize("X", [
    [['def', 1, 55], ['abc', 2, 55]],
    np.array([[10, 1, 55], [5, 2, 55]]),
    np.array([['b', 'A', 'cat'], ['a', 'B', 'cat']], dtype=object)
    ], ids=['mixed', 'numeric', 'object'])
def test_one_hot_encoder(X):
    Xtr = check_categorical_onehot(np.array(X)[:, [0]])
    assert_allclose(Xtr, [[0, 1], [1, 0]])

    Xtr = check_categorical_onehot(np.array(X)[:, [0, 1]])
    assert_allclose(Xtr, [[0, 1, 1, 0], [1, 0, 0, 1]])

    Xtr = OneHotEncoder(categories='auto').fit_transform(X)
    assert_allclose(Xtr.toarray(), [[0, 1, 1, 0,  1], [1, 0, 0, 1, 1]])


def test_one_hot_encoder_inverse():
    for sparse_ in [True, False]:
        X = [['abc', 2, 55], ['def', 1, 55], ['abc', 3, 55]]
        enc = OneHotEncoder(sparse=sparse_)
        X_tr = enc.fit_transform(X)
        exp = np.array(X, dtype=object)
        assert_array_equal(enc.inverse_transform(X_tr), exp)

        X = [[2, 55], [1, 55], [3, 55]]
        enc = OneHotEncoder(sparse=sparse_, categories='auto')
        X_tr = enc.fit_transform(X)
        exp = np.array(X)
        assert_array_equal(enc.inverse_transform(X_tr), exp)

        # with unknown categories
        X = [['abc', 2, 55], ['def', 1, 55], ['abc', 3, 55]]
        enc = OneHotEncoder(sparse=sparse_, handle_unknown='ignore',
                            categories=[['abc', 'def'], [1, 2],
                                        [54, 55, 56]])
        X_tr = enc.fit_transform(X)
        exp = np.array(X, dtype=object)
        exp[2, 1] = None
        assert_array_equal(enc.inverse_transform(X_tr), exp)

        # with an otherwise numerical output, still object if unknown
        X = [[2, 55], [1, 55], [3, 55]]
        enc = OneHotEncoder(sparse=sparse_, categories=[[1, 2], [54, 56]],
                            handle_unknown='ignore')
        X_tr = enc.fit_transform(X)
        exp = np.array(X, dtype=object)
        exp[2, 0] = None
        exp[:, 1] = None
        assert_array_equal(enc.inverse_transform(X_tr), exp)

        # incorrect shape raises
        X_tr = np.array([[0, 1, 1], [1, 0, 1]])
        msg = re.escape('Shape of the passed X data is not correct')
        assert_raises_regex(ValueError, msg, enc.inverse_transform, X_tr)


@pytest.mark.parametrize("X, cat_exp, cat_dtype", [
    ([['abc', 55], ['def', 55]], [['abc', 'def'], [55]], np.object_),
    (np.array([[1, 2], [3, 2]]), [[1, 3], [2]], np.integer),
    (np.array([['A', 'cat'], ['B', 'cat']], dtype=object),
     [['A', 'B'], ['cat']], np.object_),
    (np.array([['A', 'cat'], ['B', 'cat']]),
     [['A', 'B'], ['cat']], np.str_)
    ], ids=['mixed', 'numeric', 'object', 'string'])
def test_one_hot_encoder_categories(X, cat_exp, cat_dtype):
    # order of categories should not depend on order of samples
    for Xi in [X, X[::-1]]:
        enc = OneHotEncoder(categories='auto')
        enc.fit(Xi)
        # assert enc.categories == 'auto'
        assert isinstance(enc.categories_, list)
        for res, exp in zip(enc.categories_, cat_exp):
            assert res.tolist() == exp
            assert np.issubdtype(res.dtype, cat_dtype)


@pytest.mark.parametrize("X, X2, cats, cat_dtype", [
    (np.array([['a', 'b']], dtype=object).T,
     np.array([['a', 'd']], dtype=object).T,
     [['a', 'b', 'c']], np.object_),
    (np.array([[1, 2]], dtype='int64').T,
     np.array([[1, 4]], dtype='int64').T,
     [[1, 2, 3]], np.int64),
    (np.array([['a', 'b']], dtype=object).T,
     np.array([['a', 'd']], dtype=object).T,
     [np.array(['a', 'b', 'c'])], np.object_),
    ], ids=['object', 'numeric', 'object-string-cat'])
def test_one_hot_encoder_specified_categories(X, X2, cats, cat_dtype):
    enc = OneHotEncoder(categories=cats)
    exp = np.array([[1., 0., 0.],
                    [0., 1., 0.]])
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert list(enc.categories[0]) == list(cats[0])
    assert enc.categories_[0].tolist() == list(cats[0])
    # manually specified categories should have same dtype as
    # the data when coerced from lists
    assert enc.categories_[0].dtype == cat_dtype

    # when specifying categories manually, unknown categories should already
    # raise when fitting
    enc = OneHotEncoder(categories=cats)
    with pytest.raises(ValueError, match="Found unknown categories"):
        enc.fit(X2)
    enc = OneHotEncoder(categories=cats, handle_unknown='ignore')
    exp = np.array([[1., 0., 0.], [0., 0., 0.]])
    assert_array_equal(enc.fit(X2).transform(X2).toarray(), exp)


def test_one_hot_encoder_unsorted_categories():
    X = np.array([['a', 'b']], dtype=object).T

    enc = OneHotEncoder(categories=[['b', 'a', 'c']])
    exp = np.array([[0., 1., 0.],
                    [1., 0., 0.]])
    assert_array_equal(enc.fit(X).transform(X).toarray(), exp)
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert enc.categories_[0].tolist() == ['b', 'a', 'c']
    assert np.issubdtype(enc.categories_[0].dtype, np.object_)

    # unsorted passed categories still raise for numerical values
    X = np.array([[1, 2]]).T
    enc = OneHotEncoder(categories=[[2, 1, 3]])
    msg = 'Unsorted categories are not supported'
    with pytest.raises(ValueError, match=msg):
        enc.fit_transform(X)


def test_one_hot_encoder_specified_categories_mixed_columns():
    # multiple columns
    X = np.array([['a', 'b'], [0, 2]], dtype=object).T
    enc = OneHotEncoder(categories=[['a', 'b', 'c'], [0, 1, 2]])
    exp = np.array([[1., 0., 0., 1., 0., 0.],
                    [0., 1., 0., 0., 0., 1.]])
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert enc.categories_[0].tolist() == ['a', 'b', 'c']
    assert np.issubdtype(enc.categories_[0].dtype, np.object_)
    assert enc.categories_[1].tolist() == [0, 1, 2]
    # integer categories but from object dtype data
    assert np.issubdtype(enc.categories_[1].dtype, np.object_)


def test_one_hot_encoder_pandas():
    pd = pytest.importorskip('pandas')

    X_df = pd.DataFrame({'A': ['a', 'b'], 'B': [1, 2]})

    Xtr = check_categorical_onehot(X_df)
    assert_allclose(Xtr, [[1, 0, 1, 0], [0, 1, 0, 1]])


def test_one_hot_encoder_feature_names():
    enc = OneHotEncoder()
    X = [['Male', 1, 'girl', 2, 3],
         ['Female', 41, 'girl', 1, 10],
         ['Male', 51, 'boy', 12, 3],
         ['Male', 91, 'girl', 21, 30]]

    enc.fit(X)
    feature_names = enc.get_feature_names()
    assert isinstance(feature_names, np.ndarray)

    assert_array_equal(['x0_Female', 'x0_Male',
                        'x1_1', 'x1_41', 'x1_51', 'x1_91',
                        'x2_boy', 'x2_girl',
                        'x3_1', 'x3_2', 'x3_12', 'x3_21',
                        'x4_3',
                        'x4_10', 'x4_30'], feature_names)

    feature_names2 = enc.get_feature_names(['one', 'two',
                                            'three', 'four', 'five'])

    assert_array_equal(['one_Female', 'one_Male',
                        'two_1', 'two_41', 'two_51', 'two_91',
                        'three_boy', 'three_girl',
                        'four_1', 'four_2', 'four_12', 'four_21',
                        'five_3', 'five_10', 'five_30'], feature_names2)

    with pytest.raises(ValueError, match="input_features should have length"):
        enc.get_feature_names(['one', 'two'])


def test_one_hot_encoder_feature_names_unicode():
    enc = OneHotEncoder()
    X = np.array([[u'c❤t1', u'dat2']], dtype=object).T
    enc.fit(X)
    feature_names = enc.get_feature_names()
    assert_array_equal([u'x0_c❤t1', u'x0_dat2'], feature_names)
    feature_names = enc.get_feature_names(input_features=[u'n👍me'])
    assert_array_equal([u'n👍me_c❤t1', u'n👍me_dat2'], feature_names)


@pytest.mark.parametrize("X", [np.array([[1, np.nan]]).T,
                               np.array([['a', np.nan]], dtype=object).T],
                         ids=['numeric', 'object'])
@pytest.mark.parametrize("handle_unknown", ['error', 'ignore'])
def test_one_hot_encoder_raise_missing(X, handle_unknown):
    ohe = OneHotEncoder(categories='auto', handle_unknown=handle_unknown)

    with pytest.raises(ValueError, match="Input contains NaN"):
        ohe.fit(X)

    with pytest.raises(ValueError, match="Input contains NaN"):
        ohe.fit_transform(X)

    ohe.fit(X[:1, :])

    with pytest.raises(ValueError, match="Input contains NaN"):
        ohe.transform(X)


@pytest.mark.parametrize("X", [
    [['abc', 2, 55], ['def', 1, 55]],
    np.array([[10, 2, 55], [20, 1, 55]]),
    np.array([['a', 'B', 'cat'], ['b', 'A', 'cat']], dtype=object)
    ], ids=['mixed', 'numeric', 'object'])
def test_ordinal_encoder(X):
    enc = OrdinalEncoder()
    exp = np.array([[0, 1, 0],
                    [1, 0, 0]], dtype='int64')
    assert_array_equal(enc.fit_transform(X), exp.astype('float64'))
    enc = OrdinalEncoder(dtype='int64')
    assert_array_equal(enc.fit_transform(X), exp)


@pytest.mark.parametrize("X, X2, cats, cat_dtype", [
    (np.array([['a', 'b']], dtype=object).T,
     np.array([['a', 'd']], dtype=object).T,
     [['a', 'b', 'c']], np.object_),
    (np.array([[1, 2]], dtype='int64').T,
     np.array([[1, 4]], dtype='int64').T,
     [[1, 2, 3]], np.int64),
    (np.array([['a', 'b']], dtype=object).T,
     np.array([['a', 'd']], dtype=object).T,
     [np.array(['a', 'b', 'c'])], np.object_),
    ], ids=['object', 'numeric', 'object-string-cat'])
def test_ordinal_encoder_specified_categories(X, X2, cats, cat_dtype):
    enc = OrdinalEncoder(categories=cats)
    exp = np.array([[0.], [1.]])
    assert_array_equal(enc.fit_transform(X), exp)
    assert list(enc.categories[0]) == list(cats[0])
    assert enc.categories_[0].tolist() == list(cats[0])
    # manually specified categories should have same dtype as
    # the data when coerced from lists
    assert enc.categories_[0].dtype == cat_dtype

    # when specifying categories manually, unknown categories should already
    # raise when fitting
    enc = OrdinalEncoder(categories=cats)
    with pytest.raises(ValueError, match="Found unknown categories"):
        enc.fit(X2)


def test_ordinal_encoder_inverse():
    X = [['abc', 2, 55], ['def', 1, 55]]
    enc = OrdinalEncoder()
    X_tr = enc.fit_transform(X)
    exp = np.array(X, dtype=object)
    assert_array_equal(enc.inverse_transform(X_tr), exp)

    # incorrect shape raises
    X_tr = np.array([[0, 1, 1, 2], [1, 0, 1, 0]])
    msg = re.escape('Shape of the passed X data is not correct')
    assert_raises_regex(ValueError, msg, enc.inverse_transform, X_tr)


@pytest.mark.parametrize("X", [np.array([[1, np.nan]]).T,
                               np.array([['a', np.nan]], dtype=object).T],
                         ids=['numeric', 'object'])
def test_ordinal_encoder_raise_missing(X):
    ohe = OrdinalEncoder()

    with pytest.raises(ValueError, match="Input contains NaN"):
        ohe.fit(X)

    with pytest.raises(ValueError, match="Input contains NaN"):
        ohe.fit_transform(X)

    ohe.fit(X[:1, :])

    with pytest.raises(ValueError, match="Input contains NaN"):
        ohe.transform(X)


def test_encoder_dtypes():
    # check that dtypes are preserved when determining categories
    enc = OneHotEncoder(categories='auto')
    exp = np.array([[1., 0., 1., 0.], [0., 1., 0., 1.]], dtype='float64')

    for X in [np.array([[1, 2], [3, 4]], dtype='int64'),
              np.array([[1, 2], [3, 4]], dtype='float64'),
              np.array([['a', 'b'], ['c', 'd']]),  # string dtype
              np.array([[1, 'a'], [3, 'b']], dtype='object')]:
        enc.fit(X)
        assert all([enc.categories_[i].dtype == X.dtype for i in range(2)])
        assert_array_equal(enc.transform(X).toarray(), exp)

    X = [[1, 2], [3, 4]]
    enc.fit(X)
    assert all([np.issubdtype(enc.categories_[i].dtype, np.integer)
                for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)

    X = [[1, 'a'], [3, 'b']]
    enc.fit(X)
    assert all([enc.categories_[i].dtype == 'object' for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)


def test_encoder_dtypes_pandas():
    # check dtype (similar to test_categorical_encoder_dtypes for dataframes)
    pd = pytest.importorskip('pandas')

    enc = OneHotEncoder(categories='auto')
    exp = np.array([[1., 0., 1., 0.], [0., 1., 0., 1.]], dtype='float64')

    X = pd.DataFrame({'A': [1, 2], 'B': [3, 4]}, dtype='int64')
    enc.fit(X)
    assert all([enc.categories_[i].dtype == 'int64' for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)

    X = pd.DataFrame({'A': [1, 2], 'B': ['a', 'b']})
    enc.fit(X)
    assert all([enc.categories_[i].dtype == 'object' for i in range(2)])
    assert_array_equal(enc.transform(X).toarray(), exp)


def test_one_hot_encoder_warning():
    enc = OneHotEncoder()
    X = [['Male', 1], ['Female', 3]]
    np.testing.assert_no_warnings(enc.fit_transform, X)


def test_categorical_encoder_stub():
    from sklearn.preprocessing import CategoricalEncoder
    assert_raises(RuntimeError, CategoricalEncoder, encoding='ordinal')
