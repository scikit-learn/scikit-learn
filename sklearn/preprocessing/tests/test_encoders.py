from __future__ import division

import re

import numpy as np
from scipy import sparse
import pytest

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raises_regex
from sklearn.utils.testing import assert_allclose
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_no_warnings

from sklearn.preprocessing._encoders import _transform_selected
from sklearn.preprocessing.data import Binarizer
from sklearn.preprocessing._encoders import OneHotEncoder
from sklearn.preprocessing._encoders import OrdinalEncoder


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


def _check_transform_selected(X, X_expected, dtype, sel):
    for M in (X, sparse.csr_matrix(X)):
        Xtr = _transform_selected(M, Binarizer().transform, dtype, sel)
        assert_array_equal(toarray(Xtr), X_expected)


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
@pytest.mark.parametrize("input_dtype", [np.int32, np.float32, np.float64])
def test_transform_selected(output_dtype, input_dtype):
    X = np.asarray([[3, 2, 1], [0, 1, 1]], dtype=input_dtype)

    X_expected = np.asarray([[1, 2, 1], [0, 1, 1]], dtype=output_dtype)
    _check_transform_selected(X, X_expected, output_dtype, [0])
    _check_transform_selected(X, X_expected, output_dtype,
                              [True, False, False])

    X_expected = np.asarray([[1, 1, 1], [0, 1, 1]], dtype=output_dtype)
    _check_transform_selected(X, X_expected, output_dtype, [0, 1, 2])
    _check_transform_selected(X, X_expected, output_dtype, [True, True, True])
    _check_transform_selected(X, X_expected, output_dtype, "all")

    _check_transform_selected(X, X, output_dtype, [])
    _check_transform_selected(X, X, output_dtype, [False, False, False])


@pytest.mark.parametrize("output_dtype", [np.int32, np.float32, np.float64])
@pytest.mark.parametrize("input_dtype", [np.int32, np.float32, np.float64])
def test_transform_selected_copy_arg(output_dtype, input_dtype):
    # transformer that alters X
    def _mutating_transformer(X):
        X[0, 0] = X[0, 0] + 1
        return X

    original_X = np.asarray([[1, 2], [3, 4]], dtype=input_dtype)
    expected_Xtr = np.asarray([[2, 2], [3, 4]], dtype=output_dtype)

    X = original_X.copy()
    Xtr = _transform_selected(X, _mutating_transformer, output_dtype,
                              copy=True, selected='all')

    assert_array_equal(toarray(X), toarray(original_X))
    assert_array_equal(toarray(Xtr), expected_Xtr)


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
    enc = OneHotEncoder()
    Xtr1 = enc.fit_transform(X)

    enc = OneHotEncoder(sparse=False)
    Xtr2 = enc.fit_transform(X)

    assert_allclose(Xtr1.toarray(), Xtr2)

    assert sparse.isspmatrix_csr(Xtr1)
    return Xtr1.toarray()


def test_one_hot_encoder():
    X = [['abc', 1, 55], ['def', 2, 55]]

    Xtr = check_categorical_onehot(np.array(X)[:, [0]])
    assert_allclose(Xtr, [[1, 0], [0, 1]])

    Xtr = check_categorical_onehot(np.array(X)[:, [0, 1]])
    assert_allclose(Xtr, [[1, 0, 1, 0], [0, 1, 0, 1]])

    Xtr = OneHotEncoder().fit_transform(X)
    assert_allclose(Xtr.toarray(), [[1, 0, 1, 0,  1], [0, 1, 0, 1, 1]])


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


def test_one_hot_encoder_categories():
    X = [['abc', 1, 55], ['def', 2, 55]]

    # order of categories should not depend on order of samples
    for Xi in [X, X[::-1]]:
        enc = OneHotEncoder()
        enc.fit(Xi)
        # assert enc.categories == 'auto'
        assert isinstance(enc.categories_, list)
        cat_exp = [['abc', 'def'], [1, 2], [55]]
        for res, exp in zip(enc.categories_, cat_exp):
            assert res.tolist() == exp


def test_one_hot_encoder_specified_categories():
    X = np.array([['a', 'b']], dtype=object).T

    enc = OneHotEncoder(categories=[['a', 'b', 'c']])
    exp = np.array([[1., 0., 0.],
                    [0., 1., 0.]])
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert enc.categories[0] == ['a', 'b', 'c']
    assert enc.categories_[0].tolist() == ['a', 'b', 'c']
    assert np.issubdtype(enc.categories_[0].dtype, np.str_)

    # unsorted passed categories raises for now
    enc = OneHotEncoder(categories=[['c', 'b', 'a']])
    msg = re.escape('Unsorted categories are not yet supported')
    assert_raises_regex(ValueError, msg, enc.fit_transform, X)

    # multiple columns
    X = np.array([['a', 'b'], [0, 2]], dtype=object).T
    enc = OneHotEncoder(categories=[['a', 'b', 'c'], [0, 1, 2]])
    exp = np.array([[1., 0., 0., 1., 0., 0.],
                    [0., 1., 0., 0., 0., 1.]])
    assert_array_equal(enc.fit_transform(X).toarray(), exp)
    assert enc.categories_[0].tolist() == ['a', 'b', 'c']
    assert np.issubdtype(enc.categories_[0].dtype, np.str_)
    assert enc.categories_[1].tolist() == [0, 1, 2]
    assert np.issubdtype(enc.categories_[1].dtype, np.integer)

    # when specifying categories manually, unknown categories should already
    # raise when fitting
    X = np.array([['a', 'b', 'c']]).T
    enc = OneHotEncoder(categories=[['a', 'b']])
    assert_raises(ValueError, enc.fit, X)
    enc = OneHotEncoder(categories=[['a', 'b']], handle_unknown='ignore')
    exp = np.array([[1., 0.], [0., 1.], [0., 0.]])
    assert_array_equal(enc.fit(X).transform(X).toarray(), exp)


def test_one_hot_encoder_pandas():
    pd = pytest.importorskip('pandas')

    X_df = pd.DataFrame({'A': ['a', 'b'], 'B': [1, 2]})

    Xtr = check_categorical_onehot(X_df)
    assert_allclose(Xtr, [[1, 0, 1, 0], [0, 1, 0, 1]])


def test_ordinal_encoder():
    X = [['abc', 2, 55], ['def', 1, 55]]

    enc = OrdinalEncoder()
    exp = np.array([[0, 1, 0],
                    [1, 0, 0]], dtype='int64')
    assert_array_equal(enc.fit_transform(X), exp.astype('float64'))
    enc = OrdinalEncoder(dtype='int64')
    assert_array_equal(enc.fit_transform(X), exp)


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
