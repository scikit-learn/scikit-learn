"""
Test the pipeline module.
"""

import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import SkipTest

from sklearn.base import BaseEstimator
from sklearn.experimental import ColumnTransformer

from sklearn.preprocessing import StandardScaler, Normalizer


class Trans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return X


class SparseMatrixTrans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        n_samples = len(X)
        return sparse.eye(n_samples, n_samples).tocsr()


def test_column_transformer():
    # array
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T

    X_res_first1D = np.array([0, 1, 2])
    X_res_second1D = np.array([2, 4, 6])
    X_res_first2D = X_res_first1D.reshape(-1, 1)
    X_res_both = X_array

    # scalar
    ct = ColumnTransformer([('trans', Trans(), 0)])
    assert_array_equal(ct.fit_transform(X_array), X_res_first1D)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_first1D)

    ct = ColumnTransformer([('trans', Trans(), [0])])
    assert_array_equal(ct.fit_transform(X_array), X_res_first2D)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_first2D)

    # list
    ct = ColumnTransformer([('trans', Trans(), [0, 1])])
    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)

    ct = ColumnTransformer([('trans1', Trans(), [0]),
                            ('trans2', Trans(), [1])])
    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)

    # slice
    ct = ColumnTransformer([('trans', Trans(), slice(0, 1))])
    assert_array_equal(ct.fit_transform(X_array), X_res_first2D)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_first2D)

    ct = ColumnTransformer([('trans', Trans(), slice(0, 2))])
    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)

    # boolean mask
    ct = ColumnTransformer([('trans', Trans(), np.array([True, False]))])
    assert_array_equal(ct.fit_transform(X_array), X_res_first2D)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_first2D)

    # test with transformer_weights
    transformer_weights = {'trans1': .1, 'trans2': 10}
    both = ColumnTransformer([('trans1', Trans(), [0]),
                              ('trans2', Trans(), [1])],
                             transformer_weights=transformer_weights)
    res = np.vstack([transformer_weights['trans1'] * X_res_first1D,
                     transformer_weights['trans2'] * X_res_second1D]).T
    assert_array_equal(both.fit_transform(X_array), res)
    assert_array_equal(both.fit(X_array).transform(X_array), res)

    both = ColumnTransformer([('trans', Trans(), [0, 1])],
                             transformer_weights={'trans': .1})
    assert_array_equal(both.fit_transform(X_array), 0.1 * X_res_both)
    assert_array_equal(both.fit(X_array).transform(X_array), 0.1 * X_res_both)


def test_column_transformer_dataframe():
    # dataframe
    try:
        import pandas as pd
        X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
        X_df = pd.DataFrame(X_array, columns=['first', 'second'])
    except ImportError:
        raise SkipTest("pandas is not installed: skipping ColumnTransformer "
                       "tests for DataFrames.")

    X_res_first1D = np.array([0, 1, 2])
    X_res_first2D = X_res_first1D.reshape(-1, 1)
    X_res_both = X_array

    # String keys: label based

    # scalar
    ct = ColumnTransformer([('trans', Trans(), 'first')])
    assert_array_equal(ct.fit_transform(X_df), X_res_first1D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first1D)

    # list
    ct = ColumnTransformer([('trans', Trans(), ['first'])])
    assert_array_equal(ct.fit_transform(X_df), X_res_first2D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first2D)

    ct = ColumnTransformer([('trans', Trans(), ['first', 'second'])])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)

    ct = ColumnTransformer([('trans1', Trans(), ['first']),
                            ('trans2', Trans(), ['second'])])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)

    # slice
    ct = ColumnTransformer([('trans', Trans(), slice('first', 'second'))])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)

    # int keys: positional

    # scalar
    ct = ColumnTransformer([('trans', Trans(), 0)])
    assert_array_equal(ct.fit_transform(X_df), X_res_first1D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first1D)

    # list
    ct = ColumnTransformer([('trans', Trans(), [0])])
    assert_array_equal(ct.fit_transform(X_df), X_res_first2D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first2D)

    ct = ColumnTransformer([('trans', Trans(), [0, 1])])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)

    ct = ColumnTransformer([('trans1', Trans(), [0]),
                            ('trans2', Trans(), [1])])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)

    # slice
    ct = ColumnTransformer([('trans', Trans(), slice(0, 1))])
    assert_array_equal(ct.fit_transform(X_df), X_res_first2D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first2D)

    ct = ColumnTransformer([('trans', Trans(), slice(0, 2))])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)

    # boolean mask
    ct = ColumnTransformer([('trans', Trans(), np.array([True, False]))])
    assert_array_equal(ct.fit_transform(X_df), X_res_first2D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first2D)

    s_mask = pd.Series([True, False], index=['first', 'second'])
    ct = ColumnTransformer([('trans', Trans(), s_mask)])
    assert_array_equal(ct.fit_transform(X_df), X_res_first2D)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first2D)

    # test with transformer_weights
    transformer_weights = {'trans1': .1, 'trans2': 10}
    both = ColumnTransformer([('trans1', Trans(), ['first']),
                              ('trans2', Trans(), ['second'])],
                             transformer_weights=transformer_weights)
    res = np.vstack([transformer_weights['trans1'] * X_df['first'],
                     transformer_weights['trans2'] * X_df['second']]).T
    assert_array_equal(both.fit_transform(X_df), res)
    assert_array_equal(both.fit(X_df).transform(X_df), res)

    # test multiple columns
    both = ColumnTransformer([('trans', Trans(), ['first', 'second'])],
                             transformer_weights={'trans': .1})
    assert_array_equal(both.fit_transform(X_df), 0.1 * X_res_both)
    assert_array_equal(both.fit(X_df).transform(X_df), 0.1 * X_res_both)

    both = ColumnTransformer([('trans', Trans(), [0, 1])],
                             transformer_weights={'trans': .1})
    assert_array_equal(both.fit_transform(X_df), 0.1 * X_res_both)
    assert_array_equal(both.fit(X_df).transform(X_df), 0.1 * X_res_both)

    # ensure pandas object is passes through

    class TransAssert(BaseEstimator):

        def fit(self, X, y=None):
            return self

        def transform(self, X, y=None):
            assert_true(isinstance(X, (pd.DataFrame, pd.Series)))

    ct = ColumnTransformer([('trans', TransAssert(), 'first')])
    ct.fit_transform(X_df)
    ct = ColumnTransformer([('trans', TransAssert(), ['first', 'second'])])
    ct.fit_transform(X_df)


def test_column_transformer_sparse_array():
    X_sparse = sparse.eye(3, 2).tocsr()

    # no distinction between 1D and 2D
    X_res_first = X_sparse[:, 0].A
    X_res_both = X_sparse.A

    col_trans = ColumnTransformer([('trans', Trans(), 0)])
    assert_true(sparse.issparse(col_trans.fit_transform(X_sparse)))
    assert_array_equal(col_trans.fit_transform(X_sparse).A, X_res_first)
    assert_array_equal(col_trans.fit(X_sparse).transform(X_sparse).A,
                       X_res_first)

    col_trans = ColumnTransformer([('trans', Trans(), [0])])
    assert_true(sparse.issparse(col_trans.fit_transform(X_sparse)))
    assert_array_equal(col_trans.fit_transform(X_sparse).A, X_res_first)
    assert_array_equal(col_trans.fit(X_sparse).transform(X_sparse).A,
                       X_res_first)

    col_trans = ColumnTransformer([('trans', Trans(), [0, 1])])
    assert_true(sparse.issparse(col_trans.fit_transform(X_sparse)))
    assert_array_equal(col_trans.fit_transform(X_sparse).A, X_res_both)
    assert_array_equal(col_trans.fit(X_sparse).transform(X_sparse).A,
                       X_res_both)

    col_trans = ColumnTransformer([('trans', Trans(), slice(0, 1))])
    assert_true(sparse.issparse(col_trans.fit_transform(X_sparse)))
    assert_array_equal(col_trans.fit_transform(X_sparse).A, X_res_first)
    assert_array_equal(col_trans.fit(X_sparse).transform(X_sparse).A,
                       X_res_first)

    col_trans = ColumnTransformer([('trans', Trans(), slice(0, 2))])
    assert_true(sparse.issparse(col_trans.fit_transform(X_sparse)))
    assert_array_equal(col_trans.fit_transform(X_sparse).A, X_res_both)
    assert_array_equal(col_trans.fit(X_sparse).transform(X_sparse).A,
                       X_res_both)


def test_column_transformer_2D_array_items():
    union = ColumnTransformer(
        [("norm1", Normalizer(norm='l1'), [0, 1]),
         ("norm2", Normalizer(norm='l1'), [2, 3])])
    X = np.array([[0., 1., 2., 2.],
                  [1., 1., 0., 1.]])
    X_res = np.array([[0., 1., 0.5, 0.5],
                      [0.5, 0.5, 0., 1.]])
    assert_array_equal(union.fit_transform(X), X_res)


def test_column_transformer_sparse_stacking():
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    col_trans = ColumnTransformer([('trans1', Trans(), [0]),
                                   ('trans2', SparseMatrixTrans(), 1)])
    col_trans.fit(X_array)
    X_trans = col_trans.transform(X_array)
    assert_true(sparse.issparse(X_trans))
    assert_equal(X_trans.shape, (X_trans.shape[0], X_trans.shape[0] + 1))
    assert_array_equal(X_trans.toarray()[:, 1:], np.eye(X_trans.shape[0]))


def test_column_transformer_error_msg_1D():
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T

    col_trans = ColumnTransformer([('trans', StandardScaler(), 0)])
    assert_raise_message(ValueError, "1D data passed to a transformer",
                         col_trans.fit, X_array)
    assert_raise_message(ValueError, "1D data passed to a transformer",
                         col_trans.fit_transform, X_array)
