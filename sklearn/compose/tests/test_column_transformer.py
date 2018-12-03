"""
Test the ColumnTransformer.
"""

import numpy as np
from scipy import sparse
import pytest

from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_dict_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_allclose_dense_sparse
from sklearn.utils.testing import assert_almost_equal

from sklearn.base import BaseEstimator
from sklearn.externals import six
from sklearn.compose import ColumnTransformer, make_column_transformer
from sklearn.exceptions import NotFittedError, DataConversionWarning
from sklearn.preprocessing import StandardScaler, Normalizer, OneHotEncoder
from sklearn.feature_extraction import DictVectorizer


class Trans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        # 1D Series -> 2D DataFrame
        if hasattr(X, 'to_frame'):
            return X.to_frame()
        # 1D array -> 2D array
        if X.ndim == 1:
            return np.atleast_2d(X).T
        return X


class DoubleTrans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X):
        return 2*X


class SparseMatrixTrans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        n_samples = len(X)
        return sparse.eye(n_samples, n_samples).tocsr()


class TransNo2D(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return X


class TransRaise(BaseEstimator):

    def fit(self, X, y=None):
        raise ValueError("specific message")

    def transform(self, X, y=None):
        raise ValueError("specific message")


def test_column_transformer():
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T

    X_res_first1D = np.array([0, 1, 2])
    X_res_second1D = np.array([2, 4, 6])
    X_res_first = X_res_first1D.reshape(-1, 1)
    X_res_both = X_array

    cases = [
        # single column 1D / 2D
        (0, X_res_first),
        ([0], X_res_first),
        # list-like
        ([0, 1], X_res_both),
        (np.array([0, 1]), X_res_both),
        # slice
        (slice(0, 1), X_res_first),
        (slice(0, 2), X_res_both),
        # boolean mask
        (np.array([True, False]), X_res_first),
    ]

    for selection, res in cases:
        ct = ColumnTransformer([('trans', Trans(), selection)],
                               remainder='drop')
        assert_array_equal(ct.fit_transform(X_array), res)
        assert_array_equal(ct.fit(X_array).transform(X_array), res)

        # callable that returns any of the allowed specifiers
        ct = ColumnTransformer([('trans', Trans(), lambda x: selection)],
                               remainder='drop')
        assert_array_equal(ct.fit_transform(X_array), res)
        assert_array_equal(ct.fit(X_array).transform(X_array), res)

    ct = ColumnTransformer([('trans1', Trans(), [0]),
                            ('trans2', Trans(), [1])])
    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)
    assert len(ct.transformers_) == 2

    # test with transformer_weights
    transformer_weights = {'trans1': .1, 'trans2': 10}
    both = ColumnTransformer([('trans1', Trans(), [0]),
                              ('trans2', Trans(), [1])],
                             transformer_weights=transformer_weights)
    res = np.vstack([transformer_weights['trans1'] * X_res_first1D,
                     transformer_weights['trans2'] * X_res_second1D]).T
    assert_array_equal(both.fit_transform(X_array), res)
    assert_array_equal(both.fit(X_array).transform(X_array), res)
    assert len(both.transformers_) == 2

    both = ColumnTransformer([('trans', Trans(), [0, 1])],
                             transformer_weights={'trans': .1})
    assert_array_equal(both.fit_transform(X_array), 0.1 * X_res_both)
    assert_array_equal(both.fit(X_array).transform(X_array), 0.1 * X_res_both)
    assert len(both.transformers_) == 1


def test_column_transformer_dataframe():
    pd = pytest.importorskip('pandas')

    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_df = pd.DataFrame(X_array, columns=['first', 'second'])

    X_res_first = np.array([0, 1, 2]).reshape(-1, 1)
    X_res_both = X_array

    cases = [
        # String keys: label based

        # scalar
        ('first', X_res_first),
        # list
        (['first'], X_res_first),
        (['first', 'second'], X_res_both),
        # slice
        (slice('first', 'second'), X_res_both),

        # int keys: positional

        # scalar
        (0, X_res_first),
        # list
        ([0], X_res_first),
        ([0, 1], X_res_both),
        (np.array([0, 1]), X_res_both),
        # slice
        (slice(0, 1), X_res_first),
        (slice(0, 2), X_res_both),

        # boolean mask
        (np.array([True, False]), X_res_first),
        (pd.Series([True, False], index=['first', 'second']), X_res_first),
    ]

    for selection, res in cases:
        ct = ColumnTransformer([('trans', Trans(), selection)],
                               remainder='drop')
        assert_array_equal(ct.fit_transform(X_df), res)
        assert_array_equal(ct.fit(X_df).transform(X_df), res)

        # callable that returns any of the allowed specifiers
        ct = ColumnTransformer([('trans', Trans(), lambda X: selection)],
                               remainder='drop')
        assert_array_equal(ct.fit_transform(X_df), res)
        assert_array_equal(ct.fit(X_df).transform(X_df), res)

    ct = ColumnTransformer([('trans1', Trans(), ['first']),
                            ('trans2', Trans(), ['second'])])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] != 'remainder'

    ct = ColumnTransformer([('trans1', Trans(), [0]),
                            ('trans2', Trans(), [1])])
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] != 'remainder'

    # test with transformer_weights
    transformer_weights = {'trans1': .1, 'trans2': 10}
    both = ColumnTransformer([('trans1', Trans(), ['first']),
                              ('trans2', Trans(), ['second'])],
                             transformer_weights=transformer_weights)
    res = np.vstack([transformer_weights['trans1'] * X_df['first'],
                     transformer_weights['trans2'] * X_df['second']]).T
    assert_array_equal(both.fit_transform(X_df), res)
    assert_array_equal(both.fit(X_df).transform(X_df), res)
    assert len(both.transformers_) == 2
    assert ct.transformers_[-1][0] != 'remainder'

    # test multiple columns
    both = ColumnTransformer([('trans', Trans(), ['first', 'second'])],
                             transformer_weights={'trans': .1})
    assert_array_equal(both.fit_transform(X_df), 0.1 * X_res_both)
    assert_array_equal(both.fit(X_df).transform(X_df), 0.1 * X_res_both)
    assert len(both.transformers_) == 1
    assert ct.transformers_[-1][0] != 'remainder'

    both = ColumnTransformer([('trans', Trans(), [0, 1])],
                             transformer_weights={'trans': .1})
    assert_array_equal(both.fit_transform(X_df), 0.1 * X_res_both)
    assert_array_equal(both.fit(X_df).transform(X_df), 0.1 * X_res_both)
    assert len(both.transformers_) == 1
    assert ct.transformers_[-1][0] != 'remainder'

    # ensure pandas object is passes through

    class TransAssert(BaseEstimator):

        def fit(self, X, y=None):
            return self

        def transform(self, X, y=None):
            assert isinstance(X, (pd.DataFrame, pd.Series))
            if isinstance(X, pd.Series):
                X = X.to_frame()
            return X

    ct = ColumnTransformer([('trans', TransAssert(), 'first')],
                           remainder='drop')
    ct.fit_transform(X_df)
    ct = ColumnTransformer([('trans', TransAssert(), ['first', 'second'])])
    ct.fit_transform(X_df)

    # integer column spec + integer column names -> still use positional
    X_df2 = X_df.copy()
    X_df2.columns = [1, 0]
    ct = ColumnTransformer([('trans', Trans(), 0)], remainder='drop')
    assert_array_equal(ct.fit_transform(X_df), X_res_first)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first)

    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'drop'
    assert_array_equal(ct.transformers_[-1][2], [1])


@pytest.mark.parametrize("pandas", [True, False], ids=['pandas', 'numpy'])
@pytest.mark.parametrize("column", [[], np.array([False, False])],
                         ids=['list', 'bool'])
def test_column_transformer_empty_columns(pandas, column):
    # test case that ensures that the column transformer does also work when
    # a given transformer doesn't have any columns to work on
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_res_both = X_array

    if pandas:
        pd = pytest.importorskip('pandas')
        X = pd.DataFrame(X_array, columns=['first', 'second'])
    else:
        X = X_array

    ct = ColumnTransformer([('trans1', Trans(), [0, 1]),
                            ('trans2', Trans(), column)])
    assert_array_equal(ct.fit_transform(X), X_res_both)
    assert_array_equal(ct.fit(X).transform(X), X_res_both)
    assert len(ct.transformers_) == 2
    assert isinstance(ct.transformers_[1][1], Trans)

    ct = ColumnTransformer([('trans1', Trans(), column),
                            ('trans2', Trans(), [0, 1])])
    assert_array_equal(ct.fit_transform(X), X_res_both)
    assert_array_equal(ct.fit(X).transform(X), X_res_both)
    assert len(ct.transformers_) == 2
    assert isinstance(ct.transformers_[0][1], Trans)

    ct = ColumnTransformer([('trans', Trans(), column)],
                           remainder='passthrough')
    assert_array_equal(ct.fit_transform(X), X_res_both)
    assert_array_equal(ct.fit(X).transform(X), X_res_both)
    assert len(ct.transformers_) == 2  # including remainder
    assert isinstance(ct.transformers_[0][1], Trans)

    fixture = np.array([[], [], []])
    ct = ColumnTransformer([('trans', Trans(), column)],
                           remainder='drop')
    assert_array_equal(ct.fit_transform(X), fixture)
    assert_array_equal(ct.fit(X).transform(X), fixture)
    assert len(ct.transformers_) == 2  # including remainder
    assert isinstance(ct.transformers_[0][1], Trans)


def test_column_transformer_sparse_array():
    X_sparse = sparse.eye(3, 2).tocsr()

    # no distinction between 1D and 2D
    X_res_first = X_sparse[:, 0]
    X_res_both = X_sparse

    for col in [0, [0], slice(0, 1)]:
        for remainder, res in [('drop', X_res_first),
                               ('passthrough', X_res_both)]:
            ct = ColumnTransformer([('trans', Trans(), col)],
                                   remainder=remainder,
                                   sparse_threshold=0.8)
            assert sparse.issparse(ct.fit_transform(X_sparse))
            assert_allclose_dense_sparse(ct.fit_transform(X_sparse), res)
            assert_allclose_dense_sparse(ct.fit(X_sparse).transform(X_sparse),
                                         res)

    for col in [[0, 1], slice(0, 2)]:
        ct = ColumnTransformer([('trans', Trans(), col)],
                               sparse_threshold=0.8)
        assert sparse.issparse(ct.fit_transform(X_sparse))
        assert_allclose_dense_sparse(ct.fit_transform(X_sparse), X_res_both)
        assert_allclose_dense_sparse(ct.fit(X_sparse).transform(X_sparse),
                                     X_res_both)


def test_column_transformer_list():
    X_list = [
        [1, float('nan'), 'a'],
        [0, 0, 'b']
    ]
    expected_result = np.array([
        [1, float('nan'), 1, 0],
        [-1, 0, 0, 1],
    ])

    ct = ColumnTransformer([
        ('numerical', StandardScaler(), [0, 1]),
        ('categorical', OneHotEncoder(), [2]),
    ])

    with pytest.warns(DataConversionWarning):
        # TODO: this warning is not very useful in this case, would be good
        # to get rid of it
        assert_array_equal(ct.fit_transform(X_list), expected_result)
        assert_array_equal(ct.fit(X_list).transform(X_list), expected_result)


def test_column_transformer_sparse_stacking():
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    col_trans = ColumnTransformer([('trans1', Trans(), [0]),
                                   ('trans2', SparseMatrixTrans(), 1)],
                                  sparse_threshold=0.8)
    col_trans.fit(X_array)
    X_trans = col_trans.transform(X_array)
    assert sparse.issparse(X_trans)
    assert_equal(X_trans.shape, (X_trans.shape[0], X_trans.shape[0] + 1))
    assert_array_equal(X_trans.toarray()[:, 1:], np.eye(X_trans.shape[0]))
    assert len(col_trans.transformers_) == 2
    assert col_trans.transformers_[-1][0] != 'remainder'

    col_trans = ColumnTransformer([('trans1', Trans(), [0]),
                                   ('trans2', SparseMatrixTrans(), 1)],
                                  sparse_threshold=0.1)
    col_trans.fit(X_array)
    X_trans = col_trans.transform(X_array)
    assert not sparse.issparse(X_trans)
    assert X_trans.shape == (X_trans.shape[0], X_trans.shape[0] + 1)
    assert_array_equal(X_trans[:, 1:], np.eye(X_trans.shape[0]))


def test_column_transformer_mixed_cols_sparse():
    df = np.array([['a', 1, True],
                   ['b', 2, False]],
                  dtype='O')

    ct = make_column_transformer(
        (OneHotEncoder(), [0]),
        ('passthrough', [1, 2]),
        sparse_threshold=1.0
    )

    # this shouldn't fail, since boolean can be coerced into a numeric
    # See: https://github.com/scikit-learn/scikit-learn/issues/11912
    X_trans = ct.fit_transform(df)
    assert X_trans.getformat() == 'csr'
    assert_array_equal(X_trans.toarray(), np.array([[1, 0, 1, 1],
                                                    [0, 1, 2, 0]]))

    ct = make_column_transformer(
        (OneHotEncoder(), [0]),
        ('passthrough', [0]),
        sparse_threshold=1.0
    )
    with pytest.raises(ValueError,
                       match="For a sparse output, all columns should"):
        # this fails since strings `a` and `b` cannot be
        # coerced into a numeric.
        ct.fit_transform(df)


def test_column_transformer_sparse_threshold():
    X_array = np.array([['a', 'b'], ['A', 'B']], dtype=object).T
    # above data has sparsity of 4 / 8 = 0.5

    # apply threshold even if all sparse
    col_trans = ColumnTransformer([('trans1', OneHotEncoder(), [0]),
                                   ('trans2', OneHotEncoder(), [1])],
                                  sparse_threshold=0.2)
    res = col_trans.fit_transform(X_array)
    assert not sparse.issparse(res)
    assert not col_trans.sparse_output_

    # mixed -> sparsity of (4 + 2) / 8 = 0.75
    for thres in [0.75001, 1]:
        col_trans = ColumnTransformer(
            [('trans1', OneHotEncoder(sparse=True), [0]),
             ('trans2', OneHotEncoder(sparse=False), [1])],
            sparse_threshold=thres)
        res = col_trans.fit_transform(X_array)
        assert sparse.issparse(res)
        assert col_trans.sparse_output_

    for thres in [0.75, 0]:
        col_trans = ColumnTransformer(
            [('trans1', OneHotEncoder(sparse=True), [0]),
             ('trans2', OneHotEncoder(sparse=False), [1])],
            sparse_threshold=thres)
        res = col_trans.fit_transform(X_array)
        assert not sparse.issparse(res)
        assert not col_trans.sparse_output_

    # if nothing is sparse -> no sparse
    for thres in [0.33, 0, 1]:
        col_trans = ColumnTransformer(
            [('trans1', OneHotEncoder(sparse=False), [0]),
             ('trans2', OneHotEncoder(sparse=False), [1])],
            sparse_threshold=thres)
        res = col_trans.fit_transform(X_array)
        assert not sparse.issparse(res)
        assert not col_trans.sparse_output_


def test_column_transformer_error_msg_1D():
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T

    col_trans = ColumnTransformer([('trans', StandardScaler(), 0)])
    assert_raise_message(ValueError, "1D data passed to a transformer",
                         col_trans.fit, X_array)
    assert_raise_message(ValueError, "1D data passed to a transformer",
                         col_trans.fit_transform, X_array)

    col_trans = ColumnTransformer([('trans', TransRaise(), 0)])
    for func in [col_trans.fit, col_trans.fit_transform]:
        assert_raise_message(ValueError, "specific message", func, X_array)


def test_2D_transformer_output():
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T

    # if one transformer is dropped, test that name is still correct
    ct = ColumnTransformer([('trans1', 'drop', 0),
                            ('trans2', TransNo2D(), 1)])
    assert_raise_message(ValueError, "the 'trans2' transformer should be 2D",
                         ct.fit_transform, X_array)
    # because fit is also doing transform, this raises already on fit
    assert_raise_message(ValueError, "the 'trans2' transformer should be 2D",
                         ct.fit, X_array)


def test_2D_transformer_output_pandas():
    pd = pytest.importorskip('pandas')

    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_df = pd.DataFrame(X_array, columns=['col1', 'col2'])

    # if one transformer is dropped, test that name is still correct
    ct = ColumnTransformer([('trans1', TransNo2D(), 'col1')])
    assert_raise_message(ValueError, "the 'trans1' transformer should be 2D",
                         ct.fit_transform, X_df)
    # because fit is also doing transform, this raises already on fit
    assert_raise_message(ValueError, "the 'trans1' transformer should be 2D",
                         ct.fit, X_df)


@pytest.mark.parametrize("remainder", ['drop', 'passthrough'])
def test_column_transformer_invalid_columns(remainder):
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T

    # general invalid
    for col in [1.5, ['string', 1], slice(1, 's'), np.array([1.])]:
        ct = ColumnTransformer([('trans', Trans(), col)], remainder=remainder)
        assert_raise_message(ValueError, "No valid specification",
                             ct.fit, X_array)

    # invalid for arrays
    for col in ['string', ['string', 'other'], slice('a', 'b')]:
        ct = ColumnTransformer([('trans', Trans(), col)], remainder=remainder)
        assert_raise_message(ValueError, "Specifying the columns",
                             ct.fit, X_array)


def test_column_transformer_invalid_transformer():

    class NoTrans(BaseEstimator):
        def fit(self, X, y=None):
            return self

        def predict(self, X):
            return X

    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    ct = ColumnTransformer([('trans', NoTrans(), [0])])
    assert_raise_message(TypeError, "All estimators should implement fit",
                         ct.fit, X_array)


def test_make_column_transformer():
    scaler = StandardScaler()
    norm = Normalizer()
    ct = make_column_transformer((scaler, 'first'), (norm, ['second']))
    names, transformers, columns = zip(*ct.transformers)
    assert_equal(names, ("standardscaler", "normalizer"))
    assert_equal(transformers, (scaler, norm))
    assert_equal(columns, ('first', ['second']))

    # XXX remove in v0.22
    with pytest.warns(DeprecationWarning,
                      match='`make_column_transformer` now expects'):
        ct1 = make_column_transformer(([0], norm))
    ct2 = make_column_transformer((norm, [0]))
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    assert_almost_equal(ct1.fit_transform(X_array),
                        ct2.fit_transform(X_array))

    with pytest.warns(DeprecationWarning,
                      match='`make_column_transformer` now expects'):
        make_column_transformer(('first', 'drop'))

    with pytest.warns(DeprecationWarning,
                      match='`make_column_transformer` now expects'):
        make_column_transformer(('passthrough', 'passthrough'),
                                ('first', 'drop'))


def test_make_column_transformer_pandas():
    pd = pytest.importorskip('pandas')
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_df = pd.DataFrame(X_array, columns=['first', 'second'])
    norm = Normalizer()
    # XXX remove in v0.22
    with pytest.warns(DeprecationWarning,
                      match='`make_column_transformer` now expects'):
        ct1 = make_column_transformer((X_df.columns, norm))
    ct2 = make_column_transformer((norm, X_df.columns))
    assert_almost_equal(ct1.fit_transform(X_df),
                        ct2.fit_transform(X_df))


def test_make_column_transformer_kwargs():
    scaler = StandardScaler()
    norm = Normalizer()
    ct = make_column_transformer((scaler, 'first'), (norm, ['second']),
                                 n_jobs=3, remainder='drop',
                                 sparse_threshold=0.5)
    assert_equal(ct.transformers, make_column_transformer(
        (scaler, 'first'), (norm, ['second'])).transformers)
    assert_equal(ct.n_jobs, 3)
    assert_equal(ct.remainder, 'drop')
    assert_equal(ct.sparse_threshold, 0.5)
    # invalid keyword parameters should raise an error message
    assert_raise_message(
        TypeError,
        'Unknown keyword arguments: "transformer_weights"',
        make_column_transformer, (scaler, 'first'), (norm, ['second']),
        transformer_weights={'pca': 10, 'Transf': 1}
    )


def test_make_column_transformer_remainder_transformer():
    scaler = StandardScaler()
    norm = Normalizer()
    remainder = StandardScaler()
    ct = make_column_transformer((scaler, 'first'), (norm, ['second']),
                                 remainder=remainder)
    assert ct.remainder == remainder


def test_column_transformer_get_set_params():
    ct = ColumnTransformer([('trans1', StandardScaler(), [0]),
                            ('trans2', StandardScaler(), [1])])

    exp = {'n_jobs': None,
           'remainder': 'drop',
           'sparse_threshold': 0.3,
           'trans1': ct.transformers[0][1],
           'trans1__copy': True,
           'trans1__with_mean': True,
           'trans1__with_std': True,
           'trans2': ct.transformers[1][1],
           'trans2__copy': True,
           'trans2__with_mean': True,
           'trans2__with_std': True,
           'transformers': ct.transformers,
           'transformer_weights': None}

    assert_dict_equal(ct.get_params(), exp)

    ct.set_params(trans1__with_mean=False)
    assert_false(ct.get_params()['trans1__with_mean'])

    ct.set_params(trans1='passthrough')
    exp = {'n_jobs': None,
           'remainder': 'drop',
           'sparse_threshold': 0.3,
           'trans1': 'passthrough',
           'trans2': ct.transformers[1][1],
           'trans2__copy': True,
           'trans2__with_mean': True,
           'trans2__with_std': True,
           'transformers': ct.transformers,
           'transformer_weights': None}

    assert_dict_equal(ct.get_params(), exp)


def test_column_transformer_named_estimators():
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T
    ct = ColumnTransformer([('trans1', StandardScaler(), [0]),
                            ('trans2', StandardScaler(with_std=False), [1])])
    assert_false(hasattr(ct, 'transformers_'))
    ct.fit(X_array)
    assert hasattr(ct, 'transformers_')
    assert isinstance(ct.named_transformers_['trans1'], StandardScaler)
    assert isinstance(ct.named_transformers_.trans1, StandardScaler)
    assert isinstance(ct.named_transformers_['trans2'], StandardScaler)
    assert isinstance(ct.named_transformers_.trans2, StandardScaler)
    assert_false(ct.named_transformers_.trans2.with_std)
    # check it are fitted transformers
    assert_equal(ct.named_transformers_.trans1.mean_, 1.)


def test_column_transformer_cloning():
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T

    ct = ColumnTransformer([('trans', StandardScaler(), [0])])
    ct.fit(X_array)
    assert_false(hasattr(ct.transformers[0][1], 'mean_'))
    assert hasattr(ct.transformers_[0][1], 'mean_')

    ct = ColumnTransformer([('trans', StandardScaler(), [0])])
    ct.fit_transform(X_array)
    assert_false(hasattr(ct.transformers[0][1], 'mean_'))
    assert hasattr(ct.transformers_[0][1], 'mean_')


def test_column_transformer_get_feature_names():
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T
    ct = ColumnTransformer([('trans', Trans(), [0, 1])])
    # raise correct error when not fitted
    assert_raises(NotFittedError, ct.get_feature_names)
    # raise correct error when no feature names are available
    ct.fit(X_array)
    assert_raise_message(AttributeError,
                         "Transformer trans (type Trans) does not provide "
                         "get_feature_names", ct.get_feature_names)

    # working example
    X = np.array([[{'a': 1, 'b': 2}, {'a': 3, 'b': 4}],
                  [{'c': 5}, {'c': 6}]], dtype=object).T
    ct = ColumnTransformer(
        [('col' + str(i), DictVectorizer(), i) for i in range(2)])
    ct.fit(X)
    assert_equal(ct.get_feature_names(), ['col0__a', 'col0__b', 'col1__c'])

    # passthrough transformers not supported
    ct = ColumnTransformer([('trans', 'passthrough', [0, 1])])
    ct.fit(X)
    assert_raise_message(
        NotImplementedError, 'get_feature_names is not yet supported',
        ct.get_feature_names)

    ct = ColumnTransformer([('trans', DictVectorizer(), 0)],
                           remainder='passthrough')
    ct.fit(X)
    assert_raise_message(
        NotImplementedError, 'get_feature_names is not yet supported',
        ct.get_feature_names)

    # drop transformer
    ct = ColumnTransformer(
        [('col0', DictVectorizer(), 0), ('col1', 'drop', 1)])
    ct.fit(X)
    assert_equal(ct.get_feature_names(), ['col0__a', 'col0__b'])


def test_column_transformer_special_strings():

    # one 'drop' -> ignore
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T
    ct = ColumnTransformer(
        [('trans1', Trans(), [0]), ('trans2', 'drop', [1])])
    exp = np.array([[0.], [1.], [2.]])
    assert_array_equal(ct.fit_transform(X_array), exp)
    assert_array_equal(ct.fit(X_array).transform(X_array), exp)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] != 'remainder'

    # all 'drop' -> return shape 0 array
    ct = ColumnTransformer(
        [('trans1', 'drop', [0]), ('trans2', 'drop', [1])])
    assert_array_equal(ct.fit(X_array).transform(X_array).shape, (3, 0))
    assert_array_equal(ct.fit_transform(X_array).shape, (3, 0))
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] != 'remainder'

    # 'passthrough'
    X_array = np.array([[0., 1., 2.], [2., 4., 6.]]).T
    ct = ColumnTransformer(
        [('trans1', Trans(), [0]), ('trans2', 'passthrough', [1])])
    exp = X_array
    assert_array_equal(ct.fit_transform(X_array), exp)
    assert_array_equal(ct.fit(X_array).transform(X_array), exp)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] != 'remainder'

    # None itself / other string is not valid
    for val in [None, 'other']:
        ct = ColumnTransformer(
            [('trans1', Trans(), [0]), ('trans2', None, [1])])
        assert_raise_message(TypeError, "All estimators should implement",
                             ct.fit_transform, X_array)
        assert_raise_message(TypeError, "All estimators should implement",
                             ct.fit, X_array)


def test_column_transformer_remainder():
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T

    X_res_first = np.array([0, 1, 2]).reshape(-1, 1)
    X_res_second = np.array([2, 4, 6]).reshape(-1, 1)
    X_res_both = X_array

    # default drop
    ct = ColumnTransformer([('trans1', Trans(), [0])])
    assert_array_equal(ct.fit_transform(X_array), X_res_first)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_first)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'drop'
    assert_array_equal(ct.transformers_[-1][2], [1])

    # specify passthrough
    ct = ColumnTransformer([('trans', Trans(), [0])], remainder='passthrough')
    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'passthrough'
    assert_array_equal(ct.transformers_[-1][2], [1])

    # column order is not preserved (passed through added to end)
    ct = ColumnTransformer([('trans1', Trans(), [1])],
                           remainder='passthrough')
    assert_array_equal(ct.fit_transform(X_array), X_res_both[:, ::-1])
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both[:, ::-1])
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'passthrough'
    assert_array_equal(ct.transformers_[-1][2], [0])

    # passthrough when all actual transformers are skipped
    ct = ColumnTransformer([('trans1', 'drop', [0])],
                           remainder='passthrough')
    assert_array_equal(ct.fit_transform(X_array), X_res_second)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_second)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'passthrough'
    assert_array_equal(ct.transformers_[-1][2], [1])

    # error on invalid arg
    ct = ColumnTransformer([('trans1', Trans(), [0])], remainder=1)
    assert_raise_message(
        ValueError,
        "remainder keyword needs to be one of \'drop\', \'passthrough\', "
        "or estimator.", ct.fit, X_array)
    assert_raise_message(
        ValueError,
        "remainder keyword needs to be one of \'drop\', \'passthrough\', "
        "or estimator.", ct.fit_transform, X_array)

    # check default for make_column_transformer
    ct = make_column_transformer((Trans(), [0]))
    assert ct.remainder == 'drop'


@pytest.mark.parametrize("key", [[0], np.array([0]), slice(0, 1),
                                 np.array([True, False])])
def test_column_transformer_remainder_numpy(key):
    # test different ways that columns are specified with passthrough
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_res_both = X_array

    ct = ColumnTransformer([('trans1', Trans(), key)],
                           remainder='passthrough')
    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'passthrough'
    assert_array_equal(ct.transformers_[-1][2], [1])


@pytest.mark.parametrize(
    "key", [[0], slice(0, 1), np.array([True, False]), ['first'], 'pd-index',
            np.array(['first']), np.array(['first'], dtype=object),
            slice(None, 'first'), slice('first', 'first')])
def test_column_transformer_remainder_pandas(key):
    # test different ways that columns are specified with passthrough
    pd = pytest.importorskip('pandas')
    if isinstance(key, six.string_types) and key == 'pd-index':
        key = pd.Index(['first'])

    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_df = pd.DataFrame(X_array, columns=['first', 'second'])
    X_res_both = X_array

    ct = ColumnTransformer([('trans1', Trans(), key)],
                           remainder='passthrough')
    assert_array_equal(ct.fit_transform(X_df), X_res_both)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][1] == 'passthrough'
    assert_array_equal(ct.transformers_[-1][2], [1])


@pytest.mark.parametrize("key", [[0], np.array([0]), slice(0, 1),
                                 np.array([True, False, False])])
def test_column_transformer_remainder_transformer(key):
    X_array = np.array([[0, 1, 2],
                        [2, 4, 6],
                        [8, 6, 4]]).T
    X_res_both = X_array.copy()

    # second and third columns are doubled when remainder = DoubleTrans
    X_res_both[:, 1:3] *= 2

    ct = ColumnTransformer([('trans1', Trans(), key)],
                           remainder=DoubleTrans())

    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert isinstance(ct.transformers_[-1][1], DoubleTrans)
    assert_array_equal(ct.transformers_[-1][2], [1, 2])


def test_column_transformer_no_remaining_remainder_transformer():
    X_array = np.array([[0, 1, 2],
                        [2, 4, 6],
                        [8, 6, 4]]).T

    ct = ColumnTransformer([('trans1', Trans(), [0, 1, 2])],
                           remainder=DoubleTrans())

    assert_array_equal(ct.fit_transform(X_array), X_array)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_array)
    assert len(ct.transformers_) == 1
    assert ct.transformers_[-1][0] != 'remainder'


def test_column_transformer_drops_all_remainder_transformer():
    X_array = np.array([[0, 1, 2],
                        [2, 4, 6],
                        [8, 6, 4]]).T

    # columns are doubled when remainder = DoubleTrans
    X_res_both = 2 * X_array.copy()[:, 1:3]

    ct = ColumnTransformer([('trans1', 'drop', [0])],
                           remainder=DoubleTrans())

    assert_array_equal(ct.fit_transform(X_array), X_res_both)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_both)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert isinstance(ct.transformers_[-1][1], DoubleTrans)
    assert_array_equal(ct.transformers_[-1][2], [1, 2])


def test_column_transformer_sparse_remainder_transformer():
    X_array = np.array([[0, 1, 2],
                        [2, 4, 6],
                        [8, 6, 4]]).T

    ct = ColumnTransformer([('trans1', Trans(), [0])],
                           remainder=SparseMatrixTrans(),
                           sparse_threshold=0.8)

    X_trans = ct.fit_transform(X_array)
    assert sparse.issparse(X_trans)
    # SparseMatrixTrans creates 3 features for each column. There is
    # one column in ``transformers``, thus:
    assert X_trans.shape == (3, 3 + 1)

    exp_array = np.hstack(
        (X_array[:, 0].reshape(-1, 1), np.eye(3)))
    assert_array_equal(X_trans.toarray(), exp_array)
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert isinstance(ct.transformers_[-1][1], SparseMatrixTrans)
    assert_array_equal(ct.transformers_[-1][2], [1, 2])


def test_column_transformer_drop_all_sparse_remainder_transformer():
    X_array = np.array([[0, 1, 2],
                        [2, 4, 6],
                        [8, 6, 4]]).T
    ct = ColumnTransformer([('trans1', 'drop', [0])],
                           remainder=SparseMatrixTrans(),
                           sparse_threshold=0.8)

    X_trans = ct.fit_transform(X_array)
    assert sparse.issparse(X_trans)

    #  SparseMatrixTrans creates 3 features for each column, thus:
    assert X_trans.shape == (3, 3)
    assert_array_equal(X_trans.toarray(), np.eye(3))
    assert len(ct.transformers_) == 2
    assert ct.transformers_[-1][0] == 'remainder'
    assert isinstance(ct.transformers_[-1][1], SparseMatrixTrans)
    assert_array_equal(ct.transformers_[-1][2], [1, 2])


def test_column_transformer_get_set_params_with_remainder():
    ct = ColumnTransformer([('trans1', StandardScaler(), [0])],
                           remainder=StandardScaler())

    exp = {'n_jobs': None,
           'remainder': ct.remainder,
           'remainder__copy': True,
           'remainder__with_mean': True,
           'remainder__with_std': True,
           'sparse_threshold': 0.3,
           'trans1': ct.transformers[0][1],
           'trans1__copy': True,
           'trans1__with_mean': True,
           'trans1__with_std': True,
           'transformers': ct.transformers,
           'transformer_weights': None}

    assert ct.get_params() == exp

    ct.set_params(remainder__with_std=False)
    assert not ct.get_params()['remainder__with_std']

    ct.set_params(trans1='passthrough')
    exp = {'n_jobs': None,
           'remainder': ct.remainder,
           'remainder__copy': True,
           'remainder__with_mean': True,
           'remainder__with_std': False,
           'sparse_threshold': 0.3,
           'trans1': 'passthrough',
           'transformers': ct.transformers,
           'transformer_weights': None}

    assert ct.get_params() == exp


def test_column_transformer_no_estimators():
    X_array = np.array([[0, 1, 2],
                        [2, 4, 6],
                        [8, 6, 4]]).astype('float').T
    ct = ColumnTransformer([], remainder=StandardScaler())

    params = ct.get_params()
    assert params['remainder__with_mean']

    X_trans = ct.fit_transform(X_array)
    assert X_trans.shape == X_array.shape
    assert len(ct.transformers_) == 1
    assert ct.transformers_[-1][0] == 'remainder'
    assert ct.transformers_[-1][2] == [0, 1, 2]


def test_column_transformer_no_estimators_set_params():
    ct = ColumnTransformer([]).set_params(n_jobs=2)
    assert ct.n_jobs == 2


def test_column_transformer_callable_specifier():
    # assert that function gets the full array / dataframe
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T
    X_res_first = np.array([[0, 1, 2]]).T

    def func(X):
        assert_array_equal(X, X_array)
        return [0]

    ct = ColumnTransformer([('trans', Trans(), func)],
                           remainder='drop')
    assert_array_equal(ct.fit_transform(X_array), X_res_first)
    assert_array_equal(ct.fit(X_array).transform(X_array), X_res_first)
    assert callable(ct.transformers[0][2])
    assert ct.transformers_[0][2] == [0]

    pd = pytest.importorskip('pandas')
    X_df = pd.DataFrame(X_array, columns=['first', 'second'])

    def func(X):
        assert_array_equal(X.columns, X_df.columns)
        assert_array_equal(X.values, X_df.values)
        return ['first']

    ct = ColumnTransformer([('trans', Trans(), func)],
                           remainder='drop')
    assert_array_equal(ct.fit_transform(X_df), X_res_first)
    assert_array_equal(ct.fit(X_df).transform(X_df), X_res_first)
    assert callable(ct.transformers[0][2])
    assert ct.transformers_[0][2] == ['first']
