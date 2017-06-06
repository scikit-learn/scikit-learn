from itertools import chain, product

import numpy as np
import scipy.sparse as sp

from sklearn.base import BaseEstimator
from sklearn.feature_extraction import ColumnTransformer

from sklearn.utils.testing import assert_array_equal, assert_equal, assert_true
from sklearn.utils.validation import check_array


class Trans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        #TODO fix this in ColumnTransformer to always pass 2D data
        return check_array(np.asarray(X).reshape(-1, 1))


class SparseMatrixTrans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        n_samples = len(X)
        return sp.eye(n_samples, n_samples).tocsr()


def test_column_selection():
    # dictionary
    X_dict = {'first': [0, 1, 2],
              'second': [2, 4, 6]}
    # recarray
    X_recarray = np.recarray((3,),
                             dtype=[('first', np.int), ('second', np.int)])
    X_recarray['first'] = X_dict['first']
    X_recarray['second'] = X_dict['second']

    # array
    X_array = np.array([[0, 1, 2], [2, 4, 6]]).T

    Xs_name = [X_dict, X_recarray]
    Xs_positional = [X_array]

    # dataframe
    try:
        import pandas as pd
        X_df = pd.DataFrame(X_dict)
        Xs_name.append(X_df)
        Xs_positional.append(X_df)
    except:
        print("Pandas not found, not testing ColumnTransformer with"
              " DataFrame.")
    X_res_first = np.array(X_dict['first']).reshape(-1, 1)
    X_res_second = np.array(X_dict['second']).reshape(-1, 1)
    X_res_both = np.vstack([X_dict['first'], X_dict['second']]).T

    for X, (first, second) in chain(product(Xs_name, [('first', 'second')]),
                                    product(Xs_positional, [(0, 1)])):

        first_feat = ColumnTransformer({'trans': (Trans(), first)})
        second_feat = ColumnTransformer({'trans': (Trans(), second)})
        both = ColumnTransformer({'trans1': (Trans(), first),
                                  'trans2': (Trans(), second)})
        assert_array_equal(first_feat.fit_transform(X), X_res_first)
        assert_array_equal(second_feat.fit_transform(X), X_res_second)
        assert_array_equal(both.fit_transform(X), X_res_both)
        # fit then transform
        assert_array_equal(first_feat.fit(X).transform(X), X_res_first)
        assert_array_equal(second_feat.fit(X).transform(X), X_res_second)
        assert_array_equal(both.fit(X).transform(X), X_res_both)



    # test with transformer_weights
    transformer_weights = {'trans1': .1, 'trans2': 10}
    for X in Xs_name:
        both = ColumnTransformer({'trans1': (Trans(), 'first'),
                                  'trans2': (Trans(), 'second')},
                                 transformer_weights=transformer_weights)
        res = np.vstack([transformer_weights['trans1'] * np.array(X['first']),
                         transformer_weights['trans2'] * np.array(X['second'])]).T
        assert_array_equal(both.fit_transform(X), res)
        # fit then transform
        assert_array_equal(both.fit(X).transform(X), res)


def test_sparse_stacking():
    X_dict = {'first': [0, 1, 2],
              'second': [2, 4, 6]}
    col_trans = ColumnTransformer({'trans1': (Trans(), 'first'), 'trans2':
                                   (SparseMatrixTrans(), 'second')})
    col_trans.fit(X_dict)
    X_trans = col_trans.transform(X_dict)
    assert_true(sp.issparse(X_trans))
    assert_equal(X_trans.shape, (X_trans.shape[0], X_trans.shape[0] + 1))
    assert_array_equal(X_trans.toarray()[:, 1:], np.eye(X_trans.shape[0]))
