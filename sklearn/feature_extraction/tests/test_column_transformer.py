import numpy as np

from sklearn.base import BaseEstimator
from sklearn.feature_extraction import ColumnTransformer

from sklearn.utils.testing import assert_array_equal


class Trans(BaseEstimator):
    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        return X


def test_fields():
    # dictionary
    X_dict = {'first': [[0], [1], [2]],
              'second': [[2], [4], [6]]}
    # recarray
    X_recarray = np.recarray((3, 1),
                             dtype=[('first', np.int), ('second', np.int)])
    X_recarray['first'] = X_dict['first']
    X_recarray['second'] = X_dict['second']

    for X in [X_dict, X_recarray]:
        first_feat = ColumnTransformer({'trans': (Trans(), 'first')})
        second_feat = ColumnTransformer({'trans': (Trans(), 'second')})
        both = ColumnTransformer({'trans1': (Trans(), 'first'),
                                  'trans2': (Trans(), 'second')})
        assert_array_equal(first_feat.fit_transform(X), X['first'])
        assert_array_equal(second_feat.fit_transform(X), X['second'])
        assert_array_equal(both.fit_transform(X), np.hstack([X['first'], X['second']]))
        # fit then transform
        assert_array_equal(first_feat.fit(X).transform(X), X['first'])
        assert_array_equal(second_feat.fit(X).transform(X), X['second'])
        assert_array_equal(both.fit(X).transform(X), np.hstack([X['first'], X['second']]))
