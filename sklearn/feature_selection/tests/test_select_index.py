from scipy.sparse import csc_matrix

from sklearn.feature_selection import SelectColumnsByIndicies

data = [[0, 1, 2, 3, 4], [0, 2, 2, 3, 5], [1, 1, 2, 4, 0]]


def test_select_index():
    # Test SelectColumnsByIndicies.
    for X in [data, csc_matrix(data)]:
        X = SelectColumnsByIndicies(indicies=[0, 2, 4]).fit_transform(X)
        assert (len(data), 3) == X.shape
