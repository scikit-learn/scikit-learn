"""
Benchmarks np.bincount method vs np.mean for feature agglomeration in
../sklearn/cluster/_feature_agglomeration. Use of np.bincount provides
a significant speed up if the pooling function is np.mean.

np.bincount performs better especially as the size of X and n_clusters
increase.
"""

import numpy as np
from sklearn.cluster import FeatureAgglomeration
import time


def fit_agglomeration():
    rng = np.random.RandomState(0)
    X = rng.randn(100000, 1000)
    agglo = FeatureAgglomeration(n_clusters=5)
    agglo.fit(X)
    return X, agglo


def get_transformed_array(X, agglo, method):
    size = np.bincount(agglo.labels_)
    n_samples = X.shape[0]
    nX = []
    if len(agglo.labels_) != X.shape[1]:
        raise ValueError("X has a different number of features than "
                         "during fitting.")
    if method == "bincount":
        # a fast way to compute the mean of grouped features
        nX = np.array([np.bincount(agglo.labels_, X[i, :]) / size
                       for i in range(n_samples)])
    elif method == "np_mean":
        for l in np.unique(agglo.labels_):
            nX.append(np.mean(X[:, agglo.labels_ == l], axis=1))
        nX = np.array(nX).T
    else:
        raise ValueError("Method can have a value of 'bincount' or 'np.mean'")
    return nX


if __name__ == "__main__":
    X, agglo = fit_agglomeration()

    tick = time.time()
    result_bincount = get_transformed_array(X, agglo, "bincount")
    time_bincount = time.time() - tick

    tick = time.time()
    result_np_mean = get_transformed_array(X, agglo, "np_mean")
    time_np_mean = time.time() - tick

    print('==================')
    print('Took %s seconds using np.bincount' % (time_bincount))
    print('Took %s seconds using np.mean' % (time_np_mean))
    print('==================')
    print("np.bincount is %s times faster" % (time_np_mean/time_bincount))
