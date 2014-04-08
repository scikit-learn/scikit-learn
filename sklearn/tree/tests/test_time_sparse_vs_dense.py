"""
Testing for the tree module (sklearn.tree).
"""
import sys
sys.path.insert(0, "../../..")

import pickle

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, coo_matrix
from functools import partial
import time


from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal

from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import ExtraTreeClassifier
from sklearn.tree import ExtraTreeRegressor

from sklearn import tree
from sklearn import datasets


CLF_CRITERIONS = ("gini", "entropy")
REG_CRITERIONS = ("mse", )

CLF_TREES = {
    "DecisionTreeClassifier": DecisionTreeClassifier,
    "Presort-DecisionTreeClassifier": partial(DecisionTreeClassifier,
                                              splitter="presort-best"),
    "ExtraTreeClassifier": ExtraTreeClassifier,
}

REG_TREES = {
    "DecisionTreeRegressor": DecisionTreeRegressor,
    "Presort-DecisionTreeRegressor": partial(DecisionTreeRegressor,
                                             splitter="presort-best"),
    "ExtraTreeRegressor": ExtraTreeRegressor,
}

ALL_TREES = dict()
ALL_TREES.update(CLF_TREES)
ALL_TREES.update(REG_TREES)


X_small = [[0, 0, 4, 0, 0, 0, 1, -14, 0, -4, 0, 0, 0, 0, ],
          [0, 0, 5, 3, 0, -4, 0, 0, 1, -5, 0.2, 0, 4, 1, ],
          [-1, -1, 0, 0, -4.5, 0, 0, 2.1, 1, 0, 0, -4.5, 0, 1, ],
          [-1, -1, 0, -1.2, 0, 0, 0, 0, 0, 0, 0.2, 0, 0, 1, ],
          [-1, -1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1, ],
          [-1, -2, 0, 4, -3, 10, 4, 0, -3.2, 0, 4, 3, -4, 1, ],
          [2.11, 0, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0.5, 0, -3, 1, ],
          [2.11, 0, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0, 0, -2, 1, ],
          [2.11, 8, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0, 0, -2, 1, ],
          [2.11, 8, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0.5, 0, -1, 0, ],
          [2, 8, 5, 1, 0.5, -4, 10, 0, 1, -5, 3, 0, 2, 0, ],
          [2, 0, 1, 1, 1, -1, 1, 0, 0, -2, 3, 0, 1, 0, ],
          [2, 0, 1, 2, 3, -1, 10, 2, 0, -1, 1, 2, 2, 0, ],
          [1, 1, 0, 2, 2, -1, 1, 2, 0, -5, 1, 2, 3, 0, ],
          [3, 1, 0, 3, 0, -4, 10, 0, 1, -5, 3, 0, 3, 1, ],
          [2.11, 8, -6, -0.5, 0, 1, 0, 0, -3.2, 6, 0.5, 0, -3, 1, ],
          [2.11, 8, -6, -0.5, 0, 1, 0, 0, -3.2, 6, 1.5, 1, -1, -1, ],
          [2.11, 8, -6, -0.5, 0, 10, 0, 0, -3.2, 6, 0.5, 0, -1, -1, ],
          [2, 0, 5, 1, 0.5, -2, 10, 0, 1, -5, 3, 1, 0, -1, ],
          [2, 0, 1, 1, 1, -2, 1, 0, 0, -2, 0, 0, 0, 1, ],
          [2, 1, 1, 1, 2, -1, 10, 2, 0, -1, 0, 2, 1, 1, ],
          [1, 1, 0, 0, 1, -3, 1, 2, 0, -5, 1, 2, 1, 1, ],
          [3, 1, 0, 1, 0, -4, 1, 0, 1, -2, 0, 0, 1, 0, ]]

y_small = [1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0,
           0, 0]
y_small_reg = [1.0, 2.1, 1.2, 0.05, 10, 2.4, 3.1, 1.01, 0.01, 2.98, 3.1, 1.1,
               0.0, 1.2, 2, 11, 0, 0, 4.5, 0.201, 1.06, 0.9, 0]

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(1)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def assert_tree_enquality(d, s, message):
    from sklearn.tree._tree import TREE_LEAF
    assert_equal(d.node_count, s.node_count,
                 message + ": inequal number of node")

    assert_array_equal(d.children_right, s.children_right,
                       message + ": inequal n_node_samples")
    assert_array_equal(d.children_left, s.children_left,
                       message + ": inequal children_left")

    external = d.children_right == TREE_LEAF
    internal = np.logical_not(external)

    assert_array_equal(d.feature[internal], s.feature[internal],
                       message + ": inequal features")
    assert_array_equal(d.threshold[internal], s.threshold[internal],
                       message + ": inequal threshold")
    assert_array_equal(d.n_node_samples, s.n_node_samples,
                       message + ": inequal n_node_samples")

    assert_almost_equal(d.impurity, s.impurity,
                        err_msg=message + ": inequal impurity")

    assert_array_almost_equal(d.value[external], s.value[external],
                              err_msg=message + ": inequal value")


def time_test_random_sparse_matrix_Reg():
    #nr = 10000
    #nf = 18000
    if True:
        nr = 40
        nf = 20

        n_test = 50
        for it in range(5):
            X_ = np.random.randint(100, size=(nr, nf))
            y_ = np.random.randint(10, size=nr)
            X_test = csr_matrix(np.random.randint(100, size=(n_test, nf)))

            for j in range(X_.shape[1]):
                for i in range(X_.shape[0]):
                    r = np.random.uniform()
                    if r < 0.99:
                        X_[i, j] = 0
                    # if r < 0.2:
                    #     X_[i, j] *= -1
                    # elif r > 0.5:
                    #     X_[i, j] = 0
            s = RandomForestRegressor(n_estimators=10,
                                       max_depth=100,
                                       random_state=0).fit(csc_matrix(X_), y_)
            start = time.time()
            s.predict(X_test)
            end = time.time()
            print "sparse time : ", end - start

            start = time.time()
            s.predict(X_test.todense())
            end = time.time()
            print "dense time : ", end - start

            d = RandomForestRegressor(n_estimators=1000,
                                       max_depth=100,
                                       random_state=0).fit(X_, y_)

            start = time.time()
            d.predict(X_test)
            end = time.time()
            print "sparse time : ", end - start

            start = time.time()
            d.predict(X_test.todense())
            end = time.time()
            print "dense time : ", end - start


def time_test_random_sparse_matrix_Cls():
    #nr = 10000
    #nf = 18000
    if True:
        nr = 40
        nf = 20

        n_test = 50
        for it in range(5):
            X_ = np.random.randint(100, size=(nr, nf))
            y_ = np.random.randint(10, size=nr)
            X_test = csr_matrix(np.random.randint(100, size=(n_test, nf)))

            for j in range(X_.shape[1]):
                for i in range(X_.shape[0]):
                    r = np.random.uniform()
                    if r < 0.99:
                        X_[i, j] = 0
                    # if r < 0.2:
                    #     X_[i, j] *= -1
                    # elif r > 0.5:
                    #     X_[i, j] = 0
            s = RandomForestClassifier(n_estimators=10,
                                       max_depth=100,
                                       random_state=0).fit(csc_matrix(X_), y_)
            start = time.time()
            s.predict(X_test)
            end = time.time()
            print "sparse time : ", end - start

            start = time.time()
            s.predict(X_test)
            end = time.time()
            print "dense time : ", end - start

            d = RandomForestClassifier(n_estimators=10,
                                       max_depth=100,
                                       random_state=0).fit(X_, y_)

            start = time.time()
            d.predict(X_test)
            end = time.time()
            print "sparse time : ", end - start

            start = time.time()
            d.predict(X_test)
            end = time.time()
            print "dense time : ", end - start



if __name__ == '__main__':

    print "RandomForestRegressor"
    #time_test_random_sparse_matrix_Reg()

    print "RandomForestClassifier"
    time_test_random_sparse_matrix_Cls()
