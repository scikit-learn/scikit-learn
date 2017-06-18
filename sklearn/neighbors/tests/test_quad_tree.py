import pickle
import numpy as np
from sklearn.neighbors.quad_tree import QuadTree


def test_quadtree_similar_point():
    # Introduce a point into a quad tree where a similar point already exists.
    # Test will hang if it doesn't complete.
    Xs = []

    # check the case where points are actually different
    Xs.append(np.array([[1, 2], [3, 4]], dtype=np.float32))
    # check the case where points are the same on X axis
    Xs.append(np.array([[1.0, 2.0], [1.0, 3.0]], dtype=np.float32))
    # check the case where points are arbitrarily close on X axis
    Xs.append(np.array([[1.00001, 2.0], [1.00002, 3.0]], dtype=np.float32))
    # check the case where points are the same on Y axis
    Xs.append(np.array([[1.0, 2.0], [3.0, 2.0]], dtype=np.float32))
    # check the case where points are arbitrarily close on Y axis
    Xs.append(np.array([[1.0, 2.00001], [3.0, 2.00002]], dtype=np.float32))
    # check the case where points are arbitrarily close on both axes
    Xs.append(np.array([[1.00001, 2.00001], [1.00002, 2.00002]],
              dtype=np.float32))

    # check the case where points are arbitrarily close on both axes
    # close to machine epsilon - x axis
    Xs.append(np.array([[1, 0.0003817754041], [2, 0.0003817753750]],
              dtype=np.float32))

    # check the case where points are arbitrarily close on both axes
    # close to machine epsilon - y axis
    Xs.append(np.array([[0.0003817754041, 1.0], [0.0003817753750, 2.0]],
              dtype=np.float32))

    for X in Xs:
        tree = QuadTree(n_dimensions=2, verbose=0)
        tree.build_tree(X)
        tree.check_coherence()


def test_quad_tree_pickle():
    np.random.seed(0)
    X = np.random.random((10, 3))

    tree = QuadTree(n_dimensions=2, verbose=0)
    tree.build_tree(X)

    def check_pickle_protocol(protocol):
        s = pickle.dumps(tree, protocol=protocol)
        bt2 = pickle.loads(s)

        for x in X:
            assert tree.get_cell(x) == bt2.get_cell(x)

    for protocol in (0, 1, 2):
        yield check_pickle_protocol, protocol
