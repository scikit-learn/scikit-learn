import numpy as np

from numpy.testing import assert_array_almost_equal
from scikits.learn import neighbors, manifold
from scikits.learn.utils.fixes import product

eigen_solvers = ['auto', 'dense', 'arpack']
path_methods = ['auto', 'FW', 'D']


def test_isomap_simple_grid():
    # Isomap should preserve distances when all neighbors are used
    N_per_side = 5
    Npts = N_per_side ** 2
    n_neighbors = Npts - 1

    # grid of equidistant points in 2D, out_dim = n_dim
    X = np.array(list(product(range(N_per_side), repeat=2)))

    # distances from each point to all others
    G = neighbors.kneighbors_graph(X, n_neighbors,
                                   mode='distance').toarray()

    clf = manifold.Isomap(n_neighbors=n_neighbors, out_dim=2)

    for eigen_solver in eigen_solvers:
        for path_method in path_methods:
            clf.fit(X, eigen_solver=eigen_solver,
                    path_method=path_method)

            G_iso = neighbors.kneighbors_graph(clf.embedding_,
                                               n_neighbors,
                                               mode='distance').toarray()
            assert_array_almost_equal(G, G_iso)


if __name__ == '__main__':
    import nose
    nose.runmodule()
