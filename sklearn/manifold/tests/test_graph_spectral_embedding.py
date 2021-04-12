import unittest
import numpy as np

from sklearn.manifold import GraphSpectralEmbedding

def _er_nm(n, m, directed=False):
    A = np.zeros((n, n))
    if directed:
        idx = np.where(~np.eye(n, dtype=bool))
    else:
        idx = np.triu_indices(n, k=1)
    # get idx in 1d coordinates by ravelling
    triu = np.ravel_multi_index(idx, A.shape)
    # choose M of them
    triu = np.random.choice(triu, size=m, replace=False)
    # unravel back
    triu = np.unravel_index(triu, A.shape)
    A[triu] = 1

    if not directed:
        A = np.triu(A)
        A = A + A.T - np.diag(np.diag(A))
    return A
def _test_output_dim_directed(self, method):
    n_components = 4
    embed = method(n_components=n_components, concat=True)
    n = 10
    M = 20
    A = _er_nm(n, M, directed=True) + 5
    self.assertEqual(embed.fit_transform(A).shape, (n, 8))
    self.assertEqual(embed.latent_left_.shape, (n, 4))
    self.assertEqual(embed.latent_right_.shape, (n, 4))


def _test_output_dim(self, method, sparse=False, *args, **kwargs):
    n_components = 4
    embed = method(n_components=n_components)
    n = 10
    M = 20
    A = _er_nm(n, M) + 5
    if sparse:
        A = csr_matrix(A)
    embed.fit(A)
    self.assertEqual(embed.latent_left_.shape, (n, 4))
    self.assertTrue(embed.latent_right_ is None)

class TestAdjacencySpectralEmbed(unittest.TestCase):
    # def setUp(self):
    #     np.random.seed(9001)
    #     n = [10, 10]
    #     p = np.array([[0.9, 0.1], [0.1, 0.9]])
    #     wt = [[normal, poisson], [poisson, normal]]
    #     wtargs = [
    #         [dict(loc=3, scale=1), dict(lam=5)],
    #         [dict(lam=5), dict(loc=3, scale=1)],
    #     ]
    #     self.testgraphs = dict(
    #         Guw=sbm(n=n, p=p),
    #         Gw=sbm(n=n, p=p, wt=wt, wtargs=wtargs),
    #         Guwd=sbm(n=n, p=p, directed=True),
    #         Gwd=sbm(n=n, p=p, wt=wt, wtargs=wtargs, directed=True),
    #     )
    #     self.ase = AdjacencySpectralEmbed(n_components=2)

    def test_output_dim(self):
        _test_output_dim(self, GraphSpectralEmbedding)

    def test_output_dim_directed(self):
        _test_output_dim_directed(self, GraphSpectralEmbedding)
