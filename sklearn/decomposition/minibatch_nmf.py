import numpy as np
from scipy import sparse

from sklearn.utils import check_random_state
from sklearn.utils.extmath import row_norms, safe_sparse_dot, randomized_svd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import gen_batches
# from sklearn.utils import check_array

from sklearn.cluster.k_means_ import _k_init
from sklearn.decomposition.nmf import _special_sparse_dot
from sklearn.decomposition.nmf import norm


class MiniBatchNMF(BaseEstimator, TransformerMixin):
    """
    Mini batch non-negative matrix factorization by minimizing the
    Kullback-Leibler divergence.

    Parameters
    ----------

    n_components: int, default=10
        Number of topics of the matrix Factorization.

    batch_size: int, default=100

    r: float, default=1
        Weight parameter for the update of the W matrix

    tol: float, default=1E-3
        Tolerance for the convergence of the matrix W

    mix_iter: int, default=2

    max_iter: int, default=10

    ngram_range: tuple, default=(2, 4)

    init: str, default 'k-means++'
        Initialization method of the W matrix.

    random_state: default=None

    Attributes
    ----------

    References
    ----------
    """

    def __init__(self, n_components=10, batch_size=512,
                 r=.001, init='k-means++',
                 tol=1E-4, min_iter=2, max_iter=5, ngram_range=(2, 4),
                 add_words=False, random_state=None,
                 rescale_W=True, max_iter_e_step=20):

        self.n_components = n_components
        self.r = r
        self.batch_size = batch_size
        self.tol = tol
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.init = init
        self.add_words = add_words
        self.random_state = check_random_state(random_state)
        self.rescale_W = rescale_W
        self.max_iter_e_step = max_iter_e_step

    def _rescale_W(self, W, A, B):
        s = W.sum(axis=1, keepdims=True)
        np.divide(W, s, out=W, where=(s != 0))
        np.divide(A, s, out=A, where=(s != 0))
        return W, A, B

    def _rescale_H(self, V, H):
        epsilon = 1e-10  # in case of a document having length=0
        H *= np.maximum(epsilon, V.sum(axis=1).A)
        H /= H.sum(axis=1, keepdims=True)
        return H

    def _e_step(self, Vt, W, Ht,
                tol=1E-3, max_iter=20):
        if self.rescale_W:
            W_WT1 = W
        else:
            WT1 = np.sum(W, axis=1)
            W_WT1 = W / WT1[:, np.newaxis]
        squared_tol = tol**2
        squared_norm = 1
        for iter in range(max_iter):
            if squared_norm <= squared_tol:
                break
            Ht_W = _special_sparse_dot(Ht, W, Vt)
            Ht_W_data = Ht_W.data
            Vt_data = Vt.data
            np.divide(Vt_data, Ht_W_data, out=Ht_W_data,
                      where=(Ht_W_data != 0))
            Ht_out = Ht * safe_sparse_dot(Ht_W, W_WT1.T)
            squared_norm = np.linalg.norm(
                Ht_out - Ht) / (np.linalg.norm(Ht) + 1E-10)
            Ht[:] = Ht_out
        return Ht

    def _m_step(self, Vt, W, A, B, Ht, iter):
        Ht_W = _special_sparse_dot(Ht, W, Vt)
        Ht_W_data = Ht_W.data
        np.divide(Vt.data, Ht_W_data, out=Ht_W_data, where=(Ht_W_data != 0))
        self.rho_ = self.r ** (1 / iter)
        # self.rho_ = .98
        A *= self.rho_
        A += W * safe_sparse_dot(Ht.T, Ht_W)
        B *= self.rho_
        B += Ht.sum(axis=0).reshape(-1, 1)
        np.divide(A, B, out=W, where=(W != 0))
        if self.rescale_W:
            W, A, B = self._rescale_W(W, A, B)
        return W, A, B

    def _get_H(self, X):
        H_out = np.empty((len(X), self.n_components))
        for x, h_out in zip(X, H_out):
            h_out[:] = self.H_dict[x]
        return H_out

    def _init_vars(self, V):
        if self.init == 'k-means++':
            W = _k_init(
                V, self.n_components, row_norms(V, squared=True),
                random_state=self.random_state,
                n_local_trials=None) + .1
            W /= W.sum(axis=1, keepdims=True)
            H = np.ones((V.shape[0], self.n_components))
            H = self._rescale_H(V, H)
        elif self.init == 'random':
            W = self.random_state.gamma(
                shape=1, scale=1,
                size=(self.n_components, self.n_features_))
            W /= W.sum(axis=1, keepdims=True)
            H = np.ones((V.shape[0], self.n_components))
            H = self._rescale_H(V, H)
        elif self.init == 'nndsvd':
            eps = 1e-6
            U, S, V = randomized_svd(V, self.n_components,
                                     random_state=self.random_state)
            H, W = np.zeros(U.shape), np.zeros(V.shape)

            # The leading singular triplet is non-negative
            # so it can be used as is for initialization.
            H[:, 0] = np.sqrt(S[0]) * np.abs(U[:, 0])
            W[0, :] = np.sqrt(S[0]) * np.abs(V[0, :])

            for j in range(1, self.n_components):
                x, y = U[:, j], V[j, :]

                # extract positive and negative parts of column vectors
                x_p, y_p = np.maximum(x, 0), np.maximum(y, 0)
                x_n, y_n = np.abs(np.minimum(x, 0)), np.abs(np.minimum(y, 0))

                # and their norms
                x_p_nrm, y_p_nrm = norm(x_p), norm(y_p)
                x_n_nrm, y_n_nrm = norm(x_n), norm(y_n)

                m_p, m_n = x_p_nrm * y_p_nrm, x_n_nrm * y_n_nrm

                # choose update
                if m_p > m_n:
                    u = x_p / x_p_nrm
                    v = y_p / y_p_nrm
                    sigma = m_p
                else:
                    u = x_n / x_n_nrm
                    v = y_n / y_n_nrm
                    sigma = m_n

                lbd = np.sqrt(S[j] * sigma)
                H[:, j] = lbd * u
                W[j, :] = lbd * v

            W[W < eps] = 0
            H[H < eps] = 0
            H = np.ones((V.shape[0], self.n_components))
            H = self._rescale_H(V, H)
        else:
            raise AttributeError(
                'Initialization method %s does not exist.' % self.init)
        A = W.copy()
        B = np.ones((self.n_components, self.n_features_))
        return H, W, A, B

    def fit(self, X, y=None):
        """Fit the NMF to X.

        Parameters
        ----------
        X : string array-like, shape [n_samples, n_features]
            The data to determine the categories of each feature
        Returns
        -------
        self
        """
        n_samples, self.n_features_ = X.shape

        if sparse.issparse(X):
            H, self.W_, self.A_, self.B_ = self._init_vars(X)
            # self.rho_ = self.r**(self.batch_size / n_samples)
        # else:
            # not implemented yet

        n_batch = (n_samples - 1) // self.batch_size + 1
        self.iter = 1

        for iter in range(self.max_iter):
            for i, slice in enumerate(gen_batches(n=n_samples,
                                                  batch_size=self.batch_size)):
                if i == n_batch-1:
                    W_last = self.W_
                H[slice] = self._e_step(X[slice], self.W_, H[slice],
                                        max_iter=self.max_iter_e_step)
                self.W_, self.A_, self.B_ = self._m_step(
                    X[slice], self.W_, self.A_, self.B_, H[slice], self.iter)
                self.iter += 1
                if i == n_batch-1:
                    W_change = np.linalg.norm(
                        self.W_ - W_last) / np.linalg.norm(W_last)
            if (W_change < self.tol) and (iter >= self.min_iter - 1):
                break
        return self

    def partial_fit(self, X, y=None):
        if hasattr(self, 'iter'):
            assert X.shape[1] == self.n_features_
            n_samples, _ = X.shape

            if sparse.issparse(X):
                H = np.ones((n_samples, self.n_components))
                H = self._rescale_H(X, H)
            # else:
                # not implemented yet
        else:
            n_samples, self.n_features_ = X.shape

            if sparse.issparse(X):
                # H = np.ones((n_samples, self.n_components))
                # H = self._rescale_H(X, H)
                H, self.W_, self.A_, self.B_ = self._init_vars(X)
                self.iter = 1
                # self.rho = self.r**(self.batch_size / n_samples)
            # else:
                # not implemented yet

        for slice in gen_batches(n=n_samples, batch_size=self.batch_size):
            H[slice] = self._e_step(X[slice], self.W_, H[slice],
                                    max_iter=self.max_iter_e_step)
            self.W_, self.A_, self.B_ = self._m_step(
                X[slice], self.W_, self.A_, self.B_, H[slice], self.iter)
            self.iter += 1

    def transform(self, X):
        """Transform X using the trained matrix W.

        Parameters
        ----------
        X : array-like (str), shape [n_samples,]
            The data to encode.

        Returns
        -------
        X_new : 2-d array, shape [n_samples, n_components]
            Transformed input.
        """
        assert X.shape[1] == self.n_features_
        n_samples, _ = X.shape

        H = np.ones((n_samples, self.n_components))
        H = self._rescale_H(X, H)

        for slice in gen_batches(n=n_samples, batch_size=self.batch_size):
            H[slice] = self._e_step(X[slice], self.W_, H[slice], max_iter=50)
        return H
