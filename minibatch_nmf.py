
import numpy as np
from scipy import sparse

from sklearn.utils import check_random_state
from sklearn.utils.extmath import row_norms, safe_sparse_dot
from sklearn.base import BaseEstimator, TransformerMixin
# from sklearn.utils import check_array

from sklearn.cluster.k_means_ import _k_init
from sklearn.decomposition.nmf import _special_sparse_dot


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

    hashing: boolean, default=False
        If true, HashingVectorizer is used instead of CountVectorizer.

    hashing_n_features: int, default=2**10
        Number of features for the HashingVectorizer. Only relevant if
        hashing=True.

    hashing: boolean, default=True
        If true, the weight matrix W is rescaled at each iteration
        to have an l1 norm equal to 1 for each row.

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
                 r=.001, hashing=False,
                 hashing_n_features=2**12, init='k-means++',
                 tol=1E-4, min_iter=2, max_iter=5, ngram_range=(2, 4),
                 add_words=False, random_state=None,
                 rescale_W=True, max_iter_e_step=20):

        self.n_components = n_components
        self.r = r
        self.batch_size = batch_size
        self.tol = tol
        self.hashing = hashing
        self.hashing_n_features = hashing_n_features
        self.max_iter = max_iter
        self.min_iter = min_iter
        self.init = init
        self.add_words = add_words
        self.random_state = check_random_state(random_state)
        self.rescale_W = rescale_W
        self.max_iter_e_step = max_iter_e_step

    def _rescale_W(self, W, A, B):
        epsilon = 1E-10
        s = W.sum(axis=1, keepdims=True)
        s[s == 0] = epsilon
        W /= s
        A /= s
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
            W_WT1 = W / WT1.reshape(-1, 1)
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
        self.rho = self.r ** (1 / (iter + 1))
        A += W * safe_sparse_dot(Ht.T, Ht_W) * self.rho
        B += Ht.sum(axis=0).reshape(-1, 1) * self.rho
        np.divide(A, B, out=W, where=(B != 0))
        if self.rescale_W:
            W, A, B = self._rescale_W(A / B, A, B)
        return W, A, B

    def _get_H(self, X):
        H_out = np.empty((len(X), self.n_components))
        for x, h_out in zip(X, H_out):
            h_out[:] = self.H_dict[x]
        return H_out

    def _init_W(self, V):
        if self.init == 'k-means++':
            W = _k_init(
                V, self.n_components, row_norms(V, squared=True),
                random_state=self.random_state,
                n_local_trials=None) + .1
        elif self.init == 'random':
            W = self.random_state.gamma(
                shape=1, scale=1,
                size=(self.n_components, self.n_vocab))
        else:
            raise AttributeError(
                'Initialization method %s does not exist.' % self.init)
        W /= W.sum(axis=1, keepdims=True)
        A = np.ones((self.n_components, self.n_vocab)) * 1E-10
        B = A.copy()
        return W, A, B

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
        # needs to be changed to check is X contains strings or not
        if sparse.issparse(X):
            n_samples, self.n_vocab = X.shape
            H = np.ones((n_samples, self.n_components))
            H = self._rescale_H(X, H)
            self.W, self.A, self.B = self._init_W(X)
            # self.rho = self.r**(self.batch_size / n_samples)
        # else:
            # not implemented yet

        n_batch = (n_samples - 1) // self.batch_size + 1
        self.iter = 1

        for iter in range(self.max_iter):
            for i, (Ht, Vt) in enumerate(mini_batch(H, X, n=self.batch_size)):
                if i == n_batch-1:
                    W_last = self.W
                Ht[:] = self._e_step(Vt, self.W, Ht,
                                     max_iter=self.max_iter_e_step)
                self.W, self.A, self.B = self._m_step(Vt, self.W,
                                                      self.A, self.B, Ht,
                                                      self.iter)
                self.iter += 1
                if i == n_batch-1:
                    W_change = np.linalg.norm(
                        self.W - W_last) / np.linalg.norm(W_last)
            if (W_change < self.tol) and (iter >= self.min_iter - 1):
                break
        return self

    def partial_fit(self, X, y=None):
        if hasattr(self, 'iter'):
            assert X.shape[1] == self.n_vocab
            if sparse.issparse(X):
                n_samples, _ = X.shape
                H = np.ones((n_samples, self.n_components))
                H = self._rescale_H(X, H)
            # else:
                # not implemented yet
        else:
            if sparse.issparse(X):
                n_samples, self.n_vocab = X.shape
                H = np.ones((n_samples, self.n_components))
                H = self._rescale_H(X, H)
                self.W, self.A, self.B = self._init_W(X)
                self.iter = 1
                # self.rho = self.r**(self.batch_size / n_samples)
            # else:
                # not implemented yet

        for i, (Ht, Vt) in enumerate(mini_batch(H, X, n=self.batch_size)):
            Ht[:] = self._e_step(Vt, self.W, Ht,
                                 max_iter=self.max_iter_e_step)
            self.W, self.A, self.B = self._m_step(
                Vt, self.W, self.A, self.B, Ht, self.iter)
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
        assert X.shape[1] == self.n_vocab
        n_samples, _ = X.shape

        H = np.ones((n_samples, self.n_components))
        H = self._rescale_H(X, H)

        for Ht, Vt in mini_batch(H, X, n=self.batch_size):
            Ht[:] = self._e_step(Vt, self.W, Ht, max_iter=50)
        return H


def mini_batch(iterable1, iterable2, n=1):
    len_iter = len(iterable1)
    for idx in range(0, len_iter, n):
        this_slice = slice(idx, min(idx + n, len_iter))
        yield (iterable1[this_slice],
               iterable2[this_slice])
