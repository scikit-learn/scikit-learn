"""Restricted Boltzmann Machine
"""

# Main author: Yann N. Dauphin <dauphiya@iro.umontreal.ca>
# Author: Vlad Niculae
# Author: Gabriel Synnaeve
# License: BSD Style.

import time

import numpy as np

from ..base import BaseEstimator
from ..base import TransformerMixin
from ..externals.six.moves import xrange
from ..utils import check_arrays
from ..utils import check_random_state
from ..utils import gen_even_slices
from ..utils import issparse
from ..utils.extmath import safe_sparse_dot
from ..utils.extmath import logistic_sigmoid


class BernoulliRBM(BaseEstimator, TransformerMixin):
    """Bernoulli Restricted Boltzmann Machine (RBM).

    A Restricted Boltzmann Machine with binary visible units and
    binary hiddens. Parameters are estimated using Stochastic Maximum
    Likelihood (SML), also known as Persistent Contrastive Divergence (PCD)
    [2].

    The time complexity of this implementation is ``O(d ** 2)`` assuming
    d ~ n_features ~ n_components.

    Parameters
    ----------
    n_components : int, optional
        Number of binary hidden units.

    learning_rate : float, optional
        The learning rate for weight updates. It is *highly* recommended
        to tune this hyper-parameter. Reasonable values are in the
        10**[0., -3.] range.

    batch_size : int, optional
        Number of examples per minibatch.

    n_iter : int, optional
        Number of iterations/sweeps over the training dataset to perform
        during training.

    verbose : int, optional
        The verbosity level. Enabling it (with a non-zero value) will compute 
        the log-likelihood of each mini-batch and hence cause a runtime overhead
        in the order of 10%.

    random_state : integer or numpy.RandomState, optional
        A random number generator instance to define the state of the
        random permutations generator. If an integer is given, it fixes the
        seed. Defaults to the global numpy random number generator.

    Attributes
    ----------
    `components_` : array-like, shape (n_components, n_features), optional
        Weight matrix, where n_features in the number of visible
        units and n_components is the number of hidden units.

    `intercept_hidden_` : array-like, shape (n_components,), optional
        Biases of the hidden units.

    `intercept_visible_` : array-like, shape (n_features,), optional
        Biases of the visible units.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.neural_network import BernoulliRBM
    >>> X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
    >>> model = BernoulliRBM(n_components=2)
    >>> model.fit(X)
    BernoulliRBM(batch_size=10, learning_rate=0.1, n_components=2, n_iter=10,
           random_state=None, verbose=False)

    References
    ----------

    [1] Hinton, G. E., Osindero, S. and Teh, Y. A fast learning algorithm for
        deep belief nets. Neural Computation 18, pp 1527-1554.
        http://www.cs.toronto.edu/~hinton/absps/fastnc.pdf

    [2] Tieleman, T. Training Restricted Boltzmann Machines using
        Approximations to the Likelihood Gradient. International Conference
        on Machine Learning (ICML) 2008
    """
    def __init__(self, n_components=256, learning_rate=0.1, batch_size=10,
                 n_iter=10, verbose=False, random_state=None):
        self.n_components = n_components
        self.learning_rate = learning_rate
        self.batch_size = batch_size
        self.n_iter = n_iter
        self.verbose = verbose
        self.random_state = random_state

    def transform(self, X):
        """Compute the hidden layer activation probabilities, P(h=1|v=X).

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The data to be transformed.

        Returns
        -------
        h : array, shape (n_samples, n_components)
            Latent representations of the data.
        """
        X, = check_arrays(X, sparse_format='csr', dtype=np.float)
        return self._mean_hiddens(X)

    def _mean_hiddens(self, v):
        """Computes the probabilities P(h=1|v).

        Parameters
        ----------
        v : array-like, shape (n_samples, n_features)
            Values of the visible layer.

        Returns
        -------
        h : array-like, shape (n_samples, n_components)
            Corresponding mean field values for the hidden layer.
        """
        return logistic_sigmoid(safe_sparse_dot(v, self.components_.T)
                                + self.intercept_hidden_)

    def _sample_hiddens(self, v, rng):
        """Sample from the distribution P(h|v).

        Parameters
        ----------
        v : array-like, shape (n_samples, n_features)
            Values of the visible layer to sample from.

        rng : RandomState
            Random number generator to use.

        Returns
        -------
        h : array-like, shape (n_samples, n_components)
            Values of the hidden layer.
        """
        p = self._mean_hiddens(v)
        p[rng.uniform(size=p.shape) < p] = 1.
        return np.floor(p, p)

    def _sample_visibles(self, h, rng):
        """Sample from the distribution P(v|h).

        Parameters
        ----------
        h : array-like, shape (n_samples, n_components)
            Values of the hidden layer to sample from.

        rng : RandomState
            Random number generator to use.

        Returns
        -------
        v : array-like, shape (n_samples, n_features)
            Values of the visible layer.
        """
        p = logistic_sigmoid(np.dot(h, self.components_)
                             + self.intercept_visible_)
        p[rng.uniform(size=p.shape) < p] = 1.
        return np.floor(p, p)

    def _free_energy(self, v):
        """Computes the free energy F(v) = - log sum_h exp(-E(v,h)).

        Parameters
        ----------
        v : array-like, shape (n_samples, n_features)
            Values of the visible layer.

        Returns
        -------
        free_energy : array-like, shape (n_samples,)
            The value of the free energy.
        """
        return (- safe_sparse_dot(v, self.intercept_visible_)
                - np.log(1. + np.exp(safe_sparse_dot(v, self.components_.T)
                            + self.intercept_hidden_)).sum(axis=1))

    def gibbs(self, v):
        """Perform one Gibbs sampling step.

        Parameters
        ----------
        v : array-like, shape (n_samples, n_features)
            Values of the visible layer to start from.

        Returns
        -------
        v_new : array-like, shape (n_samples, n_features)
            Values of the visible layer after one Gibbs step.
        """
        rng = check_random_state(self.random_state)
        h_ = self._sample_hiddens(v, rng)
        v_ = self._sample_visibles(h_, rng)

        return v_

    def _fit(self, v_pos, rng):
        """Inner fit for one mini-batch.

        Adjust the parameters to maximize the likelihood of v using
        Stochastic Maximum Likelihood (SML).

        Parameters
        ----------
        v_pos : array-like, shape (n_samples, n_features)
            The data to use for training.

        rng : RandomState
            Random number generator to use for sampling.

        Returns
        -------
        pseudo_likelihood : array-like, shape (n_samples,)
            If verbose=True, pseudo-likelihood estimate for this batch.
        """
        h_pos = self._mean_hiddens(v_pos)
        v_neg = self._sample_visibles(self.h_samples_, rng)
        h_neg = self._mean_hiddens(v_neg)

        lr = float(self.learning_rate) / v_pos.shape[0]
        update = safe_sparse_dot(v_pos.T, h_pos, dense_output=True).T
        update -= np.dot(v_neg.T, h_neg).T
        self.components_ += lr * update
        self.intercept_hidden_ += lr * (h_pos.sum(axis=0) - h_neg.sum(axis=0))
        self.intercept_visible_ += lr * (np.asarray(
                                         v_pos.sum(axis=0)).squeeze() -
                                         v_neg.sum(axis=0))

        h_neg[rng.uniform(size=h_neg.shape) < h_neg] = 1.0  # sample binomial
        self.h_samples_ = np.floor(h_neg, h_neg)

        if self.verbose:
            return self.score_samples(v_pos)

    def score_samples(self, v):
        """Compute the pseudo-likelihood of v.

        Parameters
        ----------
        v : array-like, shape (n_samples, n_features)
            Values of the visible layer.

        Returns
        -------
        pseudo_likelihood : array-like, shape (n_samples,)
            Value of the pseudo-likelihood (proxy to likelihood).
        """
        rng = check_random_state(self.random_state)
        fe = self._free_energy(v)

        if issparse(v):
            v_ = v.toarray()
        else:
            v_ = v.copy()
        i_ = rng.randint(0, v.shape[1], v.shape[0])
        v_[np.arange(v.shape[0]), i_] = 1 - v_[np.arange(v.shape[0]), i_]
        fe_ = self._free_energy(v_)

        return v.shape[1] * logistic_sigmoid(fe_ - fe, log=True)

    def fit(self, X, y=None):
        """Fit the model to the data X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        self : BernoulliRBM
            The fitted model.
        """
        X, = check_arrays(X, sparse_format='csr', dtype=np.float)
        n_samples = X.shape[0]
        rng = check_random_state(self.random_state)

        self.components_ = np.asarray(
            rng.normal(0, 0.01, (self.n_components, X.shape[1])),
            order='fortran')
        self.intercept_hidden_ = np.zeros(self.n_components, )
        self.intercept_visible_ = np.zeros(X.shape[1], )
        self.h_samples_ = np.zeros((self.batch_size, self.n_components))

        n_batches = int(np.ceil(float(n_samples) / self.batch_size))
        batch_slices = list(gen_even_slices(n_batches * self.batch_size,
                                            n_batches))
        verbose = self.verbose
        for iteration in xrange(self.n_iter):
            pl = 0.
            if verbose:
                begin = time.time()

            for batch_slice in batch_slices:
                pl_batch = self._fit(X[batch_slice], rng)

                if verbose:
                    pl += pl_batch.sum()

            if verbose:
                pl /= n_samples
                end = time.time()
                print("Iteration %d, pseudo-likelihood = %.2f, time = %.2fs"
                      % (iteration, pl, end - begin))

        return self
