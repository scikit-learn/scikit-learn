""" Restricted Boltzmann Machine
"""

# Author: Yann N. Dauphin <dauphiya@iro.umontreal.ca>
# License: BSD Style.

import numpy as np

from .base import BaseEstimator, TransformerMixin
from .utils import array2d, check_random_state
from .utils.extmath import safe_sparse_dot


class RestrictedBolzmannMachine(BaseEstimator, TransformerMixin):
    """
    Restricted Boltzmann Machine (RBM)
    
    A Restricted Boltzmann Machine with binary visible units and
    binary hiddens. Parameters are estimated using Stochastic Maximum
    Likelihood (SML).
    
    The time complexity of this implementation is ``O(d ** 2)`` assuming
    d ~ n_features ~ n_components.
    
    Parameters
    ----------
    n_components : int, optional
        Number of binary hidden units
    learning_rate : float, optional
        Learning rate to use during learning. It is *highly* recommended
        to tune this hyper-parameter. Possible values are 10**[0., -3.].
    n_particles : int, optional
        Number of MCMC particles to use during learning.
    n_epochs : int, optional
        Number of epochs/sweeps over the training dataset to perform
        during training.
    verbose: bool, optional
        When True (False by default) the method outputs the progress
        of learning after each epoch.
    random_state : RandomState or an int seed (0 by default)
        A random number generator instance to define the state of the
        random permutations generator.
    
    Attributes
    ----------
    components_ : array-like, shape (n_features, n_components), optional
        Weight matrix, where n_features in the number of visible
        units and n_components is the number of hidden units.
    intercept_hidden_ : array-like, shape (n_components,), optional
        Biases of the hidden units
    intercept_visible_ : array-like, shape (n_features,), optional
        Biases of the visible units
    
    Examples
    --------
    
    >>> import numpy as np
    >>> from sklearn.rbm import RestrictedBolzmannMachine
    >>> X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
    >>> model = RestrictedBolzmannMachine(n_components=2)
    >>> model.fit(X)
    
    References
    ----------
    
    [1] Hinton, G. E., Osindero, S. and Teh, Y. A fast learning algorithm for
        deep belief nets. Neural Computation 18, pp 1527-1554.
        http://www.cs.toronto.edu/~hinton/absps/fastnc.pdf
    """
    def __init__(self, n_components=1024,
                       learning_rate=0.1,
                       n_particles=10,
                       n_epochs=10,
                       verbose=False,
                       random_state=0):
        self.n_components = n_components
        self.learning_rate = learning_rate
        self.n_particles = n_particles
        self.n_epochs = n_epochs
        self.verbose = verbose
        self.random_state = check_random_state(random_state)
    
    def _logistic_sigmoid(self, x):
        """
        Implements the logistic function.
        
        Parameters
        ----------
        x: array-like, shape (M, N)

        Returns
        -------
        x_new: array-like, shape (M, N)
        """
        return 1. / (1. + np.exp(-np.maximum(np.minimum(x, 30), -30)))
    
    def _sample_binomial(self, p):
        """
        Compute the element-wise binomial using the probabilities p.
        
        Parameters
        ----------
        x: array-like, shape (M, N)

        Notes
        -----
        This is equivalent to calling numpy.random.binomial(1, p) but is
        faster because it uses in-place operations on p.

        Returns
        -------
        x_new: array-like, shape (M, N)
        """
        p[self.random_state.uniform(size=p.shape) < p] = 1.
        
        return np.floor(p, p)
    
    def transform(self, X):
        """
        Computes the probabilities P({\bf h}_j=1|{\bf v}={\bf X}).
        
        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)

        Returns
        -------
        h: array-like, shape (n_samples, n_components)
        """
        return self.mean_h(X)
    
    def mean_h(self, v):
        """
        Computes the probabilities P({\bf h}_j=1|{\bf v}).
        
        Parameters
        ----------
        v: array-like, shape (n_samples, n_features)

        Returns
        -------
        h: array-like, shape (n_samples, n_components)
        """
        return self._logistic_sigmoid(safe_sparse_dot(v, self.components_)
            + self.intercept_hidden_)
    
    def sample_h(self, v):
        """
        Sample from the distribution P({\bf h}|{\bf v}).
        
        Parameters
        ----------
        v: array-like, shape (n_samples, n_features)
        
        Returns
        -------
        h: array-like, shape (n_samples, n_components)
        """
        return self._sample_binomial(self.mean_h(v))
    
    def mean_v(self, h):
        """
        Computes the probabilities P({\bf v}_i=1|{\bf h}).
        
        Parameters
        ----------
        h: array-like, shape (n_samples, n_components)
        
        Returns
        -------
        v: array-like, shape (n_samples, n_features)
        """
        return self._logistic_sigmoid(np.dot(h, self.components_.T)
            + self.intercept_visible_)
    
    def sample_v(self, h):
        """
        Sample from the distribution P({\bf v}|{\bf h}).
        
        Parameters
        ----------
        h: array-like, shape (n_samples, n_components)
        
        Returns
        -------
        v: array-like, shape (n_samples, n_features)
        """
        return self._sample_binomial(self.mean_v(h))
    
    def free_energy(self, v):
        """
        Computes the free energy
        \mathcal{F}({\bf v}) = - \log \sum_{\bf h} e^{-E({\bf v},{\bf h})}.
        
        Parameters
        ----------
        v: array-like, shape (n_samples, n_features)
        
        Returns
        -------
        free_energy: array-like, shape (n_samples,)
        """
        return - np.dot(v, self.intercept_visible_) - np.log(1. + np.exp(
            safe_sparse_dot(v, self.components_)+self.intercept_hidden_)).sum(1)
    
    def gibbs(self, v):
        """
        Perform one Gibbs sampling step.
        
        Parameters
        ----------
        v: array-like, shape (n_samples, n_features)
        
        Returns
        -------
        v_new: array-like, shape (n_samples, n_features)
        """
        h_ = self.sample_h(v)
        v_ = self.sample_v(h_)
        
        return v_
    
    def _fit(self, v_pos):
        """
        Adjust the parameters to maximize the likelihood of {\bf v}
        using Stochastic Maximum Likelihood (SML) [1].
        
        Parameters
        ----------
        v_pos: array-like, shape (n_samples, n_features)
        
        Returns
        -------
        pseudo_likelihood: array-like, shape (n_samples,), optional
            Pseudo Likelihood estimate for this batch.
        
        References
        ----------
        [1] Tieleman, T. Training Restricted Boltzmann Machines using
            Approximations to the Likelihood Gradient. International Conference
            on Machine Learning (ICML) 2008
        """
        h_pos = self.mean_h(v_pos)
        v_neg = self.sample_v(self.h_samples_)
        h_neg = self.mean_h(v_neg)
        
        lr = self.learning_rate / self.n_particles
        self.components_ += safe_sparse_dot(lr * v_pos.T, h_pos)
        self.components_ -= np.dot(lr * v_neg.T, h_neg)
        self.intercept_hidden_ += lr * (h_pos.sum(0) - h_neg.sum(0))
        self.intercept_visible_ += lr * (v_pos.sum(0) - v_neg.sum(0))
        
        self.h_samples_ = self._sample_binomial(h_neg)
        
        return self.pseudo_likelihood(v_pos)
    
    def pseudo_likelihood(self, v):
        """
        Compute the pseudo-likelihood of {\bf v}.
        
        Parameters
        ----------
        v: array-like, shape (n_samples, n_features)
        
        Returns
        -------
        pseudo_likelihood: array-like, shape (n_samples,)
        """
        fe = self.free_energy(v)
        
        v_ = v.copy()
        i_ = self.random_state.randint(0, v.shape[1], v.shape[0])
        v_[range(v.shape[0]), i_] = v_[range(v.shape[0]), i_] == 0
        fe_ = self.free_energy(v_)
        
        return v.shape[1] * np.log(self._logistic_sigmoid(fe_ - fe))
    
    def fit(self, X, y=None):
        """
        Fit the model to the data X.
        
        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        """
        X = array2d(X)
        
        self.components_ = np.asarray(self.random_state.normal(0, 0.01,
            (X.shape[1], self.n_components)), dtype=X.dtype)
        self.intercept_hidden_ = np.zeros(self.n_components, dtype=X.dtype)
        self.intercept_visible_ = np.zeros(X.shape[1], dtype=X.dtype)
        self.h_samples_ = np.zeros((self.n_particles, self.n_components),
            dtype=X.dtype)
        
        inds = range(X.shape[0])
        
        np.random.shuffle(inds)
        
        n_batches = int(np.ceil(len(inds) / float(self.n_particles)))
        
        for epoch in range(self.n_epochs):
            pl = 0.
            for minibatch in range(n_batches):
                pl += self._fit(X[inds[minibatch::n_batches]]).sum()
            pl /= X.shape[0]
            
            if self.verbose:
                print "Epoch %d, Pseudo-Likelihood = %.2f" % (epoch, pl)
    
    def fit_transform(self, X, y=None):
        """
        Fit the model to the data X and transform it.
        
        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.
        """
        self.fit(X, y)
        
        return self.transform(X)
