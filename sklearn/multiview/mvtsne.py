# -*- coding: utf-8 -*-
"""
MV-tSNE with log-linear opinion pooling, using GAs to optimize the weights of
each view (opinion).

It computes the multiview spectral clustering of data on a list of matrices
or distance matrices (or a mix of both), supposed to be different views of
the same data.
The function is split in two parts, one that computes P and other that
performs tsne on a given P plus an "interface" function that does everything
efficient way.
"""

___version___ = '1.0'
___author___ = 'Maria Araceli Burgueño Caballero'
___email___ = "mburgueno@uoc.edu"
___status___ = "Pre-Production"


import numpy as np
from sklearn.multiview import x2p
from sklearn.multiview import whiten as whiten_function
import itertools
import warnings
from numbers import Number
from scipy.optimize import fmin_l_bfgs_b
from scipy.special import rel_entr
from sklearn.base import BaseEstimator
from sklearn.utils import check_random_state


def find_pooling(X, is_distance, random_state, initial_dims=30, perplexity=30,
                 whiten=True, method='best'):
    """ Compute optimal pooling of P's for MV-tSNE.

    Parameters
    ----------

    X : list of matrices.
        A list of feature matrices or distance matrices.
    is_distance: array-like.
        A list or array which indicates whether a matrix with the same index
        in x is a distance matrix (true value) or not (false value).
    initial_dims : integer, default 30
        Number of dimensions to use in the reduction method.
    perplexity : integer, defuult 30
        This perplexity parameter is roughly equivalent to the optimal number
        of neighbours.
    whiten : integer, default 1.
        A boolean value indicating if the data matrices should be whitened.
    method : string
        best, log, linear (best computes both and selects the best according
        to Abbas 2009)

    Returns
    -------

    result: tuple.
        The weights and the pooled probability matrix.

    """

    nviews = len(X)
    eps = np.finfo(float).eps

    # Preprocess input matrices
    for i in np.arange(nviews):
        if not is_distance[i]:
            X[i] -= np.min(X[i])
            X[i] = X[i] / np.max(X[i])
            initial_dims = min(initial_dims, X[i].shape[1])
            if whiten:
                X[i] = whiten_function(X[i], n_comp=initial_dims)
        n = X[i].shape[0]

    # Compute the probability distribution of each input view
    # Now P is an array of nviews x n x n  (instead of n x n)
    P = np.zeros((nviews, n, n))
    for i in np.arange(nviews):
        P[i] = x2p(X[i], perplexity, 1e-05)[0]
        P[i] = 0.5 * (P[i] + P[i].T)
        P[i][P[i] < eps] = eps
        P[i] = P[i] / np.sum(P[i])

    ########################################################
    # Find optimum weight combination

    def log_linear_pooling(P, weights):
        """ Computes log.linear pooled opinion.

        Parameters
        ----------

        P : ndarray.
            Input data for the pooling (matrix).
        weights : numeric or array_like
            Weights to be used in the pooling.

        Returns
        -------

        result : tuple
            Tuple with two values:

                Pooled matrix (norm)

                Reg constant (1/sum)

        """

        P_exp_w = np.zeros(P.shape)
        if isinstance(weights, Number) or len(weights) == 1:
            # Apply function with weights value
            for i in np.arange(P.shape[0]):
                P_exp_w[i] = P[i]**weights
        else:
            # Array of weights. Check if length of P and weights are the same
            if len(P) % len(weights) != 0:
                raise Exception(
                    'Length of P and weights are not multiple. Cannot work \
                    out function')
            for i in np.arange(P.shape[0]):
                P_exp_w[i] = P[i]**weights[i]
        pooled = np.apply_over_axes(np.prod, P_exp_w, [0])[0]
        reg_const = 1 / np.sum(pooled)
        pooled *= reg_const
        return (pooled, reg_const)

    def objective_log_linear(weights):
        """General objective function for log-linear pooling
        (Abbas 2009 (9))

        Parameters
        ----------

        weights : numeric or array_like
            Input data.

        Returns
        -------

        result : float
            Log-linear pooled probability.

        """

        # Compute log-linear pooled prob with given weights
        pooling_pooled, pooling_reg_const = log_linear_pooling(P, weights)

        # Compute log-linear payoff (Abbas (9)) (here higher is worse)
        kls = np.zeros(nviews)
        pooling_pooled_p = 1.0 * pooling_pooled / np.sum(pooling_pooled)
        for i, qk in enumerate(P):
            qk = 1.0 * qk / np.sum(qk)
            vec = rel_entr(pooling_pooled_p, qk)
            kls[i] = np.sum(vec)

        payoff = np.sum(np.dot(kls, weights))

        # Introduce constraint sum(weights)=1 through a penalty
        penalty = abs(1 - np.sum(weights))
        goal = payoff + penalty
        return (-goal)

    def gradient_log_linear(weights):
        """Compute the gradient of the objective log linear.

        Parameters
        ----------

        weights : numeric or array_like
            Input weights.

        Returns
        -------

        result : ndarray
            Gradient of the objective.

        """

        pooling_pooled, pooling_reg_const = log_linear_pooling(P, weights)
        log_pooling = np.log(pooling_pooled)
        res = np.zeros(nviews)
        for i in np.arange(nviews):
            res[i] = np.sum(weights[i] * log_pooling / P[i])
        return res

    # Run the optimizer to find the best set of weights. Start point at (1/v,..
    # Using bounds to limit the weights to 0..1
    x0 = random_state.rand(nviews)  # initial guess/values
    bounds = tuple(itertools.repeat((0, 1), nviews))
    opt = fmin_l_bfgs_b(objective_log_linear, x0, fprime=gradient_log_linear,
                        bounds=bounds)

    pooled_log_p, reg_const = log_linear_pooling(P, opt[0])

    ########################################################
    # Compute KL between the different P's
    kl_ps = np.zeros((nviews, nviews), dtype='float64')
    for i in np.arange(nviews):
        for j in np.arange(nviews):
            kl_ps[i, j] = np.sum(
                P[i] * np.log((P[i] + eps) / (P[j] + eps)))

    kl_central = np.zeros(nviews)
    for i in np.arange(nviews):
        kl_central[i] = np.sum(P[i] * np.log(
            (P[i] + eps) / (pooled_log_p + eps)))

    return (pooled_log_p, opt[0], kl_ps, kl_central)


def tsne_p(P, k=2, initial_dims=30, max_iter=1000, min_cost=0,
           epoch_callback=None, epoch=100):
    """ Application of tSNE (t-distributed Stochastic Neighbor Embedding).

    Parameters
    ----------

    P: ndarray or matrix.
        Input data.
    k : integer, default 2
        The desired dimension of the resulting embedding.
    initial_dims : integer, default 30
        Number of dimensions to use in the reduction method.
    perplexity : integer, defuult 30
        This perplexity parameter is roughly equivalent to the optimal number
        of neighbours.
    max_iter : integer, default 1000
        Maximum number of iterations to perform.
    min_cost : numeric, default 0
        The minimum cost value (error) to stop iterations.
    epoch_callback : numeric, default None
        A callback function to be called after each epoch (which is a number
        of iterations controlled parameter epoch, see next).
    epoch : integer, default 100
        The number of iterations between update messages.

    Returns
    -------

    ydata: ndarray.
        Embedded space of input data.

    """

    n = P.shape[0]  # number of rows
    momentum = 0.5
    final_momentum = 0.8
    mom_switch_iter = 250

    epsilon = 500
    min_gain = 0.01
    initial_P_gain = 4
    eps = np.finfo(float).eps
    ydata = np.matrix(np.random.normal(size=(n, k)))

    # Prepare for standard tSNE
    P = P * initial_P_gain

    grads = np.zeros((ydata.shape))
    incs = np.zeros((ydata.shape))
    gains = np.ones((ydata.shape))

    cost = float("Inf")
    Q = None

    for ite in (range(1, max_iter + 1)):
        if ite % epoch == 0:
            cost = np.sum(np.sum(P * np.log((P + eps) / (Q + eps)), axis=1))

        if cost < min_cost:
            break
        if epoch_callback is not None:
            epoch_callback(ydata)

        sum_ydata = np.sum(np.square(ydata), axis=1)
        sum_ydata = np.squeeze(np.asarray(sum_ydata))
        num = 1 / (1 + sum_ydata +
                   ((-2 * np.dot(ydata, ydata.T))[:, ] + sum_ydata[:]).T)
        np.fill_diagonal(num, 0)
        Q = num / np.sum(num)
        if np.isnan(num).any():
            print("NaN in grad. descent")
        Q[Q < eps] = eps
        stiffnesses = 4 * (P - Q) * num
        for i in range(n):
            ydata_plus = (-ydata)[:, ] + ydata[i, :]
            mat_mul = np.multiply(ydata_plus.T, np.squeeze(
                np.asarray(stiffnesses[:, i]))).T
            grads[i] = np.sum(mat_mul, axis=0)

        gains = (gains + 0.2) * np.absolute(np.sign(grads) != np.sign(incs)) +\
            gains * 0.8 * np.absolute(np.sign(grads) == np.sign(incs))
        gains[gains < min_gain] = min_gain

        incs = momentum * incs - epsilon * (gains * grads)

        ydata += incs
        ydata -= np.mean(ydata, axis=0)

        if ite == mom_switch_iter:
            momentum = final_momentum
        if ite == 100:
            P /= 4
    # end ite for loop
    return ydata


class MvtSNE(BaseEstimator):
    """Multiview tSNE using an expert opinion pooling on the input probability
    matrices.

    Given a list of of input views and other parameters, mvtsne computes a
    neighbouring probability matrix for each input view, then finds the
    optimal set of weights to combine these matrices using a log-linear pool,
    and applies the pooled probability matrix as input to the standard tSNE
    procedure, where the probability matrix of the output space is adjusted
    to the pooled probability matrix using Kullback-Liebler divergence.

    Notes
    -----
    All input views must have the same number of samples (rows).

    Parameters
    ----------
    k : int, default: 2
        The desired dimension of the resulting embedding.
    initial_dims : int, default: 30
        Number of dimensions to use in the reduction method.
    perplexity : int, defuult: 30
        This perplexity parameter is roughly equivalent to the optimal number
        of neighbours.
    max_iter : int, default: 1000
        Maximum number of iterations to perform.
    min_cost : numeric, default: 0
        The minimum cost value (error) to stop iterations.
    epoch_callback : callable, default None
        A callback function to be called after each epoch (which is a number
        of iterations controlled parameter epoch, see next).
    whiten : int, default: 1
        A boolean value indicating if the data matrices should be whitened.
    epoch : int, default: 100
        The number of iterations between update messages.

    References
    ----------

        Abbas, Ali E. 2009. “A Kullback-Leibler View of Linear and Log-Linear
        Pools.” *Decision Analysis* 6 (1): 25–37. doi:10.1287/deca.1080.0133.

        Carvalho, Arthur, and Kate Larson. 2012. “A Consensual Linear Opinion
        Pool.” http://arxiv.org/abs/1204.5399.

        Van Der Maaten, Laurens, Geoffrey Hinton, and Geoffrey Hinton van der
        Maaten. 2008. “Visualizing Data using t-SNE.”
        doi:10.1007/s10479-011-0841-3.
    """

    def __init__(self, k=2, initial_dims=30, perplexity=30, max_iter=1000,
                 min_cost=0, epoch_callback=None, whiten=True, epoch=100,
                 random_state=0):
        self.k = k
        self.initial_dims = initial_dims
        self.perplexity = perplexity
        self.max_iter = max_iter
        self.min_cost = min_cost
        self.epoch_callback = epoch_callback
        self.whiten = whiten
        self.epoch = epoch
        self.random_state = random_state

    def fit(self, X, is_distance):
        """
        Computes standard tSNE algorithm to input multiview data. Return
        the weights used in the algorithm and the probabilitmatrix.

        Notes
        -----
        All input views must have the same number of samples (rows).

        Parameters
        ----------
        X : list
            A list of feature matrices or distance matrices, where each matrix
            is one of the views of the dataset.
        is_distance: array-like.
            A list or array which indicates whether a matrix with the same
            index in x is a distance matrix (true value) or not (false value).
        """
        self.fit_transform(X, is_distance)
        return self

    def fit_transform(self, X, is_distance):
        """
        Computes standard tSNE algorithm to input multiview data. Return
        the weights used in the algorithm and the probability matrix.

        Notes
        -----
        All input views must have the same number of samples (rows).

        Parameters
        ----------
        X : list
            A list of feature matrices or distance matrices, where each matrix
            is one of the views of the dataset.
        is_distance: array-like.
            A list or array which indicates whether a matrix with the same
            index in x is a distance matrix (true value) or not (false value).

        Attributes
        ----------
        embedding_: ndarray.
            Embedded space.
        weights_: ndarray.
            Ideal weights used in the embedding.

        Returns
        -------
        output : tuple.
            A tuple with two elements:

            embedding with the k-dimensional embedding of the input samples

            weights with the weights associated to each input data view.

        Raises
        ------

        ValueError: Matrices are not square matrices, k value is negative,
        data samples and is_distance parameters do not have the same length or
        scalar parameters are negative.

        Examples
        --------

        >>> import numpy as np
        >>> m = np.array([[1, 4, 7], [2, 5, 8], [3, 6, 9]])
        >>> q = np.array([[9, 6, 3], [8, 5, 2], [7, 4, 1]])
        >>> r = np.array([[2, 1, 8], [4, 5, 6], [3, 7, 9]]).T
        >>> matrices = [m, q, r]
        >>> is_distance = [False, False, False]
        >>> mvtsne = MvtSNE()
        >>> mvtsne.fit_transform(matrices, is_distance)
            (matrix([[-1347.89641563,  -415.25549328],
                     [ 1305.18939063,   398.91164491],
                     [   42.70702501,    16.34384836]]),
                     array([ 0.878037  ,  0.64703391,  0.56962457]))

        """

        if len(X) != len(is_distance):
            raise ValueError("Data samples and is_distance lengths does not "
                             "match. Data sample length: %d, is_distance "
                             "length: %d" % (len(X), len(is_distance)))

        if self.k > X[0].shape[0]:
            self.k = X[0].shape[0]
            warnings.warn("k is greater than matrix dimension. k=%d is "
                          "computed instead." % X[0].shape[0])
        elif self.k < 0:
            raise ValueError("k value must be between 0 and number of samples"
                             " of data matrix.")

        if self.initial_dims < 0 or self.min_cost < 0:
            if self.perplexity < 0 or self.max_iter < 0:
                raise ValueError("Parameters cannot be negative.")

        for i in np.arange(len(X) - 1):
            for j in np.arange(i + 1, len(X)):
                if X[i].shape[0] != X[j].shape[0]:
                    raise ValueError("Input data matrices have no same number "
                                     "of samples (rows).")

        self.random_state = check_random_state(self.random_state)

        pooling_p, weights, kl_ps, kl_central = find_pooling(
            X, is_distance=is_distance, random_state=self.random_state,
            initial_dims=self.initial_dims, perplexity=self.perplexity,
            whiten=self.whiten)
        embedding = tsne_p(pooling_p, k=self.k, initial_dims=self.initial_dims,
                           max_iter=self.max_iter, min_cost=self.min_cost,
                           epoch_callback=self.epoch_callback,
                           epoch=self.epoch)
        self.embedding_ = embedding
        self.weights_ = weights
        return (embedding, weights)

    def get_params(self, deep=True):
        return {"k": self.k, "initial_dims": self.initial_dims,
                "perplexity": self.perplexity, "max_iter": self.max_iter,
                "min_cost": self.min_cost,
                "epoch_callback": self.epoch_callback,
                "whiten": self.whiten, "epoch": self.epoch,
                "random_state": self.random_state}

    def set_params(self, **parameters):
        for parameter, value in parameters.items():
            setattr(self, parameter, value)
        return self
