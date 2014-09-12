"""

=============================================================
Online Latent Dirichlet Allocation with variational inference
=============================================================

This implementation is modified from Matthew D. Hoffman's onlineldavb code
Link: http://www.cs.princeton.edu/~mdhoffma/code/onlineldavb.tar
"""

# Author: Chyi-Kwei Yau
# Author: Matthew D. Hoffman (original onlineldavb implementation)

import numpy as np
import scipy.sparse as sp
from scipy.special import gammaln, psi

from ..base import BaseEstimator, TransformerMixin
from ..utils import (check_random_state, check_array,
                           gen_batches, gen_even_slices)
from ..externals.joblib import Parallel, delayed, cpu_count
from ..externals.six.moves import xrange


def _dirichlet_expectation(alpha):
    """
    For a vector theta ~ Dir(alpha), computes E[log(theta)] given alpha.
    """
    if (len(alpha.shape) == 1):
        return(psi(alpha) - psi(np.sum(alpha)))
    return(psi(alpha) - psi(np.sum(alpha, 1))[:, np.newaxis])


def _update_gamma(X, expElogbeta, alpha, rng, max_iters,
                  mean_change_tol, cal_delta):
    """
    E-step: update latent variable gamma
    """

    n_docs, n_vocabs = X.shape
    n_topics = expElogbeta.shape[0]

    # gamma is non-normailzed topic distribution
    gamma = rng.gamma(100., 1. / 100., (n_docs, n_topics))
    expElogtheta = np.exp(_dirichlet_expectation(gamma))
    # diff on component (only calculate it when keep_comp_change is True)
    delta_component = np.zeros(expElogbeta.shape) if cal_delta else None

    X_data = X.data
    X_indices = X.indices
    X_indptr = X.indptr

    for d in xrange(n_docs):
        ids = X_indices[X_indptr[d]:X_indptr[d + 1]]
        cnts = X_data[X_indptr[d]:X_indptr[d + 1]]
        gammad = gamma[d, :]
        expElogthetad = expElogtheta[d, :]
        expElogbetad = expElogbeta[:, ids]
        # The optimal phi_{dwk} is proportional to
        # expElogthetad_k * expElogbetad_w. phinorm is the normalizer.
        phinorm = np.dot(expElogthetad, expElogbetad) + 1e-100

        # Iterate between gamma and phi until convergence
        for it in xrange(0, max_iters):
            lastgamma = gammad
            # We represent phi implicitly to save memory and time.
            # Substituting the value of the optimal phi back into
            # the update for gamma gives this update. Cf. Lee&Seung 2001.
            gammad = alpha + expElogthetad * \
                np.dot(cnts / phinorm, expElogbetad.T)
            expElogthetad = np.exp(_dirichlet_expectation(gammad))
            phinorm = np.dot(expElogthetad, expElogbetad) + 1e-100

            meanchange = np.mean(abs(gammad - lastgamma))
            if (meanchange < mean_change_tol):
                break
        gamma[d, :] = gammad
        # Contribution of document d to the expected sufficient
        # statistics for the M step.
        if cal_delta:
            delta_component[:, ids] += np.outer(expElogthetad, cnts / phinorm)

    return (gamma, delta_component)


class OnlineLDA(BaseEstimator, TransformerMixin):

    """
    Online Latent Dirichlet Allocation implementation with variational inference

    References
    ----------
    [1] "Online Learning for Latent Dirichlet Allocation", Matthew D. Hoffman, 
        David M. Blei, Francis Bach

    [2] Matthew D. Hoffman's onlineldavb code. Link:
        http://www.cs.princeton.edu/~mdhoffma/code/onlineldavb.tar


    Parameters
    ----------
    n_topics: int, optional (default: 10)
        number of topics.

    alpha: float, optional (defalut: 0.1)
        Hyperparameter for prior on weight vectors theta. In general, it is `1 / n_topics`.

    eta: float, optional (default: 0.1)
        Hyperparameter for prior on topics beta. In general, it is `1 / n_topics`.

    kappa: float, optional (default: 0.7)
        weight for _component: it controls the weight for previous value of _component.
        The value hould be between (0.5, 1.0] to guarantee asymptotic convergence for online learning.
        It is only used in partial_fit.
        When we set kappa to 0.0 and batch_size to n_docs, then the udpate is batch VB.

    tau: float, optional (default: 1024.)
        A (positive) learning parameter that downweights early iterations. It should be greater than 1.0.
        It is only used in partial_fit

    n_docs: int, optional (default: 1e6)
        Total umber of document. It is only used in online learing(partial_fit).
        In batch learning, n_docs is X.shape[0]

    batch_size: int, optional (default: 128)
        Number of document to udpate in each EM-step

    normalize_doc: boolean, optional (default: True)
        normalize the topic distribution for transformed document or not.
        if True, sum of topic distribution for each document will be 1.0

    e_step_tol: float, optional (default: 1e-4)
        Tolerance value used in e-step break conditions.

    prex_tol: float, optional (default: 1e-1)
        Tolerance value used in preplexity break conditions.

    mean_change_tol: float, optional (default: 1e-3)
        Tolerance value used in e-step break conditions.

    n_jobs: int, optional (default: 1)
        Number of parallel jobs to run on e-step. -1 for autodetect.

    verbose : int, optional (default: 0)
        Verbosity level.

    random_state: int or RandomState instance or None, optional (default: None)
        Pseudo Random Number generator seed control.


    Attributes
    ----------
    components_: array, [n_topics, n_vocabs]
        vocab distribution for each topic. components_[i, j] represents
        vocab `j` in topic `i`

    n_iter_: int
        number of iteration

    """

    def __init__(self, n_topics=10, alpha=.1, eta=.1, kappa=.7, tau=1000.,
                 batch_size=128, n_docs=1e6, normalize_doc=True,
                 e_step_tol=1e-3, prex_tol=1e-1, mean_change_tol=1e-3,
                 n_jobs=1, verbose=0, random_state=None):
        self.n_topics = n_topics
        self.alpha = alpha
        self.eta = eta
        self.kappa = kappa
        self.tau = tau
        self.batch_size = batch_size
        self.n_docs = n_docs
        self.normalize_doc = normalize_doc
        self.e_step_tol = e_step_tol
        self.prex_tol = prex_tol
        self.mean_change_tol = mean_change_tol
        self.n_jobs = n_jobs
        self.verbose = verbose
        self.random_state = random_state

    def _init_latent_vars(self, n_vocabs):
        self.rng = check_random_state(self.random_state)
        self.n_iter_ = 1
        self.n_vocabs = n_vocabs
        init_gamma = 100.
        init_var = 1. / init_gamma
        self.components_ = self.rng.gamma(
            init_gamma, init_var, (self.n_topics, n_vocabs))

        self.Elogbeta = _dirichlet_expectation(self.components_)
        self.expElogbeta = np.exp(self.Elogbeta)

    def _e_step(self, X, cal_delta):
        """
        E-step

        set `cal_delta == True` when we need to run _m_step
        for inference, set it to False
        """

        # parell run e-step
        if self.n_jobs == -1:
            n_jobs = cpu_count()
        else:
            n_jobs = self.n_jobs

        results = Parallel(n_jobs=n_jobs, verbose=self.verbose)(
            delayed(_update_gamma)
            (X[idx_slice, :], self.expElogbeta, self.alpha,
             self.rng, 100, self.mean_change_tol, cal_delta)
            for idx_slice in gen_even_slices(X.shape[0], n_jobs))

        # merge result
        gammas, deltas = zip(*results)
        gamma = np.vstack(gammas)

        if cal_delta:
            # This step finishes computing the sufficient statistics for the
            # M step, so that
            # sstats[k, w] = \sum_d n_{dw} * phi_{dwk}
            # = \sum_d n_{dw} * exp{Elogtheta_{dk} + Elogbeta_{kw}} / phinorm_{dw}.
            delta_component = np.zeros(self.components_.shape)
            for delta in deltas:
                delta_component += delta
            delta_component *= self.expElogbeta
        else:
            delta_component = None

        return (gamma, delta_component)

    def _em_step(self, X, batch_update):
        """
        EM update for 1 iteration

        parameters
        ----------
        X:  sparse matrix

        batch_update: boolean
        update `_component` by bath VB or online VB
        """
        # e-step
        gamma, delta_component = self._e_step(X, cal_delta=True)

        # m-step
        if batch_update:
            self.components_ = self.eta + delta_component
        else:
            # online update
            rhot = np.power(self.tau + self.n_iter_, -self.kappa)
            doc_ratio = float(self.n_docs) / X.shape[0]
            self.components_ *= (1 - rhot)
            self.components_ += (rhot *
                                 (self.eta + doc_ratio * delta_component))

        self.Elogbeta = _dirichlet_expectation(self.components_)
        self.expElogbeta = np.exp(self.Elogbeta)
        self.n_iter_ += 1

        return gamma

    def _to_csr(self, X):
        """
        check & convert X to csr format
        """
        X = check_array(X, accept_sparse='csr')
        if not sp.issparse(X):
            X = sp.csr_matrix(X)

        return X

    def fit_transform(self, X, y=None, max_iters=10):
        """
        Learn a model for X and returns the transformed data

        Parameters
        ----------
        X: array or sparse matrix, shape = [n_docs, n_vocabs]
            Data matrix to be transformed by the model

        max_iters: int, (default: 10)
            Max number of iterations

        Returns
        -------
        gamma: array, [n_dics, n_topics]
            Topic distribution for each doc
        """

        """
        X = self._to_csr(X)
        n_docs, n_vocabs = X.shape

        if not hasattr(self, 'components_'):
            # initialize vocabulary & latent variables

            self.Elogbeta = _dirichlet_expectation(self.components_)
            self.expElogbeta = np.exp(self.Elogbeta)
        else:
            # make sure vacabulary size matched
            if self.components_.shape[1] != X.shape[1]:
                raise ValueError("dimension not match")

        # EM update
        gamma, delta_component = self._e_step(X)
        # self._m_step(delta_component, self.n_docs, n_docs)
        self._m_step(delta_component, self.n_docs, n_docs)

        if self.normalize_doc:
            gamma /= gamma.sum(axis=1)[:, np.newaxis]

        return gamma
        """
        X = self._to_csr(X)
        return self.fit(X, max_iters).transform(X)

    def partial_fit(self, X, y=None):
        """
        Online Learning with Min-Batch update

        Parameters
        ----------
        X: sparse matrix, shape = [n_docs, n_vocabs]
            Data matrix to be decomposed

        Returns
        -------
        self
        """

        X = self._to_csr(X)
        n_docs, n_vocabs = X.shape
        batch_size = self.batch_size

        # initialize parameters or check
        if not hasattr(self, 'components_'):
            self._init_latent_vars(n_vocabs)

        if n_vocabs != self.n_vocabs:
            raise ValueError(
                "feature dimension(vocabulary size) doesn't match.")

        for idx_slice in gen_batches(n_docs, batch_size):
            self._em_step(X[idx_slice, :], batch_update=False)

        return self

    def fit(self, X, y=None, max_iters=10):
        """
        Learn model from X. This function is for batch learning.
        So it will override the old _component variables

        Parameters
        ----------
        X: sparse matrix, shape = (n_docs, n_vocabs)
            Data matrix to be transformed by the model

        max_iters: int, (default: 10)
            Max number of iterations

        Returns
        -------
        self
        """

        X = self._to_csr(X)
        n_docs, n_vocabs = X.shape

        # initialize parameters
        self._init_latent_vars(n_vocabs)

        # change to preplexity later
        last_bound = None
        for i in xrange(max_iters):
            gamma = self._em_step(X, batch_update=True)

            # check preplexity
            bound = self.preplexity(X, gamma, sub_sampling=False)
            if self.verbose:
                print('iteration: %d, preplexity: %.4f' % (i, bound))

            if i > 0 and abs(last_bound - bound) < self.prex_tol:
                break
            last_bound = bound

        return self

    def transform(self, X):
        """
        Transform the data X according to the fitted model (run inference)

        Parameters
        ----------
        X: sparse matrix, shape = [n_docs, n_vocabs]
            Data matrix to be transformed by the model
            n_vacabs must be the same as n_vocabs in fitted model

        max_iters: int, (default: 20)
            Max number of iterations

        Returns
        -------
        data: array, [n_docs, n_topics]
            Document distribution
        """

        X = self._to_csr(X)
        n_docs, n_vocabs = X.shape

        if not hasattr(self, 'components_'):
            raise AttributeError(
                "no 'components_' attr in model. Please fit model first.")
        # make sure vocabulary size is the same in fitted model and new doc
        # matrix
        if n_vocabs != self.n_vocabs:
            raise ValueError(
                "feature dimension(vocabulary size) does not match.")

        gamma, _ = self._e_step(X, False)

        if self.normalize_doc:
            gamma /= gamma.sum(axis=1)[:, np.newaxis]

        return gamma

    def _approx_bound(self, X, gamma, sub_sampling):
        """
        calculate approximate bound for data X and topic distribution gamma


        Parameters
        ----------
        X: sparse matrix, [n_docs, n_vocabs]

        gamma: array, shape = [n_docs, n_topics]
            document distribution (can be either normalized & un-normalized)

        sub_sampling: boolean, optional, (default: False)
            Compensate for the subsampling of the population of documents
            set subsampling to `True` for online learning

        Returns
        -------
        score: float, score of gamma
        """
        X = self._to_csr(X)
        n_docs, n_topics = gamma.shape
        score = 0
        Elogtheta = _dirichlet_expectation(gamma)

        X_data = X.data
        X_indices = X.indices
        X_indptr = X.indptr

        # E[log p(docs | theta, beta)]
        for d in xrange(0, n_docs):
            ids = X_indices[X_indptr[d]:X_indptr[d + 1]]
            cnts = X_data[X_indptr[d]:X_indptr[d + 1]]
            phinorm = np.zeros(len(ids))
            for i in xrange(0, len(ids)):
                temp = Elogtheta[d, :] + self.Elogbeta[:, ids[i]]
                tmax = max(temp)
                phinorm[i] = np.log(sum(np.exp(temp - tmax))) + tmax
            score += np.sum(cnts * phinorm)

        # E[log p(theta | alpha) - log q(theta | gamma)]
        score += np.sum((self.alpha - gamma) * Elogtheta)
        score += np.sum(gammaln(gamma) - gammaln(self.alpha))
        score += sum(
            gammaln(self.alpha * self.n_topics) - gammaln(np.sum(gamma, 1)))

        # Compensate for the subsampling of the population of documents
        # E[log p(beta | eta) - log q (beta | lambda)]
        score += np.sum((self.eta - self.components_) * self.Elogbeta)
        score += np.sum(gammaln(self.components_) - gammaln(self.eta))
        score += np.sum(gammaln(self.eta * self.n_vocabs)
                        - gammaln(np.sum(self.components_, 1)))

        # Compensate for the subsampling of the population of documents
        if sub_sampling:
            doc_ratio = float(self.n_docs) / n_docs
            score *= doc_ratio

        return score

    def preplexity(self, X, gamma, sub_sampling=False):
        """
        calculate approximate bound for data X and topic distribution gamma


        Parameters
        ----------
        X: sparse matrix, [n_docs, n_vocabs]

        gamma: array, shape = [n_docs, n_topics]
            document distribution (can be either normalized & un-normalized)

        Returns
        -------
        score: float, score of gamma
        """
        X = self._to_csr(X)
        n_doc = X.shape[0]
        bound = self._approx_bound(X, gamma, sub_sampling)
        perword_bound = bound / np.sum(X.data)

        if sub_sampling:
            perword_bound = perword_bound * (float(n_doc) / self.n_docs)

        return np.exp(-1.0 * perword_bound)
