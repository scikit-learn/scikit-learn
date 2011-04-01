# Hidden Markov Models
#
# Author: Ron Weiss <ronweiss@gmail.com>

import string

import numpy as np

from .base import BaseEstimator
from .mixture import (GMM, lmvnpdf, logsum, normalize, sample_gaussian,
                 _distribute_covar_matrix_to_match_cvtype, _validate_covars)
from . import cluster
ZEROLOGPROB = -1e200


class _BaseHMM(BaseEstimator):
    """Hidden Markov Model base class.

    Representation of a hidden Markov model probability distribution.
    This class allows for easy evaluation of, sampling from, and
    maximum-likelihood estimation of the parameters of a HMM.

    See the instance documentation for details specific to a
    particular object.

    Attributes
    ----------
    n_states : int (read-only)
        Number of states in the model.
    transmat : array, shape (`n_states`, `n_states`)
        Matrix of transition probabilities between states.
    startprob : array, shape ('n_states`,)
        Initial state occupation distribution.

    Methods
    -------
    eval(X)
        Compute the log likelihood of `X` under the HMM.
    decode(X)
        Find most likely state sequence for each point in `X` using the
        Viterbi algorithm.
    rvs(n=1)
        Generate `n` samples from the HMM.
    fit(X)
        Estimate HMM parameters from `X`.
    predict(X)
        Like decode, find most likely state sequence corresponding to `X`.
    score(X)
        Compute the log likelihood of `X` under the model.

    See Also
    --------
    GMM : Gaussian mixture model
    """

    # This class implements the public interface to all HMMs that
    # derive from it, including all of the machinery for the
    # forward-backward and Viterbi algorithms.  Subclasses need only
    # implement _generate_sample_from_state(), _compute_log_likelihood(),
    # _init(), _initialize_sufficient_statistics(),
    # _accumulate_sufficient_statistics(), and _do_mstep(), all of
    # which depend on the specific emission distribution.
    #
    # Subclasses will probably also want to implement properties for
    # the emission distribution parameters to expose them publically.

    def __init__(self, n_states=1, startprob=None, transmat=None,
                 startprob_prior=None, transmat_prior=None):
        self._n_states = n_states

        if startprob is None:
            startprob = np.tile(1.0 / n_states, n_states)
        self.startprob = startprob

        if startprob_prior is None:
            startprob_prior = 1.0
        self.startprob_prior = startprob_prior

        if transmat is None:
            transmat = np.tile(1.0 / n_states, (n_states, n_states))
        self.transmat = transmat

        if transmat_prior is None:
            transmat_prior = 1.0
        self.transmat_prior = transmat_prior

    def eval(self, obs, maxrank=None, beamlogprob=-np.Inf):
        """Compute the log probability under the model and compute posteriors

        Implements rank and beam pruning in the forward-backward
        algorithm to speed up inference in large models.

        Parameters
        ----------
        obs : array_like, shape (n, n_features)
            Sequence of n_features-dimensional data points.  Each row
            corresponds to a single point in the sequence.
        maxrank : int
            Maximum rank to evaluate for rank pruning.  If not None,
            only consider the top `maxrank` states in the inner
            sum of the forward algorithm recursion.  Defaults to None
            (no rank pruning).  See The HTK Book for more details.
        beamlogprob : float
            Width of the beam-pruning beam in log-probability units.
            Defaults to -numpy.Inf (no beam pruning).  See The HTK
            Book for more details.

        Returns
        -------
        logprob : array_like, shape (n,)
            Log probabilities of the sequence `obs`
        posteriors: array_like, shape (n, n_states)
            Posterior probabilities of each state for each
            observation

        See Also
        --------
        score : Compute the log probability under the model
        decode : Find most likely state sequence corresponding to a `obs`
        """
        obs = np.asanyarray(obs)
        framelogprob = self._compute_log_likelihood(obs)
        logprob, fwdlattice = self._do_forward_pass(framelogprob, maxrank,
                                                    beamlogprob)
        bwdlattice = self._do_backward_pass(framelogprob, fwdlattice, maxrank,
                                            beamlogprob)
        gamma = fwdlattice + bwdlattice
        # gamma is guaranteed to be correctly normalized by logprob at
        # all frames, unless we do approximate inference using pruning.
        # So, we will normalize each frame explicitly in case we
        # pruned too aggressively.
        posteriors = np.exp(gamma.T - logsum(gamma, axis=1)).T
        return logprob, posteriors

    def score(self, obs, maxrank=None, beamlogprob=-np.Inf):
        """Compute the log probability under the model.

        Parameters
        ----------
        obs : array_like, shape (n, n_features)
            Sequence of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        maxrank : int
            Maximum rank to evaluate for rank pruning.  If not None,
            only consider the top `maxrank` states in the inner
            sum of the forward algorithm recursion.  Defaults to None
            (no rank pruning).  See The HTK Book for more details.
        beamlogprob : float
            Width of the beam-pruning beam in log-probability units.
            Defaults to -numpy.Inf (no beam pruning).  See The HTK
            Book for more details.

        Returns
        -------
        logprob : array_like, shape (n,)
            Log probabilities of each data point in `obs`

        See Also
        --------
        eval : Compute the log probability under the model and posteriors
        decode : Find most likely state sequence corresponding to a `obs`
        """
        obs = np.asanyarray(obs)
        framelogprob = self._compute_log_likelihood(obs)
        logprob, fwdlattice = self._do_forward_pass(framelogprob, maxrank,
                                                    beamlogprob)
        return logprob

    def decode(self, obs, maxrank=None, beamlogprob=-np.Inf):
        """Find most likely state sequence corresponding to `obs`.

        Uses the Viterbi algorithm.

        Parameters
        ----------
        obs : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        maxrank : int
            Maximum rank to evaluate for rank pruning.  If not None,
            only consider the top `maxrank` states in the inner
            sum of the forward algorithm recursion.  Defaults to None
            (no rank pruning).  See The HTK Book for more details.
        beamlogprob : float
            Width of the beam-pruning beam in log-probability units.
            Defaults to -numpy.Inf (no beam pruning).  See The HTK
            Book for more details.

        Returns
        -------
        viterbi_logprob : float
            Log probability of the maximum likelihood path through the HMM
        states : array_like, shape (n,)
            Index of the most likely states for each observation

        See Also
        --------
        eval : Compute the log probability under the model and posteriors
        score : Compute the log probability under the model
        """
        obs = np.asanyarray(obs)
        framelogprob = self._compute_log_likelihood(obs)
        logprob, state_sequence = self._do_viterbi_pass(framelogprob, maxrank,
                                                        beamlogprob)
        return logprob, state_sequence

    def predict(self, obs, **kwargs):
        """Find most likely state sequence corresponding to `obs`.

        Parameters
        ----------
        obs : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        maxrank : int
            Maximum rank to evaluate for rank pruning.  If not None,
            only consider the top `maxrank` states in the inner
            sum of the forward algorithm recursion.  Defaults to None
            (no rank pruning).  See The HTK Book for more details.
        beamlogprob : float
            Width of the beam-pruning beam in log-probability units.
            Defaults to -numpy.Inf (no beam pruning).  See The HTK
            Book for more details.

        Returns
        -------
        states : array_like, shape (n,)
            Index of the most likely states for each observation
        """
        logprob, state_sequence = self.decode(obs, **kwargs)
        return state_sequence

    def predict_proba(self, obs, **kwargs):
        """Compute the posterior probability for each state in the model

        Parameters
        ----------
        obs : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        See eval() for a list of accepted keyword arguments.

        Returns
        -------
        T : array-like, shape (n, n_states)
            Returns the probability of the sample for each state in the model.
        """
        logprob, posteriors = self.eval(obs, **kwargs)
        return posteriors

    def rvs(self, n=1):
        """Generate random samples from the model.

        Parameters
        ----------
        n : int
            Number of samples to generate.

        Returns
        -------
        obs : array_like, length `n`
            List of samples
        """

        startprob_pdf = self.startprob
        startprob_cdf = np.cumsum(startprob_pdf)
        transmat_pdf = self.transmat
        transmat_cdf = np.cumsum(transmat_pdf, 1)

        # Initial state.
        rand = np.random.rand()
        currstate = (startprob_cdf > rand).argmax()
        obs = [self._generate_sample_from_state(currstate)]

        for x in xrange(n - 1):
            rand = np.random.rand()
            currstate = (transmat_cdf[currstate] > rand).argmax()
            obs.append(self._generate_sample_from_state(currstate))

        return np.array(obs)

    def fit(self, obs, n_iter=10, thresh=1e-2, params=string.letters,
            init_params=string.letters,
            maxrank=None, beamlogprob=-np.Inf, **kwargs):
        """Estimate model parameters.

        An initialization step is performed before entering the EM
        algorithm. If you want to avoid this step, set the keyword
        argument init_params to the empty string ''. Likewise, if you
        would like just to do an initialization, call this method with
        n_iter=0.

        Parameters
        ----------
        obs : list
            List of array-like observation sequences (shape (n_i, n_features)).
        n_iter : int, optional
            Number of iterations to perform.
        thresh : float, optional
            Convergence threshold.
        params : string, optional
            Controls which parameters are updated in the training
            process.  Can contain any combination of 's' for startprob,
            't' for transmat, 'm' for means, and 'c' for covars, etc.
            Defaults to all parameters.
        init_params : string, optional
            Controls which parameters are initialized prior to
            training.  Can contain any combination of 's' for
            startprob, 't' for transmat, 'm' for means, and 'c' for
            covars, etc.  Defaults to all parameters.
        maxrank : int, optional
            Maximum rank to evaluate for rank pruning.  If not None,
            only consider the top `maxrank` states in the inner
            sum of the forward algorithm recursion.  Defaults to None
            (no rank pruning).  See "The HTK Book" for more details.
        beamlogprob : float, optional
            Width of the beam-pruning beam in log-probability units.
            Defaults to -numpy.Inf (no beam pruning).  See "The HTK
            Book" for more details.

        Notes
        -----
        In general, `logprob` should be non-decreasing unless
        aggressive pruning is used.  Decreasing `logprob` is generally
        a sign of overfitting (e.g. a covariance parameter getting too
        small).  You can fix this by getting more training data, or
        decreasing `covars_prior`.
        """
        self._init(obs, init_params)

        logprob = []
        for i in xrange(n_iter):
            # Expectation step
            stats = self._initialize_sufficient_statistics()
            curr_logprob = 0
            for seq in obs:
                framelogprob = self._compute_log_likelihood(seq)
                lpr, fwdlattice = self._do_forward_pass(framelogprob, maxrank,
                                                       beamlogprob)
                bwdlattice = self._do_backward_pass(framelogprob, fwdlattice,
                                                   maxrank, beamlogprob)
                gamma = fwdlattice + bwdlattice
                posteriors = np.exp(gamma.T - logsum(gamma, axis=1)).T
                curr_logprob += lpr
                self._accumulate_sufficient_statistics(
                    stats, seq, framelogprob, posteriors, fwdlattice,
                    bwdlattice, params)
            logprob.append(curr_logprob)

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < thresh:
                break

            # Maximization step
            self._do_mstep(stats, params, **kwargs)

        return self

    @property
    def n_states(self):
        """Number of states in the model."""
        return self._n_states

    def _get_startprob(self):
        """Mixing startprob for each state."""
        return np.exp(self._log_startprob)

    def _set_startprob(self, startprob):
        if len(startprob) != self._n_states:
            raise ValueError('startprob must have length n_states')
        if not np.allclose(np.sum(startprob), 1.0):
            raise ValueError('startprob must sum to 1.0')

        self._log_startprob = np.log(np.asanyarray(startprob).copy())

    startprob = property(_get_startprob, _set_startprob)

    def _get_transmat(self):
        """Matrix of transition probabilities."""
        return np.exp(self._log_transmat)

    def _set_transmat(self, transmat):
        if np.asanyarray(transmat).shape != (self._n_states, self._n_states):
            raise ValueError('transmat must have shape (n_states, n_states)')
        if not np.all(np.allclose(np.sum(transmat, axis=1), 1.0)):
            raise ValueError('Rows of transmat must sum to 1.0')

        self._log_transmat = np.log(np.asanyarray(transmat).copy())
        underflow_idx = np.isnan(self._log_transmat)
        self._log_transmat[underflow_idx] = -np.Inf

    transmat = property(_get_transmat, _set_transmat)

    def _do_viterbi_pass(self, framelogprob, maxrank=None,
                         beamlogprob=-np.Inf):
        nobs = len(framelogprob)
        lattice = np.zeros((nobs, self._n_states))
        traceback = np.zeros((nobs, self._n_states), dtype=np.int)

        lattice[0] = self._log_startprob + framelogprob[0]
        for n in xrange(1, nobs):
            idx = self._prune_states(lattice[n - 1], maxrank, beamlogprob)
            pr = self._log_transmat[idx].T + lattice[n - 1, idx]
            lattice[n] = np.max(pr, axis=1) + framelogprob[n]
            traceback[n] = np.argmax(pr, axis=1)
        lattice[lattice <= ZEROLOGPROB] = -np.Inf

        # Do traceback.
        reverse_state_sequence = []
        s = lattice[-1].argmax()
        logprob = lattice[-1, s]
        for frame in reversed(traceback):
            reverse_state_sequence.append(s)
            s = frame[s]

        reverse_state_sequence.reverse()
        return logprob, np.array(reverse_state_sequence)

    def _do_forward_pass(self, framelogprob, maxrank=None,
                         beamlogprob=-np.Inf):
        nobs = len(framelogprob)
        fwdlattice = np.zeros((nobs, self._n_states))

        fwdlattice[0] = self._log_startprob + framelogprob[0]
        for n in xrange(1, nobs):
            idx = self._prune_states(fwdlattice[n - 1], maxrank, beamlogprob)
            fwdlattice[n] = (logsum(self._log_transmat[idx].T
                                    + fwdlattice[n - 1, idx], axis=1)
                             + framelogprob[n])
        fwdlattice[fwdlattice <= ZEROLOGPROB] = -np.Inf

        return logsum(fwdlattice[-1]), fwdlattice

    def _do_backward_pass(self, framelogprob, fwdlattice, maxrank=None,
                          beamlogprob=-np.Inf):
        nobs = len(framelogprob)
        bwdlattice = np.zeros((nobs, self._n_states))

        for n in xrange(nobs - 1, 0, -1):
            # Do HTK style pruning (p. 137 of HTK Book version 3.4).
            # Don't bother computing backward probability if
            # fwdlattice * bwdlattice is more than a certain distance
            # from the total log likelihood.
            idx = self._prune_states(bwdlattice[n] + fwdlattice[n], None,
                                     -50)
                                     #beamlogprob)
                                     #-np.Inf)
            bwdlattice[n - 1] = logsum(self._log_transmat[:, idx] +
                                       bwdlattice[n, idx] +
                                       framelogprob[n, idx],
                                       axis=1)
        bwdlattice[bwdlattice <= ZEROLOGPROB] = -np.Inf

        return bwdlattice

    def _prune_states(self, lattice_frame, maxrank, beamlogprob):
        """ Returns indices of the active states in `lattice_frame`
        after rank and beam pruning.
        """
        # Beam pruning
        threshlogprob = logsum(lattice_frame) + beamlogprob

        # Rank pruning
        if maxrank:
            # How big should our rank pruning histogram be?
            nbins = 3 * len(lattice_frame)

            lattice_min = lattice_frame[lattice_frame > ZEROLOGPROB].min() - 1
            hst, cdf = np.histogram(lattice_frame, bins=nbins,
                                    range=(lattice_min, lattice_frame.max()))

            # Want to look at the high ranks.
            hst = hst[::-1].cumsum()
            cdf = cdf[::-1]

            rankthresh = cdf[hst >= min(maxrank, self._n_states)].max()

            # Only change the threshold if it is stricter than the beam
            # threshold.
            threshlogprob = max(threshlogprob, rankthresh)

        # Which states are active?
        state_idx, = np.nonzero(lattice_frame >= threshlogprob)
        return state_idx

    def _compute_log_likelihood(self, obs):
        pass

    def _generate_sample_from_state(self, state):
        pass

    def _init(self, obs, params):
        if 's' in params:
            self.startprob[:] = 1.0 / self._n_states
        if 't' in params:
            self.transmat[:] = 1.0 / self._n_states

    # Methods used by self.fit()

    def _initialize_sufficient_statistics(self):
        stats = {'nobs': 0,
                 'start': np.zeros(self._n_states),
                 'trans': np.zeros((self._n_states, self._n_states))}
        return stats

    def _accumulate_sufficient_statistics(self, stats, seq, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        stats['nobs'] += 1
        if 's' in params:
            stats['start'] += posteriors[0]
        if 't' in params:
            for t in xrange(len(framelogprob)):
                zeta = (fwdlattice[t - 1][:, np.newaxis] + self._log_transmat
                        + framelogprob[t] + bwdlattice[t])
                stats['trans'] += np.exp(zeta - logsum(zeta))

    def _do_mstep(self, stats, params, **kwargs):
        # Based on Huang, Acero, Hon, "Spoken Language Processing",
        # p. 443 - 445
        if 's' in params:
            self.startprob = normalize(
                np.maximum(self.startprob_prior - 1.0 + stats['start'], 1e-20))
        if 't' in params:
            self.transmat = normalize(
                np.maximum(self.transmat_prior - 1.0 + stats['trans'], 1e-20),
                axis=1)


class GaussianHMM(_BaseHMM):
    """Hidden Markov Model with Gaussian emissions

    Representation of a hidden Markov model probability distribution.
    This class allows for easy evaluation of, sampling from, and
    maximum-likelihood estimation of the parameters of a HMM.

    Attributes
    ----------
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the model.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    n_features : int (read-only)
        Dimensionality of the Gaussian emissions.
    n_states : int (read-only)
        Number of states in the model.
    transmat : array, shape (`n_states`, `n_states`)
        Matrix of transition probabilities between states.
    startprob : array, shape ('n_states`,)
        Initial state occupation distribution.
    means : array, shape (`n_states`, `n_features`)
        Mean parameters for each state.
    covars : array
        Covariance parameters for each state.  The shape depends on
        `cvtype`:
            (`n_states`,)                   if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_states`, `n_features`)           if 'diag',
            (`n_states`, `n_features`, `n_features`)  if 'full'

    Methods
    -------
    eval(X)
        Compute the log likelihood of `X` under the HMM.
    decode(X)
        Find most likely state sequence for each point in `X` using the
        Viterbi algorithm.
    rvs(n=1)
        Generate `n` samples from the HMM.
    init(X)
        Initialize HMM parameters from `X`.
    fit(X)
        Estimate HMM parameters from `X` using the Baum-Welch algorithm.
    predict(X)
        Like decode, find most likely state sequence corresponding to `X`.
    score(X)
        Compute the log likelihood of `X` under the model.

    Examples
    --------
    >>> from scikits.learn.hmm import GaussianHMM
    >>> GaussianHMM(n_states=2)
    GaussianHMM(cvtype='diag', n_states=2, means_weight=0, startprob_prior=1.0,
          startprob=array([ 0.5,  0.5]),
          transmat=array([[ 0.5,  0.5],
           [ 0.5,  0.5]]),
          transmat_prior=1.0, means_prior=None, covars_weight=1,
          covars_prior=0.01)

    See Also
    --------
    GMM : Gaussian mixture model
    """

    def __init__(self, n_states=1, cvtype='diag', startprob=None,
                 transmat=None, startprob_prior=None, transmat_prior=None,
                 means_prior=None, means_weight=0,
                 covars_prior=1e-2, covars_weight=1):
        """Create a hidden Markov model with Gaussian emissions.

        Initializes parameters such that every state has zero mean and
        identity covariance.

        Parameters
        ----------
        n_states : int
            Number of states.
        cvtype : string
            String describing the type of covariance parameters to
            use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
            Defaults to 'diag'.
        """
        super(GaussianHMM, self).__init__(n_states, startprob, transmat,
                                          startprob_prior=startprob_prior,
                                          transmat_prior=transmat_prior)

        self._cvtype = cvtype
        if not cvtype in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError('bad cvtype')

        self.means_prior = means_prior
        self.means_weight = means_weight

        self.covars_prior = covars_prior
        self.covars_weight = covars_weight

    # Read-only properties.
    @property
    def cvtype(self):
        """Covariance type of the model.

        Must be one of 'spherical', 'tied', 'diag', 'full'.
        """
        return self._cvtype

    def _get_means(self):
        """Mean parameters for each state."""
        return self._means

    def _set_means(self, means):
        means = np.asanyarray(means)
        if hasattr(self, 'n_features') and \
               means.shape != (self._n_states, self.n_features):
            raise ValueError('means must have shape (n_states, n_features)')
        self._means = means.copy()
        self.n_features = self._means.shape[1]

    means = property(_get_means, _set_means)

    def _get_covars(self):
        """Return covars as a full matrix."""
        if self.cvtype == 'full':
            return self._covars
        elif self.cvtype == 'diag':
            return [np.diag(cov) for cov in self._covars]
        elif self.cvtype == 'tied':
            return [self._covars] * self._n_states
        elif self.cvtype == 'spherical':
            return [np.eye(self.n_features) * f for f in self._covars]

    def _set_covars(self, covars):
        covars = np.asanyarray(covars)
        _validate_covars(covars, self._cvtype, self._n_states, self.n_features)
        self._covars = covars.copy()

    covars = property(_get_covars, _set_covars)

    def _compute_log_likelihood(self, obs):
        return lmvnpdf(obs, self._means, self._covars, self._cvtype)

    def _generate_sample_from_state(self, state):
        if self._cvtype == 'tied':
            cv = self._covars
        else:
            cv = self._covars[state]
        return sample_gaussian(self._means[state], cv, self._cvtype)

    def _init(self, obs, params='stmc'):
        super(GaussianHMM, self)._init(obs, params=params)

        if (hasattr(self, 'n_features')
            and self.n_features != obs[0].shape[1]):
            raise ValueError('Unexpected number of dimensions, got %s but '
                             'expected %s' % (obs[0].shape[1],
                                              self.n_features))

        self.n_features = obs[0].shape[1]

        if 'm' in params:
            self._means = cluster.KMeans(
                k=self._n_states).fit(obs[0]).cluster_centers_
        if 'c' in params:
            cv = np.cov(obs[0].T)
            if not cv.shape:
                cv.shape = (1, 1)
            self._covars = _distribute_covar_matrix_to_match_cvtype(
                cv, self._cvtype, self._n_states)

    def _initialize_sufficient_statistics(self):
        stats = super(GaussianHMM, self)._initialize_sufficient_statistics()
        stats['post'] = np.zeros(self._n_states)
        stats['obs'] = np.zeros((self._n_states, self.n_features))
        stats['obs**2'] = np.zeros((self._n_states, self.n_features))
        stats['obs*obs.T'] = np.zeros((self._n_states, self.n_features,
                                       self.n_features))
        return stats

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        super(GaussianHMM, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)

        if 'm' in params or 'c' in params:
            stats['post'] += posteriors.sum(axis=0)
            stats['obs'] += np.dot(posteriors.T, obs)

        if 'c' in params:
            if self._cvtype in ('spherical', 'diag'):
                stats['obs**2'] += np.dot(posteriors.T, obs ** 2)
            elif self._cvtype in ('tied', 'full'):
                for t, o in enumerate(obs):
                    obsobsT = np.outer(o, o)
                    for c in xrange(self._n_states):
                        stats['obs*obs.T'][c] += posteriors[t, c] * obsobsT

    def _do_mstep(self, stats, params, **kwargs):
        super(GaussianHMM, self)._do_mstep(stats, params)

        # Based on Huang, Acero, Hon, "Spoken Language Processing",
        # p. 443 - 445
        denom = stats['post'][:, np.newaxis]
        if 'm' in params:
            prior = self.means_prior
            weight = self.means_weight
            if prior is None:
                weight = 0
                prior = 0
            self._means = (weight * prior + stats['obs']) / (weight + denom)

        if 'c' in params:
            covars_prior = self.covars_prior
            covars_weight = self.covars_weight
            if covars_prior is None:
                covars_weight = 0
                covars_prior = 0

            means_prior = self.means_prior
            means_weight = self.means_weight
            if means_prior is None:
                means_weight = 0
                means_prior = 0
            meandiff = self._means - means_prior

            if self._cvtype in ('spherical', 'diag'):
                cv_num = (means_weight * (meandiff) ** 2
                          + stats['obs**2']
                          - 2 * self._means * stats['obs']
                          + self._means ** 2 * denom)
                cv_den = max(covars_weight - 1, 0) + denom
                if self._cvtype == 'spherical':
                    self._covars = (covars_prior / cv_den.mean(axis=1)
                                   + np.mean(cv_num / cv_den, axis=1))
                elif self._cvtype == 'diag':
                    self._covars = (covars_prior + cv_num) / cv_den
            elif self._cvtype in ('tied', 'full'):
                cvnum = np.empty((self._n_states, self.n_features,
                                  self.n_features))
                for c in xrange(self._n_states):
                    obsmean = np.outer(stats['obs'][c], self._means[c])

                    cvnum[c] = (means_weight * np.outer(meandiff[c],
                                                        meandiff[c])
                                + stats['obs*obs.T'][c]
                                - obsmean - obsmean.T
                                + np.outer(self._means[c], self._means[c])
                                * stats['post'][c])
                cvweight = max(covars_weight - self.n_features, 0)
                if self._cvtype == 'tied':
                    self._covars = ((covars_prior + cvnum.sum(axis=0))
                                    / (cvweight + stats['post'].sum()))
                elif self._cvtype == 'full':
                    self._covars = ((covars_prior + cvnum)
                                   / (cvweight + stats['post'][:, None, None]))


class MultinomialHMM(_BaseHMM):
    """Hidden Markov Model with multinomial (discrete) emissions

    Attributes
    ----------
    n_states : int (read-only)
        Number of states in the model.
    n_symbols : int
        Number of possible symbols emitted by the model (in the observations).
    transmat : array, shape (`n_states`, `n_states`)
        Matrix of transition probabilities between states.
    startprob : array, shape ('n_states`,)
        Initial state occupation distribution.
    emissionprob: array, shape ('n_states`, 'n_symbols`)
        Probability of emitting a given symbol when in each state.

    Methods
    -------
    eval(X)
        Compute the log likelihood of `X` under the HMM.
    decode(X)
        Find most likely state sequence for each point in `X` using the
        Viterbi algorithm.
    rvs(n=1)
        Generate `n` samples from the HMM.
    init(X)
        Initialize HMM parameters from `X`.
    fit(X)
        Estimate HMM parameters from `X` using the Baum-Welch algorithm.
    predict(X)
        Like decode, find most likely state sequence corresponding to `X`.
    score(X)
        Compute the log likelihood of `X` under the model.

    Examples
    --------
    >>> from scikits.learn.hmm import MultinomialHMM
    >>> MultinomialHMM(n_states=2)
    ... #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    MultinomialHMM(transmat=array([[ 0.5,  0.5],
           [ 0.5,  0.5]]),
            startprob_prior=1.0, n_states=2, startprob=array([ 0.5,  0.5]),
           transmat_prior=1.0)

    See Also
    --------
    GaussianHMM : HMM with Gaussian emissions
    """

    def __init__(self, n_states=1, startprob=None, transmat=None,
                 startprob_prior=None, transmat_prior=None):
        """Create a hidden Markov model with multinomial emissions.

        Parameters
        ----------
        n_states : int
            Number of states.
        """
        super(MultinomialHMM, self).__init__(n_states, startprob, transmat,
                                             startprob_prior=startprob_prior,
                                             transmat_prior=transmat_prior)

    def _get_emissionprob(self):
        """Emission probability distribution for each state."""
        return np.exp(self._log_emissionprob)

    def _set_emissionprob(self, emissionprob):
        emissionprob = np.asanyarray(emissionprob)
        if hasattr(self, 'n_symbols') and \
               emissionprob.shape != (self._n_states, self.n_symbols):
            raise ValueError('emissionprob must have shape '
                             '(n_states, n_symbols)')

        self._log_emissionprob = np.log(emissionprob)
        underflow_idx = np.isnan(self._log_emissionprob)
        self._log_emissionprob[underflow_idx] = -np.Inf
        self.n_symbols = self._log_emissionprob.shape[1]

    emissionprob = property(_get_emissionprob, _set_emissionprob)

    def _compute_log_likelihood(self, obs):
        return self._log_emissionprob[:, obs].T

    def _generate_sample_from_state(self, state):
        cdf = np.cumsum(self.emissionprob[state, :])
        rand = np.random.rand()
        symbol = (cdf > rand).argmax()
        return symbol

    def _init(self, obs, params='ste'):
        super(MultinomialHMM, self)._init(obs, params=params)

        if 'e' in params:
            emissionprob = normalize(np.random.rand(self._n_states,
                                                    self.n_symbols), 1)
            self.emissionprob = emissionprob

    def _initialize_sufficient_statistics(self):
        stats = super(MultinomialHMM, self)._initialize_sufficient_statistics()
        stats['obs'] = np.zeros((self._n_states, self.n_symbols))
        return stats

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        super(MultinomialHMM, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)
        if 'e' in params:
            for t, symbol in enumerate(obs):
                stats['obs'][:, symbol] += posteriors[t, :]

    def _do_mstep(self, stats, params, **kwargs):
        super(MultinomialHMM, self)._do_mstep(stats, params)
        if 'e' in params:
            self.emissionprob = (stats['obs']
                                 / stats['obs'].sum(1)[:, np.newaxis])


class GMMHMM(_BaseHMM):
    """Hidden Markov Model with Gaussin mixture emissions

    Attributes
    ----------
    n_states : int (read-only)
        Number of states in the model.
    transmat : array, shape (`n_states`, `n_states`)
        Matrix of transition probabilities between states.
    startprob : array, shape ('n_states`,)
        Initial state occupation distribution.
    gmms: array of GMM objects, length 'n_states`
        GMM emission distributions for each state

    Methods
    -------
    eval(X)
        Compute the log likelihood of `X` under the HMM.
    decode(X)
        Find most likely state sequence for each point in `X` using the
        Viterbi algorithm.
    rvs(n=1)
        Generate `n` samples from the HMM.
    init(X)
        Initialize HMM parameters from `X`.
    fit(X)
        Estimate HMM parameters from `X` using the Baum-Welch algorithm.
    predict(X)
        Like decode, find most likely state sequence corresponding to `X`.
    score(X)
        Compute the log likelihood of `X` under the model.

    Examples
    --------
    >>> from scikits.learn.hmm import GMMHMM
    >>> GMMHMM(n_states=2, n_mix=10, cvtype='diag')
    ... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    GMMHMM(n_mix=10, cvtype='diag', n_states=2, startprob_prior=1.0,
        startprob=array([ 0.5,  0.5]),
        transmat=array([[ 0.5,  0.5],
           [ 0.5,  0.5]]),
        transmat_prior=1.0,
        gmms=[GMM(cvtype='diag', n_states=10), GMM(cvtype='diag',
              n_states=10)])

    See Also
    --------
    GaussianHMM : HMM with Gaussian emissions
    """

    def __init__(self, n_states=1, n_mix=1, startprob=None,
                 transmat=None, startprob_prior=None, transmat_prior=None,
                 gmms=None, cvtype=None):
        """Create a hidden Markov model with GMM emissions.

        Parameters
        ----------
        n_states : int
            Number of states.
        """
        super(GMMHMM, self).__init__(n_states, startprob, transmat,
                                     startprob_prior=startprob_prior,
                                     transmat_prior=transmat_prior)

        # XXX: Hotfit for n_mix that is incompatible with the scikit's
        # BaseEstimator API
        self.n_mix = n_mix
        self.cvtype = cvtype
        if gmms is None:
            gmms = []
            for x in xrange(self.n_states):
                if cvtype is None:
                    g = GMM(n_mix)
                else:
                    g = GMM(n_mix, cvtype=cvtype)
                gmms.append(g)
        self.gmms = gmms

    def _compute_log_likelihood(self, obs):
        return np.array([g.score(obs) for g in self.gmms]).T

    def _generate_sample_from_state(self, state):
        return self.gmms[state].rvs(1).flatten()

    def _init(self, obs, params='stwmc'):
        super(GMMHMM, self)._init(obs, params=params)

        allobs = np.concatenate(obs, 0)
        for g in self.gmms:
            g.fit(allobs, n_iter=0, init_params=params)

    def _initialize_sufficient_statistics(self):
        stats = super(GMMHMM, self)._initialize_sufficient_statistics()
        stats['norm'] = [np.zeros(g.weights.shape) for g in self.gmms]
        stats['means'] = [np.zeros(np.shape(g.means)) for g in self.gmms]
        stats['covars'] = [np.zeros(np.shape(g._covars)) for g in self.gmms]
        return stats

    def _accumulate_sufficient_statistics(self, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        super(GMMHMM, self)._accumulate_sufficient_statistics(
            stats, obs, framelogprob, posteriors, fwdlattice, bwdlattice,
            params)

        for state, g in enumerate(self.gmms):
            gmm_logprob, gmm_posteriors = g.eval(obs)
            gmm_posteriors *= posteriors[:, state][:, np.newaxis]
            tmpgmm = GMM(g.n_states, cvtype=g.cvtype)
            tmpgmm.n_features = g.n_features
            tmpgmm.covars = _distribute_covar_matrix_to_match_cvtype(
                np.eye(g.n_features), g.cvtype, g.n_states)
            norm = tmpgmm._do_mstep(obs, gmm_posteriors, params)

            stats['norm'][state] += norm
            if 'm' in params:
                stats['means'][state] += tmpgmm.means * norm[:, np.newaxis]
            if 'c' in params:
                if tmpgmm.cvtype == 'tied':
                    stats['covars'][state] += tmpgmm._covars * norm.sum()
                else:
                    cvnorm = np.copy(norm)
                    shape = np.ones(tmpgmm._covars.ndim)
                    shape[0] = np.shape(tmpgmm._covars)[0]
                    cvnorm.shape = shape
                    stats['covars'][state] += tmpgmm._covars * cvnorm

    def _do_mstep(self, stats, params, covars_prior=1e-2, **kwargs):
        super(GMMHMM, self)._do_mstep(stats, params)
        # All that is left to do is to apply covars_prior to the
        # parameters updated in _accumulate_sufficient_statistics.
        for state, g in enumerate(self.gmms):
            norm = stats['norm'][state]
            if 'w' in params:
                g.weights = normalize(norm)
            if 'm' in params:
                g.means = stats['means'][state] / norm[:, np.newaxis]
            if 'c' in params:
                if g.cvtype == 'tied':
                    g.covars = ((stats['covars'][state]
                                 + covars_prior * np.eye(g.n_features))
                                / norm.sum())
                else:
                    cvnorm = np.copy(norm)
                    shape = np.ones(g._covars.ndim)
                    shape[0] = np.shape(g._covars)[0]
                    cvnorm.shape = shape
                    if g.cvtype == 'spherical' or g.cvtype == 'diag':
                        g.covars = (stats['covars'][state]
                                    + covars_prior) / cvnorm
                    elif g.cvtype == 'full':
                        eye = np.eye(g.n_features)
                        g.covars = ((stats['covars'][state]
                                     + covars_prior * eye[np.newaxis, :, :])
                                    / cvnorm)
