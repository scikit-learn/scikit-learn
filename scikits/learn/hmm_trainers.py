import string

import numpy as np

from gmm import *

class HMMTrainer(object):
    """Base class for HMM training algorithms."""

    @property
    def emission_type(self):
        pass

    def train(self, hmm, obs, niter=10, thresh=1e-2, params=string.letters,
              maxrank=None, beamlogprob=-np.Inf, **kwargs):
        """Estimate model parameters.

        Parameters
        ----------
        hmm : HMM object
            HMM to train.
        obs : list
            List of array-like observation sequences (shape (n_i, ndim)).
        niter : int
            Number of iterations to perform.
        thresh : float
            Convergence threshold.
        params : string
            Controls which parameters are updated in the training
            process.  Can contain any combination of 's' for startprob,
            't' for transmat, 'm' for means, and 'c' for covars, etc.
            Defaults to all parameters.
        maxrank : int
            Maximum rank to evaluate for rank pruning.  If not None,
            only consider the top `maxrank` states in the inner
            sum of the forward algorithm recursion.  Defaults to None
            (no rank pruning).  See "The HTK Book" for more details.
        beamlogprob : float
            Width of the beam-pruning beam in log-probability units.
            Defaults to -numpy.Inf (no beam pruning).  See The HTK
            Book for more details.

        Returns
        -------
        logprob : list
            Log probabilities of the training data after each iteration.

        Notes
        -----
        In general, `logprob` should be non-decreasing unless
        aggressive pruning is used.  Decreasing `logprob` is generally
        a sign of overfitting (e.g. a covariance parameter getting too
        small).  You can fix this by using a different trainer
        (e.g. based on model adaptation), getting more training data,
        or decreasing `covarprior`.
        """
        logprob = []
        for i in xrange(niter):
            # Expectation step
            stats = self._initialize_sufficient_statistics(hmm)
            curr_logprob = 0
            for seq in obs:
                framelogprob = hmm._compute_log_likelihood(seq)
                lpr, fwdlattice = hmm._do_forward_pass(framelogprob, maxrank,
                                                       beamlogprob)
                bwdlattice = hmm._do_backward_pass(framelogprob, fwdlattice,
                                                   maxrank, beamlogprob)
                gamma = fwdlattice + bwdlattice
                posteriors = np.exp(gamma.T - logsum(gamma, axis=1)).T
                curr_logprob += lpr
                self._accumulate_sufficient_statistics(hmm, stats, seq,
                                                       framelogprob, posteriors,
                                                       fwdlattice, bwdlattice,
                                                       params)
            logprob.append(curr_logprob)

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < thresh:
                break

            # Maximization step
            self._do_mstep(hmm, stats, params, **kwargs)

        return logprob

    def _initialize_sufficient_statistics(self, hmm):
        pass

    def _accumulate_sufficient_statistics(self, hmm, stats, seq, framelogprob, 
                                          posteriors, fwdlattice,
                                          bwdlattice, params):
        pass
    
    def _do_mstep(self, hmm, stats, params, **kwargs):
        pass


class BaseHMMBaumWelchTrainer(HMMTrainer):
    """Base class for HMM trainers.

    Uses the Baum-Welch algorithm to train the startprob and transmat
    parameters.
    """
    emission_type = None

    def _initialize_sufficient_statistics(self, hmm):
        stats = {'nobs':  0,
                 'start': np.zeros(hmm._nstates),
                 'trans': np.zeros((hmm._nstates, hmm._nstates))}
        return stats

    def _accumulate_sufficient_statistics(self, hmm, stats, seq, framelogprob, 
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        stats['nobs'] += 1
        if 's' in params:
            stats['start'] += posteriors[0]
        if 't' in params:
            for t in xrange(len(framelogprob)):
                zeta = (fwdlattice[t-1][:,np.newaxis] + hmm._log_transmat
                        + framelogprob[t] + bwdlattice[t])
                stats['trans'] += np.exp(zeta - logsum(zeta))

    def _do_mstep(self, hmm, stats, params, **kwargs):
        if 's' in params:
            hmm.startprob = stats['start'] / stats['start'].sum()
        if 't' in params:
            hmm.transmat = normalize(stats['trans'], axis=1)


class GaussianHMMBaumWelchTrainer(BaseHMMBaumWelchTrainer):
    """Baum-Welch maximum likelihood trainer for HMMs with Gaussian emissions."""
    emission_type = 'gaussian'

    def _initialize_sufficient_statistics(self, hmm):
        stats = super(GaussianHMMBaumWelchTrainer,
                      self)._initialize_sufficient_statistics(hmm)
        stats['post']      = np.zeros(hmm._nstates)
        stats['obs']       = np.zeros((hmm._nstates, hmm._ndim))
        stats['obs**2']    = np.zeros((hmm._nstates, hmm._ndim))
        stats['obs*obs.T'] = np.zeros((hmm._nstates, hmm._ndim, hmm._ndim))
        return stats

    def _accumulate_sufficient_statistics(self, hmm, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        super(GaussianHMMBaumWelchTrainer,
              self)._accumulate_sufficient_statistics(hmm, stats, obs,
                                                      framelogprob, posteriors,
                                                      fwdlattice, bwdlattice,
                                                      params)

        if 'm' in params or 'c' in params:
            stats['post'] += posteriors.sum(axis=0)
            stats['obs'] += np.dot(posteriors.T, obs)

        if 'c' in params:
            if hmm._cvtype in ('spherical', 'diag'):
                stats['obs**2'] += np.dot(posteriors.T, obs**2)
            elif hmm._cvtype in ('tied', 'full'):  
                for t, o in enumerate(obs):
                    obsobsT = np.outer(o, o)
                    for c in xrange(hmm._nstates):
                        stats['obs*obs.T'][c] += posteriors[t,c] * obsobsT
                  
    def _do_mstep(self, hmm, stats, params, covarprior=1e-2, **kwargs):
        super(GaussianHMMBaumWelchTrainer, self)._do_mstep(hmm, stats, params)

        denom = stats['post'][:,np.newaxis]
        if 'm' in params:
            hmm._means = stats['obs'] / denom

        if 'c' in params:
            if hmm._cvtype in ('spherical', 'diag'):
                cv = ((stats['obs**2']
                       - 2 * hmm._means * stats['obs']
                       + hmm._means**2 * denom
                       + covarprior)
                      / (1.0 + denom))
                if hmm._cvtype == 'spherical':
                    hmm._covars = cv.mean(axis=1)
                elif hmm._cvtype == 'diag':
                    hmm._covars = cv
            elif hmm._cvtype in ('tied', 'full'):
                cvnum = np.empty((hmm._nstates, hmm._ndim, hmm._ndim))
                cvprior = np.eye(hmm._ndim) * covarprior
                for c in xrange(hmm._nstates):
                    cvnum[c] = (stats['obs*obs.T'][c]
                                - 2 * np.outer(stats['obs'][c], hmm._means[c])
                                + np.outer(hmm._means[c] * stats['post'][c],
                                           hmm._means[c]))
                if hmm._cvtype == 'tied':
                    hmm._covars = ((cvnum.sum(0) + cvprior)
                                   / (1.0 + stats['post'].sum(0)))
                elif hmm._cvtype == 'full':
                    hmm._covars = ((cvnum + cvprior)
                                   / (1.0 + stats['post'][:,np.newaxis,np.newaxis]))


class GaussianHMMMAPTrainer(GaussianHMMBaumWelchTrainer):
    """HMM trainer based on maximum-a-posteriori (MAP) adaptation.
    """
    emission_type = 'gaussian'

    def __init__(self, startprob_prior=None, transmat_prior=None,
                 means_prior=None, means_weight=0,
                 covars_prior=None, covars_weight=0):
        """Initialize MAP trainer with the parameters for the prior
        distribution.
        """
        self.startprob_prior = startprob_prior
        self.transmat_prior = transmat_prior
        self.means_prior = means_prior
        self.means_weight = means_weight
        self.covars_prior = covars_prior
        self.covars_weight = covars_weight

    def _do_mstep(self, hmm, stats, params, **kwargs):
        # Based on Huang, Acero, Hon, "Spoken Language Processing", p. 443 - 445
        if 's' in params:
            prior = self.startprob_prior
            if prior is None:
                prior = 1.0
            hmm.startprob = normalize(np.maximum(prior - 1.0 + stats['start'],
                                                 1e-20))

        if 't' in params:
            prior = self.transmat_prior
            if prior is None:
                prior = 1.0
            hmm.transmat = normalize(np.maximum(prior - 1.0 + stats['trans'],
                                                1e-20), axis=1)

        denom = stats['post'][:,np.newaxis]
        if 'm' in params:
            prior = self.means_prior
            weight = self.means_weight
            if prior is None:
                weight = 0
                prior = 0
            hmm._means = (weight * prior + stats['obs']) / (weight + denom)

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
            meandiff = hmm._means - means_prior

            if hmm._cvtype in ('spherical', 'diag'):
                cv_num = (means_weight * (meandiff)**2
                          + stats['obs**2']
                          - 2 * hmm._means * stats['obs']
                          + hmm._means**2 * denom)
                cv_den = max(covars_weight - 1, 0) + denom
                if hmm._cvtype == 'spherical':
                    hmm._covars = (covars_prior / cv_den.mean(axis=1)
                                   + np.mean(cv_num / cv_den, axis=1))
                elif hmm._cvtype == 'diag':
                    hmm._covars = (covars_prior + cv_num) / cv_den
            elif hmm._cvtype in ('tied', 'full'):
                cvnum = np.empty((hmm._nstates, hmm._ndim, hmm._ndim))
                for c in xrange(hmm._nstates):
                    cvnum[c] = (means_weight * np.outer(meandiff[c], meandiff[c])
                                + stats['obs*obs.T'][c] 
                                - 2 * np.outer(stats['obs'][c], hmm._means[c])
                                + np.outer(hmm._means[c], hmm._means[c])
                                * stats['post'][c])
                cvweight = max(covars_weight - hmm._ndim, 0)
                if hmm._cvtype == 'tied':
                    hmm._covars = ((covars_prior + cvnum.sum(axis=0))
                                    / (cvweight + stats['post'].sum()))
                elif hmm._cvtype == 'full':
                    hmm._covars = ((covars_prior + cvnum)
                                   / (cvweight + stats['post'][:,None,None]))

 

class MultinomialHMMBaumWelchTrainer(BaseHMMBaumWelchTrainer):
    "Baum-Welch maximum likelihood trainer for HMM with multinomial emissions."
    emission_type = 'multinomial'

    def _initialize_sufficient_statistics(self, hmm):
        stats = super(MultinomialHMMBaumWelchTrainer,
                      self)._initialize_sufficient_statistics(hmm)
        stats['obs']  = np.zeros((hmm._nstates, hmm._nsymbols))
        return stats

    def _accumulate_sufficient_statistics(self, hmm, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        super(MultinomialHMMBaumWelchTrainer,
              self)._accumulate_sufficient_statistics(hmm, stats, obs,
                                                      framelogprob, posteriors,
                                                      fwdlattice, bwdlattice,
                                                      params)
        if 'e' in params:
            for t,symbol in enumerate(obs):
                stats['obs'][:,symbol] += posteriors[t,:]
                  
    def _do_mstep(self, hmm, stats, params, covarprior=1e-2, **kwargs):
        super(MultinomialHMMBaumWelchTrainer, self)._do_mstep(hmm, stats,
                                                              params)

        if 'e' in params:
            hmm.emissionprob = stats['obs'] / stats['obs'].sum(1)[:,np.newaxis]


class GMMHMMBaumWelchTrainer(BaseHMMBaumWelchTrainer):
    "Baum-Welch maximum likelihood trainer for HMM with GMM emissions."
    emission_type = 'gmm'

    def _initialize_sufficient_statistics(self, hmm):
        stats = super(GMMHMMBaumWelchTrainer,
                      self)._initialize_sufficient_statistics(hmm)
        stats['norm'] = [np.zeros(g.weights.shape) for g in hmm.gmms]
        stats['means'] = [np.zeros(np.shape(g.means)) for g in hmm.gmms]
        stats['covars'] = [np.zeros(np.shape(g._covars)) for g in hmm.gmms]
        return stats

    def _accumulate_sufficient_statistics(self, hmm, stats, obs, framelogprob,
                                          posteriors, fwdlattice, bwdlattice,
                                          params):
        super(GMMHMMBaumWelchTrainer,
              self)._accumulate_sufficient_statistics(hmm, stats, obs,
                                                      framelogprob, posteriors,
                                                      fwdlattice, bwdlattice,
                                                      params)
        for state,g in enumerate(hmm.gmms):
            gmm_logprob, gmm_posteriors = g.eval(obs)
            gmm_posteriors *= posteriors[:,state][:,np.newaxis]
            tmpgmm = GMM(g.nstates, g.ndim, cvtype=g.cvtype)
            norm = tmpgmm._do_mstep(obs, gmm_posteriors, params, min_covar=0)

            stats['norm'][state] += norm
            if 'm' in params:
                stats['means'][state] += tmpgmm.means * norm[:,np.newaxis]
            if 'c' in params:
                if tmpgmm.cvtype == 'tied':
                    stats['covars'][state] += tmpgmm._covars * norm.sum()
                else:
                    cvnorm = np.copy(norm)
                    shape = np.ones(tmpgmm._covars.ndim)
                    shape[0] = np.shape(tmpgmm._covars)[0]
                    cvnorm.shape = shape
                    stats['covars'][state] += tmpgmm._covars * cvnorm
                  
    def _do_mstep(self, hmm, stats, params, covarprior=1e-2, **kwargs):
        super(GMMHMMBaumWelchTrainer, self)._do_mstep(hmm, stats, params)
        # All we have left to do is apply covarprior to the parameters
        # we updated in _accumulate_sufficient_statistics.
        for state,g in enumerate(hmm.gmms):
            norm = stats['norm'][state]
            #print norm
            if 'w' in params:
                g.weights = normalize(norm) 
            if 'm' in params:
                g.means = stats['means'][state] / norm[:,np.newaxis]
            if 'c' in params:
                if g.cvtype == 'tied':
                    g.covars = (stats['covars'][state]
                                + covarprior * np.eye(g.ndim)) / norm.sum()
                else:
                    cvnorm = np.copy(norm)
                    shape = np.ones(g._covars.ndim)
                    shape[0] = np.shape(g._covars)[0]
                    cvnorm.shape = shape
                    if g.cvtype == 'spherical' or g.cvtype == 'diag':
                        g.covars = (stats['covars'][state]
                                    + covarprior) / cvnorm
                    elif g.cvtype == 'full':
                        g.covars = ((stats['covars'][state]
                                     + covarprior*np.eye(g.ndim)[np.newaxis,:,:])
                                    / cvnorm)
        
