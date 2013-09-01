"""Bagging meta-estimator"""

#Author: Maheshakya Wijewardena <pmaheshakya4@gmail.com>
#License: BSD 3 clause

from __future__ import division

import itertools
import numbers
import numpy as np
from warnings import warn
from abc import ABCMeta, abstractmethod
from inspect import getargspec

"""
if __name__ == "__main__" and __package__ is None:
    __package__ = "sklearn.ensemble"
"""
from ..base import ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed, cpu_count
from ..externals import six
from ..externals.six.moves import xrange
from ..metrics import r2_score
from ..utils import check_random_state, check_arrays
from ..utils.fixes import bincount, unique
from ..utils.random import sample_without_replacement

from .base import BaseEnsemble

__all__ = ["BaggingClassifier", "BaggingRegressor", "BaseBagging"]

MAX_INT = np.iinfo(np.int32).max

#Bootstrap features function is not implemented here

""" There will be many estimators with different samples for the estimator that user would want to use for bagging. Those estimators have to be distributed among jobs."""
def _partition_estimators(ensemble):
    """This private function will partition estimators among jobs"""
    
    #Computing the number of jobs required
    if ensemble.n_jobs == -1:
        n_jobs = min(cpu_count(), ensemble.n_estimators)
    else: #minimum of n_jobs and n_estimators(user specified) will be taken as the number of jobs
        n_jobs = min(ensemble.n_jobs, ensemble.n_estimators) 
    
    #Partition estimators among jobs
    n_estimators = (ensemble.n_estimators//n_jobs) * np.ones(n_jobs, dtype = np.int)
    n_estimators[:ensemble.n_estimators%n_jobs] += 1
    starts = np.cumsum(n_estimators)
    
    return n_jobs, n_estimators.tolist(), [0] + starts.tolist()

def _parallel_build_estimators(n_estimators, ensemble, X, y, sample_weight, seeds, verbose):
    """Private function used to build a batch of estimators within a job."""
    
    # Retrieve settings
    n_samples, n_features = X.shape
    max_samples = ensemble.max_samples
    
    if (not isinstance(max_samples, (np.integer, numbers.Integral)) and (0.0 < max_samples <= 1.0)):
        max_samples = int(max_samples * n_samples)


    bootstrap = ensemble.bootstrap
    support_sample_weight = ("sample_weight" in getargspec(ensemble.base_estimator.fit)[0])

    # Build estimators
    estimators = []
    estimators_samples = []
    estimators_features = []

    for i in range(n_estimators):
        if verbose > 1:
            print("building estimator %d of %d" % (i + 1, n_estimators))

        random_state = check_random_state(seeds[i])
        seed = check_random_state(random_state.randint(MAX_INT))
        estimator = ensemble._make_estimator(append=False)

        try:  # Not all estimators accept a random_state
            estimator.set_params(random_state=seed)
        except ValueError:
            pass
        
        features = sample_without_replacement(n_features, n_features, random_state=random_state)      
            

        # Draw samples, using sample weights, and then fit
        if support_sample_weight:
            if sample_weight is None:
                curr_sample_weight = np.ones((n_samples,))
            else:
                curr_sample_weight = sample_weight.copy()

            if bootstrap:
                indices = random_state.randint(0, n_samples, max_samples)
                sample_counts = bincount(indices, minlength=n_samples)
                curr_sample_weight *= sample_counts

            else:
                not_indices = sample_without_replacement(
                    n_samples,
                    n_samples - max_samples,
                    random_state=random_state)

                curr_sample_weight[not_indices] = 0

            estimator.fit(X[:, features], y, sample_weight=curr_sample_weight)
            samples = curr_sample_weight > 0.

        # Draw samples, using a mask, and then fit
        else:
            if bootstrap:
                indices = random_state.randint(0, n_samples, max_samples)
            else:
                indices = sample_without_replacement(n_samples,
                                                     max_samples,
                                                     random_state=random_state)

            sample_counts = bincount(indices, minlength=n_samples)

            estimator.fit((X[indices])[:, features], y[indices])
            samples = sample_counts > 0.

        estimators.append(estimator)
        estimators_samples.append(samples)
        estimators_features.append(features)

    return estimators, estimators_samples, estimators_features

class BaseBagging(six.with_metaclass(ABCMeta, BaseEnsemble)):
    """Base class for Bagging
     Warning: This class should not be used directly. Use derived classes instead"""  
    
    @abstractmethod
    def __init__(self, base_estimator=None, n_estimators=10, max_samples=1.0, bootstrap=True, oob_score=False, n_jobs=1, random_state=None, verbose=0):
        super(BaseBagging,self).__init__(base_estimator=base_estimator,n_estimators = n_estimators)
        self.max_samples = max_samples
        self.bootstrap = bootstrap
        self.oob_score = oob_score
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.verbose = verbose
        
    
    
    @abstractmethod
    def _set_oob_score(self, X, y):
        """Calculate out of bag predictions and score."""
        
    def _validate_y(self,y):
        return y
    
    def fit(self, X, y, sample_weight=None):
        """Builds a ensemble of estimators for Bagging from the training set (X, y).

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (integers that correspond to classes in classification, real numbers for regression).

        sample_weight : array-like, shape = [n_samples] or None
            Sample weights. If None, then samples are equally weighted.

        Returns
        -------
        self : object
            Returns self.
        """
        
        random_state = check_random_state(self.random_state)
        X, y = check_arrays(X, y) #convert data        
        
        n_samples, self.n_features_ = X.shape # Remap output
        y = self._validate_y(y)
        
        # Check parameters
        if isinstance(self.max_samples, (numbers.Integral, np.integer)):
            max_samples = self.max_samples
        else:  # "max_samples" is a float
            max_samples = int(self.max_samples * X.shape[0])

        if not (0 < max_samples <= X.shape[0]):
            raise ValueError("max_samples must be in (0, n_samples]")

        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available"
                             " if bootstrap=True")
        
        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)
        seeds = random_state.randint(MAX_INT, size=self.n_estimators)
        
        all_results = Parallel(n_jobs=n_jobs, verbose=self.verbose)(delayed(_parallel_build_estimators)(n_estimators[i], self, X, y, sample_weight,            seeds[starts[i]:starts[i + 1]],
                verbose=self.verbose) for i in range(n_jobs))

        # Reduce
        self.estimators_ = list(itertools.chain.from_iterable(t[0] for t in all_results))
        self.estimators_samples_ = list(itertools.chain.from_iterable(t[1] for t in all_results))
        self.estimators_features_ = list(itertools.chain.from_iterable(t[2] for t in all_results))

        if self.oob_score:
            self._set_oob_score(X, y)

        return self


class BaggingRegressor(BaseBagging, RegressorMixin):    
    """ A bagged regressor"""
    def __init__(self, base_estimator, n_estimators=10, max_samples=1.0, bootstrap=True, oob_score=False, n_jobs=1, random_state=None, verbose=0):
            super(BaggingRegressor, self).__init__(base_estimator, n_estimators=n_estimators, max_samples=max_samples, bootstrap=bootstrap, oob_score=oob_score, n_jobs=n_jobs, random_state=random_state, verbose=verbose)
            
    def _set_oob_score(self, X, y):
        pass


        
            
        
 

class BaggingClassifier(BaseBagging, ClassifierMixin):
    """ A bagged classifier"""    
    pass
















