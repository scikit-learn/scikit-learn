"""Bagging meta-estimator
   This module is adopted from Random forests module in Scikit-learn"""

#Author: Maheshakya Wijewardena <pmaheshakya4@gmail.com>
#License: BSD 3 clause

from __future__ import division

import itertools
import numbers
import numpy as np
from warnings import warn
from abc import ABCMeta, abstractmethod
from inspect import getargspec


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
        
        #Feature bootstraping is not implemented, therefore all features are drawed
        #For scalability purposes(add feature bootstrapping), the functionality can be addad here
        features = np.array([i for i in range(n_features)])
            

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

def _parallel_predict_proba(estimators, estimators_features, X, n_classes):
    """Private function used to compute (proba-)predictions within a job."""
    
    n_samples = X.shape[0]
    proba = np.zeros((n_samples, n_classes))

    for estimator, features in zip(estimators, estimators_features):
        try:
            proba_estimator = estimator.predict_proba(X[:, features])

            if n_classes == len(estimator.classes_):
                proba += proba_estimator

            else:
                proba[:, estimator.classes_] += \
                    proba_estimator[:, range(len(estimator.classes_))]

        except (AttributeError, NotImplementedError):
            # Resort to voting
            predictions = estimator.predict(X[:, features])

            for i in range(n_samples):
                proba[i, predictions[i]] += 1

    return proba


def _parallel_predict_log_proba(estimators, estimators_features, X, n_classes):
    """Private function used to compute log probabilities within a job."""
    n_samples = X.shape[0]
    log_proba = np.empty((n_samples, n_classes))
    log_proba.fill(-np.inf)
    all_classes = np.arange(n_classes, dtype=np.int)

    for estimator, features in zip(estimators, estimators_features):
        log_proba_estimator = estimator.predict_log_proba(X[:, features])

        if n_classes == len(estimator.classes_):
            log_proba = np.logaddexp(log_proba, log_proba_estimator)

        else:
            log_proba[:, estimator.classes_] = np.logaddexp(log_proba[:, estimator.classes_], log_proba_estimator[:, range(len(estimator.classes_))])

            missing = np.setdiff1d(all_classes, estimator.classes_)
            log_proba[:, missing] = np.logaddexp(log_proba[:, missing],
                                                 -np.inf)

    return log_proba


def _parallel_decision_function(estimators, estimators_features, X):
    """Private function used to compute decisions within a job."""
    return sum(estimator.decision_function(X[:, features])
               for estimator, features in zip(estimators,
                                              estimators_features))



def _parallel_predict_regression(estimators, estimators_features, X):
    """Private funtion which predicts with in a job in regression"""
    return sum(estimator.predict(X[:, features]) for estimator, features in zip(estimators, estimators_features))


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
    """ A bagged regressor
        
        Bagged regressor is an ensemble method which applys same estimator on random subsets of data to build
        a new estimator. Original data set is sampled randomly and they are subjected into  regression 
        seperately, then aggregated(by averaging)to form the final estimator. This ensembled estimator can be
        used to reduce the variance of other non-linear black box estimators.(eg: Decistion Tree, KNN, ect.)
        
        Parameters
        ----------
        
        base_estimator : object or None, This is Compulsory (default=None)
            The base estimator to fit on random subsets of the dataset.
            If None, a value error is raised.

        n_estimators : int, optional (default=10)
            The number of base estimators in the ensemble.
    
        max_samples : int or float, optional (default=1.0)
            The number of samples to draw from X to train each base estimator.
                - If int, then draw `max_samples` samples.
                - If float, then draw `max_samples * X.shape[0]` samples.  
    
        bootstrap : boolean, optional (default=False)
            Whether samples are drawn with replacement.
    
        oob_score : bool
            Whether to use out-of-bag samples to estimate
            the generalization error.
    
        n_jobs : int, optional (default=1)
            The number of jobs to run in parallel for both `fit` and `predict`.
            If -1, then the number of jobs is set to the number of cores.
    
        random_state : int, RandomState instance or None, optional (default=None)
            If int, random_state is the seed used by the random number generator;
            If RandomState instance, random_state is the random number generator;
            If None, the random number generator is the RandomState instance used
            by `np.random`.
    
        verbose : int, optional (default=0)
            Controls the verbosity of the building process.
    
        Attributes
        ----------
        `estimators_`: list of estimators
            The collection of fitted sub-estimators.
    
        `estimators_samples_`: list of arrays
            The subset of drawn samples (i.e., the in-bag samples) for each base
            estimator.
    
        `estimators_features_`: list of arrays
            The subset of drawn features for each base estimator.
    
        `oob_score_` : float
            Score of the training dataset obtained using an out-of-bag estimate.
    
        `oob_decision_function_` : array of shape = [n_samples, n_classes]
            Decision function computed with out-of-bag estimate on the training
            set. If n_estimators is small it might be possible that a data point
            was never left out during the bootstrap. In this case,
            `oob_decision_function_` might contain NaN. """
    
    
    def __init__(self, base_estimator, n_estimators=10, max_samples=1.0, bootstrap=True, oob_score=False, n_jobs=1, random_state=None, verbose=0): 
        if base_estimator == None:
            raise ValueError("BaseEstimator should be defined")
        super(BaggingRegressor, self).__init__(base_estimator, n_estimators=n_estimators, max_samples=max_samples, bootstrap=bootstrap, oob_score=oob_score, n_jobs=n_jobs, random_state=random_state, verbose=verbose)
           
            
    def _set_oob_score(self, X, y):
        pass
    
    def predict(self, X):
        """Predicts regression target for X.

        The predicted regression target of an input sample is computed as the
        average predicted regression targets of the estimators in the ensemble.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y: array of shape = [n_samples]
            The predicted values.
        """
        # Check data
        X, = check_arrays(X)

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(self)

        all_y = Parallel(n_jobs=n_jobs, verbose=self.verbose)(delayed(_parallel_predict_regression)(self.estimators_[starts[i]:starts[i + 1]], self.estimators_features_[starts[i]:starts[i + 1]], X) for i in range(n_jobs))

        # Reduce
        y = sum(all_y) / self.n_estimators

        return y
    
    def _set_oob_score(self, X, y):
        n_samples = y.shape[0]

        predictions = np.zeros((n_samples,))
        n_predictions = np.zeros((n_samples,))
        
        for estimator, samples, features in zip(self.estimators_,
                                                self.estimators_samples_,
                                                self.estimators_features_):
            mask = np.ones(n_samples, dtype=np.bool)
            mask[samples] = False

            predictions[mask] += estimator.predict((X[mask, :])[:, features])
            n_predictions[mask] += 1

        if (n_predictions == 0).any():
            warn("Some inputs do not have OOB scores. "
                 "This probably means too few estimators were used "
                 "to compute any reliable oob estimates.")
            n_predictions[n_predictions == 0] = 1

        predictions /= n_predictions

        self.oob_prediction_ = predictions
        self.oob_score_ = r2_score(y, predictions)
        print self.oob_score_

        
        
    
        
        
        
        
            
        
    
            
        
 

class BaggingClassifier(BaseBagging, ClassifierMixin):
    """ A bagged classifier"""    
    pass
















