"""
The :mod:`sklearn.gp_search` includes utilities to fine-tune the parameters
of an estimator through a Gaussian Process model.
"""
from __future__ import print_function

# Author: Sebastien Dubois <http://bit.ly/SebastienDubois>,
# License: BSD 3 clause

import numpy as np
from scipy.stats import norm

from .gaussian_process.gaussian_process import GaussianProcess
from .cross_validation import check_cv
from .cross_validation import _fit_and_score
from .metrics.scorer import check_scoring
from .base import is_classifier, clone


#   UTILS    #

def sample_candidates(n_candidates, param_bounds, param_isInt):

    n_parameters = param_isInt.shape[0]
    candidates = []

    for k in range(n_parameters):
        if(param_isInt[k]):
            k_sample = np.asarray(
                np.random.rand(n_candidates)
                * np.float(param_bounds[k][1]-param_bounds[k][0])
                + param_bounds[k][0],
                dtype=np.int32)
        else:
            k_sample = np.asarray(
                np.random.rand(n_candidates)
                * np.float(param_bounds[k][1]-param_bounds[k][0])
                + param_bounds[k][0])
        candidates.append(k_sample)

    candidates = np.asarray(candidates)
    candidates = candidates.T

    return compute_unique(candidates)


def compute_ei(predictions, sigma, y_best):
    ei_array = np.zeros(predictions.shape[0])
    for i in range(ei_array.shape[0]):
        z = (y_best - predictions[i]) / sigma[i]
        ei_array[i] = sigma[i] * (z * norm.cdf(z) + norm.pdf(z))

    return ei_array


def compute_unique(a):
    # keep only unique values in the ndarray a
    # http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array

    b = np.ascontiguousarray(a).view(
        np.dtype((np.void, a.dtype.itemsize * a.shape[1]))
        )
    _, idx = np.unique(b, return_index=True)
    idx = np.sort(idx)

    return a[idx]


def is_in_ndarray(item, a):
    # look for element item in ndarray a
    # returns True if item is in a, and its index

    k = 0
    idx_val = np.asarray(range(a.shape[0]))
    idxk = range(a.shape[0])
    while(k < a.shape[1]):
        idxk = (a[idxk, k] == item[k])
        if(np.sum(idxk > 0)):
            k += 1
            idx_val = idx_val[idxk]
            idxk = list(idx_val)
        else:
            return False, 0

    return True, idx_val[0]


#    GPSearchCV    #
class GPSearchCV(object):
    """
    Parameters
    ----------

    parameters : dict, parameter space on which to optimize the estimator
        The keys of the dictionnary should be the names of the parameters,
        and the values should be lists of length 2; the first element being
        the type of the parameter ('int', 'float' or 'cat' [for categorical]),
        and the second element being a list of either the bounds between which
        to search (for 'int' and 'float') or the values the parameter can take
        (for 'cat')
        Example : parameters = {'kernel' :  ['cat', ['rbf', 'poly']],
                                 'd' : ['int', [1,3]],
                                 'C' : ['float',[1,10])}

    estimator : 1) sklearn estimator or 2) callable
        1 : object type that implements the "fit" and "predict" methods,
        as a classifier or a pipeline
        2 : a function that computes the output given a dictionnary of
        parameters. The returned value should be a list of one or more
        floats if score_format == 'cv', and a float if score_format ==
        'avg'

    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    y : array-like, shape = [n_samples] or [n_samples, n_output], optional
        Target relative to X for classification or regression;
        None for unsupervised learning.

    fit_params : dict, optional
        Parameters to pass to the fit method.

    scoring : string, callable or None, optional
        A string (see sklearn's model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
        Default is None.

    cv : integer or cross-validation generator, optional
        Relevant if the estimator is an sklearn object.
        If an integer is passed, it is the number of folds.
        Specific cross-validation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects
        Default is 5.

    acquisition function : string, optional
        Function to maximize in order to choose the next parameter to test.
        - Simple : maximize the predicted output
        - UCB : maximize the upper confidence bound
        - EI : maximizes the expected improvement
        Default is 'UCB'

    n_iter : int
        Total number of iterations to perform (including n_init and
        n_final_iter).
        Default is 100.

    n_init : int, optional
        Number of random iterations to perform before the smart search.
        Default is 30.

    n_final_iter : int, optional
        Number of final iterations, ie. smart iterations but with
        acquisition_function == 'Simple'
        Default is 5.

    n_candidates : int, optional
        Number of random candidates to sample for each GP iterations
        Default is 500.

    nugget : float, optional
        The nugget to set for the Gaussian Process.
        Default is 1.e-10.


    Attributes
    ----------

    best_parameter_ : dict, the parameter set, from those tested by the
        method _fit, that maximizes the mean of the cross-validation results.

    tested_parameters_ : ndarray, the parameters tested by _fit


    Examples
    -------
    >>> from sklearn.datasets import load_digits
    >>> iris = load_digits()
    >>> X, y = iris.data, iris.target
    >>> clf = RandomForestClassifier(n_estimators=20)
    >>> parameters = {"max_depth": ['int', [3, 3]],
                    "max_features": ['int', [1, 11]],
                    "min_samples_split": ['int', [1, 11]],
                    "min_samples_leaf": ['int', [1, 11]],
                    "bootstrap": ['cat', [True, False]],
                    "criterion": ['cat', ["gini", "entropy"]]}

    >>> search = GPSearchCV(parameters,
                            estimator=clf,
                            X=X,
                            y=y,
                            n_iter=20)
    >>> search._fit()

    """

    def __init__(self,
                 parameters,
                 estimator,
                 scoring=None,
                 X=None, y=None,
                 fit_params=None,
                 refit=True,
                 cv=None,
                 acquisition_function='UCB',
                 n_iter=100,
                 n_init=10,
                 n_candidates=500,
                 gp_nugget=1.e-10,
                 verbose=True):

        self.parameters = parameters
        self.n_parameters = len(parameters)
        self.acquisition_function = acquisition_function
        self.n_iter = n_iter
        self.n_init = n_init
        self.n_candidates = n_candidates
        self.param_names = parameters.keys()
        self.param_isInt = np.array([0 if (parameters[k][0] == 'float')
                                     else 1 for k in self.param_names])
        self.param_bounds = np.zeros((self.n_parameters, 2))
        self.gp_nugget = gp_nugget
        self.verbose = verbose
        self.scoring = scoring
        self.estimator = estimator
        self.fit_params = fit_params if fit_params is not None else {}
        self.cv = cv
        self.X = X
        self.y = y

        self.best_parameter_ = None
        self.tested_parameters_ = None
        self.cv_scores_ = None

        if(callable(estimator)):
            self._callable_estimator = True
            if(verbose):
                print('Estimator is a callable and not an sklearn Estimator')
        else:
            self._callable_estimator = False

        if not self._callable_estimator:
            self.scorer_ = check_scoring(self.estimator, scoring=self.scoring)

        # init param_bounds
        for i in range(self.n_parameters):
            if(parameters[self.param_names[i]][0] == 'cat'):
                self.param_bounds[i, 0] = 0
                self.param_bounds[i, 1] = \
                    len(parameters[self.param_names[i]][1])
            else:
                self.param_bounds[i] = \
                    np.array(parameters[self.param_names[i]][1])
                if(parameters[self.param_names[i]][0] == 'int'):
                    self.param_bounds[i, 1] += 1

        if(self.verbose):
            print(self.parameters)
            print(self.param_names)
            print(self.param_isInt)
            print(self.param_bounds)

    def vector_to_dict(self, vector_parameter):
        dict_parameter = dict.fromkeys(self.param_names)
        for i in range(self.n_parameters):
            if(self.parameters[self.param_names[i]][0] == 'cat'):
                dict_parameter[self.param_names[i]] = \
                    (self.parameters[self.param_names[i]][1])[
                        int(vector_parameter[i])]
            elif(self.parameters[self.param_names[i]][0] == 'int'):
                dict_parameter[self.param_names[i]] = int(vector_parameter[i])
            else:
                dict_parameter[self.param_names[i]] = vector_parameter[i]

        return dict_parameter

    def score(self, test_parameter):
        """
        The score function to call in order to evaluate the quality
        of the parameter test_parameter

        Parameters
        ----------
        tested_parameter : dict, the parameter to test

        Returns
        -------
        score : the mean of the CV score
        """

        if not self._callable_estimator:
            cv = check_cv(self.cv, self.X, self.y,
                          classifier=is_classifier(self.estimator))
            cv_score = [_fit_and_score(
                clone(self.estimator), self.X, self.y, self.scorer_,
                train, test, False, test_parameter,
                self.fit_params, return_parameters=True)
                for train, test in cv]

            n_test_samples = 0
            score = 0
            for tmp_score, tmp_n_test_samples, _, _ in cv_score:
                tmp_score *= tmp_n_test_samples
                n_test_samples += tmp_n_test_samples
                score += tmp_score
            score /= float(n_test_samples)

        else:
            score = self.estimator(test_parameter)

        return score

    def _fit(self):
        """
        Run the hyper-parameter optimization process

        Returns
        -------
        tested_parameters_ : ndarray, the parameters tested during the process

        cv_scores_ : array of the mean CV results of the parameters tested
        """

        n_tested_parameters = 0
        tested_parameters = np.zeros((self.n_iter, self.n_parameters))
        cv_scores = np.zeros(self.n_iter)

        #  Initialize with random candidates  #
        init_candidates = sample_candidates(
            self.n_init, self.param_bounds, self.param_isInt)
        self.n_init = init_candidates.shape[0]

        for i in range(self.n_init):
            dict_candidate = self.vector_to_dict(init_candidates[i, :])
            cv_score = self.score(dict_candidate)

            if(self.verbose):
                print ('Step ' + str(i) + ' - Hyperparameter '
                       + str(dict_candidate) + ' ' + str(cv_score))

            is_in, idx = is_in_ndarray(
                init_candidates[i, :],
                tested_parameters[:n_tested_parameters, :])
            if not is_in:
                tested_parameters[n_tested_parameters, :] = \
                    init_candidates[i, :]
                cv_scores[n_tested_parameters] = cv_score
                n_tested_parameters += 1
            else:
                if(self.verbose):
                    print('Hyperparameter already tesed')
                cv_scores[idx] = (cv_scores[idx] + cv_score) / 2.

        for i in range(self.n_iter-self.n_init):

            # Model with a Gaussian Process
            gp = GaussianProcess(theta0=1. * np.ones(self.n_parameters),
                                 thetaL=0.001 * np.ones(self.n_parameters),
                                 thetaU=10. * np.ones(self.n_parameters),
                                 random_start=3,
                                 nugget=self.gp_nugget)
            gp.fit(tested_parameters[:n_tested_parameters, :],
                   cv_scores[:n_tested_parameters])

            # Sample candidates and predict their corresponding
            # acquisition values
            candidates = sample_candidates(self.n_candidates,
                                           self.param_bounds,
                                           self.param_isInt)
            if(self.acquisition_function == 'UCB'):
                predictions, MSE = gp.predict(candidates, eval_MSE=True)
                upperBound = predictions + 1.96*np.sqrt(MSE)
                best_candidate = candidates[np.argmax(upperBound)]

            elif(self.acquisition_function == 'EI'):
                predictions, MSE = gp.predict(candidates, eval_MSE=True)
                y_best = np.max(cv_scores)
                ei = compute_ei(predictions, np.sqrt(MSE), y_best)
                best_candidate = candidates[np.argmax(ei)]

            else:
                print('WARNING : acquisition_function not implemented yet : '
                      + self.acquisition_function)

            dict_candidate = self.vector_to_dict(best_candidate)
            cv_score = self.score(dict_candidate)
            if(self.verbose):
                print ('Step ' + str(i+self.n_init) + ' - Hyperparameter '
                       + str(dict_candidate) + ' ' + str(cv_score))

            is_in, idx = is_in_ndarray(
                best_candidate,
                tested_parameters[:n_tested_parameters, :])
            if not is_in:
                tested_parameters[n_tested_parameters, :] = best_candidate
                cv_scores[n_tested_parameters] = cv_score
                n_tested_parameters += 1
            else:
                if(self.verbose):
                    print('Hyperparameter already tesed')
                cv_scores[idx] = (cv_scores[idx] + cv_score) / 2.

        best_idx = np.argmax(cv_scores[:n_tested_parameters])
        vector_best_param = tested_parameters[best_idx]
        best_parameter = self.vector_to_dict(vector_best_param)

        # store
        self.best_parameter_ = best_parameter
        self.tested_parameters_ = tested_parameters[:n_tested_parameters, :]
        self.cv_scores_ = cv_scores[:n_tested_parameters]

        if(self.verbose):
            print ('\nTested ' + str(n_tested_parameters) + ' parameters')
            print ('Max cv score ' + str(cv_scores[best_idx]))
            print ('Best parameter ' + str(tested_parameters[best_idx]))
            print (best_parameter)

        return tested_parameters[:n_tested_parameters, :], \
            cv_scores[:n_tested_parameters]
