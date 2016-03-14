"""
Ensemble Selection Classifier

This module contains a Ensemble Selection Classifier for
fitted classification estimators.

"""

from collections import Counter

import operator
import numpy as np

from ..base import BaseEstimator
from ..base import ClassifierMixin
from ..utils import check_random_state
from ..metrics import accuracy_score, f1_score
from sklearn.utils.validation import check_is_fitted


def _accuracy(y, proba):
    return accuracy_score(y, np.argmax(proba, axis=1))

def _f1(y, proba):
    return f1_score(y, np.argmax(proba, axis=1))


class EnsembleSelectionClassifier(BaseEstimator, ClassifierMixin):
    """An ensemble classifier built by greedy stepwise selection.

    Parameters
    ----------
    estimators : list of (string, estimator) tuples
        The estimators from which the ensemble selection classifier is built.
        These estimators must be fitted.

    max_bag_estimators : integer, optional (default=50)
        The maximum number of estimators at each bag ensemble selection
        is terminated.

    bag_fraction : float, optional (default=0.5)
        Fraction of (post-pruning) models to randomly select for each bag
        ensemble selection.

    n_bags : integer, optional (default=20)
        Number of bag ensemble selection.

    n_best : integer, optional (default=1)
        Number of the best estimators to be added into an ensemble
        before each bag ensemble selection.

    prune_fraction : float, optional (default=0)
        Fraction of worst models to prune before each bag ensemble selection.

    score_metric : {'accuracy', 'f1'}, optional (default='accuracy')
        If 'accuracy' then use accuracy_score as the metric of ensemble
        selection;
        if 'f1' then use f1_score as the metric of ensemble selection;

        # TODO: Add more score_metric

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Attributes
    ----------
    ensemble_ : counter of classifiers
        The collection of fitted estimators and its weight.

    classes_ : array-like, shape = [n_classes_]
        The classes labels.

    n_classes_ : int
        The number of classes.

    References
    ----------
    .. [1] R. Caruana, A. Niculescu-Mizil, G. Crew, A. Ksikes, "Ensemble
           Selection from Libraries of Models", 2004.

    .. [2] R. Caruana, A. Munson, A. Niculescu-Mizil, "Getting the Most Out
           of Ensemble Selection", 2006.

    """

    _metrics = {
        'f1': _f1,
        'accuracy': _accuracy,
    }

    def __init__(self, estimators=None, max_bag_estimators=50, bag_fraction=0.5,
                 n_bags=20, n_best=1, score_metric='accuracy',
                 prune_fraction=0, random_state=None, verbose=False):
        self.estimators = estimators
        self.named_estimators = dict(estimators)
        self.max_bag_estimators = max_bag_estimators
        self.bag_fraction = bag_fraction
        self.n_bags = n_bags
        self.n_best = n_best
        self.score_metric = score_metric
        self.prune_fraction = prune_fraction
        self.random_state = random_state
        self.verbose = verbose

    def fit(self, X, y):
        """ Conduct ensemble selection on the validation set (X, y).

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Validation vectors, where n_samples is the number of samples and
            n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values.

        Returns
        -------
        self : object
        """
        if isinstance(y, np.ndarray) and len(y.shape) > 1 and y.shape[1] > 1:
            raise NotImplementedError('Multilabel and multi-output'
                                      ' classification is not supported.')

        metric_names = self._metrics.keys()
        if self.score_metric not in metric_names:
            raise ValueError('score_metric %s not in %s' %
                             (self.score_metric, metric_names))

        if self.estimators is None or len(self.estimators) == 0:
            raise AttributeError('Invalid `estimators` attribute, `estimators`'
                                 ' should be a list of'
                                 ' (string, fitted_estimator) tuples')

        self.score_metric_ = self._metrics[self.score_metric]

        # store each clf's score and predict_proba result on validation set
        self.clfs_score_ = {}
        self.clfs_proba_ = {}
        for (name, clf) in self.estimators:
            score, proba = self._get_score_of_model(X, y, clf)
            self.clfs_score_[name] = score
            self.clfs_proba_[name] = proba

        self.classes_ = np.unique(y)
        self.n_classes_ = len(self.classes_)

        self._get_final_ensemble(X, y)

    def predict(self, X):
        n_clfs = sum(self.ensemble_.values())
        ens_proba = np.zeros((X.shape[0], self.n_classes_))

        for estimator_name in self.ensemble_:
            weight = self.ensemble_[estimator_name]
            proba = self.named_estimators[estimator_name].predict_proba(X)
            ens_proba += proba * weight

        ens_proba /= n_clfs
        return np.argmax(ens_proba, axis=1)

    # TODO:
    #def predict_proba(self, X):


    def _get_score_of_model(self, X, y, estimator):
        proba = estimator.predict_proba(X)
        score = self.score_metric_(y, proba)
        return score, proba

    def _get_score_of_ensemble(self, X, y, ensemble):
        n_clfs = sum(ensemble.values())
        ens_proba = np.zeros((X.shape[0], self.n_classes_))

        for estimator_name in ensemble:
            weight = ensemble[estimator_name]
            proba = self.clfs_proba_[estimator_name]
            ens_proba += proba * weight

        ens_proba /= n_clfs
        ens_score = self.score_metric_(y, ens_proba)
        return ens_score, ens_proba

    def _build_ensemble(self, X, y, ranked_estimators_name):
        ensemble = Counter(ranked_estimators_name[:self.n_best])
        ens_score, ens_proba = self._get_score_of_ensemble(X, y, ensemble)
        n_clfs = sum(ensemble.values())

        cand_ensembles = []
        while(n_clfs < self.max_bag_estimators):
            new_ens_scores = []
            for name, new_clf in self.estimators:
                new_clf_proba = self.clfs_proba_[name]

                new_ens_proba = ((ens_proba * n_clfs + new_clf_proba)
                                 / (n_clfs + 1))
                new_ens_score = self.score_metric_(y, new_ens_proba)
                new_ens_scores.append({'score':new_ens_score,
                                       'new_estimator_name':name})

            best_ens_score = max(new_ens_scores, key=lambda x: x['score'])

            ens_score = best_ens_score['score']
            new_estimator_name = best_ens_score['new_estimator_name']
            ensemble.update({new_estimator_name: 1})

            # store current ensemble to select best later
            ens_copy = Counter(ensemble)
            cand = {'ens': ens_copy, 'score': ens_score}
            cand_ensembles.append(cand)

            n_clfs = sum(ensemble.values())

        if (n_clfs == self.max_bag_estimators):
            best_cand_ensemble = max(cand_ensembles, key = lambda x: x['score'])
            ensemble = best_cand_ensemble['ens']

        return ensemble

    def _get_final_ensemble(self, X, y):

        n_clfs = int(len(self.estimators) * (1 - self.prune_fraction))
        bag_size = int(self.bag_fraction * n_clfs)

        ranked_clfs_score = sorted(self.clfs_score_.items(),
                                   key=operator.itemgetter(1), reverse=True)
        ranked_clfs_name = map(lambda (name, clf_score): name,
                                     ranked_clfs_score)

        if (self.verbose):
            print "Ranked Model:"
            for (i, model_score) in enumerate(ranked_clfs_score):
                print "%d. " % i, model_score

        bag_ensembles = []
        rs = check_random_state(self.random_state)
        for i in xrange(self.n_bags):
            # get bag_size elements at random
            cand_indices = sorted(rs.permutation(n_clfs)[:bag_size])
            cand_ranked_clfs_name = [ranked_clfs_name[ci]
                                     for ci in cand_indices]

            bag_ensemble = self._build_ensemble(X, y, cand_ranked_clfs_name)
            bag_ensembles.append(bag_ensemble)

        # combine each ensemble built by bag ensemble selection
        self.ensemble_ = Counter()
        for bag_ensemble in bag_ensembles:
            self.ensemble_ += bag_ensemble

        # average the result of each bag ensemble selection
        for clf in self.ensemble_:
            self.ensemble_[clf] /= float(self.n_bags)

        if (self.verbose):
            print "Final Ensemble:"
            for model in self.ensemble_:
                print self.ensemble_[model], "of", model
