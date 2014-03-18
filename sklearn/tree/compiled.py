from __future__ import print_function

from ..utils import array2d
from .tree import DTYPE

import _compiled
import code_gen as cg
import numpy as np


class CompiledRegressionPredictor(object):
    """Class to construct a compiled predictor from a previously trained
    ensemble of decision trees.

    Parameters
    ----------

    clf:
      A fitted regression tree/ensemble.

    References
    ----------

    http://courses.cs.washington.edu/courses/cse501/10au/compile-machlearn.pdf
    http://crsouza.blogspot.com/2012/01/decision-trees-in-c.html
    """
    def __init__(self, clf):
        self._n_features, self._evaluator = self._build(clf)

    @classmethod
    def _build(cls, clf):
        from ..ensemble.gradient_boosting import GradientBoostingRegressor
        from ..ensemble.forest import ForestRegressor
        from .tree import DecisionTreeRegressor

        if not cls.compilable(clf):
            raise ValueError("Predictor %s cannot be compiled".format(
                clf.__class__.__name__))

        lines = None
        n_features = None
        if isinstance(clf, DecisionTreeRegressor):
            n_features = clf.n_features_
            lines = cg.code_gen_tree(tree=clf.tree_)

        if isinstance(clf, GradientBoostingRegressor):
            n_features = clf.n_features
            lines = cg.code_gen_ensemble(
                trees=[e.tree_ for e in clf.estimators_.flat],
                individual_learner_weight=clf.learning_rate,
                initial_value=clf._init_decision_function_single)

        if isinstance(clf, ForestRegressor):
            n_features = clf.n_features_
            lines = cg.code_gen_ensemble(
                trees=[e.tree_ for e in clf.estimators_],
                individual_learner_weight=1.0 / clf.n_estimators,
                initial_value=0.0)

        assert n_features is not None
        assert lines is not None

        so_f = cg.compile_code_to_object("\n".join(lines))
        return n_features, _compiled.CompiledPredictor(
            so_f, cg.EVALUATE_FN_NAME)

    @classmethod
    def compilable(cls, clf):
        """
        Verifies that the given fitted model is eligible to be compiled.

        Returns True if the model is eligible, and False otherwise.

        Parameters
        ----------

        clf:
          A fitted regression tree/ensemble.


        """
        from ..ensemble.gradient_boosting import GradientBoostingRegressor
        from ..ensemble.forest import ForestRegressor
        from .tree import DecisionTreeRegressor

        # TODO - is there an established way to check `is_fitted``?
        if isinstance(clf, DecisionTreeRegressor):
            return clf.n_outputs_ == 1 and clf.n_classes_ == 1 \
                and clf.tree_ is not None

        if isinstance(clf, GradientBoostingRegressor):
            return clf.estimators_.size and all(cls.compilable(e)
                                                for e in clf.estimators_.flat)

        if isinstance(clf, ForestRegressor):
            estimators = np.asarray(clf.estimators_)
            return estimators.size and all(cls.compilable(e)
                                           for e in estimators.flat)
        return False

    def predict(self, X):
        """Predict regression target for X.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y: array of shape = [n_samples]
            The predicted values.
        """
        if getattr(X, "dtype", None) != DTYPE or X.ndim != 2:
            X = array2d(X, dtype=DTYPE)

        # TODO - validate n_features is correct?
        n_samples, n_features = X.shape
        if self._n_features != n_features:
            raise ValueError("Number of features of the model must "
                             " match the input. Model n_features is {} and "
                             " input n_features is {}".format(
                                 self.n_features_, n_features))

        result = np.empty(n_samples, dtype=DTYPE)
        return self._evaluator.predict(X, result)
