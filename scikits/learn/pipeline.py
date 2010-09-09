"""
Pipeline: chain transforms and estimators to build a composite estimator.
"""
# Author: Edouard Duchesnay
#         Gael Varoquaux
#         Virgile Fritsch
# Licence: BSD
import copy

import numpy as np

from .base import BaseEstimator


def _set_valid_params(estimator, params):
    """ Set the params of the estimator that it accepts.

        The dictionnary params is modified: the valid parameters are
        removed.
    """
    these_params_name = set(params.keys()
                    ).intersection(estimator._get_param_names())
    these_params = dict()
    for name in these_params_name:
        these_params[name] = params.pop(name)
    estimator._set_params(**these_params)
    return params


class Pipeline(BaseEstimator):
    """ Pipeline of transforms with a final estimator 

        Sequentialy apply a list of transforms and a final estimator 
        A transform implements fit & transform methods
        A estimator implements fit & predict methods

        Example
        =======

        >>> from scikits.learn import svm, datasets
        >>> from scikits.learn.datasets import samples_generator
        >>> from scikits.learn.feature_selection import SelectKBest, f_regression
        >>> from scikits.learn.pipeline import Pipeline

        >>> # generate some data to play with
        >>> X, y = samples_generator.test_dataset_classif(k=5)

        >>> # ANOVA SVM-C
        >>> anova_filter = SelectKBest(f_regression, k=5)
        >>> clf = svm.SVC(kernel='linear')

        >>> anova_svm = Pipeline([anova_filter], clf)
        >>> _ = anova_svm.fit(X,y)

        >>> prediction = anova_svm.predict(X)
        >>> score = anova_svm.score(X)
    """

    #---------------------------------------------------------------------------
    # BaseEstimator interface
    #---------------------------------------------------------------------------

    def __init__(self, transforms=[], estimator=None, names=None):
        """
        Parameters
        ==========
        transforms: list
            List of various transform object (implementing
            fit/transform) that are chained, in the order in which
            they are chained.
        estimator: estimator object
            Object implementing fit and possibly predit or score
            that is fit with the transformed data
        names: list of names, optional
            Names that are attributed to the various transforms and
            the final estimator. This is optional if the all the
            different steps have non-overlapping parameter names.
        """
        params_names = set()
        for t in  transforms:
            assert hasattr(t, "fit") and hasattr(t, "transform"), ValueError(
                "All transforms should implement fit and transform",
                "'%s' (type %s) )" % (t, type(t))
            )
            if names is None:
                these_params = t._get_param_names()
                overlap = params_names.intersection(these_params)
                assert not overlap, \
                    ValueError(
                        "Overlapping parameters: %s, you need to "
                        "provide explicite names to the Pipeline "
                        "constructor" % list(overlap)
                    )
                params_names.update(these_params)
        assert hasattr(estimator, "fit") and hasattr(estimator, "predict"), \
            ("Predictor should implement fit and predict",
                "'%s' (type %s) )" % (estimator, type(estimator))
            )
        self.transforms = transforms
        named_steps = None
        if names is not None:
            assert len(set(names)) == len(names), ValueError(
                    'Names should be unique, %s supplied' % names
                )
            steps = copy.copy(transforms)
            steps.append(estimator)
            named_steps = dict()
            for step, name in zip(steps, names):
                named_steps[name] = step
        self.names = names
        self._named_steps = named_steps
        self.estimator = estimator


    def _get_param_names(self):
        args = list()
        if self._named_steps is None:
            for transform in self.transforms:
                args.extend(transform._get_param_names())
            args.extend(self.estimator._get_param_names())
        else:
            for name, step in self._named_steps.iteritems():
                args.extend(['%s_%s' (name, arg) 
                            for arg in step._get_param_names()])
        return args


    def _get_params(self):
        out = dict()
        if self._named_steps is None:
            for transform in self.transforms:
                out.update(transform._get_params())
            out.update(self.estimator._get_params())
        else:
            for name, step in self._named_steps.iteritems():
                for key, value in step._get_params().iteritems():
                    out['%s_%s' % (name, key)] = value
        return out


    def _set_params(self, **params):
        """ Set the parameters of the estimator.
        """
        if self._named_steps is None:
            for transform in self.transforms:
                params = _set_valid_params(transform, params)
            params = _set_valid_params(self.estimator, params)
        else:
            for name, step in self._named_steps.iteritems():
                n_chars = len(name) + 1
                these_params = dict()
                for key in params.copy():
                    if key.startswith('%s_' % name):
                        these_params[key[n_chars:]] = params.pop(key)
                step._set_params(**these_params)
        if params:
            raise ValueError('Some invalid parameters passed %s.' 
                                        % params)


    def _reinit(self):
        """Clone the estimator without estimated parameters.
        """
        return self.__class__(
                    transforms=[t._reinit() for t in self.transforms], 
                    estimator=self.estimator._reinit(),
                    names=self.names)


    def __repr__(self):
        class_name = self.__class__.__name__
        params = dict(transforms=self.transforms,
                      estimator=self.estimator,
                      names=self.names)
        if self.names is None:
            params.pop('names')

        # Do a multi-line justified repr:
        params_list = list()
        this_line_length = len(class_name)
        line_sep = ',\n' + (1+len(class_name)/2)*' '
        for i, (k, v) in enumerate(params.iteritems()):
            this_repr  = '%s=%s' % (k, repr(v))
            if i > 0: 
                if (this_line_length + len(this_repr) >= 75
                                            or '\n' in this_repr):
                    params_list.append(line_sep)
                    this_line_length += len(line_sep)
                else:
                    params_list.append(', ')
                    this_line_length += 2
            params_list.append(this_repr)
            this_line_length += len(this_repr)

        params_str = ''.join(params_list)
        return '%s(%s)' % (
                class_name,
                params_str
            )

    
    #---------------------------------------------------------------------------
    # Estimator interface
    #---------------------------------------------------------------------------

    def fit(self, X, y=None):
        Xt = X
        for transformer in self.transforms:
            Xt = transformer.fit(Xt, y).transform(Xt)
        self.estimator.fit(Xt, y)
        return self

    def predict(self, X):
        Xt = X
        for transformer in self.transforms:
            Xt = transformer.transform(Xt)
        return self.estimator.predict(Xt)

    def score(self, X, y=None):
        Xt = X
        for transformer in self.transforms:
            Xt = transformer.transform(Xt)
        return self.estimator.score(Xt, y)

