# -*- coding: utf-8 -*-

"""
The :mod:`sklearn.opencv_wrapper` module implements a thin wrapper for OpenCV classifiers listed below:
* (Boosted) random forest model
* SVM
* KNN
"""

# Author: Piotr Semenov <piotr@dropoutlabs.com>
# License: BSD 3 clause

import cv2
import numpy as np

from sklearn.base import BaseEstimator


__all__ = ['wrap']

class _Opencv_predictor_base:
    """Base class for scikit-learn compatible opencv models."""
    _positive_class_idx = 1

    def _estimate(self, features, flags=None):
        _, outputs = self.opencv_model.predict(features, flags=flags)

        return outputs.ravel()

    def predict(self, features):
        """Returns the hard prediction."""

        return self._estimate(features, flags=self._force_hard_decision)

    def decision_function(self, features):
        """Returns the soft decision."""

        return self._estimate(features, flags=self._force_soft_decision)


class _Opencv_predictor(BaseEstimator, _Opencv_predictor_base):
    """Scikit-learn compatible wrapper for opencv models."""
    __slots__ = ['opencv_model', 'classes_']
    _estimator_type = 'classifier'

    def __init__(self, model=None, class_name=None, **kwargs):
        assert model or class_name

        self.__params = kwargs
        self.__class_of_opencv_model = model.__class__.__name__ if not class_name else class_name
        self.opencv_model = model

        if not self.opencv_model:
            _, name = class_name.split('ml_', 1)
            self.opencv_model = getattr(cv2.ml, name + '_create')()
        for k, v in self.__params.items():
            if k != 'class_name':
                getattr(self.opencv_model, 'set' + k)(v)

    def fit(self, features, labels, sample_weight=None):
        _, dim = features.shape
        var_types = np.array([cv2.ml.VAR_NUMERICAL] * dim + [cv2.ml.VAR_CATEGORICAL],
                             np.uint8)
        td = cv2.ml.TrainData_create(np.array(features, np.float32),
                                     cv2.ml.ROW_SAMPLE, labels,
                                     sampleWeights=sample_weight,
                                     varType=var_types)
        self.opencv_model.train(td)
        self.classes_ = np.unique(labels)
        return self

    def get_params(self, deep=True):
        params = {'class_name': self.__class_of_opencv_model}
        if not self.__params:
            self.__params = {}

            attributes = dir(self.opencv_model)
            getter_names = [k[3:] for k in attributes if k.startswith('get')]
            setter_names = [k[3:] for k in attributes if k.startswith('set')]
            model_parameters = set(getter_names).intersection(setter_names)

            for name in model_parameters:
                self.__params.update({name: getattr(self.opencv_model, 'get' + name)()})

        params.update(self.__params)
        return params


class _Opencv_dtree(_Opencv_predictor):
    _force_hard_decision = cv2.ml.DTREES_PREDICT_MAX_VOTE
    _force_soft_decision = cv2.ml.DTREES_PREDICT_SUM


class _Opencv_svm(_Opencv_predictor):
    _force_hard_decision = False
    _force_soft_decision = True

    def __init__(self, class_name='ml_SVM', **kwargs):
        super(_Opencv_svm, self).__init__(class_name='ml_SVM', **kwargs)


class _Opencv_knn(_Opencv_predictor):
    _force_hard_decision = True
    _force_soft_decision = False

    def __init__(self, class_name='ml_KNearest', **kwargs):
        super(_Opencv_knn, self).__init__(class_name='ml_KNearest', **kwargs)

    def _estimate(self, features, flags=None):
        if flags:
            _, outputs = self.opencv_model.predict(features)
        else:
            _, _, neighbour_responses, _ = self.opencv_model.findNearest(features,
                                                                         k=self.opencv_model.getDefaultK())
            outputs = np.apply_along_axis(lambda rs: len(rs[rs == self._positive_class_idx]) / np.float32(len(rs)),
                                          1,
                                          neighbour_responses.astype(np.uint32))
        return outputs.ravel()


def wrap(model=None, class_name=None):
    """Wraps the Opencv's classifiers."""
    name = model.__class__.__name__ if model else class_name

    if name in ['ml_Boost', 'ml_RTrees']:
        wrapper = _Opencv_dtree(class_name=name)
    elif name == 'ml_SVM':
        wrapper = _Opencv_svm()
    elif name == 'ml_KNearest':
        wrapper = _Opencv_knn()
    else:
        raise AssertionError('{} is not supported.'.format(name))

    wrapper.opencv_model = model
    return wrapper
