"""
Test the opencv_wrapper module.
"""

import warnings
import numpy as np

from sklearn.base import BaseEstimator
from sklearn.opencv_wrapper import wrap

from sklearn.utils.testing import assert_true, assert_false


_supported_models = ['ml_Boost', 'ml_RTrees', 'ml_SVM', 'ml_KNearest']

def test_sklearn_base_functions():
    from sklearn.base import clone, is_classifier

    for model in _supported_models:
        wrapper = wrap(class_name=model)

        assert_true(issubclass(wrapper.__class__, BaseEstimator))
        assert_true(is_classifier(wrapper))
        with warnings.catch_warnings(record=True) as raised_warnings:
            warnings.simplefilter("always")
            try:
                clone(wrapper)
            except RuntimeError:
                assert_true(False, msg="sklearn.opencv_wrapper does not meet the sklearn.base.clone requirements.")

                assert_false(any([issubclass(w.category, DeprecationWarning) for w in raised_warnings]),
                             msg="sklearn.base.* raises the DeprecationWarnings")

def test_sklearn_model_selection():
    from sklearn import datasets
    from sklearn.metrics import classification_report
    from sklearn.model_selection import cross_val_predict

    iris = datasets.load_iris()
    samples, labels = np.array(iris.data[:, :2], np.float32), iris.target

    try:
        for model in _supported_models:
            wrapper = wrap(class_name=model)

            predictions = cross_val_predict(wrapper, samples, labels)
            report = classification_report(labels, predictions)
    except RuntimeError:
        assert_true(False, msg="Opencv_sklearn_wrapper crashes for sklearn.metrics.* and sklearn.model_selection.*")
