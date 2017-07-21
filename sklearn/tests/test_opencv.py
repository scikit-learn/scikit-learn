"""
Test the opencv_wrapper module.
"""

import warnings
import numpy as np

from sklearn.opencv_wrapper import wrap

from sklearn.base import BaseEstimator, clone, is_classifier
from sklearn import datasets
from sklearn.metrics import classification_report
from sklearn.ensemble import AdaBoostClassifier
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import cross_val_predict, train_test_split

from sklearn.utils.testing import assert_true, assert_false


_supported_models = ['ml_Boost', 'ml_RTrees', 'ml_SVM', 'ml_KNearest']


def test_sklearn_base_functions():
    for model in _supported_models:
        wrapper = wrap(class_name=model)

        assert_true(issubclass(wrapper.__class__, BaseEstimator))
        assert_true(is_classifier(wrapper))
        with warnings.catch_warnings(record=True) as raised_warnings:
            warnings.simplefilter("always")
            try:
                clone(wrapper)
            except RuntimeError:
                assert_true(False,
                            msg="sklearn.base.clone is failed.")

                assert_false(
                    any([issubclass(w.category, DeprecationWarning)
                        for w in raised_warnings]),
                    msg="sklearn.base.clone raises the DeprecationWarnings")


def test_sklearn_adaptors():
    iris = datasets.load_iris()
    samples, labels = np.array(iris.data[:, :2], np.float32), iris.target
    X_train, X_test, y_train, y_test = train_test_split(samples, labels,
                                                        test_size=0.2)

    try:
        for model in _supported_models:
            wrapper = wrap(class_name=model)
            wrapper.fit(X_train, y_train)

            clf = CalibratedClassifierCV(wrapper)
            clf.fit(X_train, y_train)
            predictions = clf.predict(X_test)

            clf = AdaBoostClassifier(wrapper, algorithm="SAMME")
            clf.fit(X_train, y_train)
            predictions = clf.predict(X_test)

    except RuntimeError:
        assert_true(False,
                    msg="sklearn.ensemble.* is failed.")


def test_sklearn_cv():
    iris = datasets.load_iris()
    samples, labels = np.array(iris.data[:, :2], np.float32), iris.target
    X_train, X_test, y_train, y_test = train_test_split(samples, labels,
                                                        test_size=0.2)

    try:
        for model in _supported_models:
            wrapper = wrap(class_name=model)
            wrapper.fit(X_train, y_train)

            predictions = cross_val_predict(wrapper, X_train, y_train, cv=5)
            report = classification_report(y_train, predictions)

    except RuntimeError:
        assert_true(False,
                    msg="sklearn.model_selection.cross_val_predict is failed.")
