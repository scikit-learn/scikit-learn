'''\
A thin wrapper for OpenCV classifiers to sklearn.
The complete list of supported OpenCV models is below:
* Boosted random forest (see cv2.ml.Boost_create function)
* Random forest (see cv2.ml.RTrees_create function)
* SVM (see cv2.ml.SVM_create function)
* KNN (see cv2.ml.KNearest_create function)

== EXAMPLE ==
    import opencv_sklearn_wrapper
    from sklearn.metrics import classification_report
    from sklearn.model_selection import class_val_predict

    model = cv2.ml.Boost_create()
    classifier = make_sklearn_wrapper(model)
    ...
    predictions = cross_val_predict(classifier, dataset, labels, cv=5)
    print(classification_report(labels, predictions))

== TESTING ==
* Use command one-liner below to run the unit tests:
    python opencv_sklearn_wrapper.py -v

* [07/14/2017] Tested for sklearn 0.18.1 + both of Python 2.7.12 and Python 3.5.2.
'''

import cv2
import math
import numpy as np

import unittest
import warnings

from sklearn.base import BaseEstimator


class Opencv_basic_predictor:
    '''Provides the basic prediction logics (valid for Opencv Boost/Random Forest/SVM).'''
    _predict__flags = None
    _predict_proba__flags = None


    def _to_posteriors(self, output):
        p = 1.0 / (1.0 + math.exp(-2.0 * output))
        return [1.0 - p, p]

    def predict(self, features):
        if self._predict__flags:
            _, predictions = self.opencv_model.predict(features, flags=self._predict__flags)
        else:
            _, predictions = self.opencv_model.predict(features)
        return predictions.ravel()

    def predict_proba(self, features):
        if self._predict_proba__flags:
            _, predictions = self.opencv_model.predict(features, flags=self._predict_proba__flags)
        else:
            _, predictions = self.opencv_model.predict(features)
        predictions = np.array(map(self._to_posteriors, predictions), np.float32)
        return predictions

    pass


class Opencv_sklearn_wrapper(BaseEstimator, Opencv_basic_predictor):
    '''Base class for all Opencv wrapper to scikit-learn.
    Any wrapper implementations for Opencv models must be derived.

    Technical Notes
    ---------------
    * Member classes_ and method predict_proba are required by sklearn.calibration.CalibratedClassifierCV.
    * This __init__ is required by sklearn.base.clone function to reconstruct the estimator by its parameters passed as
      explicit keyword arguments. As well, it has to store the input model parameters (including model_class) because
      the sklearn.base.clone checks not for parameter equivalence but for object exactness.
    '''

    __slots__ = ['opencv_model', 'classes_']
    _estimator_type = 'classifier'  # To pass sklearn.base.is_classifier check.


    def __init__(self, trained_model=None, model_class=None, **kwargs):
        if not model_class and not trained_model:
            raise RuntimeError('Cannot create a wrapper for null model.')

        self.__params = kwargs
        self.__class_of_opencv_model = trained_model.__class__.__name__ if not model_class else model_class
        self.opencv_model = trained_model

        if not self.opencv_model:
            _, name = model_class.split('ml_', 1)
            self.opencv_model = getattr(cv2.ml, name + '_create')()
        for k, v in self.__params.items():
            if k != 'model_class':
                getattr(self.opencv_model, 'set' + k)(v)
        return


    def fit(self, features, labels, sample_weight=None):
        _, dim = features.shape
        var_types = np.array([cv2.ml.VAR_NUMERICAL] * dim + [cv2.ml.VAR_CATEGORICAL],
                             np.uint8)
        td = cv2.ml.TrainData_create(np.array(features, np.float32),
                                     cv2.ml.ROW_SAMPLE, labels,
                                     sampleWeights=sample_weight,
                                     varType=var_types)
        self.opencv_model.train(td)
        self.classes_ =  np.unique(labels)
        return self


    def get_params(self, deep=True):
        params = {'model_class': self.__class_of_opencv_model}
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
    pass



class Opencv_sklearn_wrapper__dtree(Opencv_sklearn_wrapper):
    '''Opencv-to-sklearn wrapper for Random Forest/Boosted decision trees.'''
    _predict__flags = cv2.ml.DTREES_PREDICT_MAX_VOTE
    _predict_proba__flags = cv2.ml.DTREES_PREDICT_SUM
    pass


class Opencv_sklearn_wrapper__knn(Opencv_sklearn_wrapper):
    '''Opencv-to-sklearn wrapper for KNN.'''
    def __init__(self, model_class='ml_KNearest', **kwargs):
        super(Opencv_sklearn_wrapper__knn, self).__init__(model_class=model_class, **kwargs)
        return


    def predict(self, features):
        _, results, _, _ = self.opencv_model.findNearest(features, k=self.opencv_model.getDefaultK())
        return results.ravel()


    def predict_proba(self, features):
        _, results, _, _ = self.opencv_model.findNearest(features, k=self.opencv_model.getDefaultK())
        pass
    pass


def make_sklearn_wrapper(model=None, class_name=None):
    '''Constructs the wrapper for the Opencv's binary classifiers.
    Technically, this implements a simple dispatcher over the Opencv's entire classifier set.

    Parameters
    ----------
    model: the instance of any instance of Opencv's binary classifier.
    class_name: the class name of Opencv's model.
    '''
    name = model.__class__.__name__ if model else class_name
    if name in ['ml_Boost', 'ml_RTrees']:
        wrapper = Opencv_sklearn_wrapper__dtree(model_class=name)
    elif name == 'ml_SVM':
        wrapper = Opencv_sklearn_wrapper(model_class=name)
    elif name == 'ml_KNearest':
        wrapper = Opencv_sklearn_wrapper__knn()
    else:
        raise RuntimeError('make_sklearn_wrapper: both of keyword arguments "model" and "class_name" are None or \
the model is unsupported')

    wrapper.opencv_model = model
    return wrapper



class _Opencv_sklearn_wrapper__tests(unittest.TestCase):
    _supported_models = ['ml_Boost', 'ml_RTrees', 'ml_SVM', 'ml_KNearest']


    def test__sklearn_base(self):
        from sklearn.base import clone, is_classifier

        for model in self._supported_models:
            wrapper = Opencv_sklearn_wrapper(model_class=model)

            self.assertTrue(issubclass(wrapper.__class__, BaseEstimator))
            self.assertTrue(is_classifier(wrapper))
            with warnings.catch_warnings(record=True) as raised_warnings:
                warnings.simplefilter('always')

                try:
                    clone(wrapper)
                except RuntimeError:
                    self.fail('Opencv_sklearn_wrapper does not meet the sklearn.base.* requirements')

                self.assertFalse(any([issubclass(w.category, DeprecationWarning) for w in raised_warnings]),
                                 msg='sklearn.base.* raises the DeprecationWarnings')
        return


    def test_sklearn_model_selection(self):
        from sklearn import datasets
        from sklearn.metrics import classification_report
        from sklearn.model_selection import cross_val_predict

        iris = datasets.load_iris()
        samples, labels = np.array(iris.data[:, :2], np.float32), iris.target

        warnings.filterwarnings('ignore')
        try:
            for model in self._supported_models:
                wrapper = Opencv_sklearn_wrapper(model_class=model)

                predictions = cross_val_predict(wrapper, samples, labels)
                report = classification_report(labels, predictions)
        except RuntimeError:
            self.fail('Opencv_sklearn_wrapper crashes for sklearn.metrics.* and sklearn.model_selection.*')

    pass


if __name__ == '__main__':
    unittest.main()
