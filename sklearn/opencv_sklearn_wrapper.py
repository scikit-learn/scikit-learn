'''\
A thin wrapper for OpenCV classifiers to sklearn.
The complete list of supported OpenCV models is below:
* Boosted random forest (created by cv2.ml.Boost_create function)

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


class Opencv_sklearn_wrapper__base(BaseEstimator):
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


    def predict(self, features):
        raise Exception('No base predict method is defined to work for all Opencv models. Define yours own.')


    def predict_proba(self, features):
        raise Exception('No base predict_proba method is defined to work for all Opencv models. Define yours own.')
    pass


class Opencv_sklearn_wrapper__boost(Opencv_sklearn_wrapper__base):
    '''Opencv-to-sklearn wrapper for Boosted random forest.'''
    def _output_to_posteriors(self, sum):
        p = 1.0 / (1.0 + math.exp(-2.0 * sum))
        return [1.0 - p, p]


    def predict(self, features):
        _, predictions = self.opencv_model.predict(features, flags=cv2.ml.DTREES_PREDICT_MAX_VOTE)
        return np.squeeze(np.array(predictions, np.uint32))


    def predict_proba(self, features):
        _, predictions = self.opencv_model.predict(features, flags=cv2.ml.DTREES_PREDICT_SUM)
        predictions = np.array(map(self._output_to_posteriors, predictions), np.float32)
        return predictions
    pass


def make_sklearn_wrapper(model):
    '''Constructs the wrapper for the Opencv's binary classifiers.
    Technically, this implements a simple dispatcher over the Opencv's entire classifier set.

    Parameters
    ----------
    model: the instance of any Opencv's binary classifier.
    '''
    model_class = model.__class__.__name__
    if model_class == 'ml_Boost':
        wrapper = Opencv_sklearn_wrapper__boost(model_class=model_class)

    wrapper.opencv_model = model
    return wrapper



class _Opencv_sklearn_wrapper__tests(unittest.TestCase):
    _supported_models = ['ml_Boost']


    def test__sklearn_base(self):
        from sklearn.base import clone, is_classifier

        for model in self._supported_models:
            wrapper = Opencv_sklearn_wrapper__base(model_class=model)

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
    pass


if __name__ == '__main__':
    unittest.main()
