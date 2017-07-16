'''\
A thin wrapper for Opencv classifiers to scikit-learn.
The complete list of supported OpenCV models is below:
* Boosted Random Forest model created by cv2.ml.Boost_create()
* Random Forest model created by cv2.ml.RTrees_create()
* SVM model created by cv2.ml.SVM_create()
* KNN model created by cv2.ml.KNearest_create()

Note, no prior assumptions about classifier's algorithm behind is made.

== MOTIVATION ==
Want to get the mix of scikit-learn framework power (e.g. cross-validation, posterior calibration) and
the efortless deploy to C++ for the trained Opencv models.

Unfortunately, the inheritance from scikit-learn models (like KNeighborsClassifier) is not an
option due the reasons listed below:
* A huge number of methods to override.
* Method-set to override distinguishes from a model to a model.


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

* [07/14/2017] Tested for sklearn 0.18.1 + Opencv 3.1.0 + both of Python 2.7.12 and Python 3.5.2.
'''

import cv2
import numpy as np

import unittest
import warnings

from sklearn.base import BaseEstimator


class Opencv_binary_predictor:
    '''Provides the hard/soft decisions from the Opencv's binary classifiers.

    Technical Notes
    ---------------
    * Scikit-learn assumption: positive category is assumed to be at index 1 among the probabilities
    delivered from predictor for the input sample.
    Proof: https://github.com/scikit-learn/scikit-learn/blob/8570622a4450fe1ee3c683601454a7189dcccc14/sklearn/calibration.py#L293.
    '''
    _positive_class_idx = 1


    def _estimate(self, features, flags=None):
        '''Provides the basic prediction logics (valid for Opencv's (Boosted) Random Forest/SVM).'''
        if flags is None:
            raise RuntimeError('No flags provided to distinguish the hard/soft decision mode.')

        _, outputs = self.opencv_model.predict(features, flags=flags)
        return outputs.ravel()


    def predict(self, features):
        '''Returns the hard classification decision (aka class label).

        Parameters
        ----------
        features: the samples to classify. Type: numpy.ndarray of shape (n, m).
        '''
        return self._estimate(features, self._force_hard_decision__flags)


    def decision_function(self, features):
        '''Returns the soft classification decision.
         Technically, it must utilizy the classifier's "raw" output like listed below:
         * Random Forest: votes portion casted for the positive category.
         * SVM: oriented distance to the separating hyperplane.
           See http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html#sklearn.svm.SVC.decision_function.
         * Boosted: the weighted sum of the weak decisions (-1 for negative category, +1 for positive one).
         * KNN: neighbour's votes portion casted for the positive category.
         * Logistic Regression: the approximated posterior.
           See http://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html#sklearn.linear_model.LogisticRegression.decision_function.
         * Normal Bayes:

        Parameters
        ----------
        features: the samples to classify. Type: numpy.ndarray of shape (n, m).
        '''
        return self._estimate(features, self._force_soft_decision__flags)
    pass


class Opencv_sklearn_wrapper(BaseEstimator, Opencv_binary_predictor):
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
        self.classes_ = np.unique(labels)
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
    '''Wraps the Random Forest/Boosted decision trees.'''
    _force_hard_decision__flags = cv2.ml.DTREES_PREDICT_MAX_VOTE
    _force_soft_decision__flags = cv2.ml.DTREES_PREDICT_SUM
    pass


class Opencv_sklearn_wrapper__svm(Opencv_sklearn_wrapper):
    '''Wraps the SVM.'''
    _force_hard_decision__flags = False
    _force_soft_decision__flags = True


    def __init__(self, model_class='ml_SVM', **kwargs):
        super(Opencv_sklearn_wrapper__svm, self).__init__(model_class=model_class, **kwargs)
        return
    pass


class Opencv_sklearn_wrapper__knn(Opencv_sklearn_wrapper):
    '''Opencv-to-sklearn wrapper for KNN.'''
    _force_hard_decision__flags = True
    _force_soft_decision__flags = False

    def __init__(self, model_class='ml_KNearest', **kwargs):
        super(Opencv_sklearn_wrapper__knn, self).__init__(model_class=model_class, **kwargs)
        return


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
    pass


class Opencv_sklearn_wrapper__ann(Opencv_sklearn_wrapper):
    '''Opencv-to-sklearn wrapper for the neural networks.'''
    _force_hard_decision__flags = True
    _force_soft_decision__flags = False

    def __init__(self, model_class='ml_ANN_MLP', **kwargs):
        super(Opencv_sklearn_wrapper__ann, self).__init__(model_class=model_class, **kwargs)
        return


    def fit(self, features, labels, sample_weight=None):
        n = len(labels)
        m = np.zeros(n * 2, np.float32)
        index_set = np.uint(labels + np.arange(n) * 2)
        m[index_set] = 1.0
        self.opencv_model.train(features,
                                cv2.ml.ROW_SAMPLE,
                                m.reshape(-1, 2))
        return self
#        return super(Opencv_sklearn_wrapper__ann, self).fit(features,
#                                                            labels.reshape(-1, 1).astype(np.float32),
#                                                            sample_weight)

    def _estimate(self, features, flags=None):
        _, outputs = self.opencv_model.predict(features)
        print(outputs)
        # if flags:
        #     outputs = np.apply_along_axis(lambda rs: rs.argmax(-1),
        #                                   1,
        #                                   outputs)
        # else:
        #     outputs = outputs[:, self._positive_class_idx]
        # return outputs.ravel()
    pass


class Opencv_sklearn_wrapper_nbc(Opencv_sklearn_wrapper):
    '''Opencv-to-sklearn wrapper for Normal Bayesian classifier.'''
    _force_hard_decision__flags = True
    _force_soft_decision__flags = False
    pass


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
        wrapper = Opencv_sklearn_wrapper__svm(model_class=name)
    elif name == 'ml_KNearest':
        wrapper = Opencv_sklearn_wrapper__knn()
    elif name == 'ml_ANN_MLP':
        wrapper = Opencv_sklearn_wrapper__ann()
    else:
        raise RuntimeError('For make_sklearn_wrapper, both of keyword arguments "model" and "class_name" are None or \
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
