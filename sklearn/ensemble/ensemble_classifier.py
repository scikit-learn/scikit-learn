"""
Ensemble Classifier

This module contains a weighted majority rule ensemble classifier for scikit-learn classification estimators.

## Examples:

# Intializing different scikit-learn estimators for classification
clf1 = LogisticRegression()
clf2 = RandomForestClassifier()
clf3 = GaussianNB()

## Example 1: Using Majority Class Label
# initialization for using the majority class labels for prediction
eclf = EnsembleClassifier(clfs=[clf1, clf2, clf3])
eclf.fit(X_train, y_train)
eclf.fit(X_test)

## Example 2: Using Average Probabilties
# initialization for using average predict_proba scores to derive the class label.
eclf = EnsembleClassifier(clfs=[clf1, clf2, clf3], weights=[1,1,1])
eclf.fit(X_train, y_train)
eclf.fit(X_test)


"""

# Authors: Sebastian Raschka <se.raschka@gmail.com>
#
# Licence: BSD 3 clause

from ..base import BaseEstimator
from ..base import ClassifierMixin
from ..base import TransformerMixin
import numpy as np
import operator

class EnsembleClassifier(BaseEstimator, ClassifierMixin, TransformerMixin):
    """ 
    Ensemble classifier for scikit-learn estimators.
        
    Parameters
    ----------
    clfs : array-like, shape = [n_classifiers]
      A list of scikit-learn classifier objects.
      
    weights : array-like, shape = [n_classifiers], optional (default=`None`)
      If `None`, the majority rule voting will be applied to the predicted 
      class labels. If a list of weights (`float` or `int`) is provided, 
      the averaged raw probabilities (via `predict_proba`) will be used to 
      determine the most confident class label.
        
    Attributes
    ----------
    classes_ : array-like, shape = [n_class_labels, n_classifiers]
        Class labels predicted by each classifier if `weights=None`.
       
    probas_ : array, shape = [n_probabilities, n_classifiers]
        Predicted probabilities by each classifier if `weights=array-like`.
      
    """
    def __init__(self, clfs, weights=None):
        self.clfs = clfs
        self.weights = weights

    def fit(self, X, y):
        """ 
        Fits the scikit-learn estimators.
        
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.
            
        y : array-like, shape = [n_samples]
            Target values.
      
        Returns
        -------
        self : object
        
        """
        for clf in self.clfs:
            clf.fit(X, y)
            
        return self
            
    def predict(self, X):
        """
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.
        
        Returns
        ----------
        maj : array-like, shape = [n_class_labels]
            Predicted class labels by majority rule.
        
        """
        if self.weights:
            avg = self.predict_proba(X)
            maj = np.apply_along_axis(lambda x: max(enumerate(x), key=operator.itemgetter(1))[0], axis=1, arr=avg)
        
        else:
            self.classes_ = self._get_classes(X)
            maj = np.asarray([np.argmax(np.bincount(self.classes_[:,c])) for c in range(self.classes_.shape[1])]) 
        
        return maj
            
    
    def predict_proba(self, X):
        
        """
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.
        
        Returns
        ----------
        avg : array-like, shape = [n_samples, n_probabilities]
            Weighted average probability for each class per sample.
        
        """
        self.probas_ = self._get_probas(X)
        avg = np.average(self.probas_, axis=0, weights=self.weights)
        
        return avg

  
    def transform(self, X):
        """         
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples and
            n_features is the number of features.
      
        Returns
        -------
        If not `weights=None`:
          array-like = [n_classifier_results, n_class_proba, n_class]
            Class probabilties calculated by each classifier.
        
        Else:
          array-like = [n_classifier_results, n_class_label]
            Class labels predicted by each classifier.
        
        """
        if self.weights:
            return self._get_probas(X)
            
        else:
            return self._get_classes(X)  
    
    def _get_classes(self, X):
        """ Collects results from clf.predict calls. """
        return np.asarray([clf.predict(X) for clf in self.clfs])
        
    def _get_probas(self, X):
        """ Collects results from clf.predict calls. """
        return np.asarray([clf.predict_proba(X) for clf in self.clfs])
