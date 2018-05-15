# -*- coding: utf-8 -*-

"""
The :mod:`sklearn.aode` module implements Averaged one-dependance algorithm.
AODE averages over all of a small space of alternative naive-Bayes-like models 
that have weaker (and hence less detrimental) independence assumptions than naive Bayes.
"""

import numpy as np
from .base import BaseEstimator, ClassifierMixin, TransformerMixin
from .utils.validation import check_X_y, check_array, check_is_fitted
from .utils.multiclass import unique_labels
from scipy.sparse import issparse


from ._aode import AODE_helper, get_tables

class AODE(BaseEstimator, ClassifierMixin):
    """
    Averaged one-dependance estimator (AODE)

    AODE averages over all of a small space of alternative naive-Bayes-like models 
    that have weaker (and hence less detrimental) independence assumptions than naive Bayes.

    Read more in the :ref:`User Guide <>`.

    Parameters
    ----------
    m_val : float, optional (default=1.0)
        Additive (Laplace/Lidstone) smoothing parameter
        (0 for no smoothing).

    m_limit : int, optional (default=1)
        Impose a frequency limit for superParents
    """

    def __init__(self, m_val=1, m_limit=1):
        self.m_val = m_val
        self.m_limit = m_limit
    
    def fit(self, X, y):
        """Fit Averaged one-dependance estimator according to X, y
        
        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.
        
        y : array-like, shape (n_samples,)
            Target values.
        
        Returns
        -------
        self : object
        """
        X, y = check_X_y(X, y, accept_sparse='csr')
        self.classes_ = unique_labels(y)
        # Return the estimator
        
        self.X_ = X
        self.y_ = y

        if issparse(self.X_):
            self.X_ = self.X_.toarray()
        if issparse(self.y_):
            self.y_ = self.y_.toarray()
        

        self.num_attr_ = self.X_.shape[1] # Number of attributes
        self.num_samples_ = self.X_.shape[0]

                
        # Get number of unique values of x, each col
        self.num_uniq_x_vals_ = np.zeros(self.num_attr_, dtype=np.int)
        self.x_indexes_ = [{} for i in range(self.num_attr_)]
        idx = 0
        for x in range(self.num_attr_):
            idx = 0
            for i in range(self.num_samples_):
                if self.X_[i,x] not in self.x_indexes_[x]:
                    self.num_uniq_x_vals_[x] += 1
                    self.x_indexes_[x][self.X_[i,x]] = idx
                    idx += 1
        
    
        # Count Y
        self.y_vals_, y_indexes, y_counts_ = np.unique(self.y_, return_inverse=True,return_counts=True)
        self.num_classes_ = len(self.y_vals_)
        
        
        # Create empty frequency tables
        self.freq_table_ = [[[0 for k in range(self.num_uniq_x_vals_[j])] for j in range(self.num_attr_)] for i in range(self.num_classes_)]
        self.joint_freq_table_ = [[[[[0 for m in range(self.num_uniq_x_vals_[k])] for l in range(self.num_uniq_x_vals_[j])] for k in range(0,j)] for j in range(self.num_attr_)] for i in range(self.num_classes_)]

        # Fill frequency tables
        self.freq_table_, self.joint_freq_table_ = get_tables(self.X_, self.num_samples_, y_indexes, self.num_attr_, self.x_indexes_, self.num_uniq_x_vals_, self.num_classes_, self.freq_table_, self.joint_freq_table_)

        return self
    
    
    def predict(self, X):
        """
        Perform classification on an array of test vectors X.
        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
        Returns
        -------
        C : array, shape = [n_samples]
            Predicted target values for X
        """
        # Check is fit had been called
        check_is_fitted(self, ['X_', 'y_'])
        
        X = check_array(X, accept_sparse='csr')
        
        if ((self.num_attr_) != X.shape[1]):
            raise ValueError()

        if issparse(X):
            X = X.toarray()

        num_predictions = X.shape[0]
        ret = np.empty(num_predictions, dtype=self.y_.dtype)

        x_val_indexes = [x[:] for x in [[-1] * self.num_attr_] * num_predictions]
        for pred in range(num_predictions):
            for i in range(self.num_attr_):
                if (X[pred,i]) in self.x_indexes_[i]:
                    x_val_indexes[pred][i] = self.x_indexes_[i][X[pred,i]]

        aode_helper = AODE_helper(self.freq_table_, self.joint_freq_table_, self.num_uniq_x_vals_, self.num_attr_, self.num_classes_, self.num_samples_, self.m_val, self.m_limit)

        sol = aode_helper.calculate_predictions(x_val_indexes, num_predictions)

        for i in range(num_predictions):
            ret[i] = self.y_vals_[sol[i]]
        
        return ret

    def predict_proba(self, X):
        """
        Return probability estimates for the test vector X.
        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
        Returns
        -------
        C : array-like, shape = [n_samples, n_classes]
            Returns the probability of the samples for each class in
            the model. The columns correspond to the classes in sorted
            order, as they appear in the attribute `classes_`.
        """
        check_is_fitted(self, ['X_', 'y_'])
        
        X = check_array(X, accept_sparse='csr')
        
        if ((self.num_attr_) != X.shape[1]):
            raise ValueError()

        if issparse(X):
            X = X.toarray()

        num_predictions = X.shape[0]

        x_val_indexes = [x[:] for x in [[-1] * self.num_attr_] * num_predictions]
        for pred in range(num_predictions):
            for i in range(self.num_attr_):
                if (X[pred,i]) in self.x_indexes_[i]:
                    x_val_indexes[pred][i] = self.x_indexes_[i][X[pred,i]]

        aode_helper = AODE_helper(self.freq_table_, self.joint_freq_table_, self.num_uniq_x_vals_, self.num_attr_, self.num_classes_, self.num_samples_, self.m_val, self.m_limit)

        probs = aode_helper.calculate_probabilities(x_val_indexes, num_predictions)
        
        ret = np.asarray(probs, dtype=np.float64)

        return ret