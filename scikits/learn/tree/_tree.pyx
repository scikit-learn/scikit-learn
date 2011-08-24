# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
#
# Adapted for CART by Brian Holt <bdholt1@gmail.com>
#

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h":
    cdef extern double log(double x)
    cdef extern double pow(double base, double exponent)

 
"""
 Classification entropy measures
 
    From Hastie et al. Elements of Statistical Learning, 2009.
         
    If a target is a classification outcome taking on values 0,1,...,K-1
    In node m, representing a region Rm with Nm observations, let
            
       pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)
          
    be the proportion of class k observations in node m   
"""  

cdef double eval_gini(np.ndarray[np.float_t, ndim=1] pm):
    """
        
        Gini index = \sum_{k=0}^{K-1} pmk (1 - pmk)
                   = 1 - \sum_{k=0}^{K-1} pmk ** 2
            
    """
    cdef int K = pm.shape[0]
    
    cdef double sum = 0.
    cdef Py_ssize_t k    
    for k in range(K): 
        sum += pm[k]
         
    cdef double H = 1.
    for k in range(K): 
        H -=  (pm[k] / sum) * (pm[k] / sum) 
         
    return H   

cdef double eval_entropy(np.ndarray[np.float_t, ndim=1] pm):
    """
        
        Cross Entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
            
    """
    cdef int K = pm.shape[0]  
    
    cdef double sum = 0.
    cdef Py_ssize_t k    
    for k in range(K): 
        sum += pm[k]       
      
    cdef double H = 0.
    for k in range(K):
        if pm[k] > 0 :    
            H +=  (-pm[k] / sum) * log(pm[k] / sum)
         
    return H   

cdef double eval_miss(np.ndarray[np.float_t, ndim=1] pm):
    """
        
        Misclassification error = (1 - pmk)
            
    """     
    cdef int K = pm.shape[0]  
    
    cdef double sum = 0.
    cdef Py_ssize_t k    
    for k in range(K): 
        sum += pm[k]    
    
    cdef double H = 1. - (pm.max() / sum)
    
    return H
    
"""
 Regression entropy measures
 
"""      
cdef double eval_mse(np.ndarray[np.float_t, ndim=1] labels):
    """             
        MSE =  \sum_i (y_i - c0)^2  / N
    """     
    cdef int n_labels = labels.shape[0]

    cdef float c0
    cdef Py_ssize_t i = 0
    for i in range(n_labels):
        c0 += labels[i]
    c0 /= n_labels

    cdef double H = 0.
    for i in range(n_labels):
        H += pow(labels[i] - c0, 2) 
    H /= n_labels
        
    return H

cdef double call_classification_crit(int criterion, 
                                     np.ndarray[np.float_t, ndim=1] pm):
    """
    Map the criterion to the corresponding function call.
    """
    cdef double retval = 0.                        
    if criterion == 0:
        retval = eval_gini(pm)
    elif criterion == 1:
        retval = eval_entropy(pm)
    elif criterion == 2:
        retval = eval_miss(pm)        
    else:
        raise ValueError("Classification criterion %s invalid!"
                         "Must be [0, 1, 2]" % (criterion))                    
                     
    return retval
    
def _find_best_split_classification(np.ndarray[np.float_t, ndim=2] features, 
                                    np.ndarray[np.int_t, ndim=1] labels, 
                                    int crit):                                 
    """
    Find the best split for classification.  This is achieved efficiently by
    sorting the data at each dimension and then progressively adjusting
    the proportion of classes left and right of the split followed by 
    evaluation.
    """        
                                    
    cdef int n_samples = features.shape[0]
    cdef int n_features = features.shape[1]
    cdef int K = int(np.abs(labels.max())) + 1
    cdef np.ndarray[np.float_t, ndim=1] pm_left = \
        np.zeros((K,), dtype=np.float)
    cdef np.ndarray[np.float_t, ndim=1] pm_right = \
        np.zeros((K,), dtype=np.float)

    best = None
    cdef double split_error = \
        call_classification_crit(crit, np.bincount(labels).astype(np.float))
    if split_error == 0:
        return best # short circuit
    
    cdef int i = 0, j = 0
    cdef np.ndarray[np.float_t, ndim=1] features_at_i    
    cdef np.ndarray[np.int_t, ndim=1] sorted_idx
    cdef np.ndarray[np.int_t, ndim=1] sorted_labels
    cdef np.ndarray[np.float_t, ndim=1] sorted_features    
    cdef double t, e1, e2, error
    
    for i in range(n_features):
        features_at_i = features[:, i]
        sorted_idx = np.argsort(features_at_i)
        sorted_labels = labels[sorted_idx]
        sorted_features = features_at_i[sorted_idx]
        pm_left[:] = 0.
        pm_right[:] = np.bincount(labels)
        
        for j in xrange(n_samples - 1):
            pm_left[sorted_labels[j]] += 1
            pm_right[sorted_labels[j]] -= 1
            t = (sorted_features[j] + sorted_features[j + 1]) / 2. 
            e1 = (j + 1) / <double>n_samples * \
                call_classification_crit(crit, pm_left)
            e2 = (n_samples - (j + 1)) / <double>n_samples * \
                call_classification_crit(crit, pm_right)              
            error = e1 + e2
            
            if error < split_error:
                split_error = error
                best = i, t, error
                if split_error == 0:
                    return best # short circuit
    return best


cdef double call_regression_crit(int criterion, 
                                 np.ndarray[np.float_t, ndim=1] targets):
    """
    Map the criterion to the corresponding function call.
    """    
    cdef double retval = 0.                        
    if criterion == 0:
        retval = eval_mse(targets)      
    else:
        raise ValueError("Regression criterion %s invalid!"
                         "Must be [0]" % (criterion))                    
                     
    return retval


def _find_best_split_regression(np.ndarray[np.float_t, ndim=2] features, 
                                np.ndarray[np.float_t, ndim=1] targets, 
                                int crit):
    """
    Find the best split for regression.  Mirroring classification, this
    sorts the data at each dimension and then progressively evaluates
    the entropy scores at each split.  This function could probably do with
    a bit more optimising since an array slice over the labels is taken in
    the inner loop. 
    """  
    
    cdef int n_samples = features.shape[0]
    cdef int n_features = features.shape[1]

    best = None
    cdef double split_error = \
        call_regression_crit(crit, targets)
    if split_error == 0:
        return best # short circuit

    cdef int i = 0, j = 0
    cdef np.ndarray[np.float_t, ndim=1] features_at_i
    cdef np.ndarray[np.int_t, ndim=1] sorted_idx
    cdef np.ndarray[np.float_t, ndim=1] sorted_targets
    cdef np.ndarray[np.float_t, ndim=1] sorted_features
    cdef double t, e1, e2, error
          
    for i in xrange(n_features):
        features_at_i = features[:, i]
        sorted_idx = np.argsort(features_at_i)
        sorted_targets = targets[sorted_idx]
        sorted_features = features_at_i[sorted_idx]
        
        for j in xrange(n_samples - 1):
            t = (sorted_features[j] + sorted_features[j + 1]) / 2.
            e1 = (j + 1) / <double>n_samples * \
                call_regression_crit(crit, sorted_targets[0:j+1])
            e2 = (n_samples - (j + 1)) / <double>n_samples * \
                call_regression_crit(crit, sorted_targets[j + 1:n_samples])
            error = e1 + e2
            
            if error < split_error:
                split_error = error
                best = i, t, error
                if split_error == 0:
                    return best # short circuit
    return best
