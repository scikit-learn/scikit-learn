# encoding: utf-8
# cython: cdivision=True
# cython: profile=True
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

#@cython.profile(False)
cpdef double eval_gini(np.ndarray labels, 
                       np.ndarray[np.float_t, ndim=1] pm) except * :
    """
        
        Gini index = \sum_{k=0}^{K-1} pmk (1 - pmk)
                   = 1 - \sum_{k=0}^{K-1} pmk ** 2
            
    """       
    cdef int *_labels = <int *>labels.data
    cdef int n_labels = labels.shape[0]
    
    cdef double *_pm = <double *>pm.data
    cdef int K = pm.shape[0]    
    
    cdef int i = 0
    for i in range(K):
        _pm[i] = 0.
    
    cdef int value
    for i in range(n_labels):
        value = int(_labels[i])
        _pm[value] += 1. / n_labels
    
    cdef double H = 1.
    cdef Py_ssize_t k
    for k in range(K):    
        H -=  pow(_pm[k],2) 
         
    return H   

cpdef double eval_entropy(np.ndarray labels,
                          np.ndarray[np.float_t, ndim=1] pm) except * :
    """
        
        Cross Entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
            
    """
    cdef int *_labels = <int *>labels.data
    cdef int n_labels = labels.shape[0]
    
    cdef double *_pm = <double *>pm.data
    cdef int K = pm.shape[0]    
    
    cdef int i = 0
    for i in range(K):
        _pm[i] = 0.
    
    cdef int value
    for i in range(n_labels):
        value = int(_labels[i])
        _pm[value] += 1. / n_labels
        
    cdef double H = 0.
    cdef Py_ssize_t k
    for k in range(K):
        if pm[k] > 0 :    
            H +=  -pm[k] * log(pm[k])
         
    return H   

cpdef double eval_miss(np.ndarray labels,
                       np.ndarray[np.float_t, ndim=1] pm) except * :
    """
        
        Misclassification error = (1 - pmk)
            
    """
    cdef int *_labels = <int *>labels.data
    cdef int n_labels = labels.shape[0]
    
    cdef double *_pm = <double *>pm.data
    cdef int K = pm.shape[0]    
    
    cdef int i = 0
    for i in range(K):
        _pm[i] = 0.
    
    cdef int value
    for i in range(n_labels):
        value = int(_labels[i])
        _pm[value] += 1. / n_labels
        
    cdef double H = 1. - pm.max()
    
"""
 Regression entropy measures
 
"""      
cpdef double eval_mse(np.ndarray labels,
                      np.ndarray[np.float_t, ndim=1] pm) except * :
    """             
        MSE =  \sum_i (y_i - c0)^2  / N
        
        pm is a redundant argument.    
    """   
    
    cdef float c0 = np.mean(labels)
    cdef int N = int(labels.shape[0])        
                       
    cdef double H = 0.
    cdef Py_ssize_t i
    for i in range(N):
        H += pow(labels[i] - c0, 2) 
    H /= N
        
    return H
    
    