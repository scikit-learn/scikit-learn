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
import sys
from time import time

cimport numpy as np
cimport cython

cdef extern from "math.h":
    cdef extern double log(double x)
    cdef extern double pow(double base, double exponent)

"""
 Helper function
""" 
cpdef bincount_k(np.ndarray[np.int_t, ndim=1] X, int K): 
    '''
    a = bincount_k(X, K)

    Returns numpy array consisting of a bincount of X, padded to length K
    
    np.bincount is efficient for computing the occurances of labels, but
    there isn't an easy way to resize the array to size K (for K categories)

    Parameters
    ----------
    X  : array-like, shape = (n,)
    K  : integer, the number of categories

    Returns
    -------
    Padded bincount
    '''
    
    cdef np.ndarray[np.int_t, ndim=1] a = \
        np.zeros((K,), dtype=np.int)
    cdef np.ndarray[np.int_t, ndim=1] b = np.bincount(X)
    
    cdef int bins = len(b)
    assert K >= bins
    
    cdef Py_ssize_t j
    for j in range(0,bins):
        a[j] = b[j]       
    
    return a
 
"""
 Classification entropy measures
 
    From Hastie et al. Elements of Statistical Learning, 2009.
         
    If a target is a classification outcome taking on values 0,1,...,K-1
    In node m, representing a region Rm with Nm observations, let
            
       pmk = 1/ Nm \sum_{x_i in Rm} I(yi = k)
          
    be the proportion of class k observations in node m   
"""  
cpdef double eval_gini(np.ndarray[np.int_t, ndim=1] labels)  except * :
    """
        
        Gini index = \sum_{k=0}^{K-1} pmk (1 - pmk)
                   = 1 - \sum_{k=0}^{K-1} pmk ** 2
            
    """
    cdef int K = int(labels.max()) + 1
    N = float(labels.shape[0])
    
    cdef np.ndarray[np.float64_t, ndim=1] pm = \
        bincount_k(labels, K) / N
        
    cdef double H = 1.
    cdef Py_ssize_t k
    for k in range(K):    
        H -=  pow(pm[k],2) 
         
    return H   

cpdef double eval_entropy(np.ndarray[np.int_t, ndim=1] labels)  except * :
    """
        
        Cross Entropy = - \sum_{k=0}^{K-1} pmk log(pmk)
            
    """
    cdef int K = int(labels.max()) + 1
    N = float(labels.shape[0])
    
    cdef np.ndarray[np.float64_t, ndim=1] pm = \
        bincount_k(labels, K) / N
        
    cdef double H = 0.
    cdef Py_ssize_t k
    for k in range(K):
        if pm[k] > 0 :    
            H +=  -pm[k] * log(pm[k])
         
    return H   

cpdef double eval_miss(np.ndarray[np.int_t, ndim=1] labels)  except * :
    """
        
        Misclassification error = (1 - pmk)
            
    """
    cdef int K = int(labels.max()) + 1
    N = float(labels.shape[0])
    
    cdef np.ndarray[np.float64_t, ndim=1] pm = \
        bincount_k(labels, K) / N
        
    cdef double H = 1. - pm.max()
    
"""
 Regression entropy measures
 
"""      
cpdef double eval_mse(np.ndarray[np.float64_t, ndim=1] labels)  except * :
    """             
        MSE =  \sum_i (y_i - c0)^2  / N
            
    """   
    
    cdef float c0 = np.mean(labels)
    cdef int N = int(labels.shape[0])        
                       
    cdef double H = 0.
    cdef Py_ssize_t i
    for i in range(N):
        H += pow(labels[i] - c0, 2) 
    H /= N
        
    return H
    
    