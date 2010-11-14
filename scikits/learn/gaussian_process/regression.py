#!/usr/bin/python
# -*- coding: utf-8 -*-

################
# Dependencies #
################

import numpy as np

############################
# Defaut regression models #
############################

def regpoly0(x):
    """
    Zero order polynomial (constant, p = 1) regression model.
    
    regpoly0 : x --> f(x) = 1
    
    Parameters
    ----------
    x : array_like
        An array with shape (n_eval, n_features) giving the locations x at which the regression model should be evaluated.
    
    Returns
    -------
    f : array_like
        An array with shape (n_eval, p) with the values of the regression model.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    n_eval = x.shape[0]
    f = np.ones([n_eval,1])
    
    return f

def regpoly1(x):
    """
    First order polynomial (hyperplane, p = n) regression model.
    
    regpoly1 : x --> f(x) = [ x_1, ..., x_n ].T
    
    Parameters
    ----------
    x : array_like
        An array with shape (n_eval, n_features) giving the locations x at which the regression model should be evaluated.
    
    Returns
    -------
    f : array_like
        An array with shape (n_eval, p) with the values of the regression model.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    n_eval = x.shape[0]
    f = np.hstack([np.ones([n_eval,1]), x])
    
    return f

def regpoly2(x):
    """
    Second order polynomial (hyperparaboloid, p = n*(n-1)/2) regression model.
    
    regpoly2 : x --> f(x) = [ x_i*x_j,  (i,j) = 1,...,n ].T
                                            i > j
    
    Parameters
    ----------
    x : array_like
        An array with shape (n_eval, n_features) giving the locations x at which the regression model should be evaluated.
    
    Returns
    -------
    f : array_like
        An array with shape (n_eval, p) with the values of the regression model.
    """
    
    x = np.asanyarray(x, dtype=np.float)
    n_eval, n_features = x.shape
    f = np.hstack([np.ones([n_eval,1]), x])
    for  k in range(n_features):
        f = np.hstack([f, x[:,k,np.newaxis] * x[:,k:]])
    
    return f