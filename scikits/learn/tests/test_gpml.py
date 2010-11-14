"""
Testing for Gaussian Process for Machine Learning module (scikits.learn.gpml)
"""

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
                          assert_almost_equal, assert_raises, assert_

from .. import gpml, datasets, cross_val, metrics

diabetes = datasets.load_diabetes()

def test_regression_1d_x_sinx():
    """
    MLE estimation of a Gaussian Process model with a squared exponential
    correlation model (correxp2). Check random start optimization.
    
    Test the interpolating property.
    """
    
    f = lambda x: x * np.sin(x)
    X = np.array([1., 3., 5., 6., 7., 8.])
    y = f(X)
    gpm = gpml.GaussianProcessModel(theta0=1e-2, thetaL=1e-4, thetaU=1e-1, random_start=10, verbose=False).fit(X, y)
    y_pred, MSE = gpm.predict(X, eval_MSE=True)
    
    assert (np.all(np.abs((y_pred - y) / y)  < 1e-6) and np.all(MSE  < 1e-6))

def test_regression_diabetes(n_jobs = 1, verbose = 0):
    """
    MLE estimation of a Gaussian Process model with an anisotropic squared exponential
    correlation model.
    
    Test the model using cross-validation module (quite time-consuming).
    
    Poor performance: Leave-one-out estimate of explained variance is about 0.5
    at best... To be investigated!
    TODO: find a dataset that would prove GPML performance!
    """
    
    X, y = diabetes['data'], diabetes['target']
    
    gpm = gpml.GaussianProcessModel(corr=gpml.correxp2, theta0=1e-4, nugget=1e-2, verbose=False).fit(X, y)

    gpm.thetaL, gpm.thetaU = None, None
    y_pred = cross_val.cross_val_score(gpm, X, y=y, cv=cross_val.LeaveOneOut(y.size), n_jobs=n_jobs, verbose=verbose) + y
    Q2 = metrics.explained_variance(y_pred, y)
    
    assert Q2 > 0.45
