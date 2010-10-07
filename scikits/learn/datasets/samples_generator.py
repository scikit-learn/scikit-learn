"""
Generate samples of synthetic data sets.
"""

# Author: B. Thirion, G. Varoquaux, A. Gramfort, V. Michel
# License: BSD 3 clause

import numpy as np
import numpy.random as nr



def samples_classif():
    pass

######################################################################
# Generate Dataset for test
######################################################################

def test_dataset_classif(n_samples=100, n_features=100, param=[1,1],
                             k=0, seed=None):
    """
    Generate an snp matrix

    Parameters
    ----------
    n_samples : 100, int,
        the number of subjects
    n_features : 100, int,
        the number of features
    param : [1,1], list,
        parameter of a dirichlet density
        that is used to generate multinomial densities
        from which the n_featuress will be samples
    k : 0, int,
        number of informative features
    seed : None, int or np.random.RandomState
        if seed is an instance of np.random.RandomState,
        it is used to initialize the random generator

    Returns
    -------
    x : array of shape(n_samples, n_features),
        the design matrix
    y : array of shape (n_samples),
        the subject labels

    """
    assert k<=n_features, ValueError('cannot have %d informative features and'
                                   ' %d features' % (k, n_features))
    if isinstance(seed, np.random.RandomState):
        random = seed
    elif seed is not None:
        random = np.random.RandomState(seed)
    else:
        random = np.random

    x = random.randn(n_samples, n_features)
    y = np.zeros(n_samples)
    param = np.ravel(np.array(param)).astype(np.float)
    for n in range(n_samples):
        y[n] = np.nonzero(random.multinomial(1, param/param.sum()))[0]
    x[:,:k] += 3*y[:,np.newaxis]
    return x, y

def test_dataset_reg(n_samples=100, n_features=100, k=0, seed=None):
    """
    Generate an snp matrix

    Parameters
    ----------
    n_samples : 100, int,
        the number of subjects
    n_features : 100, int,
        the number of features
    k : 0, int,
        number of informative features
    seed : None, int or np.random.RandomState
        if seed is an instance of np.random.RandomState,
        it is used to initialize the random generator

    Returns
    -------
    x : array of shape(n_samples, n_features),
        the design matrix
    y : array of shape (n_samples),
        the subject data

    """
    assert k<n_features, ValueError('cannot have %d informative fetaures and'
                                   ' %d features' % (k, n_features))
    if isinstance(seed, np.random.RandomState):
        random = seed
    elif seed is not None:
        random = np.random.RandomState(seed)
    else:
        random = np.random

    x = random.randn(n_samples, n_features)
    y = random.randn(n_samples)
    x[:,:k] += y[:, np.newaxis]
    return x, y





######################################################################
# Generate Dataset for regression
######################################################################


def sparse_uncorrelated(nb_samples=100, nb_features=10):
    """
    Function creating simulated data with sparse uncorrelated design.
    (cf.Celeux et al. 2009,  Bayesian regularization in regression)
    X = NR.normal(0,1)
    Y = NR.normal(X[:,0]+2*X[:,1]-2*X[:,2]-1.5*X[:,3])
    The number of features is at least 10.

    Parameters
    ----------
    nb_samples : int
                 number of samples (default is 100).
    nb_features : int
                  number of features (default is 5).

    Returns
    -------
    X : numpy array of shape (nb_samples, nb_features) for input samples
    Y : numpy array of shape (nb_samples) for labels
    """
    X = nr.normal(loc=0, scale=1, size=(nb_samples, nb_features))
    Y = nr.normal(loc=X[:, 0] + 2 * X[:, 1] - 2 * X[:,2] - 1.5 * X[:, 3],
                  scale = np.ones(nb_samples))
    return X, Y


def friedman(nb_samples=100, nb_features=10,noise_std=1):
    """
    Function creating simulated data with non linearities 
    (cf.Friedman 1993)
    X = NR.normal(0,1)
    Y = 10*sin(X[:,0]*X[:,1]) + 20*(X[:,2]-0.5)**2 + 10*X[:,3] + 5*X[:,4]
    The number of features is at least 5.

    Parameters
    ----------
    nb_samples : int
                 number of samples (default is 100).
    nb_features : int
                  number of features (default is 10).
    noise_std : float
		std of the noise, which is added as noise_std*NR.normal(0,1)
    Returns
    -------
    X : numpy array of shape (nb_samples, nb_features) for input samples
    Y : numpy array of shape (nb_samples) for labels

    """
    X = nr.normal(loc=0, scale=1, size=(nb_samples, nb_features))
    Y = 10*np.sin(X[:,0]*X[:,1]) + 20*(X[:,2]-0.5)**2 + 10*X[:,3] + 5*X[:,4]
    Y += noise_std*nr.normal(loc=0,scale=1,size=(nb_samples))
    return X,Y
