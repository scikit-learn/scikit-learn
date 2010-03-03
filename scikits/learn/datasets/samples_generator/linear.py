import numpy as np


def sparse_uncorrelated(nb_samples=100,nb_features=10):
    """
    Function creating simulated data with sparse uncorrelated design.
    (cf.Celeux et al. 2009,  Bayesian regularization in regression)
    X = NR.normal(0,1)
    Y = NR.normal(2+X[:,2]+2*X[:,3]-2*X[:,6]-1.5*X[:,7])
    The number of features is at least 10.

    Parameters
    ----------
    nb_samples : int
                 number of samples (defaut is 100).
    nb_features : int
                  number of features (defaut is 10).
    
    Returns
    -------
    X : numpy array of shape (nb_samples,nb_features)
	simulated samples.
    Y : numpy array of shape (nb_samples)
    """
    X = nr.normal(loc=0,scale=1,size=(nb_samples,nb_features))
    Y = nr.normal(loc=2+x[:,2]+2*x[:,3]-2*x[:,6]-1.5*x[:,7],
    scale=np.ones(nb_samples))
    return X,Y