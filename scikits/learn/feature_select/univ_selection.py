import numpy as np
from scipy.stats import f_oneway




######################################################################
# Generate Dataset
######################################################################

def generate_dataset_classif(n_samples=100, n_features=100, param=[1,1],
                             k=0, seed=None):
    """
    Generate an snp matrix

    Parameters
    ----------
    n_samples : 100, int,
        the number of subjects
    n_features : 100, int,
        the number of featyres
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
    assert k<n_features, ValueError('cannot have %d informative fetaures and'
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
    x[:k] += 3*y
    return x, y.astype(np.int)

def generate_dataset_reg(n_samples=100, n_features=100, k=0, seed=None):
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
    x[:k] += 3*y
    return x, y


######################################################################
# Scoring functions
######################################################################


def quick_lm_for1Dcontrast(y, x, center=True):
    """
    Quick linear model for testing the effect of a single regressor,
    sequentially for many regressors
    This is done in 3 steps:
    1. the regressor of interest and the data are orthogonalized
    wrt constant regressors
    2. the cross correlation between data and regressors is computed
    3. it is converted to an F score then to a p-value

    Parameters
    ----------
    y : array of shape(n_samples)
        the data matrix
    x : array of shape (n_samples, n_features)
        the set of regressors sthat will tested sequentially
    center : True, bool,
        If true, x and y are centered
    
    Returns    
    -------
    pval, array of shape(m) the set of p-values
    """
    
    # orthogonalize everything wrt to confounds
    if center:
        y = y.copy() - np.mean(y)
        x = x.copy() - np.mean(x,1)
        
    # compute the correlation
    x = (x.T/np.sqrt(np.sum(x**2,1))).T
    y = (y.T/np.sqrt(np.sum(y**2,1))).T
    corr = np.dot(y,x)

    # convert to p-value
    dof = y.size-1
    F = corr**2/(1-corr**2)*dof
    pv = st.f.sf(F, 1, dof)
    return F, pv






######################################################################
# Univariate Selection
######################################################################





if __name__ == "__main__":    
    x, y = generate_dataset_classif(seed=2)
    x, y = generate_dataset_reg(seed=2)

