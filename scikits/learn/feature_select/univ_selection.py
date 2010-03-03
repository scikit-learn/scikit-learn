"""
Univariate features selection.

"""

# Author: B. Thirion, G. Varoquaux
# License: BSD 3 clause

import numpy as np
from scipy import stats 


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
    x[:,:k] += 3*y[:,np.newaxis]
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
    x[:,:k] += y[:, np.newaxis]
    return x, y


######################################################################
# Scoring functions
######################################################################

def f_classif(x, y):
    """
    Compute the Anova F-value for the provided sample

    Parameters
    ----------
    x : array of shape (n_samples, n_features)
        the set of regressors sthat will tested sequentially
    y : array of shape(n_samples)
        the data matrix
    
    Returns    
    -------
    F : array of shape (m),
        the set of F values
    pval : array of shape(m),
        the set of p-values
    """
    x = np.asanyarray(x)
    args = [x[y==k] for k in np.unique(y)]
    return stats.f_oneway(*args)
    

def f_regression(x, y, center=True):
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
    x : array of shape (n_samples, n_features)
        the set of regressors sthat will tested sequentially
    y : array of shape(n_samples)
        the data matrix

    center : True, bool,
        If true, x and y are centered
    
    Returns    
    -------
    F : array of shape (m),
        the set of F values
    pval : array of shape(m)
        the set of p-values
    """
    
    # orthogonalize everything wrt to confounds
    y = y.copy()
    x = x.copy()
    if center:
        y -= np.mean(y)
        x -= np.mean(x, 0)
        
    # compute the correlation
    x /= np.sqrt(np.sum(x**2,0))
    y /= np.sqrt(np.sum(y**2))
    corr = np.dot(y, x)

    # convert to p-value
    dof = y.size-1
    F = corr**2/(1-corr**2)*dof
    pv = stats.f.sf(F, 1, dof)
    return F, pv


######################################################################
# Selection function
######################################################################

def select_percentile(p_values, percentile):
    score = stats.scoreatpercentile(p_values, percentile)
    return (p_values < score)


######################################################################
# Univariate Selection
######################################################################

class UnivSelection(object):

    def __init__(self, estimator=None, 
                       score_func=f_regression,
                       select_func=None, select_args=()):
        """ An object to do univariate selection before using a
            classifier.

            Parameters
            -----------
            estimator: None or an estimator instance
                If an estimator is given, it is used to predict on the
                features selected.
            score_func: A callable
                The function used to score features. Should be::
                    
                    _, p_values = score_func(x, y)

                The first output argument is ignored.
            select_func: A callable
                The function used to select features. Should be::
                    
                    support = select_func(p_values, *select_args)
                If None is passed, the 10% lowest p_values are
                selected.
            select_args: A list or tuple
                The arguments passed to select_func
        """
        assert hasattr(select_args, '__iter__'), ValueError(
                "The select args should be a list-like."
            )
        assert callable(score_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func))
            )
        if select_func is None:
            if len(select_args) == 0:
                select_args = (10,)
            select_func = lambda p: select_percentile(p, *select_args)
        assert callable(select_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (select_func, type(select_func))
            )
        self.estimator = estimator
        self.score_func = score_func
        self.select_func = select_func


    #--------------------------------------------------------------------------
    # Estimator interface
    #--------------------------------------------------------------------------

    def fit(self, x, y):
        _, p_values_   = self.score_func(x, y)
        self.support_  = self.select_func(p_values_)
        self.p_values_ = p_values_
        if self.estimator is not None:
            self.estimator.fit(x[self.support_], y)
        return self


    def predict(self, x):
        support_ = self.support_
        if self.estimator is None:
            return support_
        else:
            return self.estimator.predict(x[support_])


if __name__ == "__main__":    
    x, y = generate_dataset_classif(n_samples=50, n_features=20, k=5, seed=2)
    F, pv = f_classif(x, y)
    univ_selection = UnivSelection(score_func=f_classif, select_args=(25,))
    univ_selection.fit(x, y)
    print univ_selection.support_.astype(int)
    assert np.all(univ_selection.p_values_ == pv)
    #x, y = generate_dataset_reg(n_samples=50, n_features=20, k=5, seed=2)
    #F, pv = f_regression(x, y)
    from scikits.learn import svm
    clf =  svm.SVM(kernel_type='linear')
    clf.fit(x, y)
    print clf.predict(x[:5])
    #print svm.predict(x,y,x)
