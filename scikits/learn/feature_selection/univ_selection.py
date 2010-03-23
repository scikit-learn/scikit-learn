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
    return x, y

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
    dof = y.size-2
    F = corr**2/(1-corr**2)*dof
    pv = stats.f.sf(F, 1, dof)
    return F, pv


######################################################################
# Selection function
######################################################################

def select_percentile(p_values, percentile):
    """ Select the best percentile of the p_values
    """
    alpha = stats.scoreatpercentile(p_values, percentile)
    return (p_values <= alpha)

def select_k_best(p_values, k):
    """Select the k lowest p-values
    """
    assert k<len(p_values), ValueError('cannot select %d features'
                                       ' among %d ' % (k, len(p_values)))
    #alpha = stats.scoreatpercentile(p_values, 100.*k/len(p_values))
    alpha = np.sort(p_values)[k]
    return (p_values < alpha)


def select_fpr(p_values, alpha):
    """Select the pvalues below alpha
    """
    return (p_values < alpha)

def select_fdr(p_values, alpha):
    """
    Select the p-values corresponding to an estimated false discovery rate
    of alpha
    This uses the Benjamini-Hochberg procedure
    """
    sv = np.sort(p_values)
    threshold = sv[sv < alpha*np.arange(len(p_values))].max()
    return (p_values < threshold)

def select_fwe(p_values, alpha):
    """
    Select the p-values corresponding to a corrected p-value of alpha
    """
    return (p_values<alpha/len(p_values))



######################################################################
# Univariate Selection
######################################################################

class UnivSelection(object):

    def __init__(self, estimator=None,
                       score_func=f_regression,
                       select_func=None, select_args=(10,)):
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
        if not hasattr(select_args, '__iter__'):
            select_args = list(select_args)
        assert callable(score_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (score_func, type(score_func))
            )
        if select_func is None:
            select_func = select_percentile
        assert callable(select_func), ValueError(
                "The score function should be a callable, '%s' (type %s) "
                "was passed." % (select_func, type(select_func))
            )
        self.estimator = estimator
        self.score_func = score_func
        self.select_func = select_func
        self.select_args = select_args


    #--------------------------------------------------------------------------
    # Estimator interface
    #--------------------------------------------------------------------------

    def fit(self, x, y):
        _, p_values_   = self.score_func(x, y)
        self.support_  = self.select_func(p_values_,*self.select_args)
        self.p_values_ = p_values_
        if self.estimator is not None:
            self.estimator.fit(x[:,self.support_], y)
        return self


    def predict(self, x=None):
        support_ = self.support_
        if x is None or self.estimator is None:
            return support_
        else:
            return self.estimator.predict(x[:,support_])


if __name__ == "__main__":
    x, y = generate_dataset_classif(n_samples=50, n_features=20, k=5, seed=2)
    F, pv = f_classif(x, y)
    univ_selection = UnivSelection(score_func=f_classif, select_args=(25,))
    univ_selection.fit(x, y)
    print univ_selection.support_.astype(int)

    univ_selection = UnivSelection(score_func=f_classif,
                                   select_func=select_k_best,
                                   select_args=(5,))
    univ_selection.fit(x, y)
    print univ_selection.support_.astype(int)

    univ_selection = UnivSelection(score_func=f_classif,
                                   select_func=select_fpr,
                                   select_args=(0.001,))
    univ_selection.fit(x, y)

    print univ_selection.support_.astype(int)
    univ_selection = UnivSelection(score_func=f_classif,
                                   select_func=select_fwe,
                                   select_args=(0.05,))
    univ_selection.fit(x, y)
    print univ_selection.support_.astype(int)

    univ_selection = UnivSelection(score_func=f_classif,
                                   select_func=select_fdr,
                                   select_args=(0.05,))
    univ_selection.fit(x, y)
    print univ_selection.support_.astype(int)




    assert np.all(univ_selection.p_values_ == pv)
    #x, y = generate_dataset_reg(n_samples=50, n_features=20, k=5, seed=2)
    #F, pv = f_regression(x, y)
    from scikits.learn import svm
    clf =  svm.SVM(kernel_type='linear')
    y = np.asfarray(y)
    clf.fit(x, y)
    print clf.predict(x)
    #print svm.predict(x,y,x)
