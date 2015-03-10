# Authors: Cecilia Damon <cecilia.damon@institut-hypercube.org>
"""
Mutual Information estimator based on the Renyi quadratic entropy and the Cauchy Schwartz divergence
combined with the Parzen window density estimator for continuous variable

References :
-------
    ..[1] J.C. Principe, Information Theoretic Learning: Renyi’s Entropy and Kernel 47 Perspectives, Information Science and Statistics, DOI 10.1007/978-1-4419-1570-2 2, Springer Science+Business Media, LLC 2010)
    ..[2] D. XU and D. Erdogmuns, Renyi's entropy, divergence and their nonparametric estimators
    ..[3] Peng et al: feature selection based on mutual information : criteria of max-dependency, max-relevance and min-redundancy (IEEE transactions on pattern analysis and machine intelligence, 2005)
    ..[4] Kari Torkkola, Feature Extraction by Non-Parametric Mutual Information Maximization (Journal of Machine Learning Research 3 (2003) 1415-1438, http://www.jmlr.org/papers/volume3/torkkola03a/torkkola03a.pdf)
    ..[5] Leonardo Macrini, Leonardo Gonçalves, Application of Rényi Entropy and Mutual Information of Cauchy-Schwartz in Selecting Variables

ICS : Mutual Information computed with the Cauchy-Schwartz divergence
ICS(f, g) = -log (Int f*g / sqrt (Int f^2*Int g^2)), Int means the integral

hr2 :  R�nyi quadratic entropy hr2
gw : Parzen windowing for continuous variable density estimation, i.e. sum of local interactions,
defined by the kernel K (gaussian here 2*pi*h^2 exp(-(x-xi)^2/ 2*h^2)), over all pairs of samples :
gw(x) = -log(1/n^2 sum_i K(x-xi, h) with h the window width or smoothing parameter
Rq : Can converge to the true density if n goes to infinity
h' = 1.06*s*n^(1/5) =>  estimated with sample standard deviation s or the interquartile range Iq,  h' = 0.9*min(s, Iq)*n^(1/5)
In the bivariate case h' = 0.85 min((s1^2+s2^2/2)^(-1/2), Iq1+Iq2/2 ) n^(1/6)

ICS(X,Y) = Divergence between joint density fxy and marginal densities fx and fy
ICS(X,Y) = - log (Int fxy*fx*fy / sqrt(Int fxy^2 * Int fx^2*fy^2)) = hr2(fxy*fx*fy)-0.5*hr2(fxy)-0.5*hr2(fx*fy)
with hr2(X) = -log sum f^2(x) = -log(1/n^2 sum_i sum_j K(xi-xj, 2*h'^2), n is the number of samples
"""

import numpy as np
from collections import defaultdict
from functools import reduce
import numpy.lib.arraysetops as at
from operator import mul
from joblib import Parallel, delayed

from sklearn.utils import (check_array, check_X_y)
from sklearn.neighbors import NearestNeighbors
from sklearn.feature_selection import (GenericUnivariateSelect, SelectPercentile, SelectKBest)
from sklearn.feature_selection.univariate_selection import _BaseFilter
from sklearn.utils.validation import check_is_fitted

try: import psyco; psyco.full()
except: pass


def ConvertToContinue(c, sigma = 0.01):
    '''
    Convert a discrete variable in continuous variable by applying a gaussian distribution in each point
    Parameters
    ----------
        c=discrete variable
        sigma=standard deviation of the gaussian distribution
    Returns
    -------
        newc=continuous variable
    '''
    cu = at.unique(c)
    newc = c.copy()
    newc = newc.astype(float)
    for cui in cu :
        ind = np.where(c==cui)[0]
        newc[ind] = sigma*np.random.randn(len(ind))+cui
    #MP.plot(c, '.')
    #MP.plot(newc, 'r+')
    #MP.ylim(cu[0]-0.5, cu[-1]+05)
    return newc

def KNearestNeighbors(X, k):
    '''
    Fit the KNearestNeighbors model using X as training data
    Parameters
    ----------
        X: sample, array of shape =[n_samples, n_features]
        k: number of neighbors
    returns
    -------
    if k>0 and k < n_samples return the fitted model , otherwise all the samples are considered as neighbors
    '''
    if k > 0 and k < X.shape[0]:
        neigh = NearestNeighbors(n_neighbors=k)
        neigh.fit(X)
        return neigh
    else: return range(X.shape[0])


def knn(s, neigh):
    '''
    Return the indices of the k nearest neighbors of sample s
    otherwise all the indices are considered
    '''
    if type(neigh)!=range: return neigh.kneighbors(s, return_distance=False)[0]
    return neigh

def Parallel_MI_RenyiCC_d(i, freqsi,xu,yu,xc,yc, N):
    xui = np.where(xu==i[0])[0][0]
    yui = np.where(yu==i[1])[0][0]
    hr2a = freqsi*(xc[xui]/N)*(yc[yui]/N)
    hr2b = freqsi**2
    return [hr2a, hr2b]

def Parallel_MI_RenyiCC_d_Multi(i, freqsi,u, N):
    uci = [(u[uv,1][np.where(u[uv,0]==i[uv])[0][0]])/N for uv in range(u.shape[0])]
    hr2a = freqsi*reduce(mul, uci)
    hr2b = freqsi**2
    return [hr2a, hr2b]

def Sum_Dot_Vect(uc, N):
    if len(uc)>2:
        return np.sum(((uc[0]/N)**2)*Sum_Dot_Vect(uc[1:]))
    elif len(uc)==2 :
        S =  np.sum(np.dot(np.reshape((uc[0]/N)**2,(len(uc[0]),1)), np.reshape((uc[1]/N)**2,(1,len(uc[1])))))
        return S

def Parallel_MI_RenyiCC_c_Multi(i, X, h):
    #print("sample:",i, " - neighbors:",X)
    pw = [ParzenWindow(i[v]-X[:,v], h) for v in range(len(i))]
    hr2a = reduce(mul,pw)
    w = zip(*(i-X).T)
    hr2b = ParzenWindow(w, h, len(i))
    return [hr2a, hr2b] + pw

def Parallel_MI_RenyiCC_cd_Multi(xyui, X, hx, neigh, k):
    nxyui = xyui.shape[0]
    #pwa = np.array([[ParzenWindow(xyui[i,j]-X[:,j], hx) for i in range(nxyui)] for j in range(X.shape[1])])#nbfeatures*nbsamples
    pwa = np.array([[ParzenWindow(xyui[i,j]-X[knn(xyui[i,:], neigh),j], hx) for i in range(nxyui)]\
                                                                        for j in range(X.shape[1])])#nbfeatures*nbsamples
    hr2a = nxyui*np.sum(reduce(mul,pwa))
    #if X.shape[1]>1: hr2b = np.sum(ParzenWindow(zip(*(i-xyui).T),hx, X.shape[1]) for i in xyui)
    #else : hr2b = np.sum(ParzenWindow(i-xyui,hx) for i in xyui)
    subneigh = KNearestNeighbors(xyui, k)
    if X.shape[1]>1: hr2b = np.sum(ParzenWindow(zip(*(i-xyui[knn(i, subneigh),:]).T),hx, X.shape[1]) for i in xyui)
    else : hr2b = np.sum(ParzenWindow(i-xyui[knn(i, subneigh),:],hx) for i in xyui)
    nxyu = nxyui**2
    return [hr2a, hr2b, nxyu]

def Parallel_MI_RenyiCC_cd_Multi_hr2c_dim0(i, X, hx):
    return [ParzenWindow(i[j]-X[:,j], hx) for j in range(X.shape[1])]

def Parallel_MI_RenyiCC_cd_Multi_hr2c_dim1(Xj, hx, neigh):
    #return np.sum(ParzenWindow(xi-Xj, hx) for xi in Xj)
    return np.sum(ParzenWindow(xi-Xj[knn(xi, neigh),:], hx) for xi in Xj)

def check_array_type(X, y, type):
    if y is not None:
        check_X_y(X, y)
        if type in ['cd', 'd']:
            yu, yc = at.unique(y, return_counts=True)
            if sum(yc)!=len(y): np.testing.assert_equal(y.dtype, np.int)
            elif y.dtype!=np.int: y = np.asarray(y, np.int)
        if type == 'c': np.testing.assert_equal(y.dtype, np.float)
        if type in ['c','d']:
            X = np.hstack((X,y.reshape((len(y),1))))
            y = None
    else :
        check_array(X)
        if type == 'c': np.testing.assert_equal(X.dtype, np.float)
        if type == 'd': np.testing.assert_equal(X.dtype, np.int)
        if type=='cd' : raise ValueError("y has to be defined for type cd")
    return X, y


def MI_RenyiCC_Multi(X,y=None, k=0, type='c', njobs=4):
    """
    Mutual Information estimator based on the Renyi quadratic entropy and the Cauchy Schwartz divergence
    Parameters
    ----------
        X = data of shape = [n_samples, n_features]
        type = type of the computation according to the variable types
             'd' for discrete variables ,'c' for continuous variables (by default) and 'cd' for estimating MI of
             continuous variables with a discrete target y
        y = discrete target (for classification study), array of shape(n_samples)
        k = the number of neighbors to considered for the parzen window esimation (0 by default means that we considered
            all the samples)
        njobs = number of parallel job for computation (4 by default)

    Returns
    -------
        MI_QRCS = Mutual Information score , i.e. equal to 0 if variables in X are independant

    """
    N = X.shape[0]
    X, y = check_array_type(X, y, type)
    if type == 'd':
        u = np.array([at.unique(x, return_counts=True) for x in X.T])
        freqs = DiscDensity(zip(*X.T), N)
        hr2c = Sum_Dot_Vect(u[:,1], N)
        hr2 = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_d_Multi)(i, freqs[i],u, N) for i in freqs)
        s = np.sum(np.array(hr2),0)
        hr2a = s[0]
        hr2b = s[1]
        #print("hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c)
    elif type == 'c' :
        neigh = KNearestNeighbors(X, k)
        iqrx = [np.subtract(*np.percentile(x, [75, 25])) for x in X.T]
        varx = [np.var(x) for x in X.T]
        h = 0.85*min( 1/np.sqrt(np.mean(varx)), np.mean(iqrx))*N**(-1/6)
        #hr2 = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_c_Multi)(i, X, h) for i in zip(*X.T))
        hr2 = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_c_Multi)(i, X[knn(i, neigh),:], h) for i in zip(*X.T))
        s = np.sum(np.array(hr2),0)
        hr2a = s[0]
        hr2b = s[1]
        pw = [s[i] for i in range(2,len(s))]
        hr2c = (1/N**4)*reduce(mul, pw)
        hr2a = (1/N**3)*hr2a
        hr2b = (1/N**2)*hr2b
        #print("hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c)
    elif type == "cd":
        yu, yc = at.unique(y, return_counts=True)
        #hr2y = -np.log(np.sum((yc/N)**2))
        if X.shape[1]==1: hx = 0.9*min(np.std(X),np.subtract(*np.percentile(X, [75, 25])))*N**(-1/5)**2
        else :
            iqrx = [np.subtract(*np.percentile(x, [75, 25])) for x in X.T]
            varx = [np.var(x) for x in X.T]
            hx = 0.85*min( 1/np.sqrt(np.mean(varx)), np.mean(iqrx))*N**(-1/6)
        xyu = defaultdict(list)
        z = zip(*np.hstack((X,np.reshape(y,(N,1)))).T)
        for i in z: xyu[int(i[-1:][0])].append(i[:-1])
        neigh = KNearestNeighbors(X, k)
        hr2 = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_cd_Multi)(np.array(xyu[yui]), X, hx, neigh, k) for yui in yu)
        s = np.sum(np.array(hr2),0)
        hr2a = s[0]
        hr2b = s[1]
        nxyu = s[2]
        #Parallelize loop according to the biggest dimension between the number of samples and the number of features
        #Notes : by using the knn, the two parallelize estimations are different since the first compute the knn of the
        #ith sample through all the features dimension while the 2nd compute knn of the ith sample along each feature
        #dimension
        if X.shape[0]>X.shape[1]:
            #hr2cp = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_cd_Multi_hr2c_dim0)(i, X, hx) for i in range(N))
            hr2cp = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_cd_Multi_hr2c_dim0)\
                                                            (i, X[knn(i, neigh),:], hx) for i in X)
            hr2cp = reduce(mul,np.sum(np.array(hr2cp),0))
            hr2c = (1/N**4)*nxyu*hr2cp
        else:
            hr2cp = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_cd_Multi_hr2c_dim1)(X[:,j], hx, neigh) for j in range(X.shape[1]))
            hr2cp = reduce(mul,hr2cp)
            hr2c = (1/N**4)*nxyu*hr2cp
        #print("hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c)
        hr2a = (1/N**3)*hr2a
        hr2b = (1/N**2)*hr2b
        #hr2x = -np.log((1/N**2)*hr2x)

    hr2a = max(10**(-100), hr2a)
    hr2b = max(10**(-100), hr2b)
    hr2c = max(10**(-100), hr2c)
    #print("hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c)
    lhr2a = -np.log(hr2a)
    lhr2b = -np.log(hr2b)
    lhr2c = -np.log(hr2c)
    MI_QRCS =  lhr2a - 0.5*lhr2b - 0.5*lhr2c
    return MI_QRCS


def MI_RenyiCC(x, y, type, njobs=4):
    """
    Mutual Information estimator based on the Renyi quadratic entropy and the Cauchy Schwartz divergence
    Compute Renyi Quadratic Entropies hr2(p(x,y)*p(x)*p(y)), hr2 p(x,y) and hr2 p(x)p(y) for all types of variables couple
    Parameters
    ----------
        x, y = two variables
        type = type of the computation according to the variable types
             'dd' for 2 discret variables ,'cc' for 2 continue variables or 'cd' for 2 mixed variables
        njobs = number of parallel job for computation (4 by default)
    Returns :
        MI_QRCS = hr2(p(x,y)*p(x)*p(y))-1/2hr2 p(x,y) - 1/2hr2 p(x)p(y) , i.e. equal to 0 if x and y are independant

    Notes
    -----
    MI_RenyiCC_Multi may be used for bivariate variable => could be removed
    """
    N = len(x)
    if type == 'dd':
        xu, xc = at.unique(x, return_counts=True)
        yu, yc = at.unique(y, return_counts=True)
        hr2x = -np.log(np.sum((xc/N)**2))
        hr2y = -np.log(np.sum((yc/N)**2))
        freqs = DiscDensity(zip(x,y), N)
        hr2c = np.sum(np.dot(np.reshape((yc/N)**2,(len(yc),1)), np.reshape((xc/N)**2,(1,len(xc)))))
        hr2 = Parallel(n_jobs=njobs)(delayed(Parallel_MI_RenyiCC_d)(i, freqs[i],xu, yu,xc,yc, N) for i in freqs)
        s = np.sum(np.array(hr2),0)
        hr2a = s[0]
        hr2b = s[1]
        #print("hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c)
    elif type == 'cc' :
        hr2a = 0; hr2b=0; hr2c = 0
        iqrx = np.subtract(*np.percentile(x, [75, 25]))
        iqry = np.subtract(*np.percentile(y, [75, 25]))
        h = 0.85*min(1/np.sqrt((np.var(x)+np.var(y))/2),(iqrx+iqry)/2)*N**(-1/6)
        hr2x = 0; hr2y = 0
        pwX = 0; pwY = 0
        for i in zip(x,y) :
            hr2x += ParzenWindow(i[0]-x, 0.9*min(np.std(x),iqrx)*N**(-1/5)**2)
            hr2y += ParzenWindow(i[1]-y, 0.9*min(np.std(x),iqrx)*N**(-1/5)**2)
            pwx = ParzenWindow(i[0]-x, h)
            pwy = ParzenWindow(i[1]-y, h)
            hr2a += pwx*pwy
            w = zip(i[0]-x,i[1]-y)
            hr2b += ParzenWindow(w, h, 2)
            pwX += pwx; pwY += pwy
        hr2c += (1/N**4)*(pwX * pwY)
        hr2a = (1/N**3)*hr2a
        hr2b = (1/N**2)*hr2b
        #print("-hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c,"-pw:",[pwX,pwY])
        hr2x = -np.log((1/N**2)*hr2x)
        hr2y = -np.log((1/N**2)*hr2y)
    elif type == "cd":
        yu, yc = at.unique(y, return_counts=True)
        hr2y = -np.log(np.sum((yc/N)**2))
        xyu = defaultdict(list)
        iqrx = np.subtract(*np.percentile(x, [75, 25]))
        hx = 0.9*min(np.std(x),iqrx)*N**(-1/5)**2
        hr2x = 0; hr2a = 0; hr2b = 0; hr2c =0
        nxyu = 0
        for i in zip(x,y):
            xyu[i[1]].append(i[0])
            hr2x += ParzenWindow(i[0]-x, hx)
        for yui in yu:
            nxyui = len(xyu[yui])
            varxyui = np.var(xyu[yui])
            iqrxyui = np.subtract(*np.percentile(xyu[yui], [75, 25]))
            h = 0.85*min(1/np.sqrt((np.var(x)+varxyui)/2),(iqrx+iqrxyui)/2)*N**(-1/6)
            hr2a += nxyui*np.sum(ParzenWindow(j-x, hx) for j in xyu[yui])
            hr2b += np.sum(ParzenWindow(j-xyu[yui],hx) for j in xyu[yui])
            nxyu += nxyui**2
        hr2c = (1/N**4)*nxyu*np.sum(ParzenWindow(xi-x, hx) for xi in x)
        #print("hr2a:",hr2a,"-hr2b:",hr2b,"-hr2c:",hr2c)
        hr2a = (1/N**3)*hr2a
        hr2b = (1/N**2)*hr2b
        hr2x = -np.log((1/N**2)*hr2x)

    lhr2a = -np.log(hr2a)
    lhr2b = -np.log(hr2b)
    lhr2c = -np.log(hr2c)
    MI_QRCS =  lhr2a - 0.5*lhr2b - 0.5*lhr2c
    return MI_QRCS

def DiscDensity(var, N):
    """
    Estimation of the frequencies of data values
    Parameters
    ----------
        var = bivariate variable (x,y)
        N = number of data samples
    Returns
    -------
        freqs = Frequency of each bivariate data and each univariate data
    """
    freqs = defaultdict(int)
    for i in var: freqs[i] += 1
    for k in freqs.keys() :freqs[k] /= N
    return freqs

def ParzenWindow(w, h, d=1):
    """""""""""""""""""""""
    Average of the gaussian window functions centered on each data point of w for marginal or joint density
    estimation at a particular point x
    Parameters
    ----------
        w = vector of distances from a point x to the other points
        h = window width
        d = length or the variable dimension (1 by default if unit variable for marginal density and 2 if bivariate
            variable for joint density.)
    Returns
    -------
        pw = Estimation of the parzen window function for the density estimation f(w)
    """""""""""""""""""""""
    if d>1:
        pw = np.sum(np.prod(GaussianWindow(list(w),h),1))
    else: pw = np.sum(GaussianWindow(w,h))#np.sum(np.exp(-w**2/(2*phi))/den)
    return pw

def GaussianWindow(w,h):
    """
    Gaussian kernel function with a variance of 2h^2
    :param w: vector of distances from a point x to the other points
    :param h: Window width
    :return: VAlue of the gaussian function
    """
    phi = 2*h**2
    den = (2*np.pi*phi)**(1/2)
    return den*np.exp(-np.power(w,2)/(2*phi))

class SelectKFirst(_BaseFilter):
    """Select features according to the k first indices in ranking array.
    If k<0 : return the k last indices corresponding to the k first in a backward features selection
    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y, and returning a pair of arrays
        (scores, ranking).

    k : int or "all", optional, default=10
        Number of top features to select.
        The "all" option bypasses selection, for use in a parameter search.

    Attributes
    ----------
    scores_ : array-like, shape=(n_features,)
        Scores of features.

    ranking_ : array-like, shape=(n_features,)
        ranking of feature score_func function.
    """

    def __init__(self, score_func, k=10):
        super(SelectKFirst, self).__init__(score_func)
        self.k = k

    def _check_params(self, X, y):
        if not (self.k == "all" or -X.shape[1] <= self.k <= X.shape[1]):
            raise ValueError("k should be >=0, <= n_features to select the k"
                             "first features or >=-n_features, < 0 to select the k"
                             "last features; got %r."
                             "Use k='all' to return all features."
                             % self.k)

    def _get_support_mask(self):
        check_is_fitted(self, 'scores_', 'ranking_')

        if self.k == 'all':
            return np.ones(self.ranking_.shape, dtype=bool)
        elif self.k == 0:
            return np.zeros(self.ranking_.shape, dtype=bool)
        else:
            mask = np.zeros(self.scores_.shape, dtype=bool)
            if self.k<0 : mask[self.ranking_[self.k:]] = 1
            else : mask[self.ranking_[:self.k]] = 1
            return mask

class  _MI_Filter(GenericUnivariateSelect):
    """Initialize the univariate feature selection based on MI.

    Parameters
    ----------
    score_func : callable
        Function taking two arrays X and y and optional parameters, and returning a
        scores and a ranking array. Ranking is obtained if the score_func function
         integrates a features selection procedure as forward or backward algorithms,
         otherwise,ranking array is simply the features range corresponding to the scores array.
     """

    _selection_modes = {'percentile':   SelectPercentile,
                        'k_best':       SelectKBest,
                         'k_first':     SelectKFirst}

    def __init__(self, score_func, mode='percentile', param=1e-5, **params):
        super(_MI_Filter, self).__init__(score_func, mode, param)
        self.params = params

    def fit(self, X, y):
        """Run score function on (X, y) with optional parameters and get the appropriate features.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            The training input samples.

        y : array-like, shape = [n_samples]
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        X, y = check_X_y(X, y, ['csr', 'csc', 'coo'])

        if not callable(self.score_func):
            raise TypeError("The score function should be a callable, %s (%s) "
                            "was passed."
                            % (self.score_func, type(self.score_func)))

        self._check_params(X, y)
        self.scores_, self.ranking_ = self.score_func(X, y,**self.params)
        self.scores_ = np.asarray(self.scores_)
        self.ranking_ = np.asarray(self.ranking_)
        return self

    def _get_support_mask(self):
        check_is_fitted(self, 'scores_', 'ranking_')
        selector = self._make_selector()
        selector.scores_ = self.scores_
        selector.ranking_ = self.ranking_
        return selector._get_support_mask()


def univariate_f_MI(X, y, **params): #BaseEstimator
    """Compute bivariate MI for each output/feature combination.

    This score can be used to select the n_features features with the
    highest values according to MI from X relative to the output.

    Recall that the MI measures dependence between variables, so using
    this function "weeds out" the features that are the
    most likely to be independent of the output and therefore irrelevant for
    classification/regression purposes.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = (n_samples, n_features_in)
        Sample vectors.

    y : array-like, shape = (n_samples,)
        Target vector (class or continuous labels).

    Returns
    ----------
    scores : array, shape = (n_features,)
        MI of each feature.


    """
    p = X.shape[1]
    scores = np.zeros(p)
    for i in range(p):
        scores[i] = MI_RenyiCC_Multi(X[:,[i]],y, **params)
    return scores, np.arange(p)

def multivariate_forward_f_MI(X, y, **params):
    """Compute a multivariate forward (greedy) features selection based on MI scores
    At each step, we add the most important feature x maximizing the relevance with the output and the
    set of previous selected features
    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = (n_samples, n_features_in)
        Sample vectors.

    y : array-like, shape = (n_samples,)
        Target vector (class or continuous labels).

    Returns
    -------
    scores : array, shape = (n_features,); scores used for features selection at each step
    S : array, shape = (n_features,); Features Ranking.
    """
    p = X.shape[1]
    scores = np.zeros(p)
    ranking = []; R = np.arange(p)
    scores = []
    while len(ranking)!=p:
        if len(ranking)>0 and len(ranking)<p-1:
            rel_scores = np.array([MI_RenyiCC_Multi(np.hstack((X[:,ranking],x.reshape((len(x),1)))), y, **params) for x in X[:,R].T])
        elif len(ranking)==0 :
            rel_scores, r = univariate_f_MI(X, y, **params)
        else:
            rel_scores = MI_RenyiCC_Multi(X, y, **params)
        pos = np.argmax(rel_scores)
        m = R[pos]
        sc = np.max(rel_scores)
        ranking.append(m)
        R = np.delete(R, pos)
        scores.append(sc)
    return np.array(scores), np.array(ranking)

def univariate_forward_f_MI(X, y, **params):
    """Compute a forward (greedy) features selection based on MI scores
    At each step, we add the feature x maximizing the relevance with the output (i.e. max I(x,y))
    and minimizing the redundancy with the previous selected features
    References
    ----------
        ... [1] F.H.Long, H.C.Peng, and C.Ding, Feature selection based on mutual
                information: criteria of max-dependency ,max-relevance, and minredundancy,
                IEEE Transactions on Pattern Analysis and Machine Intelligence,
                vol. 27, no. 8, pp. 1226–1238, August 2005.
        ... [2]  R. Battiti, Using mutual information for selection features in
                 supervised neural net learning,
                IEEE Transactions on Neural Networks, vol. 5, no. 4, pp. 537–550, 1994

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = (n_samples, n_features_in)
        Sample vectors.

    y : array-like, shape = (n_samples,)
        Target vector (class or continuous labels).

    Returns
    ----------
    scores : array, shape = (n_features,); scores used for features selection at each step
    S : array, shape = (n_features,); Features Ranking.
    """
    p = X.shape[1]
    scores = np.zeros(p)
    ranking = []; R = np.arange(p)
    scores = []
    nparams = params.copy()
    if nparams['type'] == ['cd']: nparams['type']='c'
    rel_scores, r = univariate_f_MI(X, y, **params)
    print("rel_scores:", rel_scores)
    '''
    while len(ranking)!=p:
        if len(ranking)>0 and len(ranking)<p-1:
            red_scores = []
            for si in ranking:
                temp, r = univariate_f_MI(X[:,R] , X[:,si], **nparams)
                red_scores.append(temp)
                print("red with selected feature:",si ,":", temp)
            sum_red_scores = np.sum(np.array(red_scores),0)/len(ranking)
            print("R:",R,"-rel_scores-red_scores:", rel_scores[R]- sum_red_scores)
            pos = np.argmax(rel_scores[R]- sum_red_scores)
            m = R[pos]
            sc = np.max(rel_scores[R]- sum_red_scores)
        else :
            pos = np.argmax(rel_scores[R])
            m = R[np.argmax(rel_scores[R])]
            sc = np.max(rel_scores[R])
        ranking.append(m)
        R = np.delete(R, pos)
        scores.append(sc)
    '''
    old_red_scores = []; old_ranking=[]; old_R = []
    while len(ranking)!=p:
        if len(ranking)>0 and len(ranking)<p-1:
            if len(old_red_scores) == 0 : ranking_temp = ranking
            else : ranking_temp = list(set(ranking)- set(old_ranking))
            red_scores = []
            if len(old_ranking)>0:
                [red_scores.append(np.delete(old_red_scores[i], np.where(np.in1d(old_R[i], ranking[i+1:], assume_unique=True))[0])) for i in range(len(old_ranking))]
            for si in ranking_temp:
                temp, r = univariate_f_MI(X[:,R] , X[:,si], **nparams)
                red_scores.append(temp)
                old_red_scores.append(temp)
                old_R.append(R)
                old_ranking.append(si)
            sum_red_scores = np.sum(np.array(red_scores),0)/len(ranking)
            pos = np.argmax(rel_scores[R]- sum_red_scores)
            m = R[pos]
            sc = np.max(rel_scores[R]- sum_red_scores)
        else :
            pos = np.argmax(rel_scores[R])
            m = R[np.argmax(rel_scores[R])]
            sc = np.max(rel_scores[R])
        ranking.append(m)
        ranking_temp = ranking
        R = np.delete(R, pos)
        scores.append(sc)
    return np.array(scores), np.array(ranking)

def multivariate_backward_f_MI(X, y, **params):
    """Compute a multivariate backward features selection based on MI scores
    At each step, we remove the least important feature x , i.e.
    maximizing the relevance with the output and the set of remaining features.
    The backward algo has the advantage that interactions among all remaining features are considered.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = (n_samples, n_features_in)
        Sample vectors.

    y : array-like, shape = (n_samples,)
        Target vector (class or continuous labels).

    Returns
    -------
    scores : array, shape = (n_features,); scores used for features selection at each step
    S : array, shape = (n_features,); Features Ranking.
    """
    p = X.shape[1]
    scores = np.zeros(p)
    ranking = []; R = np.arange(p)
    scores = []
    while len(ranking)!=p:
        if len(R)>1:
            rel_scores = np.array([MI_RenyiCC_Multi(np.delete(X,ranking+[i],1), y, **params) for i in R ])
        else:
            rel_scores = MI_RenyiCC_Multi(X, y, **params)
        pos = np.argmax(rel_scores)
        m = R[pos]
        sc = np.max(rel_scores)
        ranking.append(m)
        R = np.delete(R, pos)
        scores.append(sc)
    return np.array(scores), np.array(ranking)

def get_support(scores, k):
    mask = np.zeros(scores.shape, dtype=bool)
    mask[np.argsort(scores, kind="mergesort")[-k:]] = 1
    return mask

