"""
Mutual Information estimator based on the Renyi quadratic entropy and the Cauchy Schwartz divergence
combined with the Parzen window density estimator for continuous variable

References :
-P. Operacional, R�nyi entropy and Cauchy-Schwartz mutual information applied ti mifs-u variable selection
algorithm : a comparative study. 2011. ISSN
-D. XU and D. Erdogmuns, Renyi's entropy, divergence and their nonparametric estimators

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

# Authors: Cecilia Damon <cecilia.damon@institut-hypercube.org>
import numpy as np
from numpy.linalg import pinv
from collections import defaultdict
from functools import reduce
import numpy.lib.arraysetops as at
from itertools import starmap

from joblib import Parallel, delayed
try: import psyco; psyco.full()
except: pass


def Density(var, N):
    d = defaultdict(int)
    for i in var:
        d[i] += 1
        for j in i : d[i[j]] += 1
        np.vstack((x,y)).T
    res=-sum(map(logPos,freqs))

def ConvertToContinue(c, sigma = 0.01):
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

def Parallel_MI_RenyiCC_d(i, freqsi,xu,yu,xc,yc, N):
    xui = np.where(xu==i[0])[0][0]
    yui = np.where(yu==i[1])[0][0]
    hr2a = freqsi*(xc[xui]/N)*(yc[yui]/N)
    hr2b = freqsi**2
    return [hr2a, hr2b]

def MI_RenyiCC(x, y, type, njobs=4):
    """
    Mutual Information estimator based on the Renyi quadratic entropy and the Cauchy Schwartz divergence
    Compute Renyi Quadratic Entropies hr2(p(x,y)*p(x)*p(y)), hr2 p(x,y) and hr2 p(x)p(y) for all types of variables couple
    INPUT :
        x, y = two variables
        type = type of the computation according to the variable types
             'dd' for 2 discret variables ,'cc' for 2 continue variables or 'cd' for 2 mixed variables
        njobs = number of parallel job for computation (4 by default)
    OUTPUT :
        MI_QRCS = hr2(p(x,y)*p(x)*p(y))-1/2hr2 p(x,y) - 1/2hr2 p(x)p(y) , i.e. equal to 0 if x and y are independant
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
    INPUT :
        var = bivariate variable (x,y)
        N = number of data samples
    OUTPUT :
        freqs = Frequency of each bivariate data and each univariate data
    """
    freqs = defaultdict(int)
    for i in var: freqs[i] += 1
    for k in freqs.keys() :freqs[k] /= N#**2
    return freqs

def ParzenWindow(w, h, d=1):
    """""""""""""""""""""""
    Average of the gaussian window functions centered on each data point of w for marginal or joint density
    estimation at a particular point x
    INPUT :
        w = vector of distances from a point x to the other points
        h = window width
        d = length or the variable dimension (1 by default if unit variable for marginal density and 2 if bivariate
            variable for joint density.)
    OUTPUT :
        pw = Estimation of the parzen window function for the density estimation f(w)
    """""""""""""""""""""""
    if d>1: pw = sum(starmap(np.dot,GaussianWindow(list(w),h)))
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

def main():
    import matplotlib.pylab as MP
    np.random.seed(4)
    x= np.random.randint(0,3,100)
    #y= np.exp(x+10)  + np.random.randint(0,10, 100)
    y= np.random.randint(0,10, 100)
    y2 = x
    print("MI dd independante:", MI_RenyiCC(x, y, "dd"))
    print("MI dd dependante:", MI_RenyiCC(x, y2, "dd"))

    x= np.random.randn(100)
    y= np.random.randn(100)
    y2 = x
    print("MI cc independante:", MI_RenyiCC(x, y, "cc"))
    print("MI cc dependante:", MI_RenyiCC(x, y2, "cc"))

    x= 5*np.random.randn(100)+5
    y= np.random.randint(0,10, 100)
    y2 = np.round(x)
    print("MI cd independante:", MI_RenyiCC(x, y, "cd"))
    print("MI cd dependante:", MI_RenyiCC(x, y2, "cd"))
    yc = ConvertToContinue(y,0.01)
    y2c = ConvertToContinue(y2, 0.01)
    print("MI cd convert in cc independante:", MI_RenyiCC(x, yc, "cc"))
    print("MI cd convert in cc dependante:", MI_RenyiCC(x, y2c, "cc"))
    # Mutual information between two correlated gaussian variables
    # Entropy of a 2-dimensional gaussian variable
    n = 500#00
    rng = np.random.RandomState(0)
    #P = np.random.randn(2, 2)
    P = np.array([[1, 0], [0.5, 1]])
    C = np.dot(P, P.T)
    U = rng.randn(2, n)
    Z = np.dot(P, U).T
    X = Z[:, 0]
    X = X.reshape(len(X), 1)
    Y = Z[:, 1]
    Y = Y.reshape(len(Y), 1)
    MP.plot(X,Y,'.')
    print("MI cc 1 :", MI_RenyiCC(X, Y, "cc"))

    # Test that our estimators are well-behaved with regards to
    # degenerate solutions
    rng = np.random.RandomState(0)
    x = rng.randn(n)
    X = np.c_[x, x]
    MP.plot(x,x,'.')
    print("MI cc 2 :", MI_RenyiCC(x, x, "cc"))

    # Mutual information between two correlated gaussian variables
    # Entropy of a 2-dimensional gaussian variable
    rng = np.random.RandomState(0)
    #P = np.random.randn(2, 2)
    P = np.array([[1, 0], [.9, .1]])
    C = np.dot(P, P.T)
    U = rng.randn(2, n)
    Z = np.dot(P, U).T
    X = Z[:, 0]
    X = X.reshape(len(X), 1)
    Y = Z[:, 1]
    Y = Y.reshape(len(Y), 1)
    MP.plot(X.ravel(), Y.ravel(),'.')
    print("MI cc 3 :", MI_RenyiCC(X.ravel(), Y.ravel(), "cc"))

if __name__ == "__main__":
    main()

