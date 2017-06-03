# coding=utf8
"""
BoostCluster

Notes
-----
References:
[1] Yi Liu, Rong Jin, and Anil K. Jain. BoostCluster: Boosting Clustering by Pairwise Constraints.

"""

# Author: Christophe Saint-Jean
# Licence: BSD
import scipy as sc
import numpy as np

from ..base import ClusterMixin


class BoostCluster(ClusterMixin):
    def __init__(self, A, s, T, spectral = True, **kwargs):
        self._A = A
        self._spectral = spectral
        self._s = s
        self._T = T
        self._A_params = kwargs

    def fit(self,X,y=None,Sp,Sm):
        self._K =  sc.zeros_like(Sp)
        for t in range(self._T):
            p = Sp*np.exp(-K)
            q = Sm*np.exp(K)
            sum_p = np.sum(p)
            sum_q = np.sum(q)
            T = (p/sum_p) - (q/sum_q)
            l,v = sc.linalg.eigh(T,self._s)
            P = np.
            Xp = np.dot(P.T,X)
            G = self._A.fit_predict(Xp,self._A_params)
            Delta = np.from_function(lambda i,j: G[i] == G[j], Sp.shape, dtype = np.uint8)
            alpha = 0.5 * (np.log()+np.log()-np.log()-np.log())
            K = K + alpha*Delta
        if self._spectral is True:
            self._A.fit_predict(K,self._A_params)
        else:
            Xp = sc.linalg.eigh(K,self._s+1)
            self._A.fit_predict(Xp,self._A_params)
