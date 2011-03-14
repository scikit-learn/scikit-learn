# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 23:17:39 2011

@author: vene
"""
from __future__ import division
import numpy as np

class Bunch:
     def __init__(self, **kwds):
         self.__dict__.update(kwds)
         
def dict_argmax(dictionary):
    val = lambda x: x[1]
    return max(dictionary.items(), key=val)[0]
    
class CRO():
    """
    Closeness to Rank One Hierarchical Clustering
    
    Model that clusters the columns of a matrix into a given number of clusters
    by joining at each step the clusters with the largest CRO value.
    
    Parameters
    ----------
        n_comp: int or None
            Target number of components (clusters) to extract.
            Defaults to 1
        
        epsilon: double or None
            Padding value for output matrices. The value influences sparsity
            and NMF convergence time, but does not influence the performance
            of the initialization.
        
    Attributes
    ----------
        components_, data_:
            Output matrices to be used for NMF initialization
        clusters: 
            List of clusters extracted, each one having:
                size: int, the number of columns in the cluster
                data: array, the submatrix corresponding to the cluster
                svd: tuple (u, s, v), the rank-one approximation of the data
    
    Examples
    --------
    The example in the paper outputs the given result
    >>> X = [[1, 2, 0, 3, 1],
    ...      [0, 0, 1, 0, 0],
    ...      [0, 0, 1, 0, 0],
    ...      [2, 4, 2, 6, 3],
    ...      [3, 6, 4, 9, 4],
    ...      [0, 0, 2, 0, 0]]
    >>> X = np.array(X)
    >>> model = CRO(1)
    >>> model.fit(X)
    (0, 1)
    (0, 2)
    (0, 2)
    (0, 1)
    >>> model.clusters[0].data
    array([[1, 2, 3, 1, 0],
           [0, 0, 0, 0, 1],
           [0, 0, 0, 0, 1],
           [2, 4, 6, 3, 2],
           [3, 6, 9, 4, 4],
           [0, 0, 0, 0, 2]])
           
    Notes
    -----
    See the paper "A method of initialization for nonnegative matrix
    factorization" by Yong-Deok Kim and Seungjin Choi, available at:
    http://www.postech.ac.kr/~seungjin/publications/icassp07_ydkim.pdf
    
    """
    def __init__(self, n_comp = 1, epsilon = 1e-5):
        """
        Initializes the CRO model
        """
        self.n_comp = n_comp
        self.clusters = []
        self.epsilon = epsilon
        
    def fit(self, X):
        """
        Clusters the matrix X
        """
        
        # Each column is an individual clusters
        n_samples, n_features = X.shape
        for col in X.T:
            norm = np.linalg.norm(col)
            svd = (col / norm, norm, np.ones(1))
            self.clusters.append(Bunch(size = 1,
                                       data = col,
                                       svd = svd))
        
        for step in xrange(n_features, self.n_comp, -1):
            cros = {}
            for i in xrange(step - 1):
                for j in xrange(i + 1, step):
                    cro = self.pairwise_cro(self.clusters[i], self.clusters[j])
                    cros[i, j] = cro 
            
            pair = dict_argmax(cros)
            print pair   
            t, s = self.clusters[pair[0]], self.clusters[pair[1]]
            self.clusters[pair[0]] = self._merge(t, s)
            del self.clusters[pair[1]]
        
        self.data_ = np.zeros((n_samples, 0))
        self.components_ = np.zeros((self.n_comp, n_features))
        j = 0
        for i, cl in enumerate(self.clusters):
            self.data_ = np.c_[self.data_, cl.svd[1] * cl.svd[0]]
            self.components_[i, j:j+cl.size] = cl.svd[2]
            j += cl.size
        
        self.data_[self.data_ == 0] += self.epsilon
        self.components_[self.components_ == 0] += self.epsilon
        
        
    def _merge(self, target, source):
        """
        Merges two clusters and updates the rank-one approximation
        """
        
        size = target.size + source.size
        data = np.c_[target.data, source.data]
        L = np.c_[target.svd[0] * target.svd[1],
                  source.svd[0] * source.svd[1]]
        R = np.r_[np.c_[target.svd[2], np.zeros(np.shape(target.svd[2]))],
                  np.c_[np.zeros(np.shape(source.svd[2])), source.svd[2]]]
        _, S, V = np.linalg.svd(np.dot(L.T, L))
        S = np.sqrt(S) + 1e-8
        assert (S != 0).all()
        U = np.atleast_2d(np.dot(np.dot(L, V), np.diag(1 / S))) #WORKS, more precision
        #U = np.atleast_2d(np.dot(np.dot(L, V), 1 / S)  #SHOULD WORK
        
        svd = U[:, 0], S[0], np.dot(R, V[0])
        return Bunch(size=size, data=data, svd=svd)
        
    def pairwise_cro(self, u, v):
        """
        CRO between two clusters
        """
        #assert(u.size == v.size == 1)
        #X = np.c_[u.data, v.data]
        #sigma = max(np.linalg.eigvals(np.dot(X.T, X)))
        #return sigma / np.linalg.norm(X) ** 2
        result = self._merge(u, v)
        return (result.svd[1] / np.linalg.norm(result.data)) ** 2      
