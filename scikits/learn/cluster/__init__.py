"""
Clustering algorithms
"""

from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import mean_shift, MeanShift, estimate_bandwidth
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .k_means_ import k_means, KMeans
from .som_ import SelfOrganizingMap

import numpy as np

def pseudo_F(X, labels, centroids):
    '''
    The pseudo F statistic :

    pseudo F = [( [(T - PG)/(G - 1)])/( [(PG)/(n - G)])] 

    The pseudo F statistic was suggested by Calinski and Harabasz (1974)

    Calinski, T. and J. Harabasz. 1974. 
    A dendrite method for cluster analysis. Commun. Stat. 3: 1-27.
    http://dx.doi.org/10.1080/03610927408827101
    '''
    mean = np.mean(X,axis=0) 
    B = np.sum([ (c - mean)**2 for c in centroids])
    W = np.sum([ (x-centroids[labels[i]])**2 
                 for i, x in enumerate(X)])
    c = len(centroids)
    n = len(X)
    return (B /(c-1))/(W/ (n-c))
    
