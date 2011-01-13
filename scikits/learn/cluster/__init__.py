"""
Clustering algorithms
"""

from .spectral import spectral_clustering, SpectralClustering
from .mean_shift_ import mean_shift, MeanShift, estimate_bandwidth
from .affinity_propagation_ import affinity_propagation, AffinityPropagation
from .k_means_ import k_means, KMeans
from .som_ import SelfOrganizingMap

import numpy as np

def calinski_index(X,labels,centroids):
    mean = np.mean(X,axis=0) 
    B = np.sum([ (c - mean)**2 for c in centroids])
    W = np.sum([ (x-centroids[labels[i]])**2 
                 for i,x in enumerate(X)])
    c = len(centroids)
    n = len(X)
    return (B /(c-1))/(W/ (n-c))
    
