
"""
Dimensionality reduction with geodesic distances
"""

import os
import numpy
import numpy.random
import numpy.linalg
import math

#from scikits.optimization import *

from .distances import kneigh as distances_kneigh
from .distances import numpy_floyd
from .euclidian_mds import mds as euclidian_mds
#from .cca_function import CostFunction as CCA_CostFunction
#from .cost_function import CostFunction as RobustCostFunction
#from .NLM import NLM_CostFunction

#from .multiresolution_dimensionality_reduction import \
#    optimize_cost_function as multiresolution_optimize_cost_function
#from .cca_multiresolution_dimensionality_reduction import \
#    optimize_cost_function as cca_multiresolution_optimize_cost_function
#from .robust_dimensionality_reduction import optimize_cost_function 
#    as robust_dimensionality_optimize_cost_function

def reduct(reduction, function, samples, **kwargs):
    """
    Data reduction with geodesic distance approximation

    Parameters
    ----------
    reduction:
      The reduction technique

    samples : matrix
      The points to consider.

    temp_file : string
      name of a file for caching the distance matrix

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this function will be passed to its constructor.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not given.
    """
    if 'temp_file' in kwargs and os.path.exists(kwargs['temp_file']):
        dists = numpy.fromfile(kwargs['temp_file'])
        size = int(math.sqrt(dists.shape[0]))
        dists.shape = (size, size)
    else:
        if 'neigh' in kwargs:
            neighborer = kwargs['neigh'](samples, **kwargs)
        else:
            neighborer = distances_kneigh(samples, kwargs.get('n_neighbors', 9))

        dists = populateDistanceMatrixFromneighbors(samples, neighborer)
        numpy_floyd(dists)
        if 'temp_file' in kwargs:
            dists.tofile(kwargs['temp_file'])
        del neighborer

    return reduction(dists, function, **kwargs)

def populateDistanceMatrixFromneighbors(points, neighborer):
    """
    Creates a matrix with infinite value safe for points that are neighbors
    """
    distances = numpy.ones((points.shape[0], points.shape[0]),
        dtype = numpy.float)
    distances *= 1e30000
    for indice in xrange(0, len(points)):
        neighborList = neighborer[indice]
        for element in neighborList:
            distances[indice, element] = math.sqrt(
                numpy.sum((points[indice] - points[element])**2))
            distances[element, indice] = math.sqrt(
                numpy.sum((points[indice] - points[element])**2))

    return distances

class Isomap(object):
    """
    Isomap embedding object

    Parameters
    ----------
    temp_file : string
      name of a file for caching the distance matrix

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this function will be passed to its constructor.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not given.

    Attributes
    ----------
    embedding_ : array_like
        Embedding of the learning data
    
    X_ : array_like
        Original data that is embedded

    See Also
    --------

   
    Notes
    -----
    
    .. [1] Tenenbaum, J. B., de Silva, V. and Langford, J. C.,
           "A Global Geometric Framework for Nonlinear Dimensionality 
           Reduction,"
           Science, 290(5500), pp. 2319-2323, 2000
    
    Examples
    --------  
    
    """
    def __init__(self, **embedded_opts):
        self.__embedded_opts = embedded_opts

    def transform(self, X):
        """
        Parameters
        ----------
        X : array_like
        The learning dataset
        
        Returns
        -------
        Self
        """
        def function(*args, **kwargs):
            return None
        self.X_ = X
        self.embedding_, reduced_parameter_set = reduct(euclidian_mds, 
            function, X, **self.__embedded_opts)
        return self

def ccaCompression(samples, nb_coords, **kwargs):
    """
    CCA compression :
      - samples is an array with the samples for the compression
      - nb_coords is the number of coordinates that must be retained
      - temp_file is a temporary file used for caching the distance matrix
      - neigh is the neighboring class that will be used
      - neighbors is the number of k-neighbors if the K-neighborhood is used
      - window_size is the window size to use
      - max_dist is the maximum distance to preserve
    """
    return reduct(cca_multiresolution_optimize_cost_function,
        CCA_CostFunction, samples, nb_coords, **kwargs)

def robustCompression(samples, nb_coords, **kwargs):
    """
    Robust compression :
      - samples is an array with the samples for the compression
      - nb_coords is the number of coordinates that must be retained
      - temp_file is a temporary file used for caching the distance matrix
      - neigh is the neighboring class that will be used
      - neighbors is the number of k-neighbors if the K-neighborhood is used
      - window_size is the window size to use
    """
    return reduct(robust_dimensionality_optimize_cost_function,
        RobustCostFunction, samples, nb_coords, **kwargs)

def robustMultiresolutionCompression(samples, nb_coords, **kwargs):
    """
    Robust multiresolution compression :
      - samples is an array with the samples for the compression
      - nb_coords is the number of coordinates that must be retained
      - temp_file is a temporary file used for caching the distance matrix
      - neigh is the neighboring class that will be used
      - neighbors is the number of k-neighbors if the K-neighborhood is used
      - window_size is the window size to use
    """
    return reduct(multiresolution_optimize_cost_function, RobustCostFunction,
        samples, nb_coords, **kwargs)

def geodesicNLM(samples, nb_coords, **kwargs):
    """
    Data reduction with NonLinear Mapping algorithm
    (JR. J. Sammon. A nonlinear mapping for data structure analysis.
    IEEE Transactions on Computers, C-18(No. 5):401--409, May 1969):
      - samples is an array with the samples for the compression
      - nb_coords is the number of coordinates that must be retained
    Geodesic distances are used here.
    """
    return reduct(optimize_cost_function,
        NLM_CostFunction, samples, nb_coords, **kwargs)

