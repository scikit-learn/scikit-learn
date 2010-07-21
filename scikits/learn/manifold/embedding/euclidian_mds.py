
import numpy
import math

from .tools import dist2hd
#from .dimensionality_reduction import optimize_cost_function
#from .NLM import NLM_CostFunction

def reduct(reduction, function, samples, nb_coords, **kwargs):
    """
    Data reduction with euclidian distance approximation:
      - reduction is the algorithm to use
      - function is the function to optimize
      - samples is an array with the samples for the compression
      - nb_coords is the number of coordinates that must be retained
    """
    distances = dist2hd(samples, samples)
    return reduction(distances, function, nb_coords, **kwargs)

def mds(distances, function, nb_coords = None, **kargs):
    """
    Computes a new set of coordinates based on the distance matrix passed as a
    parameter, in fact it is a classical MDS
    """
    if nb_coords is None:
      raise RuntimeError("You must indicate the number of coordinates 'nb_coords'")
    square_distances = -distances ** 2 / 2.
    correlations = square_distances + numpy.mean(
        square_distances) - numpy.mean(square_distances, axis=0) - numpy.mean(
        square_distances, axis=1)[numpy.newaxis].T
    (u, s, vh) = numpy.linalg.svd(correlations)
    return u[:, :nb_coords] * numpy.sqrt(s[:nb_coords]), {
        'scaling' : u[:nb_coords], 'eigen_vectors' : u[:, :nb_coords]}

def NLM(samples, nb_coords, **kargs):
    """
    Data reduction with NonLinear Mapping algorithm
    (JR. J. Sammon. A nonlinear mapping for data structure analysis.
    IEEE Transactions on Computers, C-18(No. 5):401--409, May 1969):
      - samples is an array with the samples for the compression
      - nb_coords is the number of coordinates that must be retained
    """
    return reduct(optimize_cost_function, NLM_CostFunction, samples, nb_coords,
       **kwargs)
