
"""
Dimensionality reduction with similarities
"""

import numpy
import numpy.random
import numpy.linalg
import math

__all__ = ['LLE', 'laplacianEigenmap', 'diffusionMap', ]

from similarities import LLE

from .similarities import laplacian_maps, sparse_heat_kernel, \
    normalized_heat_kernel
from .tools import centered_normalized

def laplacianEigenmap(samples, nb_coords, **kwargs):
    """
    Computes the Laplacian eigenmap coordinates for a set of points
    Parameters:
      - samples are the samples that will be reduced
      - nb_coords is the number of coordinates in the manifold
      - parameter is the temperature of the heat kernel
      - neigh is the neighboorer used (optional, default KNeighboor)
      - neighboor is the number of neighboors (optional, default 9)
    """
    return laplacian_maps(samples, nb_coords, method=sparse_heat_kernel,
        **kwargs)

def diffusionMap(samples, nb_coords, **kwargs):
    """
    Computes the diffusion map coordinates for a set of points
    Parameters:
      - samples are the samples that will be reduced
      - nb_coords is the number of coordinates in the manifold
      - parameter is the temperature of the heat kernel
      - neigh is the neighboorer used (optional, default KNeighboor)
      - neighboor is the number of neighboors (optional, default 9)
    """
    return laplacian_maps(centered_normalized(samples), nb_coords, 
        method=normalized_heat_kernel, **kwargs)
