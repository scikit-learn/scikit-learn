
"""
Dimensionality reduction with similarities
"""

# Matthieu Brucher
# Last Change : 2008-04-15 10:32

import numpy
import numpy.random
import numpy.linalg
import math

__all__ = ['LLE', 'laplacianEigenmap', 'diffusionMap', ]

from similarities import LLE

import similarities
import tools

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
  return similarities.laplacian_map(samples, nb_coords, method=similarities.sparse_heat_kernel, **kwargs)

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
  return similarities.laplacian_map(tools.centered_normalized(samples), nb_coords, method=similarities.normalized_heat_kernel, **kwargs)
