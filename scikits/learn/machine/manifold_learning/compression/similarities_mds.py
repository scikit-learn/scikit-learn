
"""
Dimensionality reduction with similarities
"""

# Matthieu Brucher
# Last Change : 2007-12-20 15:12

import numpy
import numpy.random
import numpy.linalg
import math

__all__ = ['lle', 'laplacian_eigenmap', 'diffusion_map', ]

from similarities import lle

import similarities
import tools

def laplacian_eigenmap(samples, nbCoords, **kwargs):
  """
  Computes the Laplacian eigenmap coordinates for a set of points
  Parameters:
  - samples are the samples that will be reduced
  - nbCoords is the number of coordinates in the manifold
  - parameter is the temperature of the heat kernel
  - neigh is the neighboorer used (optional, default KNeighboor)
  - neighboor is the number of neighboors (optional, default 9)
  """
  return similarities.laplacian_maps(samples, nbCoords, method=similarities.sparse_heat_kernel, **kwargs)

def diffusion_map(samples, nbCoords, **kwargs):
  """
  Computes the diffusion map coordinates for a set of points
  Parameters:
  - samples are the samples that will be reduced
  - nbCoords is the number of coordinates in the manifold
  - parameter is the temperature of the heat kernel
  - neigh is the neighboorer used (optional, default KNeighboor)
  - neighboor is the number of neighboors (optional, default 9)
  """
  return similarities.laplacian_maps(tools.centered_normalized(samples), nbCoords, method=similarities.normalized_heat_kernel, **kwargs)
