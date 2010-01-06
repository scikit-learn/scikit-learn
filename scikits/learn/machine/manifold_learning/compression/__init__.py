
# Matthieu Brucher
# Last Change : 2008-01-14 09:42

"""
Compression module
"""

from euclidian_mds import *
from geodesic_mds import *
from similarities_mds import *
from pca import *
from similarities import hessian_map

__all__ = ['isomap', 'isomapCompression', 'multiIsomapCompression', 'ccaCompression', 'robustCompression', 'robustMultiresolutionCompression', 'geodesicNLM', 'geodesicRNLM', 'robustGeodesicRNLMCompression',
           'PCA',
           'lle', 'laplacian_eigenmap', 'diffusion_map',
           'hessian_map',
           ]
