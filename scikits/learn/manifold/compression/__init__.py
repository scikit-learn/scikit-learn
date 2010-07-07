
# Matthieu Brucher
# Last Change : 2008-04-15 10:33

"""
Compression module
"""

from euclidian_mds import *
from geodesic_mds import *
from similarities_mds import *
from pca import *
from similarities import hessianMap

__all__ = ['isomap', 'isomapCompression', 'multiIsomapCompression', 'ccaCompression', 'robustCompression', 'robustMultiresolutionCompression', 'geodesicNLM',
           'PCA',
           'LLE', 'laplacianEigenmap', 'diffusionMap',
           'hessianMap',
           ]
