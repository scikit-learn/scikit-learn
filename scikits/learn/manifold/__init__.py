# -*- coding: utf-8 -*-

"""
Manifold Learning Module
"""

from .embedding.similarities import LLE, HessianMap
from .embedding.similarities_mds import LaplacianEigenmap, DiffusionMap

from .mapping.barycenter import Barycenter
