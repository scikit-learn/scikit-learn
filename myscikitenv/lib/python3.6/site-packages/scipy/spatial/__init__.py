"""
=============================================================
Spatial algorithms and data structures (:mod:`scipy.spatial`)
=============================================================

.. currentmodule:: scipy.spatial

Spatial Transformations
=======================

These are contained in the `scipy.spatial.transform` submodule.

Nearest-neighbor Queries
========================
.. autosummary::
   :toctree: generated/

   KDTree      -- class for efficient nearest-neighbor queries
   cKDTree     -- class for efficient nearest-neighbor queries (faster impl.)
   Rectangle

Distance metrics are contained in the :mod:`scipy.spatial.distance` submodule.

Delaunay Triangulation, Convex Hulls and Voronoi Diagrams
=========================================================

.. autosummary::
   :toctree: generated/

   Delaunay    -- compute Delaunay triangulation of input points
   ConvexHull  -- compute a convex hull for input points
   Voronoi     -- compute a Voronoi diagram hull from input points
   SphericalVoronoi -- compute a Voronoi diagram from input points on the surface of a sphere
   HalfspaceIntersection -- compute the intersection points of input halfspaces

Plotting Helpers
================

.. autosummary::
   :toctree: generated/

   delaunay_plot_2d     -- plot 2-D triangulation
   convex_hull_plot_2d  -- plot 2-D convex hull
   voronoi_plot_2d      -- plot 2-D voronoi diagram

.. seealso:: :ref:`Tutorial <qhulltutorial>`


Simplex representation
======================
The simplices (triangles, tetrahedra, ...) appearing in the Delaunay
tessellation (N-dim simplices), convex hull facets, and Voronoi ridges
(N-1 dim simplices) are represented in the following scheme::

    tess = Delaunay(points)
    hull = ConvexHull(points)
    voro = Voronoi(points)

    # coordinates of the j-th vertex of the i-th simplex
    tess.points[tess.simplices[i, j], :]        # tessellation element
    hull.points[hull.simplices[i, j], :]        # convex hull facet
    voro.vertices[voro.ridge_vertices[i, j], :] # ridge between Voronoi cells

For Delaunay triangulations and convex hulls, the neighborhood
structure of the simplices satisfies the condition:

    ``tess.neighbors[i,j]`` is the neighboring simplex of the i-th
    simplex, opposite to the j-vertex. It is -1 in case of no
    neighbor.

Convex hull facets also define a hyperplane equation::

    (hull.equations[i,:-1] * coord).sum() + hull.equations[i,-1] == 0

Similar hyperplane equations for the Delaunay triangulation correspond
to the convex hull facets on the corresponding N+1 dimensional
paraboloid.

The Delaunay triangulation objects offer a method for locating the
simplex containing a given point, and barycentric coordinate
computations.

Functions
---------

.. autosummary::
   :toctree: generated/

   tsearch
   distance_matrix
   minkowski_distance
   minkowski_distance_p
   procrustes

"""

from __future__ import division, print_function, absolute_import

from .kdtree import *
from .ckdtree import *
from .qhull import *
from ._spherical_voronoi import SphericalVoronoi
from ._plotutils import *
from ._procrustes import procrustes

__all__ = [s for s in dir() if not s.startswith('_')]
__all__ += ['distance', 'transform']

from . import distance, transform

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
