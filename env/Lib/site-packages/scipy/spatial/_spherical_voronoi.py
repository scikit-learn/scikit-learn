"""
Spherical Voronoi Code

.. versionadded:: 0.18.0

"""
#
# Copyright (C)  Tyler Reddy, Ross Hemsley, Edd Edmondson,
#                Nikolai Nowaczyk, Joe Pitt-Francis, 2015.
#
# Distributed under the same BSD license as SciPy.
#

import warnings
import numpy as np
import scipy
from . import _voronoi
from scipy.spatial import cKDTree

__all__ = ['SphericalVoronoi']


class SphericalVoronoi:
    """ Voronoi diagrams on the surface of a sphere.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, ndim)
        Coordinates of points from which to construct a spherical
        Voronoi diagram.
    radius : float, optional
        Radius of the sphere (Default: 1)
    center : ndarray of floats, shape (ndim,)
        Center of sphere (Default: origin)
    threshold : float
        Threshold for detecting duplicate points and
        mismatches between points and sphere parameters.
        (Default: 1e-06)

    Attributes
    ----------
    points : double array of shape (npoints, ndim)
        the points in `ndim` dimensions to generate the Voronoi diagram from
    radius : double
        radius of the sphere
    center : double array of shape (ndim,)
        center of the sphere
    vertices : double array of shape (nvertices, ndim)
        Voronoi vertices corresponding to points
    regions : list of list of integers of shape (npoints, _ )
        the n-th entry is a list consisting of the indices
        of the vertices belonging to the n-th point in points

    Raises
    ------
    ValueError
        If there are duplicates in `points`.
        If the provided `radius` is not consistent with `points`.

    Notes
    -----
    The spherical Voronoi diagram algorithm proceeds as follows. The Convex
    Hull of the input points (generators) is calculated, and is equivalent to
    their Delaunay triangulation on the surface of the sphere [Caroli]_.
    The Convex Hull neighbour information is then used to
    order the Voronoi region vertices around each generator. The latter
    approach is substantially less sensitive to floating point issues than
    angle-based methods of Voronoi region vertex sorting.

    Empirical assessment of spherical Voronoi algorithm performance suggests
    quadratic time complexity (loglinear is optimal, but algorithms are more
    challenging to implement).

    References
    ----------
    .. [Caroli] Caroli et al. Robust and Efficient Delaunay triangulations of
                points on or close to a sphere. Research Report RR-7004, 2009.

    See Also
    --------
    Voronoi : Conventional Voronoi diagrams in N dimensions.

    Examples
    --------
    Do some imports and take some points on a cube:

    >>> from matplotlib import colors
    >>> from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    >>> import matplotlib.pyplot as plt
    >>> from scipy.spatial import SphericalVoronoi
    >>> from mpl_toolkits.mplot3d import proj3d
    >>> # set input data
    >>> points = np.array([[0, 0, 1], [0, 0, -1], [1, 0, 0],
    ...                    [0, 1, 0], [0, -1, 0], [-1, 0, 0], ])

    Calculate the spherical Voronoi diagram:

    >>> radius = 1
    >>> center = np.array([0, 0, 0])
    >>> sv = SphericalVoronoi(points, radius, center)

    Generate plot:

    >>> # sort vertices (optional, helpful for plotting)
    >>> sv.sort_vertices_of_regions()
    >>> fig = plt.figure()
    >>> ax = fig.add_subplot(111, projection='3d')
    >>> # plot the unit sphere for reference (optional)
    >>> u = np.linspace(0, 2 * np.pi, 100)
    >>> v = np.linspace(0, np.pi, 100)
    >>> x = np.outer(np.cos(u), np.sin(v))
    >>> y = np.outer(np.sin(u), np.sin(v))
    >>> z = np.outer(np.ones(np.size(u)), np.cos(v))
    >>> ax.plot_surface(x, y, z, color='y', alpha=0.1)
    >>> # plot generator points
    >>> ax.scatter(points[:, 0], points[:, 1], points[:, 2], c='b')
    >>> # plot Voronoi vertices
    >>> ax.scatter(sv.vertices[:, 0], sv.vertices[:, 1], sv.vertices[:, 2],
    ...                    c='g')
    >>> # indicate Voronoi regions (as Euclidean polygons)
    >>> for region in sv.regions:
    ...    random_color = colors.rgb2hex(np.random.rand(3))
    ...    polygon = Poly3DCollection([sv.vertices[region]], alpha=1.0)
    ...    polygon.set_color(random_color)
    ...    ax.add_collection3d(polygon)
    >>> plt.show()

    """
    def __init__(self, points, radius=1, center=None, threshold=1e-06):

        if radius is None:
            radius = 1.
            warnings.warn('`radius` is `None`. '
                          'This will raise an error in a future version. '
                          'Please provide a floating point number '
                          '(i.e. `radius=1`).',
                          DeprecationWarning)

        self.points = points
        self.radius = radius
        self.threshold = threshold
        self._dim = len(points[0])
        if center is None:
            self.center = np.zeros(self._dim)
        else:
            self.center = np.array(center)

        # test degenerate input
        self._rank = np.linalg.matrix_rank(self.points - self.center,
                                           tol=self.threshold * self.radius)
        if self._rank <= 1:
            raise ValueError("Rank of input points must be at least 2")

        if cKDTree(self.points).query_pairs(self.threshold * self.radius):
            raise ValueError("Duplicate generators present.")

        radii = np.linalg.norm(self.points - self.center, axis=1)
        max_discrepancy = np.abs(radii - self.radius).max()
        if max_discrepancy >= self.threshold * self.radius:
            raise ValueError("Radius inconsistent with generators.")
        self.vertices = None
        self.regions = None
        self._tri = None
        self._calc_vertices_regions()

    def _handle_geodesic_input(self):

        # center the points
        centered = self.points - self.center

        # calculate an orthogonal transformation using SVD
        _, _, vh = np.linalg.svd(centered)

        # calculate the north and south poles in this basis
        poles = [[0, 0, self.radius], [0, 0, -self.radius]] @ vh

        # project points into inverse basis (such that z-components are zero)
        circle = centered @ vh.T[:, :2]

        # simplicial neighbors are adjacent on the circle
        angles = np.arctan2(circle[:, 1], circle[:, 0])
        indices = np.argsort(angles)

        # Voronoi vertices lie halfway between neighboring pairs
        vertices = centered[indices] + centered[np.roll(indices, 1)]
        vertices /= np.linalg.norm(vertices, axis=1)[:, np.newaxis]
        vertices *= self.radius

        # north and south poles are also Voronoi vertices
        vertices = np.concatenate((vertices, poles))

        # each region contains two vertices from the plane and the north and
        # south poles
        invf = np.argsort(indices)
        invb = np.argsort(np.roll(indices, 1))

        n = len(self.points)
        regions = np.vstack([invf,            # forward neighbor
                             [n] * n,         # north pole
                             invb,            # backward neighbor
                             [n + 1] * n]).T  # south pole

        self.regions = [list(region) for region in regions]
        self.vertices = vertices + self.center

    def _calc_vertices_regions(self):
        """
        Calculates the Voronoi vertices and regions of the generators stored
        in self.points. The vertices will be stored in self.vertices and the
        regions in self.regions.

        This algorithm was discussed at PyData London 2015 by
        Tyler Reddy, Ross Hemsley and Nikolai Nowaczyk
        """
        if self._dim == 3 and self._rank == 2:
            self._handle_geodesic_input()
            return

        # get Convex Hull
        self._tri = scipy.spatial.ConvexHull(self.points)
        # get circumcenters of Convex Hull triangles from facet equations
        # for 3D input circumcenters will have shape: (2N-4, 3)
        self.vertices = self.radius * self._tri.equations[:, :-1] + self.center
        # calculate regions from triangulation
        # for 3D input simplex_indices will have shape: (2N-4,)
        simplex_indices = np.arange(self._tri.simplices.shape[0])
        # for 3D input tri_indices will have shape: (6N-12,)
        tri_indices = np.column_stack([simplex_indices] * self._dim).ravel()
        # for 3D input point_indices will have shape: (6N-12,)
        point_indices = self._tri.simplices.ravel()
        # for 3D input indices will have shape: (6N-12,)
        indices = np.argsort(point_indices, kind='mergesort')
        # for 3D input flattened_groups will have shape: (6N-12,)
        flattened_groups = tri_indices[indices].astype(np.intp)
        # intervals will have shape: (N+1,)
        intervals = np.cumsum(np.bincount(point_indices + 1))

        # split flattened groups to get nested list of unsorted regions
        groups = [list(flattened_groups[intervals[i]:intervals[i + 1]])
                  for i in range(len(intervals) - 1)]
        self.regions = groups

    def sort_vertices_of_regions(self):
        """Sort indices of the vertices to be (counter-)clockwise ordered.

        Raises
        ------
        TypeError
            If the points are not three-dimensional.

        Notes
        -----
        For each region in regions, it sorts the indices of the Voronoi
        vertices such that the resulting points are in a clockwise or
        counterclockwise order around the generator point.

        This is done as follows: Recall that the n-th region in regions
        surrounds the n-th generator in points and that the k-th
        Voronoi vertex in vertices is the circumcenter of the k-th triangle
        in _tri.simplices.  For each region n, we choose the first triangle
        (=Voronoi vertex) in _tri.simplices and a vertex of that triangle
        not equal to the center n. These determine a unique neighbor of that
        triangle, which is then chosen as the second triangle. The second
        triangle will have a unique vertex not equal to the current vertex or
        the center. This determines a unique neighbor of the second triangle,
        which is then chosen as the third triangle and so forth. We proceed
        through all the triangles (=Voronoi vertices) belonging to the
        generator in points and obtain a sorted version of the vertices
        of its surrounding region.
        """
        if self._dim != 3:
            raise TypeError("Only supported for three-dimensional point sets")
        if self._rank == 2:
            return  # regions are sorted by construction
        _voronoi.sort_vertices_of_regions(self._tri.simplices, self.regions)
