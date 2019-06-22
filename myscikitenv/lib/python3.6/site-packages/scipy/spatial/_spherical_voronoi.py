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

import numpy as np
import scipy
import itertools
from . import _voronoi
from scipy.spatial.distance import pdist

__all__ = ['SphericalVoronoi']

def sphere_check(points, radius, center):
    """ Determines distance of generators from theoretical sphere
    surface.

    """
    actual_squared_radii = (((points[...,0] - center[0]) ** 2) +
                            ((points[...,1] - center[1]) ** 2) +
                            ((points[...,2] - center[2]) ** 2))
    max_discrepancy = (np.sqrt(actual_squared_radii) - radius).max()
    return abs(max_discrepancy)

def calc_circumcenters(tetrahedrons):
    """ Calculates the cirumcenters of the circumspheres of tetrahedrons.

    An implementation based on
    http://mathworld.wolfram.com/Circumsphere.html

    Parameters
    ----------
    tetrahedrons : an array of shape (N, 4, 3)
        consisting of N tetrahedrons defined by 4 points in 3D

    Returns
    ----------
    circumcenters : an array of shape (N, 3)
        consisting of the N circumcenters of the tetrahedrons in 3D

    """

    num = tetrahedrons.shape[0]
    a = np.concatenate((tetrahedrons, np.ones((num, 4, 1))), axis=2)

    sums = np.sum(tetrahedrons ** 2, axis=2)
    d = np.concatenate((sums[:, :, np.newaxis], a), axis=2)

    dx = np.delete(d, 1, axis=2)
    dy = np.delete(d, 2, axis=2)
    dz = np.delete(d, 3, axis=2)

    dx = np.linalg.det(dx)
    dy = -np.linalg.det(dy)
    dz = np.linalg.det(dz)
    a = np.linalg.det(a)

    nominator = np.vstack((dx, dy, dz))
    denominator = 2*a
    return (nominator / denominator).T


def project_to_sphere(points, center, radius):
    """
    Projects the elements of points onto the sphere defined
    by center and radius.

    Parameters
    ----------
    points : array of floats of shape (npoints, ndim)
             consisting of the points in a space of dimension ndim
    center : array of floats of shape (ndim,)
            the center of the sphere to project on
    radius : float
            the radius of the sphere to project on

    returns: array of floats of shape (npoints, ndim)
            the points projected onto the sphere
    """

    lengths = scipy.spatial.distance.cdist(points, np.array([center]))
    return (points - center) / lengths * radius + center


class SphericalVoronoi:
    """ Voronoi diagrams on the surface of a sphere.

    .. versionadded:: 0.18.0

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, 3)
        Coordinates of points from which to construct a spherical
        Voronoi diagram.
    radius : float, optional
        Radius of the sphere (Default: 1)
    center : ndarray of floats, shape (3,)
        Center of sphere (Default: origin)
    threshold : float
        Threshold for detecting duplicate points and
        mismatches between points and sphere parameters.
        (Default: 1e-06)

    Attributes
    ----------
    points : double array of shape (npoints, 3)
        the points in 3D to generate the Voronoi diagram from
    radius : double
        radius of the sphere
        Default: None (forces estimation, which is less precise)
    center : double array of shape (3,)
        center of the sphere
        Default: None (assumes sphere is centered at origin)
    vertices : double array of shape (nvertices, 3)
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
    A 3D Delaunay tetrahedralization is obtained by including the origin of
    the coordinate system as the fourth vertex of each simplex of the Convex
    Hull. The circumcenters of all tetrahedra in the system are calculated and
    projected to the surface of the sphere, producing the Voronoi vertices.
    The Delaunay tetrahedralization neighbour information is then used to
    order the Voronoi region vertices around each generator. The latter
    approach is substantially less sensitive to floating point issues than
    angle-based methods of Voronoi region vertex sorting.

    The surface area of spherical polygons is calculated by decomposing them
    into triangles and using L'Huilier's Theorem to calculate the spherical
    excess of each triangle [Weisstein]_. The sum of the spherical excesses is
    multiplied by the square of the sphere radius to obtain the surface area
    of the spherical polygon. For nearly-degenerate spherical polygons an area
    of approximately 0 is returned by default, rather than attempting the
    unstable calculation.

    Empirical assessment of spherical Voronoi algorithm performance suggests
    quadratic time complexity (loglinear is optimal, but algorithms are more
    challenging to implement). The reconstitution of the surface area of the
    sphere, measured as the sum of the surface areas of all Voronoi regions,
    is closest to 100 % for larger (>> 10) numbers of generators.

    References
    ----------
    .. [Caroli] Caroli et al. Robust and Efficient Delaunay triangulations of
                points on or close to a sphere. Research Report RR-7004, 2009.
    .. [Weisstein] "L'Huilier's Theorem." From MathWorld -- A Wolfram Web
                Resource. http://mathworld.wolfram.com/LHuiliersTheorem.html

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

    def __init__(self, points, radius=None, center=None, threshold=1e-06):
        self.points = points
        if np.any(center):
            self.center = center
        else:
            self.center = np.zeros(3)
        if radius:
            self.radius = radius
        else:
            self.radius = 1

        if pdist(self.points).min() <= threshold * self.radius:
            raise ValueError("Duplicate generators present.")

        max_discrepancy = sphere_check(self.points,
                                       self.radius,
                                       self.center)
        if max_discrepancy >= threshold * self.radius:
            raise ValueError("Radius inconsistent with generators.")
        self.vertices = None
        self.regions = None
        self._tri = None
        self._calc_vertices_regions()

    def _calc_vertices_regions(self):
        """
        Calculates the Voronoi vertices and regions of the generators stored
        in self.points. The vertices will be stored in self.vertices and the
        regions in self.regions.

        This algorithm was discussed at PyData London 2015 by
        Tyler Reddy, Ross Hemsley and Nikolai Nowaczyk
        """

        # perform 3D Delaunay triangulation on data set
        # (here ConvexHull can also be used, and is faster)
        self._tri = scipy.spatial.ConvexHull(self.points)

        # add the center to each of the simplices in tri to get the same
        # tetrahedrons we'd have gotten from Delaunay tetrahedralization
        # tetrahedrons will have shape: (2N-4, 4, 3)
        tetrahedrons = self._tri.points[self._tri.simplices]
        tetrahedrons = np.insert(
            tetrahedrons,
            3,
            np.array([self.center]),
            axis=1
        )

        # produce circumcenters of tetrahedrons from 3D Delaunay
        # circumcenters will have shape: (2N-4, 3)
        circumcenters = calc_circumcenters(tetrahedrons)

        # project tetrahedron circumcenters to the surface of the sphere
        # self.vertices will have shape: (2N-4, 3)
        self.vertices = project_to_sphere(
            circumcenters,
            self.center,
            self.radius
        )

        # calculate regions from triangulation
        # simplex_indices will have shape: (2N-4,)
        simplex_indices = np.arange(self._tri.simplices.shape[0])
        # tri_indices will have shape: (6N-12,)
        tri_indices = np.column_stack([simplex_indices, simplex_indices,
            simplex_indices]).ravel()
        # point_indices will have shape: (6N-12,)
        point_indices = self._tri.simplices.ravel()

        # array_associations will have shape: (6N-12, 2)
        array_associations = np.dstack((point_indices, tri_indices))[0]
        array_associations = array_associations[np.lexsort((
                                                array_associations[...,1],
                                                array_associations[...,0]))]
        array_associations = array_associations.astype(np.intp)

        # group by generator indices to produce
        # unsorted regions in nested list
        groups = [list(list(zip(*list(g)))[1])
                  for k, g in itertools.groupby(array_associations,
                                                lambda t: t[0])]

        self.regions = groups

    def sort_vertices_of_regions(self):
        """Sort indices of the vertices to be (counter-)clockwise ordered.

        Notes
        -----
        For each region in regions, it sorts the indices of the Voronoi
        vertices such that the resulting points are in a clockwise or
        counterclockwise order around the generator point.

        This is done as follows: Recall that the n-th region in regions
        surrounds the n-th generator in points and that the k-th
        Voronoi vertex in vertices is the projected circumcenter of the
        tetrahedron obtained by the k-th triangle in _tri.simplices (and the
        origin). For each region n, we choose the first triangle (=Voronoi
        vertex) in _tri.simplices and a vertex of that triangle not equal to
        the center n. These determine a unique neighbor of that triangle,
        which is then chosen as the second triangle. The second triangle
        will have a unique vertex not equal to the current vertex or the
        center. This determines a unique neighbor of the second triangle,
        which is then chosen as the third triangle and so forth. We proceed
        through all the triangles (=Voronoi vertices) belonging to the
        generator in points and obtain a sorted version of the vertices
        of its surrounding region.
        """

        _voronoi.sort_vertices_of_regions(self._tri.simplices, self.regions)
