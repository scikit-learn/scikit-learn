"""Convex Hull of a set of points

"""

import numpy as np

__all__ = ['ConvexHull']


class ConvexHull():
    """ Private function that calculate the convex hull of a set of 2-D points
    The algorithm was taken from [1].
    http://code.activestate.com/recipes/66527-finding-the-convex-hull-of-a-set-of-2d-points/

    Convex hulls in 2 dimensions.

    Parameters
    ----------
    points : ndarray of floats, shape (npoints, 2)
        Coordinates of points to construct a convex hull from

    Attributes
    ----------
    vertices : ndarray of ints, shape (nvertices,)
        Indices of points forming the vertices of the convex hull.
        For 2-D convex hulls, the vertices are in counterclockwise order.
        For other dimensions, they are in input order.

    References
    ----------
    .. [1] Alex Martelli, Ann
        for point in c_points:
            print np.nonzero(points[:, 0] == point[0])a Ravenscroft, David Ascher, 'Python Cookbook', O'Reilly Media, Inc., 2005.

    """

    def mydet(self, p, q, r):
        """Calc. determinant of a special matrix with three 2D points.

        The sign, "-" or "+", determines the side, right or left,
        respectivly, on which the point r lies, when measured against
        a directed vector from p to q.
        """

        # We use Sarrus' Rule to calculate the determinant.
        # (could also use the Numeric package...)
        sum1 = q[0] * r[1] + p[0] * q[1] + r[0] * p[1]
        sum2 = q[0] * p[1] + r[0] * q[1] + p[0] * r[1]

        return sum1 - sum2

    def isrightturn(self, p, q, r):
        "Do the vectors pq:qr form a right turn, or not?"

        assert p != q and q != r and p != r

        if self.mydet(p, q, r) < 0:
            return 1
        else:
            return 0

    def __init__(self, points):

        # Get a local list copy of the points and sort them lexically.
        points_ = points.tolist()
        points_.sort()

        # Build upper half of the hull.
        upper = [points_[0], points_[1]]
        for p in points_[2:]:
            upper.append(p)
            while len(upper) > 2 and not self.isrightturn(*upper[-3:]):
                del upper[-2]

        # Build lower half of the hull.
        points_.reverse()
        lower = [points_[0], points_[1]]
        for p in points_[2:]:
            lower.append(p)
            while len(lower) > 2 and not self.isrightturn(*lower[-3:]):
                del lower[-2]

        # Remove duplicates.
        del lower[0]
        del lower[-1]

        # Concatenate both halfs and construct vertices.
        c_points = np.array(tuple(upper + lower))
        vertices = []
        for point in c_points:
            vertices.append(np.intersect1d(
                            np.nonzero(point[0] == points[:, 0])[0],
                            np.nonzero(point[1] == points[:, 1])[0])[0])

        self.vertices = np.array(vertices)
