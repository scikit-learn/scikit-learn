from __future__ import division
import math
import numpy as np
from scipy import spatial

from .._shared.utils import get_bound_method_class, safe_as_int


def _to_ndimage_mode(mode):
    """Convert from `numpy.pad` mode name to the corresponding ndimage mode."""
    mode_translation_dict = dict(edge='nearest', symmetric='reflect',
                                 reflect='mirror')
    if mode in mode_translation_dict:
        mode = mode_translation_dict[mode]
    return mode


def _center_and_normalize_points(points):
    """Center and normalize image points.

    The points are transformed in a two-step procedure that is expressed
    as a transformation matrix. The matrix of the resulting points is usually
    better conditioned than the matrix of the original points.

    Center the image points, such that the new coordinate system has its
    origin at the centroid of the image points.

    Normalize the image points, such that the mean distance from the points
    to the origin of the coordinate system is sqrt(2).

    Parameters
    ----------
    points : (N, 2) array
        The coordinates of the image points.

    Returns
    -------
    matrix : (3, 3) array
        The transformation matrix to obtain the new points.
    new_points : (N, 2) array
        The transformed image points.

    References
    ----------
    .. [1] Hartley, Richard I. "In defense of the eight-point algorithm."
           Pattern Analysis and Machine Intelligence, IEEE Transactions on 19.6
           (1997): 580-593.

    """

    centroid = np.mean(points, axis=0)

    rms = math.sqrt(np.sum((points - centroid) ** 2) / points.shape[0])

    norm_factor = math.sqrt(2) / rms

    matrix = np.array([[norm_factor, 0, -norm_factor * centroid[0]],
                       [0, norm_factor, -norm_factor * centroid[1]],
                       [0, 0, 1]])

    pointsh = np.row_stack([points.T, np.ones((points.shape[0]),)])

    new_pointsh = np.dot(matrix, pointsh).T

    new_points = new_pointsh[:, :2]
    new_points[:, 0] /= new_pointsh[:, 2]
    new_points[:, 1] /= new_pointsh[:, 2]

    return matrix, new_points


def _umeyama(src, dst, estimate_scale):
    """Estimate N-D similarity transformation with or without scaling.

    Parameters
    ----------
    src : (M, N) array
        Source coordinates.
    dst : (M, N) array
        Destination coordinates.
    estimate_scale : bool
        Whether to estimate scaling factor.

    Returns
    -------
    T : (N + 1, N + 1)
        The homogeneous similarity transformation matrix. The matrix contains
        NaN values only if the problem is not well-conditioned.

    References
    ----------
    .. [1] "Least-squares estimation of transformation parameters between two
            point patterns", Shinji Umeyama, PAMI 1991, DOI: 10.1109/34.88573

    """

    num = src.shape[0]
    dim = src.shape[1]

    # Compute mean of src and dst.
    src_mean = src.mean(axis=0)
    dst_mean = dst.mean(axis=0)

    # Subtract mean from src and dst.
    src_demean = src - src_mean
    dst_demean = dst - dst_mean

    # Eq. (38).
    A = np.dot(dst_demean.T, src_demean) / num

    # Eq. (39).
    d = np.ones((dim,), dtype=np.double)
    if np.linalg.det(A) < 0:
        d[dim - 1] = -1

    T = np.eye(dim + 1, dtype=np.double)

    U, S, V = np.linalg.svd(A)

    # Eq. (40) and (43).
    rank = np.linalg.matrix_rank(A)
    if rank == 0:
        return np.nan * T
    elif rank == dim - 1:
        if np.linalg.det(U) * np.linalg.det(V) > 0:
            T[:dim, :dim] = np.dot(U, V)
        else:
            s = d[dim - 1]
            d[dim - 1] = -1
            T[:dim, :dim] = np.dot(U, np.dot(np.diag(d), V))
            d[dim - 1] = s
    else:
        T[:dim, :dim] = np.dot(U, np.dot(np.diag(d), V))

    if estimate_scale:
        # Eq. (41) and (42).
        scale = 1.0 / src_demean.var(axis=0).sum() * np.dot(S, d)
    else:
        scale = 1.0

    T[:dim, dim] = dst_mean - scale * np.dot(T[:dim, :dim], src_mean.T)
    T[:dim, :dim] *= scale

    return T


class GeometricTransform(object):
    """Base class for geometric transformations.

    """
    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Source coordinates.

        Returns
        -------
        coords : (N, 2) array
            Destination coordinates.

        """
        raise NotImplementedError()

    def inverse(self, coords):
        """Apply inverse transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Destination coordinates.

        Returns
        -------
        coords : (N, 2) array
            Source coordinates.

        """
        raise NotImplementedError()

    def residuals(self, src, dst):
        """Determine residuals of transformed destination coordinates.

        For each transformed source coordinate the euclidean distance to the
        respective destination coordinate is determined.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        residuals : (N, ) array
            Residual for coordinate.

        """
        return np.sqrt(np.sum((self(src) - dst)**2, axis=1))

    def __add__(self, other):
        """Combine this transformation with another.

        """
        raise NotImplementedError()


class FundamentalMatrixTransform(GeometricTransform):
    """Fundamental matrix transformation.

    The fundamental matrix relates corresponding points between a pair of
    uncalibrated images. The matrix transforms homogeneous image points in one
    image to epipolar lines in the other image.

    The fundamental matrix is only defined for a pair of moving images. In the
    case of pure rotation or planar scenes, the homography describes the
    geometric relation between two images (`ProjectiveTransform`). If the
    intrinsic calibration of the images is known, the essential matrix describes
    the metric relation between the two images (`EssentialMatrixTransform`).

    References
    ----------
    .. [1] Hartley, Richard, and Andrew Zisserman. Multiple view geometry in
           computer vision. Cambridge university press, 2003.

    Parameters
    ----------
    matrix : (3, 3) array, optional
        Fundamental matrix.

    Attributes
    ----------
    params : (3, 3) array
        Fundamental matrix.

    """

    def __init__(self, matrix=None):
        if matrix is None:
            # default to an identity transform
            matrix = np.eye(3)
        if matrix.shape != (3, 3):
            raise ValueError("Invalid shape of transformation matrix")
        self.params = matrix

    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Source coordinates.

        Returns
        -------
        coords : (N, 3) array
            Epipolar lines in the destination image.

        """
        coords_homogeneous = np.column_stack([coords, np.ones(coords.shape[0])])
        return np.dot(coords_homogeneous, self.params.T)

    def inverse(self, coords):
        """Apply inverse transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Destination coordinates.

        Returns
        -------
        coords : (N, 3) array
            Epipolar lines in the source image.

        """
        coords_homogeneous = np.column_stack([coords, np.ones(coords.shape[0])])
        return np.dot(coords_homogeneous, self.params)

    def _setup_constraint_matrix(self, src, dst):
        """Setup and solve the homogeneous epipolar constraint matrix::

            dst' * F * src = 0.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        F_normalized : (3, 3) array
            The normalized solution to the homogeneous system. If the system
            is not well-conditioned, this matrix contains NaNs.
        src_matrix : (3, 3) array
            The transformation matrix to obtain the normalized source
            coordinates.
        dst_matrix : (3, 3) array
            The transformation matrix to obtain the normalized destination
            coordinates.

        """
        assert src.shape == dst.shape
        assert src.shape[0] >= 8

        # Center and normalize image points for better numerical stability.
        try:
            src_matrix, src = _center_and_normalize_points(src)
            dst_matrix, dst = _center_and_normalize_points(dst)
        except ZeroDivisionError:
            self.params = np.full((3, 3), np.nan)
            return 3 * [np.full((3, 3), np.nan)]

        # Setup homogeneous linear equation as dst' * F * src = 0.
        A = np.ones((src.shape[0], 9))
        A[:, :2] = src
        A[:, :3] *= dst[:, 0, np.newaxis]
        A[:, 3:5] = src
        A[:, 3:6] *= dst[:, 1, np.newaxis]
        A[:, 6:8] = src

        # Solve for the nullspace of the constraint matrix.
        _, _, V = np.linalg.svd(A)
        F_normalized = V[-1, :].reshape(3, 3)

        return F_normalized, src_matrix, dst_matrix

    def estimate(self, src, dst):
        """Estimate fundamental matrix using 8-point algorithm.

        The 8-point algorithm requires at least 8 corresponding point pairs for
        a well-conditioned solution, otherwise the over-determined solution is
        estimated.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        F_normalized, src_matrix, dst_matrix = \
            self._setup_constraint_matrix(src, dst)

        # Enforcing the internal constraint that two singular values must be
        # non-zero and one must be zero.
        U, S, V = np.linalg.svd(F_normalized)
        S[2] = 0
        F = np.dot(U, np.dot(np.diag(S), V))

        self.params = np.dot(dst_matrix.T, np.dot(F, src_matrix))

        return True

    def residuals(self, src, dst):
        """Compute the Sampson distance.

        The Sampson distance is the first approximation to the geometric error.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        residuals : (N, ) array
            Sampson distance.

        """
        src_homogeneous = np.column_stack([src, np.ones(src.shape[0])])
        dst_homogeneous = np.column_stack([dst, np.ones(dst.shape[0])])

        F_src = np.dot(self.params, src_homogeneous.T)
        Ft_dst = np.dot(self.params.T, dst_homogeneous.T)

        dst_F_src = np.sum(dst_homogeneous * F_src.T, axis=1)

        return np.abs(dst_F_src) / np.sqrt(F_src[0] ** 2 + F_src[1] ** 2
                                           + Ft_dst[0] ** 2 + Ft_dst[1] ** 2)


class EssentialMatrixTransform(FundamentalMatrixTransform):
    """Essential matrix transformation.

    The essential matrix relates corresponding points between a pair of
    calibrated images. The matrix transforms normalized, homogeneous image
    points in one image to epipolar lines in the other image.

    The essential matrix is only defined for a pair of moving images capturing a
    non-planar scene. In the case of pure rotation or planar scenes, the
    homography describes the geometric relation between two images
    (`ProjectiveTransform`). If the intrinsic calibration of the images is
    unknown, the fundamental matrix describes the projective relation between
    the two images (`FundamentalMatrixTransform`).

    References
    ----------
    .. [1] Hartley, Richard, and Andrew Zisserman. Multiple view geometry in
           computer vision. Cambridge university press, 2003.

    Parameters
    ----------
    rotation : (3, 3) array, optional
        Rotation matrix of the relative camera motion.
    translation : (3, 1) array, optional
        Translation vector of the relative camera motion. The vector must
        have unit length.
    matrix : (3, 3) array, optional
        Essential matrix.

    Attributes
    ----------
    params : (3, 3) array
        Essential matrix.

    """

    def __init__(self, rotation=None, translation=None, matrix=None):
        if rotation is not None:
            if translation is None:
                raise ValueError("Both rotation and translation required")
            if rotation.shape != (3, 3):
                raise ValueError("Invalid shape of rotation matrix")
            if abs(np.linalg.det(rotation) - 1) > 1e-6:
                raise ValueError("Rotation matrix must have unit determinant")
            if translation.size != 3:
                raise ValueError("Invalid shape of translation vector")
            if abs(np.linalg.norm(translation) - 1) > 1e-6:
                raise ValueError("Translation vector must have unit length")
            # Matrix representation of the cross product for t.
            t_x = np.array([0, -translation[2], translation[1],
                            translation[2], 0, -translation[0],
                            -translation[1], translation[0], 0]).reshape(3, 3)
            self.params = np.dot(t_x, rotation)
        elif matrix is not None:
            if matrix.shape != (3, 3):
                raise ValueError("Invalid shape of transformation matrix")
            self.params = matrix
        else:
            # default to an identity transform
            self.params = np.eye(3)

    def estimate(self, src, dst):
        """Estimate essential matrix using 8-point algorithm.

        The 8-point algorithm requires at least 8 corresponding point pairs for
        a well-conditioned solution, otherwise the over-determined solution is
        estimated.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        E_normalized, src_matrix, dst_matrix = \
            self._setup_constraint_matrix(src, dst)

        # Enforcing the internal constraint that two singular values must be
        # equal and one must be zero.
        U, S, V = np.linalg.svd(E_normalized)
        S[0] = (S[0] + S[1]) / 2.0
        S[1] = S[0]
        S[2] = 0
        E = np.dot(U, np.dot(np.diag(S), V))

        self.params = np.dot(dst_matrix.T, np.dot(E, src_matrix))

        return True


class ProjectiveTransform(GeometricTransform):
    """Projective transformation.

    Apply a projective transformation (homography) on coordinates.

    For each homogeneous coordinate :math:`\mathbf{x} = [x, y, 1]^T`, its
    target position is calculated by multiplying with the given matrix,
    :math:`H`, to give :math:`H \mathbf{x}`::

      [[a0 a1 a2]
       [b0 b1 b2]
       [c0 c1 1 ]].

    E.g., to rotate by theta degrees clockwise, the matrix should be::

      [[cos(theta) -sin(theta) 0]
       [sin(theta)  cos(theta) 0]
       [0            0         1]]

    or, to translate x by 10 and y by 20::

      [[1 0 10]
       [0 1 20]
       [0 0 1 ]].

    Parameters
    ----------
    matrix : (3, 3) array, optional
        Homogeneous transformation matrix.

    Attributes
    ----------
    params : (3, 3) array
        Homogeneous transformation matrix.

    """

    _coeffs = range(8)

    def __init__(self, matrix=None):
        if matrix is None:
            # default to an identity transform
            matrix = np.eye(3)
        if matrix.shape != (3, 3):
            raise ValueError("invalid shape of transformation matrix")
        self.params = matrix

    @property
    def _inv_matrix(self):
        return np.linalg.inv(self.params)

    def _apply_mat(self, coords, matrix):
        coords = np.array(coords, copy=False, ndmin=2)

        x, y = np.transpose(coords)
        src = np.vstack((x, y, np.ones_like(x)))
        dst = np.dot(src.transpose(), matrix.transpose())

        # rescale to homogeneous coordinates
        dst[:, 0] /= dst[:, 2]
        dst[:, 1] /= dst[:, 2]

        return dst[:, :2]

    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Source coordinates.

        Returns
        -------
        coords : (N, 2) array
            Destination coordinates.

        """
        return self._apply_mat(coords, self.params)

    def inverse(self, coords):
        """Apply inverse transformation.

        Parameters
        ----------
        coords : (N, 2) array
            Destination coordinates.

        Returns
        -------
        coords : (N, 2) array
            Source coordinates.

        """
        return self._apply_mat(coords, self._inv_matrix)

    def estimate(self, src, dst):
        """Estimate the transformation from a set of corresponding points.

        You can determine the over-, well- and under-determined parameters
        with the total least-squares method.

        Number of source and destination coordinates must match.

        The transformation is defined as::

            X = (a0*x + a1*y + a2) / (c0*x + c1*y + 1)
            Y = (b0*x + b1*y + b2) / (c0*x + c1*y + 1)

        These equations can be transformed to the following form::

            0 = a0*x + a1*y + a2 - c0*x*X - c1*y*X - X
            0 = b0*x + b1*y + b2 - c0*x*Y - c1*y*Y - Y

        which exist for each set of corresponding points, so we have a set of
        N * 2 equations. The coefficients appear linearly so we can write
        A x = 0, where::

            A   = [[x y 1 0 0 0 -x*X -y*X -X]
                   [0 0 0 x y 1 -x*Y -y*Y -Y]
                    ...
                    ...
                  ]
            x.T = [a0 a1 a2 b0 b1 b2 c0 c1 c3]

        In case of total least-squares the solution of this homogeneous system
        of equations is the right singular vector of A which corresponds to the
        smallest singular value normed by the coefficient c3.

        In case of the affine transformation the coefficients c0 and c1 are 0.
        Thus the system of equations is::

            A   = [[x y 1 0 0 0 -X]
                   [0 0 0 x y 1 -Y]
                    ...
                    ...
                  ]
            x.T = [a0 a1 a2 b0 b1 b2 c3]

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        try:
            src_matrix, src = _center_and_normalize_points(src)
            dst_matrix, dst = _center_and_normalize_points(dst)
        except ZeroDivisionError:
            self.params = np.nan * np.empty((3, 3))
            return False

        xs = src[:, 0]
        ys = src[:, 1]
        xd = dst[:, 0]
        yd = dst[:, 1]
        rows = src.shape[0]

        # params: a0, a1, a2, b0, b1, b2, c0, c1
        A = np.zeros((rows * 2, 9))
        A[:rows, 0] = xs
        A[:rows, 1] = ys
        A[:rows, 2] = 1
        A[:rows, 6] = - xd * xs
        A[:rows, 7] = - xd * ys
        A[rows:, 3] = xs
        A[rows:, 4] = ys
        A[rows:, 5] = 1
        A[rows:, 6] = - yd * xs
        A[rows:, 7] = - yd * ys
        A[:rows, 8] = xd
        A[rows:, 8] = yd

        # Select relevant columns, depending on params
        A = A[:, list(self._coeffs) + [8]]

        _, _, V = np.linalg.svd(A)

        H = np.zeros((3, 3))
        # solution is right singular vector that corresponds to smallest
        # singular value
        H.flat[list(self._coeffs) + [8]] = - V[-1, :-1] / V[-1, -1]
        H[2, 2] = 1

        # De-center and de-normalize
        H = np.dot(np.linalg.inv(dst_matrix), np.dot(H, src_matrix))

        self.params = H

        return True

    def __add__(self, other):
        """Combine this transformation with another.

        """
        if isinstance(other, ProjectiveTransform):
            # combination of the same types result in a transformation of this
            # type again, otherwise use general projective transformation
            if type(self) == type(other):
                tform = self.__class__
            else:
                tform = ProjectiveTransform
            return tform(other.params.dot(self.params))
        elif (hasattr(other, '__name__')
                and other.__name__ == 'inverse'
                and hasattr(get_bound_method_class(other), '_inv_matrix')):
            return ProjectiveTransform(other.__self__._inv_matrix.dot(self.params))
        else:
            raise TypeError("Cannot combine transformations of differing "
                            "types.")


class AffineTransform(ProjectiveTransform):
    """2D affine transformation.

    Has the following form::

        X = a0*x + a1*y + a2 =
          = sx*x*cos(rotation) - sy*y*sin(rotation + shear) + a2

        Y = b0*x + b1*y + b2 =
          = sx*x*sin(rotation) + sy*y*cos(rotation + shear) + b2

    where ``sx`` and ``sy`` are scale factors in the x and y directions,
    and the homogeneous transformation matrix is::

        [[a0  a1  a2]
         [b0  b1  b2]
         [0   0    1]]

    Parameters
    ----------
    matrix : (3, 3) array, optional
        Homogeneous transformation matrix.
    scale : (sx, sy) as array, list or tuple, optional
        Scale factors.
    rotation : float, optional
        Rotation angle in counter-clockwise direction as radians.
    shear : float, optional
        Shear angle in counter-clockwise direction as radians.
    translation : (tx, ty) as array, list or tuple, optional
        Translation parameters.

    Attributes
    ----------
    params : (3, 3) array
        Homogeneous transformation matrix.

    """

    _coeffs = range(6)

    def __init__(self, matrix=None, scale=None, rotation=None, shear=None,
                 translation=None):
        params = any(param is not None
                     for param in (scale, rotation, shear, translation))

        if params and matrix is not None:
            raise ValueError("You cannot specify the transformation matrix and"
                             " the implicit parameters at the same time.")
        elif matrix is not None:
            if matrix.shape != (3, 3):
                raise ValueError("Invalid shape of transformation matrix.")
            self.params = matrix
        elif params:
            if scale is None:
                scale = (1, 1)
            if rotation is None:
                rotation = 0
            if shear is None:
                shear = 0
            if translation is None:
                translation = (0, 0)

            sx, sy = scale
            self.params = np.array([
                [sx * math.cos(rotation), -sy * math.sin(rotation + shear), 0],
                [sx * math.sin(rotation),  sy * math.cos(rotation + shear), 0],
                [                      0,                                0, 1]
            ])
            self.params[0:2, 2] = translation
        else:
            # default to an identity transform
            self.params = np.eye(3)

    @property
    def scale(self):
        sx = math.sqrt(self.params[0, 0] ** 2 + self.params[1, 0] ** 2)
        sy = math.sqrt(self.params[0, 1] ** 2 + self.params[1, 1] ** 2)
        return sx, sy

    @property
    def rotation(self):
        return math.atan2(self.params[1, 0], self.params[0, 0])

    @property
    def shear(self):
        beta = math.atan2(- self.params[0, 1], self.params[1, 1])
        return beta - self.rotation

    @property
    def translation(self):
        return self.params[0:2, 2]


class PiecewiseAffineTransform(GeometricTransform):
    """2D piecewise affine transformation.

    Control points are used to define the mapping. The transform is based on
    a Delaunay triangulation of the points to form a mesh. Each triangle is
    used to find a local affine transform.

    Attributes
    ----------
    affines : list of AffineTransform objects
        Affine transformations for each triangle in the mesh.
    inverse_affines : list of AffineTransform objects
        Inverse affine transformations for each triangle in the mesh.

    """

    def __init__(self):
        self._tesselation = None
        self._inverse_tesselation = None
        self.affines = None
        self.inverse_affines = None

    def estimate(self, src, dst):
        """Estimate the transformation from a set of corresponding points.

        Number of source and destination coordinates must match.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        # forward piecewise affine
        # triangulate input positions into mesh
        self._tesselation = spatial.Delaunay(src)
        # find affine mapping from source positions to destination
        self.affines = []
        for tri in self._tesselation.vertices:
            affine = AffineTransform()
            affine.estimate(src[tri, :], dst[tri, :])
            self.affines.append(affine)

        # inverse piecewise affine
        # triangulate input positions into mesh
        self._inverse_tesselation = spatial.Delaunay(dst)
        # find affine mapping from source positions to destination
        self.inverse_affines = []
        for tri in self._inverse_tesselation.vertices:
            affine = AffineTransform()
            affine.estimate(dst[tri, :], src[tri, :])
            self.inverse_affines.append(affine)

        return True

    def __call__(self, coords):
        """Apply forward transformation.

        Coordinates outside of the mesh will be set to `- 1`.

        Parameters
        ----------
        coords : (N, 2) array
            Source coordinates.

        Returns
        -------
        coords : (N, 2) array
            Transformed coordinates.

        """

        out = np.empty_like(coords, np.double)

        # determine triangle index for each coordinate
        simplex = self._tesselation.find_simplex(coords)

        # coordinates outside of mesh
        out[simplex == -1, :] = -1

        for index in range(len(self._tesselation.vertices)):
            # affine transform for triangle
            affine = self.affines[index]
            # all coordinates within triangle
            index_mask = simplex == index

            out[index_mask, :] = affine(coords[index_mask, :])

        return out

    def inverse(self, coords):
        """Apply inverse transformation.

        Coordinates outside of the mesh will be set to `- 1`.

        Parameters
        ----------
        coords : (N, 2) array
            Source coordinates.

        Returns
        -------
        coords : (N, 2) array
            Transformed coordinates.

        """

        out = np.empty_like(coords, np.double)

        # determine triangle index for each coordinate
        simplex = self._inverse_tesselation.find_simplex(coords)

        # coordinates outside of mesh
        out[simplex == -1, :] = -1

        for index in range(len(self._inverse_tesselation.vertices)):
            # affine transform for triangle
            affine = self.inverse_affines[index]
            # all coordinates within triangle
            index_mask = simplex == index

            out[index_mask, :] = affine(coords[index_mask, :])

        return out


class EuclideanTransform(ProjectiveTransform):
    """2D Euclidean transformation.

    Has the following form::

        X = a0 * x - b0 * y + a1 =
          = x * cos(rotation) - y * sin(rotation) + a1

        Y = b0 * x + a0 * y + b1 =
          = x * sin(rotation) + y * cos(rotation) + b1

    where the homogeneous transformation matrix is::

        [[a0  b0  a1]
         [b0  a0  b1]
         [0   0    1]]

    The Euclidean transformation is a rigid transformation with rotation and
    translation parameters. The similarity transformation extends the Euclidean
    transformation with a single scaling factor.

    Parameters
    ----------
    matrix : (3, 3) array, optional
        Homogeneous transformation matrix.
    rotation : float, optional
        Rotation angle in counter-clockwise direction as radians.
    translation : (tx, ty) as array, list or tuple, optional
        x, y translation parameters.

    Attributes
    ----------
    params : (3, 3) array
        Homogeneous transformation matrix.

    """

    def __init__(self, matrix=None, rotation=None, translation=None):
        params = any(param is not None
                     for param in (rotation, translation))

        if params and matrix is not None:
            raise ValueError("You cannot specify the transformation matrix and"
                             " the implicit parameters at the same time.")
        elif matrix is not None:
            if matrix.shape != (3, 3):
                raise ValueError("Invalid shape of transformation matrix.")
            self.params = matrix
        elif params:
            if rotation is None:
                rotation = 0
            if translation is None:
                translation = (0, 0)

            self.params = np.array([
                [math.cos(rotation), - math.sin(rotation), 0],
                [math.sin(rotation),   math.cos(rotation), 0],
                [                 0,                    0, 1]
            ])
            self.params[0:2, 2] = translation
        else:
            # default to an identity transform
            self.params = np.eye(3)

    def estimate(self, src, dst):
        """Estimate the transformation from a set of corresponding points.

        You can determine the over-, well- and under-determined parameters
        with the total least-squares method.

        Number of source and destination coordinates must match.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        self.params = _umeyama(src, dst, False)

        return True

    @property
    def rotation(self):
        return math.atan2(self.params[1, 0], self.params[1, 1])

    @property
    def translation(self):
        return self.params[0:2, 2]


class SimilarityTransform(EuclideanTransform):
    """2D similarity transformation.

    Has the following form::

        X = a0 * x - b0 * y + a1 =
          = s * x * cos(rotation) - s * y * sin(rotation) + a1

        Y = b0 * x + a0 * y + b1 =
          = s * x * sin(rotation) + s * y * cos(rotation) + b1

    where ``s`` is a scale factor and the homogeneous transformation matrix is::

        [[a0  b0  a1]
         [b0  a0  b1]
         [0   0    1]]

    The similarity transformation extends the Euclidean transformation with a
    single scaling factor in addition to the rotation and translation
    parameters.

    Parameters
    ----------
    matrix : (3, 3) array, optional
        Homogeneous transformation matrix.
    scale : float, optional
        Scale factor.
    rotation : float, optional
        Rotation angle in counter-clockwise direction as radians.
    translation : (tx, ty) as array, list or tuple, optional
        x, y translation parameters.

    Attributes
    ----------
    params : (3, 3) array
        Homogeneous transformation matrix.

    """

    def __init__(self, matrix=None, scale=None, rotation=None,
                 translation=None):
        params = any(param is not None
                     for param in (scale, rotation, translation))

        if params and matrix is not None:
            raise ValueError("You cannot specify the transformation matrix and"
                             " the implicit parameters at the same time.")
        elif matrix is not None:
            if matrix.shape != (3, 3):
                raise ValueError("Invalid shape of transformation matrix.")
            self.params = matrix
        elif params:
            if scale is None:
                scale = 1
            if rotation is None:
                rotation = 0
            if translation is None:
                translation = (0, 0)

            self.params = np.array([
                [math.cos(rotation), - math.sin(rotation), 0],
                [math.sin(rotation),   math.cos(rotation), 0],
                [                 0,                    0, 1]
            ])
            self.params[0:2, 0:2] *= scale
            self.params[0:2, 2] = translation
        else:
            # default to an identity transform
            self.params = np.eye(3)

    def estimate(self, src, dst):
        """Estimate the transformation from a set of corresponding points.

        You can determine the over-, well- and under-determined parameters
        with the total least-squares method.

        Number of source and destination coordinates must match.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """

        self.params = _umeyama(src, dst, True)

        return True

    @property
    def scale(self):
        if abs(math.cos(self.rotation)) < np.spacing(1):
            # sin(self.rotation) == 1
            scale = self.params[1, 0]
        else:
            scale = self.params[0, 0] / math.cos(self.rotation)
        return scale


class PolynomialTransform(GeometricTransform):
    """2D polynomial transformation.

    Has the following form::

        X = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i ))
        Y = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i ))

    Parameters
    ----------
    params : (2, N) array, optional
        Polynomial coefficients where `N * 2 = (order + 1) * (order + 2)`. So,
        a_ji is defined in `params[0, :]` and b_ji in `params[1, :]`.

    Attributes
    ----------
    params : (2, N) array
        Polynomial coefficients where `N * 2 = (order + 1) * (order + 2)`. So,
        a_ji is defined in `params[0, :]` and b_ji in `params[1, :]`.

    """

    def __init__(self, params=None):
        if params is None:
            # default to transformation which preserves original coordinates
            params = np.array([[0, 1, 0], [0, 0, 1]])
        if params.shape[0] != 2:
            raise ValueError("invalid shape of transformation parameters")
        self.params = params

    def estimate(self, src, dst, order=2):
        """Estimate the transformation from a set of corresponding points.

        You can determine the over-, well- and under-determined parameters
        with the total least-squares method.

        Number of source and destination coordinates must match.

        The transformation is defined as::

            X = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i ))
            Y = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i ))

        These equations can be transformed to the following form::

            0 = sum[j=0:order]( sum[i=0:j]( a_ji * x**(j - i) * y**i )) - X
            0 = sum[j=0:order]( sum[i=0:j]( b_ji * x**(j - i) * y**i )) - Y

        which exist for each set of corresponding points, so we have a set of
        N * 2 equations. The coefficients appear linearly so we can write
        A x = 0, where::

            A   = [[1 x y x**2 x*y y**2 ... 0 ...             0 -X]
                   [0 ...                 0 1 x y x**2 x*y y**2 -Y]
                    ...
                    ...
                  ]
            x.T = [a00 a10 a11 a20 a21 a22 ... ann
                   b00 b10 b11 b20 b21 b22 ... bnn c3]

        In case of total least-squares the solution of this homogeneous system
        of equations is the right singular vector of A which corresponds to the
        smallest singular value normed by the coefficient c3.

        Parameters
        ----------
        src : (N, 2) array
            Source coordinates.
        dst : (N, 2) array
            Destination coordinates.
        order : int, optional
            Polynomial order (number of coefficients is order + 1).

        Returns
        -------
        success : bool
            True, if model estimation succeeds.

        """
        xs = src[:, 0]
        ys = src[:, 1]
        xd = dst[:, 0]
        yd = dst[:, 1]
        rows = src.shape[0]

        # number of unknown polynomial coefficients
        order = safe_as_int(order)
        u = (order + 1) * (order + 2)

        A = np.zeros((rows * 2, u + 1))
        pidx = 0
        for j in range(order + 1):
            for i in range(j + 1):
                A[:rows, pidx] = xs ** (j - i) * ys ** i
                A[rows:, pidx + u // 2] = xs ** (j - i) * ys ** i
                pidx += 1

        A[:rows, -1] = xd
        A[rows:, -1] = yd

        _, _, V = np.linalg.svd(A)

        # solution is right singular vector that corresponds to smallest
        # singular value
        params = - V[-1, :-1] / V[-1, -1]

        self.params = params.reshape((2, u // 2))

        return True

    def __call__(self, coords):
        """Apply forward transformation.

        Parameters
        ----------
        coords : (N, 2) array
            source coordinates

        Returns
        -------
        coords : (N, 2) array
            Transformed coordinates.

        """
        x = coords[:, 0]
        y = coords[:, 1]
        u = len(self.params.ravel())
        # number of coefficients -> u = (order + 1) * (order + 2)
        order = int((- 3 + math.sqrt(9 - 4 * (2 - u))) / 2)
        dst = np.zeros(coords.shape)

        pidx = 0
        for j in range(order + 1):
            for i in range(j + 1):
                dst[:, 0] += self.params[0, pidx] * x ** (j - i) * y ** i
                dst[:, 1] += self.params[1, pidx] * x ** (j - i) * y ** i
                pidx += 1

        return dst

    def inverse(self, coords):
        raise Exception(
            'There is no explicit way to do the inverse polynomial '
            'transformation. Instead, estimate the inverse transformation '
            'parameters by exchanging source and destination coordinates,'
            'then apply the forward transformation.')


TRANSFORMS = {
    'euclidean': EuclideanTransform,
    'similarity': SimilarityTransform,
    'affine': AffineTransform,
    'piecewise-affine': PiecewiseAffineTransform,
    'projective': ProjectiveTransform,
    'fundamental': FundamentalMatrixTransform,
    'essential': EssentialMatrixTransform,
    'polynomial': PolynomialTransform,
}


def estimate_transform(ttype, src, dst, **kwargs):
    """Estimate 2D geometric transformation parameters.

    You can determine the over-, well- and under-determined parameters
    with the total least-squares method.

    Number of source and destination coordinates must match.

    Parameters
    ----------
    ttype : {'euclidean', similarity', 'affine', 'piecewise-affine', \
             'projective', 'polynomial'}
        Type of transform.
    kwargs : array or int
        Function parameters (src, dst, n, angle)::

            NAME / TTYPE        FUNCTION PARAMETERS
            'euclidean'         `src, `dst`
            'similarity'        `src, `dst`
            'affine'            `src, `dst`
            'piecewise-affine'  `src, `dst`
            'projective'        `src, `dst`
            'polynomial'        `src, `dst`, `order` (polynomial order,
                                                      default order is 2)

        Also see examples below.

    Returns
    -------
    tform : :class:`GeometricTransform`
        Transform object containing the transformation parameters and providing
        access to forward and inverse transformation functions.

    Examples
    --------
    >>> import numpy as np
    >>> from skimage import transform as tf

    >>> # estimate transformation parameters
    >>> src = np.array([0, 0, 10, 10]).reshape((2, 2))
    >>> dst = np.array([12, 14, 1, -20]).reshape((2, 2))

    >>> tform = tf.estimate_transform('similarity', src, dst)

    >>> np.allclose(tform.inverse(tform(src)), src)
    True

    >>> # warp image using the estimated transformation
    >>> from skimage import data
    >>> image = data.camera()

    >>> warp(image, inverse_map=tform.inverse) # doctest: +SKIP

    >>> # create transformation with explicit parameters
    >>> tform2 = tf.SimilarityTransform(scale=1.1, rotation=1,
    ...     translation=(10, 20))

    >>> # unite transformations, applied in order from left to right
    >>> tform3 = tform + tform2
    >>> np.allclose(tform3(src), tform2(tform(src)))
    True

    """
    ttype = ttype.lower()
    if ttype not in TRANSFORMS:
        raise ValueError('the transformation type \'%s\' is not'
                         'implemented' % ttype)

    tform = TRANSFORMS[ttype]()
    tform.estimate(src, dst, **kwargs)

    return tform


def matrix_transform(coords, matrix):
    """Apply 2D matrix transform.

    Parameters
    ----------
    coords : (N, 2) array
        x, y coordinates to transform
    matrix : (3, 3) array
        Homogeneous transformation matrix.

    Returns
    -------
    coords : (N, 2) array
        Transformed coordinates.

    """
    return ProjectiveTransform(matrix)(coords)
