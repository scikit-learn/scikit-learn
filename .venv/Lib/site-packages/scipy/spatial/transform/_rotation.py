from __future__ import annotations

from collections.abc import Iterable, Iterator
from types import EllipsisType, ModuleType, NotImplementedType

import numpy as np

import scipy.spatial.transform._rotation_cy as cython_backend
import scipy.spatial.transform._rotation_xp as xp_backend
from scipy.spatial.transform._rotation_groups import create_group
from scipy._lib._array_api import (
    array_namespace,
    Array,
    is_numpy,
    ArrayLike,
    is_lazy_array,
    xp_capabilities,
    xp_promote,
)
from scipy._lib.array_api_compat import device as xp_device
import scipy._lib.array_api_extra as xpx
from scipy._lib._util import _transition_to_rng, broadcastable

backend_registry = {array_namespace(np.empty(0)): cython_backend}


def select_backend(xp: ModuleType, cython_compatible: bool):
    """Select the backend for the given array library.

    We need this selection function because the Cython backend for numpy does not
    support quaternions of arbitrary dimensions. We therefore only use the Array API
    backend for numpy if we are dealing with rotations of more than one leading
    dimension.
    """
    if is_numpy(xp) and not cython_compatible:
        return xp_backend
    return backend_registry.get(xp, xp_backend)


@xp_capabilities()
def _promote(*args: tuple[ArrayLike, ...], xp: ModuleType) -> Array:
    """Promote arrays to float64 for numpy, else according to the Array API spec.

    The return array dtype follows the following rules:
    - If quat is an ArrayLike or NumPy array, we always promote to float64
    - If quat is an Array from frameworks other than NumPy, we preserve the precision
      of the input array dtype.

    The first rule is required by the cython backend signatures that expect
    cython.double views. The second rule is necessary to promote non-floating arrays
    to the correct type in frameworks that may not support double precision (e.g.
    jax by default).
    """
    if is_numpy(xp):
        args += (np.empty(0, dtype=np.float64),)  # Force float64 conversion
        out = xp_promote(*args, force_floating=True, xp=xp)
        if len(args) == 2:  # One argument was passed  + the added empty array
            return out[0]
        return out[:-1]
    return xp_promote(*args, force_floating=True, xp=xp)


class Rotation:
    """Rotation in 3 dimensions.

    This class provides an interface to initialize from and represent rotations
    with:

    - Quaternions
    - Rotation Matrices
    - Rotation Vectors
    - Modified Rodrigues Parameters
    - Euler Angles
    - Davenport Angles (Generalized Euler Angles)

    The following operations on rotations are supported:

    - Application on vectors
    - Rotation Composition
    - Rotation Inversion
    - Rotation Indexing

    A `Rotation` instance can contain a single rotation transform or rotations of
    multiple leading dimensions. E.g., it is possible to have an N-dimensional array of
    (N, M, K) rotations. When applied to other rotations or vectors, standard
    broadcasting rules apply.

    Indexing within a rotation is supported to access a subset of the rotations stored
    in a `Rotation` instance.

    To create `Rotation` objects use ``from_...`` methods (see examples below).
    ``Rotation(...)`` is not supposed to be instantiated directly.

    Attributes
    ----------
    single

    Methods
    -------
    __len__
    from_quat
    from_matrix
    from_rotvec
    from_mrp
    from_euler
    from_davenport
    as_quat
    as_matrix
    as_rotvec
    as_mrp
    as_euler
    as_davenport
    concatenate
    apply
    __mul__
    __pow__
    inv
    magnitude
    approx_equal
    mean
    reduce
    create_group
    __getitem__
    identity
    random
    align_vectors

    See Also
    --------
    Slerp

    Notes
    -----
    .. versionadded:: 1.2.0

    Examples
    --------
    >>> from scipy.spatial.transform import Rotation as R
    >>> import numpy as np

    A `Rotation` instance can be initialized in any of the above formats and
    converted to any of the others. The underlying object is independent of the
    representation used for initialization.

    Consider a counter-clockwise rotation of 90 degrees about the z-axis. This
    corresponds to the following quaternion (in scalar-last format):

    >>> r = R.from_quat([0, 0, np.sin(np.pi/4), np.cos(np.pi/4)])

    The rotation can be expressed in any of the other formats:

    >>> r.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
    [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
    [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The same rotation can be initialized using a rotation matrix:

    >>> r = R.from_matrix([[0, -1, 0],
    ...                    [1, 0, 0],
    ...                    [0, 0, 1]])

    Representation in other formats:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The rotation vector corresponding to this rotation is given by:

    >>> r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))

    Representation in other formats:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_euler('zyx', degrees=True)
    array([90.,  0.,  0.])

    The ``from_euler`` method is quite flexible in the range of input formats
    it supports. Here we initialize a single rotation about a single axis:

    >>> r = R.from_euler('z', 90, degrees=True)

    Again, the object is representation independent and can be converted to any
    other format:

    >>> r.as_quat()
    array([0.        , 0.        , 0.70710678, 0.70710678])
    >>> r.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> r.as_rotvec()
    array([0.        , 0.        , 1.57079633])

    It is also possible to initialize multiple rotations in a single instance
    using any of the ``from_...`` functions. Here we initialize a stack of 3
    rotations using the ``from_euler`` method:

    >>> r = R.from_euler('zyx', [
    ... [90, 0, 0],
    ... [0, 45, 0],
    ... [45, 60, 30]], degrees=True)

    The other representations also now return a stack of 3 rotations. For
    example:

    >>> r.as_quat()
    array([[0.        , 0.        , 0.70710678, 0.70710678],
           [0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])

    Applying the above rotations onto a vector:

    >>> v = [1, 2, 3]
    >>> r.apply(v)
    array([[-2.        ,  1.        ,  3.        ],
           [ 2.82842712,  2.        ,  1.41421356],
           [ 2.24452282,  0.78093109,  2.89002836]])

    A `Rotation` instance can be indexed and sliced as if it were an ND array:

    >>> r.as_quat()
    array([[0.        , 0.        , 0.70710678, 0.70710678],
           [0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])
    >>> p = r[0]
    >>> p.as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
    >>> q = r[1:3]
    >>> q.as_quat()
    array([[0.        , 0.38268343, 0.        , 0.92387953],
           [0.39190384, 0.36042341, 0.43967974, 0.72331741]])

    In fact it can be converted to numpy.array:

    >>> r_array = np.asarray(r)
    >>> r_array.shape
    (3,)
    >>> r_array[0].as_matrix()
    array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
           [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
           [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])

    Multiple rotations can be composed using the ``*`` operator:

    >>> r1 = R.from_euler('z', 90, degrees=True)
    >>> r2 = R.from_rotvec([np.pi/4, 0, 0])
    >>> v = [1, 2, 3]
    >>> r2.apply(r1.apply(v))
    array([-2.        , -1.41421356,  2.82842712])
    >>> r3 = r2 * r1 # Note the order
    >>> r3.apply(v)
    array([-2.        , -1.41421356,  2.82842712])

    A rotation can be composed with itself using the ``**`` operator:

    >>> p = R.from_rotvec([1, 0, 0])
    >>> q = p ** 2
    >>> q.as_rotvec()
    array([2., 0., 0.])

    Finally, it is also possible to invert rotations:

    >>> r1 = R.from_euler('z', [[90], [45]], degrees=True)
    >>> r2 = r1.inv()
    >>> r2.as_euler('zyx', degrees=True)
    array([[-90.,   0.,   0.],
           [-45.,   0.,   0.]])

    The following function can be used to plot rotations with Matplotlib by
    showing how they transform the standard x, y, z coordinate axes:

    >>> import matplotlib.pyplot as plt

    >>> def plot_rotated_axes(ax, r, name=None, offset=(0, 0, 0), scale=1):
    ...     colors = ("#FF6666", "#005533", "#1199EE")  # Colorblind-safe RGB
    ...     loc = np.array([offset, offset])
    ...     for i, (axis, c) in enumerate(zip((ax.xaxis, ax.yaxis, ax.zaxis),
    ...                                       colors)):
    ...         axlabel = axis.axis_name
    ...         axis.set_label_text(axlabel)
    ...         axis.label.set_color(c)
    ...         axis.line.set_color(c)
    ...         axis.set_tick_params(colors=c)
    ...         line = np.zeros((2, 3))
    ...         line[1, i] = scale
    ...         line_rot = r.apply(line)
    ...         line_plot = line_rot + loc
    ...         ax.plot(line_plot[:, 0], line_plot[:, 1], line_plot[:, 2], c)
    ...         text_loc = line[1]*1.2
    ...         text_loc_rot = r.apply(text_loc)
    ...         text_plot = text_loc_rot + loc[0]
    ...         ax.text(*text_plot, axlabel.upper(), color=c,
    ...                 va="center", ha="center")
    ...     ax.text(*offset, name, color="k", va="center", ha="center",
    ...             bbox={"fc": "w", "alpha": 0.8, "boxstyle": "circle"})

    Create three rotations - the identity and two Euler rotations using
    intrinsic and extrinsic conventions:

    >>> r0 = R.identity()
    >>> r1 = R.from_euler("ZYX", [90, -30, 0], degrees=True)  # intrinsic
    >>> r2 = R.from_euler("zyx", [90, -30, 0], degrees=True)  # extrinsic

    Add all three rotations to a single plot:

    >>> ax = plt.figure().add_subplot(projection="3d", proj_type="ortho")
    >>> plot_rotated_axes(ax, r0, name="r0", offset=(0, 0, 0))
    >>> plot_rotated_axes(ax, r1, name="r1", offset=(3, 0, 0))
    >>> plot_rotated_axes(ax, r2, name="r2", offset=(6, 0, 0))
    >>> _ = ax.annotate(
    ...     "r0: Identity Rotation\\n"
    ...     "r1: Intrinsic Euler Rotation (ZYX)\\n"
    ...     "r2: Extrinsic Euler Rotation (zyx)",
    ...     xy=(0.6, 0.7), xycoords="axes fraction", ha="left"
    ... )
    >>> ax.set(xlim=(-1.25, 7.25), ylim=(-1.25, 1.25), zlim=(-1.25, 1.25))
    >>> ax.set(xticks=range(-1, 8), yticks=[-1, 0, 1], zticks=[-1, 0, 1])
    >>> ax.set_aspect("equal", adjustable="box")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.tight_layout()

    Show the plot:

    >>> plt.show()

    These examples serve as an overview into the `Rotation` class and highlight
    major functionalities. For more thorough examples of the range of input and
    output formats supported, consult the individual method's examples.

    """

    def __init__(
        self,
        quat: ArrayLike,
        normalize: bool = True,
        copy: bool = True,
        scalar_first: bool = False,
    ):
        xp = array_namespace(quat)
        self._xp = xp
        quat = _promote(quat, xp=xp)
        if quat.shape[-1] != 4:
            raise ValueError(
                f"Expected `quat` to have shape (..., 4), got {quat.shape}."
            )
        # Single NumPy quats or list of quats are accelerated by the cython backend.
        # This backend needs inputs with fixed ndim, so we always expand to 2D and
        # select the 0th element if quat was single to get the correct shape. For other
        # frameworks and quaternion tensors we use the generic array API backend.
        self._single = quat.ndim == 1 and is_numpy(xp)
        if self._single:
            quat = xpx.atleast_nd(quat, ndim=2, xp=xp)
        self._backend = select_backend(xp, cython_compatible=quat.ndim < 3)
        self._quat: Array = self._backend.from_quat(
            quat, normalize=normalize, copy=copy, scalar_first=scalar_first
        )

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_quat(quat: ArrayLike, *, scalar_first: bool = False) -> Rotation:
        """Initialize from quaternions.

        Rotations in 3 dimensions can be represented using unit norm
        quaternions [1]_.

        The 4 components of a quaternion are divided into a scalar part ``w``
        and a vector part ``(x, y, z)`` and can be expressed from the angle
        ``theta`` and the axis ``n`` of a rotation as follows::

            w = cos(theta / 2)
            x = sin(theta / 2) * n_x
            y = sin(theta / 2) * n_y
            z = sin(theta / 2) * n_z

        There are 2 conventions to order the components in a quaternion:

        - scalar-first order -- ``(w, x, y, z)``
        - scalar-last order -- ``(x, y, z, w)``

        The choice is controlled by `scalar_first` argument.
        By default, it is False and the scalar-last order is assumed.

        Advanced users may be interested in the "double cover" of 3D space by
        the quaternion representation [2]_. As of version 1.11.0, the
        following subset (and only this subset) of operations on a `Rotation`
        ``r`` corresponding to a quaternion ``q`` are guaranteed to preserve
        the double cover property: ``r = Rotation.from_quat(q)``,
        ``r.as_quat(canonical=False)``, ``r.inv()``, and composition using the
        ``*`` operator such as ``r*r``.

        Parameters
        ----------
        quat : array_like, shape (..., 4)
            Each row is a (possibly non-unit norm) quaternion representing an
            active rotation. Each quaternion will be normalized to unit norm.
        scalar_first : bool, optional
            Whether the scalar component goes first or last.
            Default is False, i.e. the scalar-last order is assumed.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input quaternions.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
        .. [2] Hanson, Andrew J. "Visualizing quaternions."
            Morgan Kaufmann Publishers Inc., San Francisco, CA. 2006.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        A rotation can be initialized from a quaternion with the scalar-last
        (default) or scalar-first component order as shown below:

        >>> r = R.from_quat([0, 0, 0, 1])
        >>> r.as_matrix()
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])
        >>> r = R.from_quat([1, 0, 0, 0], scalar_first=True)
        >>> r.as_matrix()
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])

        It is possible to initialize multiple rotations in a single object by
        passing an N-dimensional array:

        >>> r = R.from_quat([[
        ... [1, 0, 0, 0],
        ... [0, 0, 0, 1]
        ... ]])
        >>> r.as_quat()
        array([[[1., 0., 0., 0.],
                [0., 0., 0., 1.]]])
        >>> r.as_quat().shape
        (1, 2, 4)

        It is also possible to have a stack of a single rotation:

        >>> r = R.from_quat([[0, 0, 0, 1]])
        >>> r.as_quat()
        array([[0., 0., 0., 1.]])
        >>> r.as_quat().shape
        (1, 4)

        Quaternions are normalized before initialization.

        >>> r = R.from_quat([0, 0, 1, 1])
        >>> r.as_quat()
        array([0.        , 0.        , 0.70710678, 0.70710678])
        """
        return Rotation(quat, normalize=True, scalar_first=scalar_first)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_matrix(matrix: ArrayLike, *, assume_valid: bool = False) -> Rotation:
        """Initialize from rotation matrix.

        Rotations in 3 dimensions can be represented with 3 x 3 orthogonal
        matrices [1]_. If the input is not orthogonal, an approximation is
        created by orthogonalizing the input matrix using the method described
        in [2]_, and then converting the orthogonal rotation matrices to
        quaternions using the algorithm described in [3]_. Matrices must be
        right-handed.

        Parameters
        ----------
        matrix : array_like, shape (..., 3, 3)
            A single matrix or an ND array of matrices, where the last two dimensions
            contain the rotation matrices.
        assume_valid : bool, optional
            Must be False unless users can guarantee the input is a valid rotation
            matrix, i.e. it is orthogonal, rows and columns have unit norm and the
            determinant is 1. Setting this to True without ensuring these properties
            is unsafe and will silently lead to incorrect results. If True,
            normalization steps are skipped, which can improve runtime performance.
            Default is False.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by the rotation
            matrices.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions
        .. [2] https://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
        .. [3] F. Landis Markley, "Unit Quaternion from Rotation Matrix",
               Journal of guidance, control, and dynamics vol. 31.2, pp.
               440-442, 2008.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Initialize a single rotation:

        >>> r = R.from_matrix([
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]])
        >>> r.single
        True
        >>> r.as_matrix().shape
        (3, 3)

        Initialize multiple rotations in a single object:

        >>> r = R.from_matrix([
        ... [
        ...     [0, -1, 0],
        ...     [1, 0, 0],
        ...     [0, 0, 1],
        ... ],
        ... [
        ...     [1, 0, 0],
        ...     [0, 0, -1],
        ...     [0, 1, 0],
        ... ]])
        >>> r.as_matrix().shape
        (2, 3, 3)
        >>> r.single
        False
        >>> len(r)
        2

        If input matrices are not special orthogonal (orthogonal with
        determinant equal to +1), then a special orthogonal estimate is stored:

        >>> a = np.array([
        ... [0, -0.5, 0],
        ... [0.5, 0, 0],
        ... [0, 0, 0.5]])
        >>> np.linalg.det(a)
        0.125
        >>> r = R.from_matrix(a)
        >>> matrix = r.as_matrix()
        >>> matrix
        array([[ 0., -1.,  0.],
               [ 1.,  0.,  0.],
               [ 0.,  0.,  1.]])
        >>> np.linalg.det(matrix)
        1.0

        It is also possible to have a stack containing a single rotation:

        >>> r = R.from_matrix([[
        ... [0, -1, 0],
        ... [1, 0, 0],
        ... [0, 0, 1]]])
        >>> r.as_matrix()
        array([[[ 0., -1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0.,  1.]]])
        >>> r.as_matrix().shape
        (1, 3, 3)

        We can also create an N-dimensional array of rotations:

        >>> r = R.from_matrix(np.tile(np.eye(3), (2, 3, 1, 1)))
        >>> r.shape
        (2, 3)

        Notes
        -----
        This function was called from_dcm before.

        .. versionadded:: 1.4.0
        """
        xp = array_namespace(matrix)
        matrix = _promote(matrix, xp=xp)
        # Resulting quat will have 1 less dimension than matrix
        backend = select_backend(xp, cython_compatible=matrix.ndim < 4)
        quat = backend.from_matrix(matrix, assume_valid=assume_valid)
        return Rotation._from_raw_quat(quat, xp=xp, backend=backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_rotvec(rotvec: ArrayLike, degrees: bool = False) -> Rotation:
        """Initialize from rotation vectors.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotation and whose norm gives the angle of rotation [1]_.

        Parameters
        ----------
        rotvec : array_like, shape (..., 3)
            A single vector or an ND array of vectors, where the last dimension
            contains the rotation vectors.
        degrees : bool, optional
            If True, then the given magnitudes are assumed to be in degrees.
            Default is False.

            .. versionadded:: 1.7.0

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input rotation
            vectors.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Initialize a single rotation:

        >>> r = R.from_rotvec(np.pi/2 * np.array([0, 0, 1]))
        >>> r.as_rotvec()
        array([0.        , 0.        , 1.57079633])
        >>> r.as_rotvec().shape
        (3,)

        Initialize a rotation in degrees, and view it in degrees:

        >>> r = R.from_rotvec(45 * np.array([0, 1, 0]), degrees=True)
        >>> r.as_rotvec(degrees=True)
        array([ 0., 45.,  0.])

        Initialize multiple rotations in one object:

        >>> r = R.from_rotvec([
        ... [0, 0, np.pi/2],
        ... [np.pi/2, 0, 0]])
        >>> r.as_rotvec()
        array([[0.        , 0.        , 1.57079633],
               [1.57079633, 0.        , 0.        ]])
        >>> r.as_rotvec().shape
        (2, 3)

        It is also possible to have a stack of a single rotation:

        >>> r = R.from_rotvec([[0, 0, np.pi/2]])
        >>> r.as_rotvec().shape
        (1, 3)

        """
        xp = array_namespace(rotvec)
        rotvec = _promote(rotvec, xp=xp)
        backend = select_backend(xp, cython_compatible=rotvec.ndim < 3)
        quat = backend.from_rotvec(rotvec, degrees=degrees)
        return Rotation._from_raw_quat(quat, xp=xp, backend=backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_euler(seq: str, angles: ArrayLike, degrees: bool = False) -> Rotation:
        """Initialize from Euler angles.

        Rotations in 3-D can be represented by a sequence of 3
        rotations around a sequence of axes. In theory, any three axes spanning
        the 3-D Euclidean space are enough. In practice, the axes of rotation are
        chosen to be the basis vectors.

        The three rotations can either be in a global frame of reference
        (extrinsic) or in a body centred frame of reference (intrinsic), which
        is attached to, and moves with, the object under rotation [1]_.

        Parameters
        ----------
        seq : string
            Specifies sequence of axes for rotations. Up to 3 characters
            belonging to the set {'X', 'Y', 'Z'} for intrinsic rotations, or
            {'x', 'y', 'z'} for extrinsic rotations. Extrinsic and intrinsic
            rotations cannot be mixed in one function call.
        angles : float or array_like, shape (...,  [1 or 2 or 3])
            Euler angles specified in radians (`degrees` is False) or degrees
            (`degrees` is True).
            Each character in `seq` defines one axis around which `angles` turns.
            The resulting rotation has the shape np.atleast_1d(angles).shape[:-1].
            Dimensionless angles are thus only valid for single character `seq`.

        degrees : bool, optional
            If True, then the given angles are assumed to be in degrees.
            Default is False.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotation represented by the sequence of
            rotations around given axes with given angles.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Initialize a single rotation along a single axis:

        >>> r = R.from_euler('x', 90, degrees=True)
        >>> r.as_quat().shape
        (4,)

        Initialize a single rotation with a given axis sequence:

        >>> r = R.from_euler('zyx', [90, 45, 30], degrees=True)
        >>> r.as_quat().shape
        (4,)

        Initialize a stack with a single rotation around a single axis:

        >>> r = R.from_euler('x', [[90]], degrees=True)
        >>> r.as_quat().shape
        (1, 4)

        Initialize a stack with a single rotation with an axis sequence:

        >>> r = R.from_euler('zyx', [[90, 45, 30]], degrees=True)
        >>> r.as_quat().shape
        (1, 4)

        Initialize multiple elementary rotations in one object:

        >>> r = R.from_euler('x', [[90], [45], [30]], degrees=True)
        >>> r.as_quat().shape
        (3, 4)

        Initialize multiple rotations in one object:

        >>> r = R.from_euler('zyx', [[90, 45, 30], [35, 45, 90]], degrees=True)
        >>> r.as_quat().shape
        (2, 4)

        """
        xp = array_namespace(angles)
        angles = _promote(angles, xp=xp)
        backend = select_backend(xp, cython_compatible=angles.ndim < 3)
        quat = backend.from_euler(seq, angles, degrees=degrees)
        return Rotation._from_raw_quat(quat, xp=xp, backend=backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_davenport(
        axes: ArrayLike,
        order: str,
        angles: ArrayLike | float,
        degrees: bool = False,
    ) -> Rotation:
        """Initialize from Davenport angles.

        Rotations in 3-D can be represented by a sequence of 3
        rotations around a sequence of axes.

        The three rotations can either be in a global frame of reference
        (extrinsic) or in a body centred frame of reference (intrinsic), which
        is attached to, and moves with, the object under rotation [1]_.

        For both Euler angles and Davenport angles, consecutive axes must
        be are orthogonal (``axis2`` is orthogonal to both ``axis1`` and
        ``axis3``). For Euler angles, there is an additional relationship
        between ``axis1`` or ``axis3``, with two possibilities:

            - ``axis1`` and ``axis3`` are also orthogonal (asymmetric sequence)
            - ``axis1 == axis3`` (symmetric sequence)

        For Davenport angles, this last relationship is relaxed [2]_, and only
        the consecutive orthogonal axes requirement is maintained.

        Parameters
        ----------
        axes : array_like, shape (3,) or (..., [1 or 2 or 3], 3)
            Axis of rotation, if one dimensional. If two or more dimensional, describes
            the sequence of axes for rotations, where each axes[..., i, :] is the ith
            axis. If more than one axis is given, then the second axis must be
            orthogonal to both the first and third axes.
        order : string
            If it is equal to 'e' or 'extrinsic', the sequence will be
            extrinsic. If it is equal to 'i' or 'intrinsic', sequence
            will be treated as intrinsic.
        angles : float or array_like, shape (..., [1 or 2 or 3])
            Angles specified in radians (`degrees` is False) or degrees
            (`degrees` is True).
            Each angle i in the last dimension of `angles` turns around the corresponding
            axis axis[..., i, :]. The resulting rotation has the shape
            np.broadcast_shapes(np.atleast_2d(axes).shape[:-2], np.atleast_1d(angles).shape[:-1])
            Dimensionless angles are thus only valid for a single axis.

        degrees : bool, optional
            If True, then the given angles are assumed to be in degrees.
            Default is False.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotation represented by the sequence of
            rotations around given axes with given angles.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations
        .. [2] Shuster, Malcolm & Markley, Landis. (2003). Generalization of
               the Euler Angles. Journal of the Astronautical Sciences. 51. 123-132.
               10.1007/BF03546304.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Davenport angles are a generalization of Euler angles, when we use the
        canonical basis axes:

        >>> ex = [1, 0, 0]
        >>> ey = [0, 1, 0]
        >>> ez = [0, 0, 1]

        Initialize a single rotation with a given axis sequence:

        >>> axes = [ez, ey, ex]
        >>> r = R.from_davenport(axes, 'extrinsic', [90, 0, 0], degrees=True)
        >>> r.as_quat().shape
        (4,)

        It is equivalent to Euler angles in this case:

        >>> r.as_euler('zyx', degrees=True)
        array([90.,  0., -0.])

        Initialize multiple rotations in one object:

        >>> r = R.from_davenport(axes, 'extrinsic', [[90, 45, 30], [35, 45, 90]], degrees=True)
        >>> r.as_quat().shape
        (2, 4)

        Using only one or two axes is also possible:

        >>> r = R.from_davenport([ez, ex], 'extrinsic', [[90, 45], [35, 45]], degrees=True)
        >>> r.as_quat().shape
        (2, 4)

        Non-canonical axes are possible, and they do not need to be normalized,
        as long as consecutive axes are orthogonal:

        >>> e1 = [2, 0, 0]
        >>> e2 = [0, 1, 0]
        >>> e3 = [1, 0, 1]
        >>> axes = [e1, e2, e3]
        >>> r = R.from_davenport(axes, 'extrinsic', [90, 45, 30], degrees=True)
        >>> r.as_quat()
        [ 0.701057,  0.430459, -0.092296,  0.560986]
        """  # noqa: E501
        xp = array_namespace(axes)
        axes, angles = _promote(axes, angles, xp=xp)
        cython_compatible = axes.ndim < 3 and angles.ndim < 2
        backend = select_backend(xp, cython_compatible=cython_compatible)
        quat = backend.from_davenport(axes, order, angles, degrees)
        return Rotation._from_raw_quat(quat, xp=xp, backend=backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_mrp(mrp: ArrayLike) -> Rotation:
        """Initialize from Modified Rodrigues Parameters (MRPs).

        MRPs are a 3 dimensional vector co-directional to the axis of rotation and whose
        magnitude is equal to ``tan(theta / 4)``, where ``theta`` is the angle of
        rotation (in radians) [1]_.

        MRPs have a singularity at 360 degrees which can be avoided by ensuring the
        angle of rotation does not exceed 180 degrees, i.e. switching the direction of
        the rotation when it is past 180 degrees.

        Parameters
        ----------
        mrp : array_like, shape (..., 3)
            A single vector or an ND array of vectors, where the last dimension
            contains the rotation parameters.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the rotations represented by input MRPs.

        References
        ----------
        .. [1] Shuster, M. D. "A Survey of Attitude Representations",
               The Journal of Astronautical Sciences, Vol. 41, No.4, 1993,
               pp. 475-476

        Notes
        -----

        .. versionadded:: 1.6.0

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Initialize a single rotation:

        >>> r = R.from_mrp([0, 0, 1])
        >>> r.as_euler('xyz', degrees=True)
        array([0.        , 0.        , 180.      ])
        >>> r.as_euler('xyz').shape
        (3,)

        Initialize multiple rotations in one object:

        >>> r = R.from_mrp([
        ... [0, 0, 1],
        ... [1, 0, 0]])
        >>> r.as_euler('xyz', degrees=True)
        array([[0.        , 0.        , 180.      ],
               [180.0     , 0.        , 0.        ]])
        >>> r.as_euler('xyz').shape
        (2, 3)

        It is also possible to have a stack of a single rotation:

        >>> r = R.from_mrp([[0, 0, np.pi/2]])
        >>> r.as_euler('xyz').shape
        (1, 3)

        """
        xp = array_namespace(mrp)
        mrp = _promote(mrp, xp=xp)
        backend = select_backend(xp, cython_compatible=mrp.ndim < 3)
        quat = backend.from_mrp(mrp)
        return Rotation._from_raw_quat(quat, xp=xp, backend=backend)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_quat(self, canonical: bool = False, *, scalar_first: bool = False) -> Array:
        """Represent as quaternions.

        Rotations in 3 dimensions can be represented using unit norm
        quaternions [1]_.

        The 4 components of a quaternion are divided into a scalar part ``w``
        and a vector part ``(x, y, z)`` and can be expressed from the angle
        ``theta`` and the axis ``n`` of a rotation as follows::

            w = cos(theta / 2)
            x = sin(theta / 2) * n_x
            y = sin(theta / 2) * n_y
            z = sin(theta / 2) * n_z

        There are 2 conventions to order the components in a quaternion:

        - scalar-first order -- ``(w, x, y, z)``
        - scalar-last order -- ``(x, y, z, w)``

        The choice is controlled by `scalar_first` argument.
        By default, it is False and the scalar-last order is used.

        The mapping from quaternions to rotations is
        two-to-one, i.e. quaternions ``q`` and ``-q``, where ``-q`` simply
        reverses the sign of each component, represent the same spatial
        rotation.

        Parameters
        ----------
        canonical : `bool`, default False
            Whether to map the redundant double cover of rotation space to a
            unique "canonical" single cover. If True, then the quaternion is
            chosen from {q, -q} such that the w term is positive. If the w term
            is 0, then the quaternion is chosen such that the first nonzero
            term of the x, y, and z terms is positive.
        scalar_first : bool, optional
            Whether the scalar component goes first or last.
            Default is False, i.e. the scalar-last order is used.

        Returns
        -------
        quat : `numpy.ndarray`, shape (..., 4)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        A rotation can be represented as a quaternion with either scalar-last
        (default) or scalar-first component order.
        This is shown for a single rotation:

        >>> r = R.from_matrix(np.eye(3))
        >>> r.as_quat()
        array([0., 0., 0., 1.])
        >>> r.as_quat(scalar_first=True)
        array([1., 0., 0., 0.])

        The resulting shape of the quaternion is always the shape of the Rotation
        object with an added last dimension of size 4. E.g. when the `Rotation` object
        contains an N-dimensional array (N, M, K) of rotations, the result will be a
        4-dimensional array:

        >>> r = R.from_rotvec(np.ones((2, 3, 4, 3)))
        >>> r.as_quat().shape
        (2, 3, 4, 4)

        Quaternions can be mapped from a redundant double cover of the
        rotation space to a canonical representation with a positive w term.

        >>> r = R.from_quat([0, 0, 0, -1])
        >>> r.as_quat()
        array([0. , 0. , 0. , -1.])
        >>> r.as_quat(canonical=True)
        array([0. , 0. , 0. , 1.])
        """
        quat = self._backend.as_quat(
            self._quat, canonical=canonical, scalar_first=scalar_first
        )
        if self._single:
            return quat[0, ...]
        return quat

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_matrix(self) -> Array:
        """Represent as rotation matrix.

        3D rotations can be represented using rotation matrices, which
        are 3 x 3 real orthogonal matrices with determinant equal to +1 [1]_.

        Returns
        -------
        matrix : ndarray, shape (..., 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Rotation_matrix#In_three_dimensions

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_matrix()
        array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
               [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
               [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
        >>> r.as_matrix().shape
        (3, 3)

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[1, 1, 0, 0]])
        >>> r.as_matrix()
        array([[[ 0.,  1.,  0.],
                [ 1.,  0.,  0.],
                [ 0.,  0., -1.]]])
        >>> r.as_matrix().shape
        (1, 3, 3)

        Represent multiple rotations:

        >>> r = R.from_rotvec([[np.pi/2, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_matrix()
        array([[[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  2.22044605e-16, -1.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00,  2.22044605e-16]],
               [[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
                [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]]])
        >>> r.as_matrix().shape
        (2, 3, 3)

        Notes
        -----
        This function was called as_dcm before.

        .. versionadded:: 1.4.0
        """
        matrix = self._backend.as_matrix(self._quat)
        if self._single:
            return matrix[0, ...]
        return matrix

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_rotvec(self, degrees: bool = False) -> Array:
        """Represent as rotation vectors.

        A rotation vector is a 3 dimensional vector which is co-directional to
        the axis of rotation and whose norm gives the angle of rotation [1]_.

        Parameters
        ----------
        degrees : boolean, optional
            Returned magnitudes are in degrees if this flag is True, else they are
            in radians. Default is False.

            .. versionadded:: 1.7.0

        Returns
        -------
        rotvec : ndarray, shape (..., 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation#Rotation_vector

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Represent a single rotation:

        >>> r = R.from_euler('z', 90, degrees=True)
        >>> r.as_rotvec()
        array([0.        , 0.        , 1.57079633])
        >>> r.as_rotvec().shape
        (3,)

        Represent a rotation in degrees:

        >>> r = R.from_euler('YX', (-90, -90), degrees=True)
        >>> s = r.as_rotvec(degrees=True)
        >>> s
        array([-69.2820323, -69.2820323, -69.2820323])
        >>> np.linalg.norm(s)
        120.00000000000001

        Represent a stack with a single rotation:

        >>> r = R.from_quat([[0, 0, 1, 1]])
        >>> r.as_rotvec()
        array([[0.        , 0.        , 1.57079633]])
        >>> r.as_rotvec().shape
        (1, 3)

        Represent multiple rotations in a single object:

        >>> r = R.from_quat([[0, 0, 1, 1], [1, 1, 0, 1]])
        >>> r.as_rotvec()
        array([[0.        , 0.        , 1.57079633],
               [1.35102172, 1.35102172, 0.        ]])
        >>> r.as_rotvec().shape
        (2, 3)

        """
        rotvec = self._backend.as_rotvec(self._quat, degrees=degrees)
        if self._single:
            return rotvec[0, ...]
        return rotvec

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_euler(
        self, seq: str, degrees: bool = False, *, suppress_warnings: bool = False
    ) -> Array:
        """Represent as Euler angles.

        Any orientation can be expressed as a composition of 3 elementary
        rotations. Once the axis sequence has been chosen, Euler angles define
        the angle of rotation around each respective axis [1]_.

        The algorithm from [2]_ has been used to calculate Euler angles for the
        rotation about a given sequence of axes.

        Euler angles suffer from the problem of gimbal lock [3]_, where the
        representation loses a degree of freedom and it is not possible to
        determine the first and third angles uniquely. In this case,
        a warning is raised (unless the ``suppress_warnings`` option is used),
        and the third angle is set to zero. Note however that the returned
        angles still represent the correct rotation.

        Parameters
        ----------
        seq : string, length 3
            3 characters belonging to the set {'X', 'Y', 'Z'} for intrinsic
            rotations, or {'x', 'y', 'z'} for extrinsic rotations [1]_.
            Adjacent axes cannot be the same.
            Extrinsic and intrinsic rotations cannot be mixed in one function
            call.
        degrees : boolean, optional
            Returned angles are in degrees if this flag is True, else they are
            in radians. Default is False.
        suppress_warnings : boolean, optional
            Disable warnings about gimbal lock. Default is False.

        Returns
        -------
        angles : ndarray, shape (..., 3)
            Shape depends on shape of inputs used to initialize object.
            The returned angles are in the range:

            - First angle belongs to [-180, 180] degrees (both inclusive)
            - Third angle belongs to [-180, 180] degrees (both inclusive)
            - Second angle belongs to:

                - [-90, 90] degrees if all axes are different (like xyz)
                - [0, 180] degrees if first and third axes are the same
                  (like zxz)

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Euler_angles#Definition_by_intrinsic_rotations
        .. [2] Bernardes E, Viollet S (2022) Quaternion to Euler angles
               conversion: A direct, general and computationally efficient
               method. PLoS ONE 17(11): e0276302.
               https://doi.org/10.1371/journal.pone.0276302
        .. [3] https://en.wikipedia.org/wiki/Gimbal_lock#In_applied_mathematics

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_euler('zxy', degrees=True)
        array([90.,  0.,  0.])
        >>> r.as_euler('zxy', degrees=True).shape
        (3,)

        Represent a stack of single rotation:

        >>> r = R.from_rotvec([[0, 0, np.pi/2]])
        >>> r.as_euler('zxy', degrees=True)
        array([[90.,  0.,  0.]])
        >>> r.as_euler('zxy', degrees=True).shape
        (1, 3)

        Represent multiple rotations in a single object:

        >>> r = R.from_rotvec([
        ... [0, 0, np.pi/2],
        ... [0, -np.pi/3, 0],
        ... [np.pi/4, 0, 0]])
        >>> r.as_euler('zxy', degrees=True)
        array([[ 90.,   0.,   0.],
               [  0.,   0., -60.],
               [  0.,  45.,   0.]])
        >>> r.as_euler('zxy', degrees=True).shape
        (3, 3)

        """
        euler = self._backend.as_euler(
            self._quat, seq, degrees=degrees, suppress_warnings=suppress_warnings
        )
        if self._single:
            return euler[0, ...]
        return euler

    @xp_capabilities(
        skip_backends=[
            ("dask.array", "missing linalg.cross/det functions and .mT attribute"),
            ("cupy", "missing .mT attribute in cupy<14.*"),
        ]
    )
    def as_davenport(
        self,
        axes: ArrayLike,
        order: str,
        degrees: bool = False,
        *,
        suppress_warnings: bool = False,
    ) -> Array:
        """Represent as Davenport angles.

        Any orientation can be expressed as a composition of 3 elementary
        rotations.

        For both Euler angles and Davenport angles, consecutive axes must
        be are orthogonal (``axis2`` is orthogonal to both ``axis1`` and
        ``axis3``). For Euler angles, there is an additional relationship
        between ``axis1`` or ``axis3``, with two possibilities:

            - ``axis1`` and ``axis3`` are also orthogonal (asymmetric sequence)
            - ``axis1 == axis3`` (symmetric sequence)

        For Davenport angles, this last relationship is relaxed [1]_, and only
        the consecutive orthogonal axes requirement is maintained.

        A slightly modified version of the algorithm from [2]_ has been used to
        calculate Davenport angles for the rotation about a given sequence of
        axes.

        Davenport angles, just like Euler angles, suffer from the problem of
        gimbal lock [3]_, where the representation loses a degree of freedom
        and it is not possible to determine the first and third angles
        uniquely. In this case, a warning is raised (unless the
        ``suppress_warnings`` option is used), and the third angle is set
        to zero. Note however that the returned angles still represent the
        correct rotation.

        Parameters
        ----------
        axes : array_like, shape (..., [1 or 2 or 3], 3) or (..., 3)
            Axis of rotation, if one dimensional. If N dimensional, describes the
            sequence of axes for rotations, where each axes[..., i, :] is the ith
            axis. If more than one axis is given, then the second axis must be
            orthogonal to both the first and third axes.
        order : string
            If it belongs to the set {'e', 'extrinsic'}, the sequence will be
            extrinsic. If it belongs to the set {'i', 'intrinsic'}, sequence
            will be treated as intrinsic.
        degrees : boolean, optional
            Returned angles are in degrees if this flag is True, else they are
            in radians. Default is False.
        suppress_warnings : boolean, optional
            Disable warnings about gimbal lock. Default is False.

        Returns
        -------
        angles : ndarray, shape (..., 3)
            Shape depends on shape of inputs used to initialize object.
            The returned angles are in the range:

            - First angle belongs to [-180, 180] degrees (both inclusive)
            - Third angle belongs to [-180, 180] degrees (both inclusive)
            - Second angle belongs to a set of size 180 degrees,
              given by: ``[-abs(lambda), 180 - abs(lambda)]``, where ``lambda``
              is the angle between the first and third axes.

        References
        ----------
        .. [1] Shuster, Malcolm & Markley, Landis. (2003). Generalization of
               the Euler Angles. Journal of the Astronautical Sciences. 51. 123-132.
               10.1007/BF03546304.
        .. [2] Bernardes E, Viollet S (2022) Quaternion to Euler angles
               conversion: A direct, general and computationally efficient method.
               PLoS ONE 17(11): e0276302. 10.1371/journal.pone.0276302
        .. [3] https://en.wikipedia.org/wiki/Gimbal_lock#In_applied_mathematics

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Davenport angles are a generalization of Euler angles, when we use the
        canonical basis axes:

        >>> ex = [1, 0, 0]
        >>> ey = [0, 1, 0]
        >>> ez = [0, 0, 1]

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True)
        array([90.,  0.,  0.])
        >>> r.as_euler('zxy', degrees=True)
        array([90.,  0.,  0.])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True).shape
        (3,)

        Represent a stack of single rotation:

        >>> r = R.from_rotvec([[0, 0, np.pi/2]])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True)
        array([[90.,  0.,  0.]])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True).shape
        (1, 3)

        Represent multiple rotations in a single object:

        >>> r = R.from_rotvec([
        ... [0, 0, 90],
        ... [45, 0, 0]], degrees=True)
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True)
        array([[90.,  0.,  0.],
               [ 0., 45.,  0.]])
        >>> r.as_davenport([ez, ex, ey], 'extrinsic', degrees=True).shape
        (2, 3)
        """
        axes = self._xp.asarray(
            axes, dtype=self._quat.dtype, device=xp_device(self._quat)
        )
        davenport = self._backend.as_davenport(
            self._quat, axes, order, degrees, suppress_warnings=suppress_warnings
        )
        if self._single:
            return davenport[0, ...]
        return davenport

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_mrp(self) -> Array:
        """Represent as Modified Rodrigues Parameters (MRPs).

        MRPs are a 3 dimensional vector co-directional to the axis of rotation and whose
        magnitude is equal to ``tan(theta / 4)``, where ``theta`` is the angle of
        rotation (in radians) [1]_.

        MRPs have a singularity at 360 degrees which can be avoided by ensuring the
        angle of rotation does not exceed 180 degrees, i.e. switching the direction of
        the rotation when it is past 180 degrees. This function will always return MRPs
        corresponding to a rotation of less than or equal to 180 degrees.

        Returns
        -------
        mrps : ndarray, shape (..., 3)
            Shape depends on shape of inputs used for initialization.

        References
        ----------
        .. [1] Shuster, M. D. "A Survey of Attitude Representations",
               The Journal of Astronautical Sciences, Vol. 41, No.4, 1993,
               pp. 475-476

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Represent a single rotation:

        >>> r = R.from_rotvec([0, 0, np.pi])
        >>> r.as_mrp()
        array([0.        , 0.        , 1.         ])
        >>> r.as_mrp().shape
        (3,)

        Represent a stack with a single rotation:

        >>> r = R.from_euler('xyz', [[180, 0, 0]], degrees=True)
        >>> r.as_mrp()
        array([[1.       , 0.        , 0.         ]])
        >>> r.as_mrp().shape
        (1, 3)

        Represent multiple rotations:

        >>> r = R.from_rotvec([[np.pi/2, 0, 0], [0, 0, np.pi/2]])
        >>> r.as_mrp()
        array([[0.41421356, 0.        , 0.        ],
               [0.        , 0.        , 0.41421356]])
        >>> r.as_mrp().shape
        (2, 3)

        Notes
        -----

        .. versionadded:: 1.6.0
        """
        mrp = self._backend.as_mrp(self._quat)
        if self._single:
            return mrp[0, ...]
        return mrp

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def concatenate(rotations: Rotation | Iterable[Rotation]) -> Rotation:
        """Concatenate a sequence of `Rotation` objects into a single object.

        This is useful if you want to, for example, take the mean of a set of
        rotations and need to pack them into a single object to do so.

        Parameters
        ----------
        rotations : sequence of `Rotation` objects
            The rotations to concatenate. If a single `Rotation` object is
            passed in, a copy is returned.

        Returns
        -------
        concatenated : `Rotation` instance
            The concatenated rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> r1 = R.from_rotvec([0, 0, 1])
        >>> r2 = R.from_rotvec([0, 0, 2])
        >>> rc = R.concatenate([r1, r2])
        >>> rc.as_rotvec()
        array([[0., 0., 1.],
               [0., 0., 2.]])
        >>> rc.mean().as_rotvec()
        array([0., 0., 1.5])

        Concatenation of a split rotation recovers the original object.

        >>> rs = [r for r in rc]
        >>> R.concatenate(rs).as_rotvec()
        array([[0., 0., 1.],
               [0., 0., 2.]])

        Note that it may be simpler to create the desired rotations by passing
        in a single list of the data during initialization, rather then by
        concatenating:

        >>> R.from_rotvec([[0, 0, 1], [0, 0, 2]]).as_rotvec()
        array([[0., 0., 1.],
               [0., 0., 2.]])

        Notes
        -----
        .. versionadded:: 1.8.0
        """
        if isinstance(rotations, Rotation):
            return Rotation(rotations.as_quat(), normalize=False, copy=True)
        if not all(isinstance(x, Rotation) for x in rotations):
            raise TypeError("input must contain Rotation objects only")

        xp = array_namespace(rotations[0].as_quat())
        quats = xp.concat(
            [xpx.atleast_nd(x.as_quat(), ndim=2, xp=xp) for x in rotations]
        )
        return Rotation._from_raw_quat(quats, xp=xp)

    @xp_capabilities(
        skip_backends=[
            ("dask.array", "missing linalg.cross/det functions and .mT attribute"),
            ("cupy", "missing .mT attribute in cupy<14.*"),
        ]
    )
    def apply(self, vectors: ArrayLike, inverse: bool = False) -> Array:
        """Apply this rotation to a set of vectors.

        If the original frame rotates to the final frame by this rotation, then
        its application to a vector can be seen in two ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical rotation of a vector being glued to the original
              frame as it rotates. In this case the vector components are
              expressed in the original frame before and after the rotation.

        In terms of rotation matrices, this application is the same as
        ``self.as_matrix() @ vectors``.

        Parameters
        ----------
        vectors : array_like, shape (..., 3)
            Each `vectors[..., :]` represents a vector in 3D space. The shape of
            rotations and shape of vectors given must follow standard numpy
            broadcasting rules: either one of them equals unity or they both
            equal each other.
        inverse : boolean, optional
            If True then the inverse of the rotation(s) is applied to the input
            vectors. Default is False.

        Returns
        -------
        rotated_vectors : ndarray, shape (..., 3)
            Result of applying rotation on input vectors.
            Shape is determined according to numpy broadcasting rules. I.e., the result
            will have the shape `np.broadcast_shapes(r.shape, v.shape[:-1]) + (3,)`

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Single rotation applied on a single vector:

        >>> vector = np.array([1, 0, 0])
        >>> r = R.from_rotvec([0, 0, np.pi/2])
        >>> r.as_matrix()
        array([[ 2.22044605e-16, -1.00000000e+00,  0.00000000e+00],
               [ 1.00000000e+00,  2.22044605e-16,  0.00000000e+00],
               [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
        >>> r.apply(vector)
        array([2.22044605e-16, 1.00000000e+00, 0.00000000e+00])
        >>> r.apply(vector).shape
        (3,)

        Single rotation applied on multiple vectors:

        >>> vectors = np.array([
        ... [1, 0, 0],
        ... [1, 2, 3]])
        >>> r = R.from_rotvec([0, 0, np.pi/4])
        >>> r.as_matrix()
        array([[ 0.70710678, -0.70710678,  0.        ],
               [ 0.70710678,  0.70710678,  0.        ],
               [ 0.        ,  0.        ,  1.        ]])
        >>> r.apply(vectors)
        array([[ 0.70710678,  0.70710678,  0.        ],
               [-0.70710678,  2.12132034,  3.        ]])
        >>> r.apply(vectors).shape
        (2, 3)

        Multiple rotations on a single vector:

        >>> r = R.from_rotvec([[0, 0, np.pi/4], [np.pi/2, 0, 0]])
        >>> vector = np.array([1,2,3])
        >>> r.as_matrix()
        array([[[ 7.07106781e-01, -7.07106781e-01,  0.00000000e+00],
                [ 7.07106781e-01,  7.07106781e-01,  0.00000000e+00],
                [ 0.00000000e+00,  0.00000000e+00,  1.00000000e+00]],
               [[ 1.00000000e+00,  0.00000000e+00,  0.00000000e+00],
                [ 0.00000000e+00,  2.22044605e-16, -1.00000000e+00],
                [ 0.00000000e+00,  1.00000000e+00,  2.22044605e-16]]])
        >>> r.apply(vector)
        array([[-0.70710678,  2.12132034,  3.        ],
               [ 1.        , -3.        ,  2.        ]])
        >>> r.apply(vector).shape
        (2, 3)

        Multiple rotations on multiple vectors. Each rotation is applied on the
        corresponding vector:

        >>> r = R.from_euler('zxy', [
        ... [0, 0, 90],
        ... [45, 30, 60]], degrees=True)
        >>> vectors = [
        ... [1, 2, 3],
        ... [1, 0, -1]]
        >>> r.apply(vectors)
        array([[ 3.        ,  2.        , -1.        ],
               [-0.09026039,  1.11237244, -0.86860844]])
        >>> r.apply(vectors).shape
        (2, 3)

        Broadcasting rules apply:

        >>> r = R.from_rotvec(np.tile([0, 0, np.pi/4], (5, 1, 4, 1)))
        >>> vectors = np.ones((3, 4, 3))
        >>> r.shape, vectors.shape
        ((5, 1, 4), (3, 4, 3))
        >>> r.apply(vectors).shape
        (5, 3, 4, 3)

        It is also possible to apply the inverse rotation:

        >>> r = R.from_euler('zxy', [
        ... [0, 0, 90],
        ... [45, 30, 60]], degrees=True)
        >>> vectors = [
        ... [1, 2, 3],
        ... [1, 0, -1]]
        >>> r.apply(vectors, inverse=True)
        array([[-3.        ,  2.        ,  1.        ],
               [ 1.09533535, -0.8365163 ,  0.3169873 ]])

        """
        vectors = self._xp.asarray(
            vectors, device=xp_device(self._quat), dtype=self._quat.dtype
        )
        single_vector = vectors.ndim == 1
        # Numpy optimization: The Cython backend typing requires us to have fixed
        # dimensions, so for the Numpy case we always broadcast the vector to 2D.
        if vectors.shape[-1] != 3:
            raise ValueError(f"Expected input of shape (..., 3), got {vectors.shape}.")
        if is_numpy(self._xp):
            vectors = xpx.atleast_nd(vectors, ndim=2, xp=self._xp)
        cython_compatible = self._quat.ndim < 3 and vectors.ndim < 3
        backend = select_backend(self._xp, cython_compatible=cython_compatible)
        result = backend.apply(self._quat, vectors, inverse=inverse)
        if self._single and single_vector:
            return result[0, ...]
        return result

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def __mul__(self, other: Rotation) -> Rotation | NotImplementedType:
        """Compose this rotation with the other.

        If `p` and `q` are two rotations, then the composition of 'q followed
        by p' is equivalent to `p * q`. In terms of rotation matrices,
        the composition can be expressed as
        ``p.as_matrix() @ q.as_matrix()``.

        Parameters
        ----------
        other : `Rotation` instance
            Object containing the rotations to be composed with this one. Note
            that rotation compositions are not commutative, so ``p * q`` is
            generally different from ``q * p``.

        Returns
        -------
        composition : `Rotation` instance
            This function supports composition of multiple rotations at a time.
            Composition follows standard numpy broadcasting rules. The resulting
            `Rotation` object will have the shape
            `np.broadcast_shapes(p.shape, q.shape)`. In dimensions with size > 1,
            rotations are composed with matching indices. In dimensions with only
            one rotation, the single rotation is composed with each rotation in the
            other object.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Composition of two single rotations:

        >>> p = R.from_quat([0, 0, 1, 1])
        >>> q = R.from_quat([1, 0, 0, 1])
        >>> p.as_matrix()
        array([[ 0., -1.,  0.],
               [ 1.,  0.,  0.],
               [ 0.,  0.,  1.]])
        >>> q.as_matrix()
        array([[ 1.,  0.,  0.],
               [ 0.,  0., -1.],
               [ 0.,  1.,  0.]])
        >>> r = p * q
        >>> r.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])

        Composition of two objects containing equal number of rotations:

        >>> p = R.from_quat([[0, 0, 1, 1], [1, 0, 0, 1]])
        >>> q = R.from_rotvec([[np.pi/4, 0, 0], [-np.pi/4, 0, np.pi/4]])
        >>> p.as_quat()
        array([[0.        , 0.        , 0.70710678, 0.70710678],
               [0.70710678, 0.        , 0.        , 0.70710678]])
        >>> q.as_quat()
        array([[ 0.38268343,  0.        ,  0.        ,  0.92387953],
               [-0.37282173,  0.        ,  0.37282173,  0.84971049]])
        >>> r = p * q
        >>> r.as_quat()
        array([[ 0.27059805,  0.27059805,  0.65328148,  0.65328148],
               [ 0.33721128, -0.26362477,  0.26362477,  0.86446082]])

        Broadcasting rules apply:
        >>> p = R.from_quat(np.tile(np.array([0, 0, 1, 1]), (5, 1, 1)))
        >>> q = R.from_quat(np.tile(np.array([1, 0, 0, 1]), (1, 6, 1)))
        >>> p.shape, q.shape
        ((5, 1), (1, 6))
        >>> r = p * q
        >>> r.shape
        (5, 6)
        """
        # Check that other is a Rotation object. We want to return NotImplemented
        # instead of raising an error to allow other types to implement __rmul__.
        # Python will then automatically try to delegate the multiplication to the
        # other type.
        # See https://github.com/scipy/scipy/issues/21541
        if not isinstance(other, Rotation):
            return NotImplemented
        if not broadcastable(self._quat.shape, other._quat.shape):
            raise ValueError(
                f"Cannot broadcast {self._quat.shape[:-1]} rotations in "
                f"first to {other._quat.shape[:-1]} rotations in second object."
            )
        cython_compatible = self._quat.ndim < 3 and other._quat.ndim < 3
        backend = select_backend(self._xp, cython_compatible=cython_compatible)
        quat = backend.compose_quat(self._quat, other._quat)
        if self._single and other._single:
            quat = quat[0]
        return Rotation(quat, normalize=True, copy=False)

    @xp_capabilities(
        skip_backends=[("dask.array", "cannot handle zero-length rotations")]
    )
    def __pow__(self, n: float | Array, modulus: None = None) -> Rotation:
        """Compose this rotation with itself `n` times.

        Composition of a rotation ``p`` with itself can be extended to
        non-integer ``n`` by considering the power ``n`` to be a scale factor
        applied to the angle of rotation about the rotation's fixed axis. The
        expression ``q = p ** n`` can also be expressed as
        ``q = Rotation.from_rotvec(n * p.as_rotvec())``.

        If ``n`` is negative, then the rotation is inverted before the power
        is applied. In other words, ``p ** -abs(n) == p.inv() ** abs(n)``.

        Parameters
        ----------
        n : float | Array
            The number of times to compose the rotation with itself. If `n` is
            an array, then it must be 0d or 1d with shape (1,).
        modulus : None
            This overridden argument is not applicable to Rotations and must be
            ``None``.

        Returns
        -------
        power : `Rotation` instance
            The resulting rotation will be of the same shape as the original rotation
            object. Each element of the output is the corresponding element of the
            input rotation raised to the power of ``n``.

        Notes
        -----
        For example, a power of 2 will double the angle of rotation, and a
        power of 0.5 will halve the angle. There are three notable cases: if
        ``n == 1`` then the original rotation is returned, if ``n == 0``
        then the identity rotation is returned, and if ``n == -1`` then
        ``p.inv()`` is returned.

        Note that fractional powers ``n`` which effectively take a root of
        rotation, do so using the shortest path smallest representation of that
        angle (the principal root). This means that powers of ``n`` and ``1/n``
        are not necessarily inverses of each other. For example, a 0.5 power of
        a +240 degree rotation will be calculated as the 0.5 power of a -120
        degree rotation, with the result being a rotation of -60 rather than
        +120 degrees.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Raising a rotation to a power:

        >>> p = R.from_rotvec([1, 0, 0])
        >>> q = p ** 2
        >>> q.as_rotvec()
        array([2., 0., 0.])
        >>> r = p ** 0.5
        >>> r.as_rotvec()
        array([0.5, 0., 0.])

        Inverse powers do not necessarily cancel out:

        >>> p = R.from_rotvec([0, 0, 120], degrees=True)
        >>> ((p ** 2) ** 0.5).as_rotvec(degrees=True)
        array([  -0.,   -0., -60.])

        """
        if modulus is not None:
            raise NotImplementedError("modulus not supported")
        quat = self._backend.pow(self._quat, n)
        if self._single:
            quat = quat[0]
        return Rotation._from_raw_quat(quat, xp=self._xp, backend=self._backend)

    @xp_capabilities(
        skip_backends=[("dask.array", "cannot handle zero-length rotations")]
    )
    def inv(self) -> Rotation:
        """Invert this rotation.

        Composition of a rotation with its inverse results in an identity
        transformation.

        Returns
        -------
        inverse : `Rotation` instance
            Object containing inverse of the rotations in the current instance.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Inverting a single rotation:

        >>> p = R.from_euler('z', 45, degrees=True)
        >>> q = p.inv()
        >>> q.as_euler('zyx', degrees=True)
        array([-45.,   0.,   0.])

        Inverting multiple rotations:

        >>> p = R.from_rotvec([[0, 0, np.pi/3], [-np.pi/4, 0, 0]])
        >>> q = p.inv()
        >>> q.as_rotvec()
        array([[-0.        , -0.        , -1.04719755],
               [ 0.78539816, -0.        , -0.        ]])

        """
        q_inv = self._backend.inv(self._quat)
        if self._single:
            q_inv = q_inv[0, ...]
        return Rotation._from_raw_quat(q_inv, xp=self._xp, backend=self._backend)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def magnitude(self) -> Array:
        """Get the magnitude(s) of the rotation(s).

        Returns
        -------
        magnitude : ndarray or float
            Angle(s) in radians, float if object contains a single rotation
            and ndarray if object contains ND rotations. The magnitude
            will always be in the range [0, pi].

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np
        >>> r = R.from_quat(np.eye(4))
        >>> r.as_quat()
        array([[ 1., 0., 0., 0.],
               [ 0., 1., 0., 0.],
               [ 0., 0., 1., 0.],
               [ 0., 0., 0., 1.]])
        >>> r.magnitude()
        array([3.14159265, 3.14159265, 3.14159265, 0.        ])

        Magnitude of a single rotation:

        >>> r[0].magnitude()
        3.141592653589793
        """
        magnitude = self._backend.magnitude(self._quat)
        if self._single:
            # Special handling for numpy and single rotations. self._single is only set
            # if xp is numpy. We therefore know that magnitude is a numpy array and
            # that it contains a single element. Previously this code returned a Python
            # float in that case. Here we return a numpy float64 scalar. All other
            # Array API libraries return 0d arrays instead.
            # See https://github.com/scipy/scipy/pull/23198#issuecomment-3003757848
            return magnitude[0]
        return magnitude

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def approx_equal(
        self, other: Rotation, atol: float | None = None, degrees: bool = False
    ) -> Array:
        """Determine if another rotation is approximately equal to this one.

        Equality is measured by calculating the smallest angle between the
        rotations, and checking to see if it is smaller than `atol`.

        Parameters
        ----------
        other : `Rotation` instance
            Object containing the rotations to measure against this one.
        atol : float, optional
            The absolute angular tolerance, below which the rotations are
            considered equal. If not given, then set to 1e-8 radians by
            default.
        degrees : bool, optional
            If True and `atol` is given, then `atol` is measured in degrees. If
            False (default), then atol is measured in radians.

        Returns
        -------
        approx_equal : ndarray or bool
            Whether the rotations are approximately equal, bool if object
            contains a single rotation and ndarray if object contains multiple
            rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np
        >>> p = R.from_quat([0, 0, 0, 1])
        >>> q = R.from_quat(np.eye(4))
        >>> p.approx_equal(q)
        array([False, False, False, True])

        Approximate equality for a single rotation:

        >>> p.approx_equal(q[0])
        False
        """
        cython_compatible = self._quat.ndim < 3 and other._quat.ndim < 3
        backend = select_backend(self._xp, cython_compatible=cython_compatible)
        return backend.approx_equal(self._quat, other._quat, atol=atol, degrees=degrees)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def mean(
        self,
        weights: ArrayLike | None = None,
        axis: None | int | tuple[int, ...] = None,
    ) -> Rotation:
        """Get the mean of the rotations.

        The mean used is the chordal L2 mean (also called the projected or
        induced arithmetic mean) [1]_. If ``A`` is a set of rotation matrices,
        then the mean ``M`` is the rotation matrix that minimizes the
        following loss function:

        .. math::

            L(M) = \\sum_{i = 1}^{n} w_i \\lVert \\mathbf{A}_i -
            \\mathbf{M} \\rVert^2 ,

        where :math:`w_i`'s are the `weights` corresponding to each matrix.

        Parameters
        ----------
        weights : array_like shape (..., N), optional
            Weights describing the relative importance of the rotations. If
            None (default), then all values in `weights` are assumed to be
            equal. If given, the shape of `weights` must be broadcastable to
            the rotation shape. Weights must be non-negative.
        axis : None, int, or tuple of ints, optional
            Axis or axes along which the means are computed. The default is to
            compute the mean of all rotations.

        Returns
        -------
        mean : `Rotation` instance
            Single rotation containing the mean of the rotations in the
            current instance.

        References
        ----------
        .. [1] Hartley, Richard, et al.,
                "Rotation Averaging", International Journal of Computer Vision
                103, 2013, pp. 267-305.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> r = R.from_euler('zyx', [[0, 0, 0],
        ...                          [1, 0, 0],
        ...                          [0, 1, 0],
        ...                          [0, 0, 1]], degrees=True)
        >>> r.mean().as_euler('zyx', degrees=True)
        array([0.24945696, 0.25054542, 0.24945696])
        """
        mean = self._backend.mean(self._quat, weights=weights, axis=axis)
        return Rotation._from_raw_quat(mean, xp=self._xp, backend=self._backend)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def reduce(
        self,
        left: Rotation | None = None,
        right: Rotation | None = None,
        return_indices: bool = False,
    ) -> Rotation | tuple[Rotation, Array, Array]:
        """Reduce this rotation with the provided rotation groups.

        Reduction of a rotation ``p`` is a transformation of the form
        ``q = l * p * r``, where ``l`` and ``r`` are chosen from `left` and
        `right` respectively, such that rotation ``q`` has the smallest
        magnitude.

        If `left` and `right` are rotation groups representing symmetries of
        two objects rotated by ``p``, then ``q`` is the rotation of the
        smallest magnitude to align these objects considering their symmetries.

        Parameters
        ----------
        left : `Rotation` instance, optional
            Object containing the left rotation(s). Default value (None)
            corresponds to the identity rotation.
        right : `Rotation` instance, optional
            Object containing the right rotation(s). Default value (None)
            corresponds to the identity rotation.
        return_indices : bool, optional
            Whether to return the indices of the rotations from `left` and
            `right` used for reduction.

        Returns
        -------
        reduced : `Rotation` instance
            Object containing reduced rotations.
        left_best, right_best: integer ndarray
            Indices of elements from `left` and `right` used for reduction.
        """
        left = left.as_quat() if left is not None else None
        right = right.as_quat() if right is not None else None
        reduced, left_idx, right_idx = self._backend.reduce(
            self._quat, left=left, right=right
        )
        if self._single:
            reduced = reduced[0, ...]
        rot = Rotation._from_raw_quat(reduced, xp=self._xp, backend=self._backend)
        if return_indices:
            left_idx = left_idx if left is not None else None
            right_idx = right_idx if right is not None else None
            return rot, left_idx, right_idx
        return rot

    @classmethod
    def create_group(cls, group: str, axis: str = "Z") -> Rotation:
        """Create a 3D rotation group.

        Parameters
        ----------
        group : string
            The name of the group. Must be one of 'I', 'O', 'T', 'Dn', 'Cn',
            where `n` is a positive integer. The groups are:

                * I: Icosahedral group
                * O: Octahedral group
                * T: Tetrahedral group
                * D: Dicyclic group
                * C: Cyclic group

        axis : integer
            The cyclic rotation axis. Must be one of ['X', 'Y', 'Z'] (or
            lowercase). Default is 'Z'. Ignored for groups 'I', 'O', and 'T'.

        Returns
        -------
        rotation : `Rotation` instance
            Object containing the elements of the rotation group.

        Notes
        -----
        This method generates rotation groups only. The full 3-dimensional
        point groups [PointGroups]_ also contain reflections.

        References
        ----------
        .. [PointGroups] `Point groups
           <https://en.wikipedia.org/wiki/Point_groups_in_three_dimensions>`_
           on Wikipedia.
        """
        # TODO: We defer the implementation of groups for arbitrary Array API frameworks
        # to the follow-up PR that adds general Array API support for Rotations.
        return create_group(cls, group, axis=axis)

    @xp_capabilities(
        jax_jit=False,
        skip_backends=[("dask.array", "cannot handle zero-length rotations")],
    )
    def __getitem__(self, indexer: int | slice | EllipsisType | None) -> Rotation:
        """Extract rotation(s) at given index(es) from object.

        Create a new `Rotation` instance containing a subset of rotations
        stored in this object.

        Parameters
        ----------
        indexer : index, slice, or index array
            Specifies which rotation(s) to extract. A single indexer must be
            specified, i.e. as if indexing a 1 dimensional array or list.

        Returns
        -------
        rotation : `Rotation` instance
            Contains
                - a single rotation, if `indexer` is a single index
                - a stack of rotation(s), if `indexer` is a slice, or and index
                  array.

        Raises
        ------
        TypeError if the instance was created as a single rotation.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R
        >>> rs = R.from_quat([
        ... [1, 1, 0, 0],
        ... [0, 1, 0, 1],
        ... [1, 1, -1, 0]])  # These quats are normalized
        >>> rs.as_quat()
        array([[ 0.70710678,  0.70710678,  0.        ,  0.        ],
               [ 0.        ,  0.70710678,  0.        ,  0.70710678],
               [ 0.57735027,  0.57735027, -0.57735027,  0.        ]])

        Indexing using a single index:

        >>> a = rs[0]
        >>> a.as_quat()
        array([0.70710678, 0.70710678, 0.        , 0.        ])

        Array slicing:

        >>> b = rs[1:3]
        >>> b.as_quat()
        array([[ 0.        ,  0.70710678,  0.        ,  0.70710678],
               [ 0.57735027,  0.57735027, -0.57735027,  0.        ]])

        List comprehension to split each rotation into its own object:

        >>> c = [r for r in rs]
        >>> print([r.as_quat() for r in c])
        [array([ 0.70710678,  0.70710678,  0.        ,  0.        ]),
         array([ 0.        ,  0.70710678,  0.        ,  0.70710678]),
         array([ 0.57735027,  0.57735027, -0.57735027,  0.        ])]

        Concatenation of split rotations will recover the original object:

        >>> R.concatenate([a, b]).as_quat()
        array([[ 0.70710678,  0.70710678,  0.        ,  0.        ],
               [ 0.        ,  0.70710678,  0.        ,  0.70710678],
               [ 0.57735027,  0.57735027, -0.57735027,  0.        ]])
        """
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation is not subscriptable.")
        is_array = isinstance(indexer, type(self._quat))
        # Masking is only specified in the Array API when the array is the sole index
        # TODO: This special case handling is mainly a result of Array API limitations.
        # Ideally we would get rid of them altogether and converge to [indexer, ...]
        # indexing.
        if is_array and indexer.dtype == self._xp.bool:
            return Rotation(self._quat[indexer], normalize=False)
        if is_array and self._xp.isdtype(indexer.dtype, "integral"):
            # xp.take is implementation-defined for zero-dim arrays, hence we raise
            # pre-emptively to have consistent behavior across frameworks.
            if self._quat.shape[0] == 0:
                raise IndexError("cannot do a non-empty take from an empty axes.")
            return Rotation(self._xp.take(self._quat, indexer, axis=0), normalize=False)
        return Rotation(self._quat[indexer, ...], normalize=False)

    @xp_capabilities(
        jax_jit=False,
        skip_backends=[("dask.array", "cannot handle zero-length rotations")],
    )
    def __setitem__(self, indexer: int | slice | EllipsisType | None, value: Rotation):
        """Set rotation(s) at given index(es) from object.

        Parameters
        ----------
        indexer : index, slice, or index array
            Specifies which rotation(s) to replace. A single indexer must be
            specified, i.e. as if indexing a 1 dimensional array or list.

        value : `Rotation` instance
            The rotations to set.

        Raises
        ------
        TypeError if the instance was created as a single rotation.

        Notes
        -----

        .. versionadded:: 1.8.0
        """
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation is not subscriptable.")

        if not isinstance(value, Rotation):
            raise TypeError("value must be a Rotation object")

        self._quat = self._backend.setitem(self._quat, value.as_quat(), indexer)

    @staticmethod
    def identity(
        num: int | None = None, *, shape: int | tuple[int, ...] | None = None
    ) -> Rotation:
        """Get identity rotation(s).

        Composition with the identity rotation has no effect.

        Parameters
        ----------
        num : int or None, optional
            Number of identity rotations to generate. If None (default), then a
            single rotation is generated.
        shape : int or tuple of ints, optional
            Shape of identity rotations to generate. If specified, `num` must
            be None.

        Returns
        -------
        identity : Rotation object
            The identity rotation.
        """
        # TODO: We should move to one single way of specifying the output shape and
        # deprecate `num`.
        if num is not None and shape is not None:
            raise ValueError("Only one of `num` or `shape` can be specified.")
        quat = cython_backend.identity(num, shape=shape)
        return Rotation._from_raw_quat(quat, xp=array_namespace(quat))

    @staticmethod
    @_transition_to_rng("random_state", position_num=2)
    def random(
        num: int | None = None,
        rng: np.random.Generator | None = None,
        *,
        shape: tuple[int, ...] | None = None,
    ) -> Rotation:
        r"""Generate rotations that are uniformly distributed on a sphere.

        Formally, the rotations follow the Haar-uniform distribution over the SO(3)
        group.

        Parameters
        ----------
        num : int or None, optional
            Number of random rotations to generate. If None (default), then a
            single rotation is generated.
        rng : `numpy.random.Generator`, optional
            Pseudorandom number generator state. When `rng` is None, a new
            `numpy.random.Generator` is created using entropy from the
            operating system. Types other than `numpy.random.Generator` are
            passed to `numpy.random.default_rng` to instantiate a `Generator`.
        shape : tuple of ints, optional
            Shape of random rotations to generate. If specified, `num` must be None.

        Returns
        -------
        random_rotation : `Rotation` instance
            Contains a single rotation if `num` is None. Otherwise contains a
            stack of `num` rotations.

        Notes
        -----
        This function is optimized for efficiently sampling random rotation
        matrices in three dimensions. For generating random rotation matrices
        in higher dimensions, see `scipy.stats.special_ortho_group`.

        Examples
        --------
        >>> from scipy.spatial.transform import Rotation as R

        Sample a single rotation:

        >>> R.random().as_euler('zxy', degrees=True)
        array([-110.5976185 ,   55.32758512,   76.3289269 ])  # random

        Sample a stack of rotations:

        >>> R.random(5).as_euler('zxy', degrees=True)
        array([[-110.5976185 ,   55.32758512,   76.3289269 ],  # random
               [ -91.59132005,  -14.3629884 ,  -93.91933182],
               [  25.23835501,   45.02035145, -121.67867086],
               [ -51.51414184,  -15.29022692, -172.46870023],
               [ -81.63376847,  -27.39521579,    2.60408416]])

        See Also
        --------
        scipy.stats.special_ortho_group

        """
        # TODO: We should move to one single way of specifying the output shape and
        # deprecate `num`.
        if num is not None and shape is not None:
            raise ValueError("Only one of `num` or `shape` can be specified.")
        sample = cython_backend.random(num, rng, shape=shape)
        return Rotation(sample, normalize=True, copy=False)

    @staticmethod
    @xp_capabilities(
        skip_backends=[
            ("dask.array", "missing linalg.cross/det functions and .mT attribute"),
            ("cupy", "missing .mT attribute in cupy<14.*"),
        ]
    )
    def align_vectors(
        a: ArrayLike,
        b: ArrayLike,
        weights: ArrayLike | None = None,
        return_sensitivity: bool = False,
    ) -> tuple[Rotation, float] | tuple[Rotation, float, Array]:
        """Estimate a rotation to optimally align two sets of vectors.

        Find a rotation between frames A and B which best aligns a set of
        vectors `a` and `b` observed in these frames. The following loss
        function is minimized to solve for the rotation matrix
        :math:`C`:

        .. math::

            L(C) = \\frac{1}{2} \\sum_{i = 1}^{n} w_i \\lVert \\mathbf{a}_i -
            C \\mathbf{b}_i \\rVert^2 ,

        where :math:`w_i`'s are the `weights` corresponding to each vector.

        The rotation is estimated with Kabsch algorithm [1]_, and solves what
        is known as the "pointing problem", or "Wahba's problem" [2]_.

        Note that the length of each vector in this formulation acts as an
        implicit weight. So for use cases where all vectors need to be
        weighted equally, you should normalize them to unit length prior to
        calling this method.

        There are two special cases. The first is if a single vector is given
        for `a` and `b`, in which the shortest distance rotation that aligns
        `b` to `a` is returned.

        The second is when one of the weights is infinity. In this case, the
        shortest distance rotation between the primary infinite weight vectors
        is calculated as above. Then, the rotation about the aligned primary
        vectors is calculated such that the secondary vectors are optimally
        aligned per the above loss function. The result is the composition
        of these two rotations. The result via this process is the same as the
        Kabsch algorithm as the corresponding weight approaches infinity in
        the limit. For a single secondary vector this is known as the
        "align-constrain" algorithm [3]_.

        For both special cases (single vectors or an infinite weight), the
        sensitivity matrix does not have physical meaning and an error will be
        raised if it is requested. For an infinite weight, the primary vectors
        act as a constraint with perfect alignment, so their contribution to
        `rssd` will be forced to 0 even if they are of different lengths.

        Parameters
        ----------
        a : array_like, shape (3,) or (N, 3)
            Vector components observed in initial frame A. Each row of `a`
            denotes a vector.
        b : array_like, shape (3,) or (N, 3)
            Vector components observed in another frame B. Each row of `b`
            denotes a vector.
        weights : array_like shape (N,), optional
            Weights describing the relative importance of the vector
            observations. If None (default), then all values in `weights` are
            assumed to be 1. One and only one weight may be infinity, and
            weights must be positive.
        return_sensitivity : bool, optional
            Whether to return the sensitivity matrix. See Notes for details.
            Default is False.

        Returns
        -------
        rotation : `Rotation` instance
            Best estimate of the rotation that transforms `b` to `a`.
        rssd : float
            Stands for "root sum squared distance". Square root of the weighted
            sum of the squared distances between the given sets of vectors
            after alignment. It is equal to ``sqrt(2 * minimum_loss)``, where
            ``minimum_loss`` is the loss function evaluated for the found
            optimal rotation.
            Note that the result will also be weighted by the vectors'
            magnitudes, so perfectly aligned vector pairs will have nonzero
            `rssd` if they are not of the same length. This can be avoided by
            normalizing them to unit length prior to calling this method,
            though note that doing this will change the resulting rotation.
        sensitivity_matrix : ndarray, shape (3, 3)
            Sensitivity matrix of the estimated rotation estimate as explained
            in Notes. Returned only when `return_sensitivity` is True. Not
            valid if aligning a single pair of vectors or if there is an
            infinite weight, in which cases an error will be raised.

        Notes
        -----
        The sensitivity matrix gives the sensitivity of the estimated rotation
        to small perturbations of the vector measurements. Specifically we
        consider the rotation estimate error as a small rotation vector of
        frame A. The sensitivity matrix is proportional to the covariance of
        this rotation vector assuming that the vectors in `a` was measured with
        errors significantly less than their lengths. To get the true
        covariance matrix, the returned sensitivity matrix must be multiplied
        by harmonic mean [4]_ of variance in each observation. Note that
        `weights` are supposed to be inversely proportional to the observation
        variances to get consistent results. For example, if all vectors are
        measured with the same accuracy of 0.01 (`weights` must be all equal),
        then you should multiple the sensitivity matrix by 0.01**2 to get the
        covariance.

        Refer to [5]_ for more rigorous discussion of the covariance
        estimation. See [6]_ for more discussion of the pointing problem and
        minimal proper pointing.

        This function does not support broadcasting or ND arrays with N > 2.

        References
        ----------
        .. [1] https://en.wikipedia.org/wiki/Kabsch_algorithm
        .. [2] https://en.wikipedia.org/wiki/Wahba%27s_problem
        .. [3] Magner, Robert,
                "Extending target tracking capabilities through trajectory and
                momentum setpoint optimization." Small Satellite Conference,
                2018.
        .. [4] https://en.wikipedia.org/wiki/Harmonic_mean
        .. [5] F. Landis Markley,
                "Attitude determination using vector observations: a fast
                optimal matrix algorithm", Journal of Astronautical Sciences,
                Vol. 41, No.2, 1993, pp. 261-280.
        .. [6] Bar-Itzhack, Itzhack Y., Daniel Hershkowitz, and Leiba Rodman,
                "Pointing in Real Euclidean Space", Journal of Guidance,
                Control, and Dynamics, Vol. 20, No. 5, 1997, pp. 916-922.

        Examples
        --------
        >>> import numpy as np
        >>> from scipy.spatial.transform import Rotation as R

        Here we run the baseline Kabsch algorithm to best align two sets of
        vectors, where there is noise on the last two vector measurements of
        the ``b`` set:

        >>> a = [[0, 1, 0], [0, 1, 1], [0, 1, 1]]
        >>> b = [[1, 0, 0], [1, 1.1, 0], [1, 0.9, 0]]
        >>> rot, rssd, sens = R.align_vectors(a, b, return_sensitivity=True)
        >>> rot.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])

        When we apply the rotation to ``b``, we get vectors close to ``a``:

        >>> rot.apply(b)
        array([[0. , 1. , 0. ],
               [0. , 1. , 1.1],
               [0. , 1. , 0.9]])

        The error for the first vector is 0, and for the last two the error is
        magnitude 0.1. The `rssd` is the square root of the sum of the
        weighted squared errors, and the default weights are all 1, so in this
        case the `rssd` is calculated as
        ``sqrt(1 * 0**2 + 1 * 0.1**2 + 1 * (-0.1)**2) = 0.141421356237308``

        >>> a - rot.apply(b)
        array([[ 0., 0.,  0. ],
               [ 0., 0., -0.1],
               [ 0., 0.,  0.1]])
        >>> np.sqrt(np.sum(np.ones(3) @ (a - rot.apply(b))**2))
        0.141421356237308
        >>> rssd
        0.141421356237308

        The sensitivity matrix for this example is as follows:

        >>> sens
        array([[0.2, 0. , 0.],
               [0. , 1.5, 1.],
               [0. , 1. , 1.]])

        Special case 1: Find a minimum rotation between single vectors:

        >>> a = [1, 0, 0]
        >>> b = [0, 1, 0]
        >>> rot, _ = R.align_vectors(a, b)
        >>> rot.as_matrix()
        array([[0., 1., 0.],
               [-1., 0., 0.],
               [0., 0., 1.]])
        >>> rot.apply(b)
        array([1., 0., 0.])

        Special case 2: One infinite weight. Here we find a rotation between
        primary and secondary vectors that can align exactly:

        >>> a = [[0, 1, 0], [0, 1, 1]]
        >>> b = [[1, 0, 0], [1, 1, 0]]
        >>> rot, _ = R.align_vectors(a, b, weights=[np.inf, 1])
        >>> rot.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])
        >>> rot.apply(b)
        array([[0., 1., 0.],
               [0., 1., 1.]])

        Here the secondary vectors must be best-fit:

        >>> a = [[0, 1, 0], [0, 1, 1]]
        >>> b = [[1, 0, 0], [1, 2, 0]]
        >>> rot, _ = R.align_vectors(a, b, weights=[np.inf, 1])
        >>> rot.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])
        >>> rot.apply(b)
        array([[0., 1., 0.],
               [0., 1., 2.]])
        """
        xp = array_namespace(a)
        a, b, weights = _promote(a, b, weights, xp=xp)
        cython_compatible = (
            (a.ndim < 3) & (b.ndim < 3) & (weights is None or weights.ndim < 2)
        )
        backend = select_backend(xp, cython_compatible=cython_compatible)
        q, rssd, sensitivity = backend.align_vectors(a, b, weights, return_sensitivity)
        if return_sensitivity:
            return Rotation._from_raw_quat(q, xp=xp, backend=backend), rssd, sensitivity
        return Rotation._from_raw_quat(q, xp=xp, backend=backend), rssd

    def __getstate__(self) -> tuple[Array, bool]:
        return (self._quat, self._single)

    def __setstate__(self, state: tuple[Array, bool]):
        quat, single = state
        xp = array_namespace(quat)
        self._xp = xp
        self._quat = xp.asarray(quat, copy=True)
        self._backend = select_backend(xp, cython_compatible=self._quat.ndim < 3)
        self._single = single

    @property
    def single(self) -> bool:
        """Whether this instance represents a single rotation."""
        return self._single or self._quat.ndim == 1

    @property
    def shape(self) -> tuple[int, ...]:
        """The shape of the rotation's leading dimensions."""
        if self._single:
            return ()
        return self._quat.shape[:-1]

    def __bool__(self) -> bool:
        """Comply with Python convention for objects to be True.

        Required because `Rotation.__len__()` is defined and not always truthy.
        """
        return True

    def __len__(self) -> int:
        """Number of rotations contained in this object.

        Multiple rotations can be stored in a single instance.

        Returns
        -------
        length : int
            Number of rotations stored in object.

        Raises
        ------
        TypeError if the instance was created as a single rotation.
        """
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation has no len().")
        return self._quat.shape[0]

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def __repr__(self) -> str:
        m = f"{self.as_matrix()!r}".splitlines()
        # bump indent (+21 characters)
        m[1:] = [" " * 21 + m[i] for i in range(1, len(m))]
        return "Rotation.from_matrix(" + "\n".join(m) + ")"

    @xp_capabilities()
    def __iter__(self) -> Iterator[Rotation]:
        """Iterate over rotations."""
        if self._single or self._quat.ndim == 1:
            raise TypeError("Single rotation is not iterable.")
        # We return a generator that yields a new Rotation object for each rotation
        # in the current object. We cannot rely on the default implementation using
        # __getitem__ because jax will not raise an IndexError for out-of-bounds
        # indices.
        for i in range(self._quat.shape[0]):
            yield Rotation(self._quat[i, ...], normalize=False, copy=False)

    @staticmethod
    def _from_raw_quat(
        quat: Array, xp: ModuleType, backend: ModuleType | None = None
    ) -> Rotation:
        """Create a Rotation skipping all sanitization steps.

        This method is intended for internal, performant creation of Rotations with
        quaternions that are guaranteed to be valid.
        """
        rot = Rotation.__new__(Rotation)
        rot._single = quat.ndim == 1 and is_numpy(xp)
        if rot._single:
            quat = xpx.atleast_nd(quat, ndim=2, xp=xp)
        rot._quat = quat
        rot._xp = xp
        if backend is None:
            backend = select_backend(xp, cython_compatible=quat.ndim < 3)
        rot._backend = backend
        return rot


class Slerp:
    """Spherical Linear Interpolation of Rotations.

    The interpolation between consecutive rotations is performed as a rotation
    around a fixed axis with a constant angular velocity [1]_. This ensures
    that the interpolated rotations follow the shortest path between initial
    and final orientations.

    Parameters
    ----------
    times : array_like, shape (N,)
        Times of the known rotations. At least 2 times must be specified.
    rotations : `Rotation` instance
        Rotations to perform the interpolation between. Must contain N
        rotations.

    Methods
    -------
    __call__

    See Also
    --------
    Rotation

    Notes
    -----
    This class only supports interpolation of rotations with a single leading
    dimension.

    .. versionadded:: 1.2.0

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Slerp#Quaternion_Slerp

    Examples
    --------
    >>> from scipy.spatial.transform import Rotation as R
    >>> from scipy.spatial.transform import Slerp

    Setup the fixed keyframe rotations and times:

    >>> key_rots = R.random(5, random_state=2342345)
    >>> key_times = [0, 1, 2, 3, 4]

    Create the interpolator object:

    >>> slerp = Slerp(key_times, key_rots)

    Interpolate the rotations at the given times:

    >>> times = [0, 0.5, 0.25, 1, 1.5, 2, 2.75, 3, 3.25, 3.60, 4]
    >>> interp_rots = slerp(times)

    The keyframe rotations expressed as Euler angles:

    >>> key_rots.as_euler('xyz', degrees=True)
    array([[ 14.31443779, -27.50095894,  -3.7275787 ],
           [ -1.79924227, -24.69421529, 164.57701743],
           [146.15020772,  43.22849451, -31.34891088],
           [ 46.39959442,  11.62126073, -45.99719267],
           [-88.94647804, -49.64400082, -65.80546984]])

    The interpolated rotations expressed as Euler angles. These agree with the
    keyframe rotations at both endpoints of the range of keyframe times.

    >>> interp_rots.as_euler('xyz', degrees=True)
    array([[  14.31443779,  -27.50095894,   -3.7275787 ],
           [   4.74588574,  -32.44683966,   81.25139984],
           [  10.71094749,  -31.56690154,   38.06896408],
           [  -1.79924227,  -24.69421529,  164.57701743],
           [  11.72796022,   51.64207311, -171.7374683 ],
           [ 146.15020772,   43.22849451,  -31.34891088],
           [  68.10921869,   20.67625074,  -48.74886034],
           [  46.39959442,   11.62126073,  -45.99719267],
           [  12.35552615,    4.21525086,  -64.89288124],
           [ -30.08117143,  -19.90769513,  -78.98121326],
           [ -88.94647804,  -49.64400082,  -65.80546984]])

    """

    @xp_capabilities(
        jax_jit=False, skip_backends=[("dask.array", "missing linalg.cross function")]
    )
    def __init__(self, times: ArrayLike, rotations: Rotation):
        if not isinstance(rotations, Rotation):
            raise TypeError("`rotations` must be a `Rotation` instance.")

        if rotations.single or len(rotations) <= 1:
            raise ValueError("`rotations` must be a sequence of at least 2 rotations.")

        q = rotations.as_quat()
        if q.ndim > 2:
            raise ValueError(
                "Rotations with more than 1 leading dimension are not supported."
            )

        xp = array_namespace(q)
        times = xp.asarray(times, device=xp_device(q), dtype=q.dtype)
        if times.ndim != 1:
            raise ValueError(
                "Expected times to be specified in a 1 dimensional array, got "
                f"{times.ndim} dimensions."
            )

        if times.shape[0] != len(rotations):
            raise ValueError(
                "Expected number of rotations to be equal to number of timestamps "
                f"given, got {len(rotations)} rotations and {times.shape[0]} "
                "timestamps."
            )
        self.times = times
        self.timedelta = xp.diff(times)

        # We cannot check for values for lazy backends, so we cannot raise an error on
        # timedelta < 0 for lazy backends. Instead, we set timedelta to nans
        neg_mask = xp.any(self.timedelta <= 0)
        if is_lazy_array(neg_mask):
            self.timedelta = xp.where(neg_mask, xp.nan, self.timedelta)
            self.times = xp.where(neg_mask, xp.nan, self.times)
        elif xp.any(neg_mask):
            raise ValueError("Times must be in strictly increasing order.")

        self.rotations = rotations[:-1]
        self.rotvecs = (self.rotations.inv() * rotations[1:]).as_rotvec()

    @xp_capabilities()
    def __call__(self, times: ArrayLike) -> Rotation:
        """Interpolate rotations.

        Compute the interpolated rotations at the given `times`.

        Parameters
        ----------
        times : array_like
            Times to compute the interpolations at. Can be a scalar or
            1-dimensional.

        Returns
        -------
        interpolated_rotation : `Rotation` instance
            Object containing the rotations computed at given `times`.

        """
        xp = array_namespace(self.times)
        device = xp_device(self.times)
        # Clearly differentiate from self.times property
        compute_times = xp.asarray(times, device=device, dtype=self.times.dtype)
        if compute_times.ndim > 1:
            raise ValueError("`times` must be at most 1-dimensional.")

        single_time = compute_times.ndim == 0
        compute_times = xpx.atleast_nd(compute_times, ndim=1, xp=xp)

        # side = 'left' (default) excludes t_min.
        ind = xp.searchsorted(self.times, compute_times) - 1
        # Include t_min. Without this step, index for t_min equals -1
        ind = xpx.at(ind, compute_times == self.times[0]).set(0)
        # We cannot error out on invalid indices for jit compiled code. To not produce
        # an index error, we set the index to 0 in case it is out of bounds, and later
        # set the result to nan.
        invalid_ind = (ind < 0) | (ind > len(self.rotations) - 1)
        if is_lazy_array(invalid_ind):
            ind = xpx.at(ind, invalid_ind).set(0)
        elif xp.any(invalid_ind):
            raise ValueError(
                f"Interpolation times must be within the range [{self.times[0]}, "
                f"{self.times[-1]}], both inclusive."
            )

        alpha = (compute_times - self.times[ind]) / self.timedelta[ind]
        alpha = xpx.at(alpha, invalid_ind).set(xp.nan)

        # The array API does not support integer arrays + ellipsis indexing and won't
        # stabilize this feature due to blockers in PyTorch. Therefore we need to
        # construct the index for the last dimension manually.
        # See https://github.com/scipy/scipy/pull/23249#discussion_r2198363047
        result = self.rotations[ind] * Rotation.from_rotvec(
            self.rotvecs[ind[:, None], xp.arange(3, device=device)] * alpha[:, None]
        )

        if single_time:
            result = result[0]

        return result
