from __future__ import annotations

from collections.abc import Iterable, Iterator
from types import EllipsisType, ModuleType, NotImplementedType
from collections.abc import Callable

import numpy as np

from scipy._lib._array_api import (
    array_namespace,
    is_numpy,
    ArrayLike,
    Array,
    xp_capabilities,
)
from scipy.spatial.transform import Rotation
from scipy.spatial.transform._rotation import _promote
import scipy.spatial.transform._rigid_transform_cy as cython_backend
import scipy.spatial.transform._rigid_transform_xp as xp_backend
import scipy._lib.array_api_extra as xpx
from scipy._lib.array_api_compat import device
from scipy._lib._array_api import xp_promote
from scipy._lib._util import broadcastable


__all__ = ["RigidTransform"]

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
def normalize_dual_quaternion(dual_quat: ArrayLike) -> Array:
    """Normalize dual quaternion."""
    xp = array_namespace(dual_quat)
    dual_quat = _promote(dual_quat, xp=xp)
    single = dual_quat.ndim == 1 and is_numpy(xp)
    if single:
        dual_quat = xpx.atleast_nd(dual_quat, ndim=2, xp=xp)
    cython_compatible = dual_quat.ndim < 3
    dq = select_backend(xp, cython_compatible).normalize_dual_quaternion(dual_quat)
    if single:
        return dq[0]
    return dq


class RigidTransform:
    """Rigid transform in 3 dimensions.

    This class provides an interface to initialize from and represent rigid
    transforms (rotation and translation) in 3D space. In different fields,
    this type of transform may be referred to as "*pose*" (especially in
    robotics), "*extrinsic parameters*", or the "*model matrix*" (especially in
    computer graphics), but the core concept is the same: a rotation and
    translation describing the orientation of one 3D coordinate frame relative
    to another. Mathematically, these transforms belong to the Special
    Euclidean group SE(3), which encodes rotation (SO(3)) plus translation.

    The following operations on rigid transforms are supported:

    - Application on vectors
    - Transformation composition
    - Transformation inversion
    - Transformation indexing

    Note that coordinate systems must be right-handed. Because of this, this
    class more precisely represents *proper* rigid transformations in SE(3)
    rather than rigid transforms in E(3) more generally [1]_.

    Indexing within a transform is supported since multiple transforms can be
    stored within a single `RigidTransform` instance.

    Multiple transforms can be stored in a single `RigidTransform` object, which can be
    initialized using N-dimensional arrays and supports broadcasting for all
    operations.

    To create `RigidTransform` objects use ``from_...`` methods (see examples
    below). ``RigidTransform(...)`` is not supposed to be instantiated directly.

    For rigorous introductions to rigid transforms, see [2]_, [3]_, and [4]_.

    Attributes
    ----------
    single
    rotation
    translation

    Methods
    -------
    __len__
    __getitem__
    __mul__
    __pow__
    from_matrix
    from_rotation
    from_translation
    from_components
    from_exp_coords
    from_dual_quat
    as_matrix
    as_components
    as_exp_coords
    as_dual_quat
    concatenate
    mean
    apply
    inv
    identity

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Rigid_transformation
    .. [2] https://motion.cs.illinois.edu/RoboticSystems/CoordinateTransformations.html
    .. [3] https://www.brainvoyager.com/bv/doc/UsersGuide/CoordsAndTransforms/SpatialTransformationMatrices.html
    .. [4] Kevin M. Lynch and Frank C. Park, "Modern Robotics: Mechanics,
           Planning, and Control" Chapter 3.3, 2017, Cambridge University Press.
           https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf#page=107.31
    .. [5] Paul Furgale, "Representing Robot Pose: The good, the bad, and the
           ugly", June 9, 2014.
           https://rpg.ifi.uzh.ch/docs/teaching/2024/FurgaleTutorial.pdf

    Notes
    -----
    .. versionadded:: 1.16.0

    Examples
    --------
    A `RigidTransform` instance can be initialized in any of the above formats
    and converted to any of the others. The underlying object is independent of
    the representation used for initialization.

    **Notation Conventions and Composition**

    The notation here largely follows the convention defined in [5]_. When we
    name transforms, we read the subscripts from right to left. So ``tf_A_B``
    represents a transform A <- B and can be interpreted as:

    - the coordinates and orientation of B relative to A
    - the transformation of points from B to A
    - the pose of B described in A's coordinate system

    .. parsed-literal::
        :class: highlight-none

        tf_A_B
           ^ ^
           | |
           | --- from B
           |
           ----- to A

    When composing transforms, the order is important. Transforms are not
    commutative, so in general ``tf_A_B * tf_B_C`` is not the same as
    ``tf_B_C * tf_A_B``. Transforms are composed and applied to vectors
    right-to-left. So ``(tf_A_B * tf_B_C).apply(p_C)`` is the same as
    ``tf_A_B.apply(tf_B_C.apply(p_C))``.

    When composed, transforms should be ordered such that the multiplication
    operator is surrounded by a single frame, so the frame "cancels out" and
    the outside frames are left. In the example below, B cancels out and the
    outside frames A and C are left. Or to put it another way, A <- C is the
    same as A <- B <- C.

    .. parsed-literal::
        :class: highlight-none

                      ----------- B cancels out
                      |      |
                      v      v
        tf_A_C = tf_A_B * tf_B_C
                    ^          ^
                    |          |
                    ------------ to A, from C are left

    When we notate vectors, we write the subscript of the frame that the vector
    is defined in. So ``p_B`` means the point ``p`` defined in frame B. To
    transform this point from frame B to coordinates in frame A, we apply the
    transform ``tf_A_B`` to the vector, lining things up such that the notated
    B frames are next to each other and "cancel out".

    .. parsed-literal::
        :class: highlight-none

                   ------------ B cancels out
                   |         |
                   v         v
        p_A = tf_A_B.apply(p_B)
                 ^
                 |
                 -------------- A is left

    **Visualization**

    >>> from scipy.spatial.transform import RigidTransform as Tf
    >>> from scipy.spatial.transform import Rotation as R
    >>> import numpy as np

    The following function can be used to plot transforms with Matplotlib
    by showing how they transform the standard x, y, z coordinate axes:

    >>> import matplotlib.pyplot as plt
    >>> colors = ("#FF6666", "#005533", "#1199EE")  # Colorblind-safe RGB
    >>> def plot_transformed_axes(ax, tf, name=None, scale=1):
    ...     r = tf.rotation
    ...     t = tf.translation
    ...     loc = np.array([t, t])
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
    ...         text_plot = text_loc_rot + t
    ...         ax.text(*text_plot, axlabel.upper(), color=c,
    ...                 va="center", ha="center")
    ...     ax.text(*tf.translation, name, color="k", va="center", ha="center",
    ...             bbox={"fc": "w", "alpha": 0.8, "boxstyle": "circle"})

    **Defining Frames**

    Let's work through an example.

    First, define the "world frame" A, also called the "base frame".
    All frames are the identity transform from their own perspective.

    >>> tf_A = Tf.identity()

    We will visualize a new frame B in A's coordinate system. So we need to
    define the transform that converts coordinates from frame B to frame A
    (A <- B).

    Physically, let's imagine constructing B from A by:

    1) Rotating A by +90 degrees around its x-axis.
    2) Translating the rotated frame 2 units in A's -x direction.

    From A's perspective, B is at [-2, 0, 0] and rotated +90 degrees about the
    x-axis, which is exactly the transform A <- B.

    >>> t_A_B = np.array([-2, 0, 0])
    >>> r_A_B = R.from_euler('xyz', [90, 0, 0], degrees=True)
    >>> tf_A_B = Tf.from_components(t_A_B, r_A_B)

    Let's plot these frames.

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_A, name="tfA")     # A plotted in A
    >>> plot_transformed_axes(ax, tf_A_B, name="tfAB")  # B plotted in A
    >>> ax.set_title("A, B frames with respect to A")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()

    Now let's visualize a new frame C in B's coordinate system.
    Let's imagine constructing C from B by:

    1) Translating B by 2 units in its +z direction.
    2) Rotating B by +30 degrees around its z-axis.

    >>> t_B_C = np.array([0, 0, 2])
    >>> r_B_C = R.from_euler('xyz', [0, 0, 30], degrees=True)
    >>> tf_B_C = Tf.from_components(t_B_C, r_B_C)

    In order to plot these frames from a consistent perspective, we need to
    calculate the transform between A and C. Note that we do not make this
    transform directly, but instead compose intermediate transforms that let us
    get from C to A:

    >>> tf_A_C = tf_A_B * tf_B_C  # A <- B <- C

    Now we can plot these three frames from A's perspective.

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_A, name="tfA")     # A plotted in A
    >>> plot_transformed_axes(ax, tf_A_B, name="tfAB")  # B plotted in A
    >>> plot_transformed_axes(ax, tf_A_C, name="tfAC")  # C plotted in A
    >>> ax.set_title("A, B, C frames with respect to A")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()

    **Transforming Vectors**

    Let's transform a vector from A, to B and C. To do this, we will first
    invert the transforms we already have from B and C, to A.

    >>> tf_B_A = tf_A_B.inv()  # B <- A
    >>> tf_C_A = tf_A_C.inv()  # C <- A

    Now we can define a point in A and use the above transforms to get its
    coordinates in B and C:

    >>> p1_A = np.array([1, 0, 0])  # +1 in x_A direction
    >>> p1_B = tf_B_A.apply(p1_A)
    >>> p1_C = tf_C_A.apply(p1_A)
    >>> print(p1_A)  # Original point 1 in A
    [1 0 0]
    >>> print(p1_B)  # Point 1 in B
    [3. 0. 0.]
    >>> print(p1_C)  # Point 1 in C
    [ 2.59807621 -1.5       -2.        ]

    We can also do the reverse. We define a point in C and transform it to A:

    >>> p2_C = np.array([0, 1, 0])  # +1 in y_C direction
    >>> p2_A = tf_A_C.apply(p2_C)
    >>> print(p2_C)  # Original point 2 in C
    [0 1 0]
    >>> print(p2_A)  # Point 2 in A
    [-2.5       -2.         0.8660254]

    Plot the frames with respect to A again, but also plot these two points:

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_A, name="tfA")     # A plotted in A
    >>> plot_transformed_axes(ax, tf_A_B, name="tfAB")  # B plotted in A
    >>> plot_transformed_axes(ax, tf_A_C, name="tfAC")  # C plotted in A
    >>> ax.scatter(p1_A[0], p1_A[1], p1_A[2], color=colors[0])  # +1 x_A
    >>> ax.scatter(p2_A[0], p2_A[1], p2_A[2], color=colors[1])  # +1 y_C
    >>> ax.set_title("A, B, C frames and points with respect to A")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()

    **Switching Base Frames**

    Up to this point, we have been visualizing frames from A's perspective.
    Let's use the transforms we defined to visualize the frames from C's
    perspective.

    Now C is the "base frame" or "world frame". All frames are the identity
    transform from their own perspective.

    >>> tf_C = Tf.identity()

    We've already defined the transform C <- A, and can obtain C <- B by
    inverting the existing transform B <- C.

    >>> tf_C_B = tf_B_C.inv()  # C <- B

    This lets us plot everything from C's perspective:

    >>> fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    >>> plot_transformed_axes(ax, tf_C, name="tfC")     # C plotted in C
    >>> plot_transformed_axes(ax, tf_C_B, name="tfCB")  # B plotted in C
    >>> plot_transformed_axes(ax, tf_C_A, name="tfCA")  # A plotted in C
    >>> ax.scatter(p1_C[0], p1_C[1], p1_C[2], color=colors[0])
    >>> ax.scatter(p2_C[0], p2_C[1], p2_C[2], color=colors[1])
    >>> ax.set_title("A, B, C frames and points with respect to C")
    >>> ax.set_aspect("equal")
    >>> ax.figure.set_size_inches(6, 5)
    >>> plt.show()
    """

    def __init__(self, matrix: ArrayLike, normalize: bool = True, copy: bool = True):
        """Initialize from a 4x4 transformation matrix.

        Rotations are not meant to be initialized directly. Please use one of
        the `from_...` methods instead.

        Parameters
        ----------
        matrix : array_like, shape (..., 4, 4)
            A single transformation matrix or a stack of transformation
            matrices.
        normalize : bool, optional
            If True, orthonormalize the rotation matrix using singular value
            decomposition. If False, the rotation matrix is not checked for
            orthogonality or right-handedness.
        copy : bool, optional
            If True, copy the input matrix. If False, a reference to the input
            matrix is used. If normalize is True, the input matrix is always
            copied regardless of the value of copy.

        Returns
        -------
        transform : `RigidTransform` instance
        """
        xp = array_namespace(matrix)
        self._xp = xp
        matrix = _promote(matrix, xp=xp)
        if matrix.shape[-2:] != (4, 4):
            raise ValueError(
                f"Expected `matrix` to have shape (..., 4, 4), got {matrix.shape}."
            )
        # We only need the _single flag for the cython backend. The Array API backend
        # uses broadcasting by default and hence returns the correct shape without
        # additional _single logic
        self._single = matrix.ndim == 2 and is_numpy(xp)
        if self._single:
            matrix = xpx.atleast_nd(matrix, ndim=3, xp=xp)

        self._backend = select_backend(xp, matrix.ndim < 4)
        self._matrix = self._backend.from_matrix(matrix, normalize, copy)

    def __repr__(self):
        m = f"{self.as_matrix()!r}".splitlines()
        # bump indent (+27 characters)
        m[1:] = [" " * 27 + m[i] for i in range(1, len(m))]
        return "RigidTransform.from_matrix(" + "\n".join(m) + ")"

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_matrix(matrix: ArrayLike) -> RigidTransform:
        """Initialize from a 4x4 transformation matrix.

        Parameters
        ----------
        matrix : array_like, shape (..., 4, 4)
            Transformation matrices. Each matrix[..., :, :] represents a 4x4
            transformation matrix.

        Returns
        -------
        transform : `RigidTransform` instance

        Notes
        -----
        4x4 rigid transformation matrices are of the form::

            [       tx]
            [   R   ty]
            [       tz]
            [ 0 0 0  1]

        where ``R`` is a 3x3 rotation matrix and ``t = [tx, ty, tz]`` is a 3x1
        translation vector. As rotation matrices must be proper orthogonal, the
        rotation component is orthonormalized using singular value decomposition
        before initialization.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        Creating a transform from a single matrix:

        >>> m = np.array([[0, 1, 0, 2],
        ...               [0, 0, 1, 3],
        ...               [1, 0, 0, 4],
        ...               [0, 0, 0, 1]])
        >>> tf = Tf.from_matrix(m)
        >>> tf.as_matrix()
        array([[0., 1., 0., 2.],
               [0., 0., 1., 3.],
               [1., 0., 0., 4.],
               [0., 0., 0., 1.]])
        >>> tf.single
        True

        Creating a transform from an N-dimensional array of matrices:

        >>> m = np.tile(np.eye(4), (2, 5, 1, 1))  # Shape (2, 5, 4, 4)
        >>> tf = Tf.from_matrix(m)
        >>> tf.shape
        (2, 5)
        >>> tf.single
        False
        >>> len(tf)
        2

        Matrices with a rotation component that is not proper orthogonal are
        orthogonalized using singular value decomposition before initialization:

        >>> tf = Tf.from_matrix(np.diag([2, 2, 2, 1]))
        >>> tf.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])
        """
        return RigidTransform(matrix, normalize=True, copy=True)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_rotation(rotation: Rotation) -> RigidTransform:
        """Initialize from a rotation, without a translation.

        When applying this transform to a vector ``v``, the result is the
        same as if the rotation was applied to the vector.
        ``Tf.from_rotation(r).apply(v) == r.apply(v)``

        Parameters
        ----------
        rotation : `Rotation` instance
            A single rotation or a rotation with N leading dimensions.

        Returns
        -------
        transform : `RigidTransform` instance

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating a transform from a single rotation:

        >>> r = R.from_euler("ZYX", [90, 30, 0], degrees=True)
        >>> r.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> tf = Tf.from_rotation(r)
        >>> tf.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> tf.single
        True

        The upper 3x3 submatrix of the transformation matrix is the rotation
        matrix:

        >>> np.allclose(tf.as_matrix()[:3, :3], r.as_matrix(), atol=1e-12)
        True

        Creating multiple transforms from a rotation with N leading dimensions:

        >>> r = R.from_euler("ZYX", [[90, 30, 0], [45, 30, 60]], degrees=True)
        >>> r.apply([1, 0, 0])
        array([[0.        , 0.8660254 , -0.5       ],
               [0.61237244, 0.61237244, -0.5       ]])
        >>> tf = Tf.from_rotation(r)
        >>> tf.apply([1, 0, 0])
        array([[0.        , 0.8660254 , -0.5       ],
               [0.61237244, 0.61237244, -0.5       ]])
        >>> tf.single
        False
        >>> len(tf)
        2
        """
        if not isinstance(rotation, Rotation):
            raise TypeError(
                "Expected `rotation` to be a `Rotation` instance, "
                f"got {type(rotation)}."
            )
        quat = rotation.as_quat()
        xp = array_namespace(quat)
        backend = select_backend(xp, quat.ndim < 3)
        matrix = backend.from_rotation(quat)
        return RigidTransform._from_raw_matrix(matrix, xp, backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_translation(translation: ArrayLike) -> RigidTransform:
        """Initialize from a translation numpy array, without a rotation.

        When applying this transform to a vector ``v``, the result is the same
        as if the translation and vector were added together. If ``t`` is the
        displacement vector of the translation, then:

        ``Tf.from_translation(t).apply(v) == t + v``

        Parameters
        ----------
        translation : array_like, shape (..., 3)
            Translation vectors. Each translation[..., :] represents a 3D
            translation vector.

        Returns
        -------
        transform : `RigidTransform` instance

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        Creating a transform from a single translation vector:

        >>> t = np.array([2, 3, 4])
        >>> t + np.array([1, 0, 0])
        array([3, 3, 4])
        >>> tf = Tf.from_translation(t)
        >>> tf.apply([1, 0, 0])
        array([3., 3., 4.])
        >>> tf.single
        True

        The top 3x1 points in the rightmost column of the transformation matrix
        is the translation vector:

        >>> tf.as_matrix()
        array([[1., 0., 0., 2.],
               [0., 1., 0., 3.],
               [0., 0., 1., 4.],
               [0., 0., 0., 1.]])
        >>> np.allclose(tf.as_matrix()[:3, 3], t)
        True

        Creating multiple transforms from an N-dimensional array of translation
        vectors:

        >>> t = np.array([[2, 3, 4], [1, 0, 0]])
        >>> t + np.array([1, 0, 0])
        array([[3, 3, 4],
               [2, 0, 0]])
        >>> tf = Tf.from_translation(t)
        >>> tf.apply([1, 0, 0])
        array([[3., 3., 4.],
               [2., 0., 0.]])
        >>> np.allclose(tf.as_matrix()[:, :3, 3], t)
        True
        >>> tf.single
        False
        >>> len(tf)
        2
        """
        xp = array_namespace(translation)
        translation = _promote(translation, xp=xp)
        backend = select_backend(xp, translation.ndim < 3)
        matrix = backend.from_translation(translation)
        return RigidTransform._from_raw_matrix(matrix, xp, backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_components(translation: ArrayLike, rotation: Rotation) -> RigidTransform:
        """Initialize a rigid transform from translation and rotation
        components.

        When creating a rigid transform from a translation and rotation, the
        translation is applied after the rotation, such that
        ``tf = Tf.from_components(translation, rotation)``
        is equivalent to
        ``tf = Tf.from_translation(translation) * Tf.from_rotation(rotation)``.

        When applying a transform to a vector ``v``, the result is the
        same as if the transform was applied to the vector in the
        following way: ``tf.apply(v) == translation + rotation.apply(v)``

        Parameters
        ----------
        translation : array_like, shape (..., 3)
            Translation vectors. Each translation[..., :] represents a 3D
            translation vector.
        rotation : `Rotation` instance
            Rotation objects. The shape must be broadcastable with the
            translation shape.

        Returns
        -------
        transform : `RigidTransform` instance
            Rigid transform objects. The shape is determined by broadcasting
            the translation and rotation shapes together.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating from a single rotation and translation:

        >>> t = np.array([2, 3, 4])
        >>> r = R.from_euler("ZYX", [90, 30, 0], degrees=True)
        >>> r.as_matrix()
        array([[ 0.       , -1.,  0.        ],
               [ 0.8660254,  0.,  0.5       ],
               [-0.5      ,  0.,  0.8660254 ]])
        >>> tf = Tf.from_components(t, r)
        >>> tf.rotation.as_matrix()
        array([[ 0.       , -1.,  0.        ],
               [ 0.8660254,  0.,  0.5       ],
               [-0.5      ,  0.,  0.8660254 ]])
        >>> tf.translation
        array([2., 3., 4.])
        >>> tf.single
        True

        When applying a transform to a vector ``v``, the result is the same as
        if the transform was applied to the vector in the following way:
        ``tf.apply(v) == translation + rotation.apply(v)``

        >>> r.apply([1, 0, 0])
        array([0.       , 0.8660254, -0.5     ])
        >>> t + r.apply([1, 0, 0])
        array([2.       , 3.8660254,  3.5     ])
        >>> tf.apply([1, 0, 0])
        array([2.       , 3.8660254,  3.5     ])
        """
        rotation_tf = RigidTransform.from_rotation(rotation)
        return RigidTransform.from_translation(translation) * rotation_tf

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_exp_coords(exp_coords: ArrayLike) -> RigidTransform:
        r"""Initialize from exponential coordinates of transform.

        This implements the exponential map that converts 6-dimensional real
        vectors to SE(3).

        An exponential coordinate vector consists of 6 elements
        ``[rx, ry, rz, vx, vy, vz]``. The first 3 encode rotation (and form a
        rotation vector used in `Rotation.from_rotvec`) and the last 3 encode
        translation (and form a translation vector for pure translations).
        The exponential mapping can be expressed as matrix exponential
        ``T = exp(tau)``, where ``T`` is a 4x4 matrix representing a rigid
        transform and ``tau`` is a 4x4 matrix formed from the elements of an
        exponential coordinate vector::

            tau = [  0 -rz  ry vx]
                  [ rz   0 -rx vy]
                  [-ry  rx   0 vz]
                  [  0   0   0  1]

        Parameters
        ----------
        exp_coords : array_like, shape (..., 6)
            Exponential coordinate vectors. Each exp_coords[..., :] represents
            a 6D exponential coordinate vector with the expected order of
            components ``[rx, ry, rz, vx, vy, vz]``. The first 3 components
            encode rotation and the last 3 encode translation.

        Returns
        -------
        transform : `RigidTransform` instance
            Rigid transform objects with the same leading dimensions as the input.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        Creating from a single 6d vector of exponential coordinates:

        >>> tf = Tf.from_exp_coords([
        ...     -2.01041204, -0.52983629, 0.65773501,
        ...     0.10386614, 0.05855009, 0.54959179])
        >>> tf.as_matrix()
        array([[0.76406621, 0.10504613, -0.63652819, -0.10209961],
               [0.59956454, -0.47987325, 0.64050295, 0.40158789],
               [-0.2381705, -0.87102639, -0.42963687, 0.19637636],
               [0., 0., 0., 1.]])
        >>> tf.single
        True

        A vector of zeros represents the identity transform:

        >>> tf = Tf.from_exp_coords(np.zeros(6))
        >>> tf.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        The last three numbers encode translation. If the first three numbers
        are zero, the last three components can be interpreted as the
        translation:

        >>> tf_trans = Tf.from_exp_coords([0, 0, 0, 4.3, -2, 3.4])
        >>> tf_trans.translation
        array([4.3, -2., 3.4])

        The first three numbers encode rotation as a rotation vector:

        >>> tf_rot = Tf.from_exp_coords([0.5, 0.3, 0.1, 0, 0, 0])
        >>> tf_rot.rotation.as_rotvec()
        array([0.5, 0.3, 0.1])

        Combining translation and rotation preserves the rotation vector,
        but changes the last three components as they encode translation and
        rotation:

        >>> (tf_trans * tf_rot).as_exp_coords()
        array([0.5, 0.3, 0.1, 3.64305882, -1.25879559, 4.46109265])
        """
        xp = array_namespace(exp_coords)
        exp_coords = xp_promote(exp_coords, force_floating=True, xp=xp)
        if exp_coords.shape[-1] != 6:
            raise ValueError(
                f"Expected `exp_coords` to have shape (..., 6), got {exp_coords.shape}."
            )
        backend = select_backend(xp, cython_compatible=exp_coords.ndim < 3)
        matrix = backend.from_exp_coords(exp_coords)
        return RigidTransform._from_raw_matrix(matrix, xp, backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def from_dual_quat(
        dual_quat: ArrayLike, *, scalar_first: bool = False
    ) -> RigidTransform:
        """Initialize from a unit dual quaternion.

        Unit dual quaternions encode orientation in a real unit quaternion
        and translation in a dual quaternion. There is a double cover, i.e.,
        the unit dual quaternions q and -q represent the same transform.

        Unit dual quaternions must have a real quaternion with unit norm and
        a dual quaternion that is orthogonal to the real quaternion to satisfy
        the unit norm constraint. This function will enforce both properties
        through normalization.

        Parameters
        ----------
        dual_quat : array_like, shape (..., 8)
            Unit dual quaternions. Each dual_quat[..., :] represents a unit
            dual quaternion. The real part is stored in the first four
            components and the dual part in the last four components.
        scalar_first : bool, optional
            Whether the scalar component goes first or last in the two
            individual quaternions that represent the real and the dual part.
            Default is False, i.e. the scalar-last order is used.

        Returns
        -------
        transform : `RigidTransform` instance
            Rigid transform objects with the same leading dimensions as the input.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        Creating from a single unit dual quaternion:

        >>> tf = Tf.from_dual_quat([
        ...     0.0617101, -0.06483886, 0.31432811, 0.94508498,
        ...     0.04985168, -0.26119618, 0.1691491, -0.07743254])
        >>> tf.as_matrix()
        array([[0.79398752, -0.60213598, -0.08376202, 0.24605262],
               [0.58613113, 0.79477941, -0.15740392, -0.4932833],
               [0.16135089, 0.07588122, 0.98397557, 0.34262676],
               [0., 0., 0., 1.]])
        >>> tf.single
        True
        """
        xp = array_namespace(dual_quat)
        dual_quat = _promote(dual_quat, xp=xp)
        backend = select_backend(xp, dual_quat.ndim < 3)
        matrix = backend.from_dual_quat(dual_quat, scalar_first=scalar_first)
        return RigidTransform._from_raw_matrix(matrix, xp, backend)

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def identity(
        num: int | None = None, *, shape: int | tuple[int, ...] | None = None
    ) -> RigidTransform:
        """Initialize an identity transform.

        Composition with the identity transform has no effect, and
        applying the identity transform to a vector has no effect.

        Parameters
        ----------
        num : int, optional
            Number of identity transforms to generate. If None (default),
            then a single transform is generated.
        shape : int or tuple of ints, optional
            Shape of the identity transforms. If specified, `num` must
            be None.

        Returns
        -------
        transform : `RigidTransform` instance
            The identity transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Creating a single identity transform:

        >>> tf = Tf.identity()
        >>> tf.as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])
        >>> tf.single
        True

        The identity transform can be applied to a vector without effect:

        >>> tf.apply([1, 2, 3])
        array([1., 2., 3.])

        The identity transform when composed with another transform has no
        effect:

        >>> rng = np.random.default_rng(123)
        >>> t = rng.random(3)
        >>> r = R.random(rng=rng)
        >>> tf = Tf.from_components(t, r)
        >>> np.allclose((Tf.identity() * tf).as_matrix(),
        ...             tf.as_matrix(), atol=1e-12)
        True

        Multiple identity transforms can be generated at once:

        >>> tf = Tf.identity(2)
        >>> tf.as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        >>> tf.single
        False
        >>> len(tf)
        2
        """
        if num is not None and shape is not None:
            raise ValueError("Only one of `num` and `shape` can be specified.")
        if num is None and shape is None:
            shape = ()
        elif num is not None:
            shape = (num,)
        elif isinstance(shape, int):
            shape = (shape,)
        elif not isinstance(shape, tuple):
            raise ValueError("`shape` must be an int or a tuple of ints or None.")
        matrix = np.tile(np.eye(4), shape + (1, 1))
        # No need for a backend call here since identity is easy to construct and we are
        # currently not offering a backend-specific identity matrix
        return RigidTransform._from_raw_matrix(matrix, array_namespace(matrix))

    @staticmethod
    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def concatenate(
        transforms: RigidTransform | Iterable[RigidTransform],
    ) -> RigidTransform:
        """
        Concatenate a sequence of `RigidTransform` objects into a
        single object.

        Parameters
        ----------
        transforms : sequence of `RigidTransform`
            If a single `RigidTransform` instance is passed in, a copy of
            it is returned.

        Returns
        -------
        transform : `RigidTransform` instance
            The concatenated transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> tf1 = Tf.from_translation([1, 0, 0])
        >>> tf2 = Tf.from_translation([[2, 0, 0], [3, 0, 0]])
        >>> Tf.concatenate([tf1, tf2]).translation
        array([[1., 0., 0.],
               [2., 0., 0.],
               [3., 0., 0.]])
        """
        if isinstance(transforms, RigidTransform):
            return RigidTransform._from_raw_matrix(
                transforms.as_matrix(), transforms._xp
            )
        if not all(isinstance(x, RigidTransform) for x in transforms):
            raise TypeError("input must contain RigidTransform objects only")

        xp = array_namespace(transforms[0].as_matrix())
        matrix = xp.concat(
            [xpx.atleast_nd(x.as_matrix(), ndim=3, xp=xp) for x in transforms]
        )
        return RigidTransform._from_raw_matrix(matrix, xp, None)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def mean(self,
        weights: ArrayLike | None = None,
        axis: None | int | tuple[int, ...] = None
    ) -> RigidTransform:
        """Get the mean of the transforms.

        The mean of a set of transforms is the same as the mean of its
        rotation and translation components.

        The mean used for the rotation component is the chordal L2 mean (also
        called the projected or induced arithmetic mean) [1]_. If ``A`` is a
        set of rotation matrices, then the mean ``M`` is the rotation matrix
        that minimizes the following loss function:

        .. math::

            L(M) = \\sum_{i = 1}^{n} w_i \\lVert \\mathbf{A}_i -
            \\mathbf{M} \\rVert^2 ,

        where :math:`w_i`'s are the `weights` corresponding to each matrix.

        Parameters
        ----------
        weights : array_like shape (..., N), optional
            Weights describing the relative importance of the transforms. If
            None (default), then all values in `weights` are assumed to be
            equal. If given, the shape of `weights` must be broadcastable to
            the transform shape. Weights must be non-negative.
        axis : None, int, or tuple of ints, optional
            Axis or axes along which the means are computed. The default is to
            compute the mean of all transforms.

        Returns
        -------
        mean : `RigidTransform` instance
            Single transform containing the mean of the transforms in the
            current instance.

        References
        ----------
        .. [1] Hartley, Richard, et al.,
                "Rotation Averaging", International Journal of Computer Vision
                103, 2013, pp. 267-305.

        Examples
        --------
        >>> import numpy as np
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> rng = np.random.default_rng(seed=123)

        The mean of a set of transforms is the same as the mean of the
        translation and rotation components:

        >>> t = rng.random((4, 3))
        >>> r = R.random(4, rng=rng)
        >>> tf = Tf.from_components(t, r)
        >>> tf.mean().as_matrix()
        array([[ 0.61593485, -0.74508342,  0.25588075,  0.66999034],
               [-0.59353615, -0.65246765, -0.47116962,  0.25481794],
               [ 0.51801458,  0.13833531, -0.84411151,  0.52429339],
               [0., 0., 0., 1.]])
        >>> Tf.from_components(t.mean(axis=0), r.mean()).as_matrix()
        array([[ 0.61593485, -0.74508342,  0.25588075,  0.66999034],
               [-0.59353615, -0.65246765, -0.47116962,  0.25481794],
               [ 0.51801458,  0.13833531, -0.84411151,  0.52429339],
               [0., 0., 0., 1.]])
        """
        mean = self._backend.mean(self._matrix, weights=weights, axis=axis)
        return RigidTransform._from_raw_matrix(mean, xp=self._xp,
                                               backend=self._backend)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_matrix(self) -> Array:
        """Return a copy of the matrix representation of the transform.

        4x4 rigid transformation matrices are of the form::

            [       tx]
            [   R   ty]
            [       tz]
            [ 0 0 0  1]

        where ``R`` is a 3x3 orthonormal rotation matrix and
        ``t = [tx, ty, tz]`` is a 3x1 translation vector.

        Returns
        -------
        matrix : numpy.ndarray, shape (..., 4, 4)
            Transformation matrices with the same leading dimensions as the transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        A transformation matrix is a 4x4 matrix formed from a 3x3 rotation
        matrix and a 3x1 translation vector:

        >>> t = np.array([2, 3, 4])
        >>> r = R.from_matrix([[0, 0, 1],
        ...                    [1, 0, 0],
        ...                    [0, 1, 0]])
        >>> tf = Tf.from_components(t, r)
        >>> tf.as_matrix()
        array([[ 0., 0., 1., 2.],
               [ 1., 0., 0., 3.],
               [ 0., 1., 0., 4.],
               [ 0., 0., 0., 1.]])

        >>> Tf.identity(2).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]],
               [[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        """
        matrix = self._xp.asarray(self._matrix, copy=True)
        if self._single:
            return matrix[0, ...]
        return matrix

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_components(self) -> tuple[Array, Rotation]:
        """Return the translation and rotation components of the transform,
        where the rotation is applied first, followed by the translation.

        4x4 rigid transformation matrices are of the form::

            [       tx]
            [   R   ty]
            [       tz]
            [ 0 0 0  1]

        Where ``R`` is a 3x3 orthonormal rotation matrix and
        ``t = [tx, ty, tz]`` is a 3x1 translation vector. This function
        returns the rotation corresponding to this rotation matrix
        ``r = Rotation.from_matrix(R)`` and the translation vector ``t``.

        When applying a transform ``tf`` to a vector ``v``, the result is the same
        as if the rotation and translation components were applied to the vector
        with the following operation:
        ``tf.apply(v) == translation + rotation.apply(v)``.

        Returns
        -------
        translation : numpy.ndarray, shape (..., 3)
            The translation of the transform.
        rotation : `Rotation` instance
            The rotation of the transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Recover the rotation and translation from a transform:

        >>> t = np.array([2, 3, 4])
        >>> r = R.from_matrix([[0, 0, 1],
        ...                    [1, 0, 0],
        ...                    [0, 1, 0]])
        >>> tf = Tf.from_components(t, r)
        >>> t, r = tf.as_components()
        >>> t
        array([2., 3., 4.])
        >>> r.as_matrix()
        array([[0., 0., 1.],
               [1., 0., 0.],
               [0., 1., 0.]])

        The transform applied to a vector is equivalent to the rotation applied
        to the vector, followed by the translation:

        >>> t + r.apply([1, 0, 0])
        array([2., 4., 4.])
        >>> tf.apply([1, 0, 0])
        array([2., 4., 4.])
        """
        return self.translation, self.rotation

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_exp_coords(self) -> Array:
        """Return the exponential coordinates of the transform.

        This implements the logarithmic map that converts SE(3) to 6-dimensional
        real vectors.

        This is an inverse of `from_exp_coords` where details on the mapping can
        be found.

        Returns
        -------
        exp_coords : numpy.ndarray, shape (..., 6)
            Exponential coordinate vectors with the same leading dimensions as the
            transform. The first three components define the rotation and the last
            three components define the translation.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        Get exponential coordinates of the identity matrix:

        >>> Tf.identity().as_exp_coords()
        array([0., 0., 0., 0., 0., 0.])
        """
        exp_coords = self._backend.as_exp_coords(self._matrix)
        if self._single:
            exp_coords = exp_coords[0, ...]
        return exp_coords

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def as_dual_quat(self, *, scalar_first: bool = False) -> Array:
        """Return the dual quaternion representation of the transform.

        Unit dual quaternions encode orientation in a real unit quaternion
        and translation in a dual quaternion. There is a double cover, i.e.,
        the unit dual quaternions q and -q represent the same transform.

        Parameters
        ----------
        scalar_first : bool, optional
            Whether the scalar component goes first or last in the two
            individual quaternions that represent the real and the dual part.
            Default is False, i.e. the scalar-last order is used.

        Returns
        -------
        dual_quat : numpy.ndarray, shape (..., 8)
            Unit dual quaternion vectors with the same leading dimensions as the
            transform. The real part is stored in the first four components and the
            dual part in the last four components.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        Get identity dual quaternion (we use scalar-last by default):

        >>> Tf.identity().as_dual_quat()
        array([0., 0., 0., 1., 0., 0., 0., 0.])

        When we want to use the scalar-first convention, we use the argument:

        >>> Tf.identity().as_dual_quat(scalar_first=True)
        array([1., 0., 0., 0., 0., 0., 0., 0.])
        """
        dual_quat = self._backend.as_dual_quat(self._matrix, scalar_first=scalar_first)
        if self._single:
            dual_quat = dual_quat[0]
        return dual_quat

    def __len__(self) -> int:
        """Return the length of the leading transform dimension.

        A transform can store an N-dimensional array of transforms. The length is the
        size of the first dimension of this array. If the transform is a single
        transform, the length is not defined and an error is raised.

        Returns
        -------
        length : int
            The number of transforms in this object.

        Raises
        ------
        TypeError
            If the transform is a single transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np
        >>> tf = Tf.identity(3)
        >>> len(tf)
        3

        An N-dimensional array of transforms returns its first dimension size:

        >>> t = np.ones((5, 2, 3, 1, 3))
        >>> tf = Tf.from_translation(t)
        >>> len(tf)
        5

        A single transform has no length:

        >>> tf = Tf.from_translation([1, 0, 0])
        >>> len(tf)  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        TypeError: Single transform has no len().
        """
        # We don't use self._single here because we also want to raise an error for
        # Array API backends that call len() on single RigidTransform objects.
        if self.single:
            raise TypeError("Single transform has no len")
        return self._matrix.shape[0]

    @xp_capabilities(
        jax_jit=False,
        skip_backends=[("dask.array", "cannot handle zero-length rigid transforms")],
    )
    def __getitem__(
        self, indexer: int | slice | EllipsisType | None | ArrayLike
    ) -> RigidTransform:
        """Extract transform(s) at given index(es) from this object.

        Creates a new `RigidTransform` instance containing a subset of
        transforms stored in this object.

        Parameters
        ----------
        indexer : int or slice or array_like
            Specifies which transform(s) to extract. A single indexer must be
            specified, i.e. as if indexing a 1 dimensional array or list.

        Returns
        -------
        transform : `RigidTransform` instance
            Contains
                - a single transform, if `indexer` is a single index
                - a stack of transform(s), if `indexer` is a slice, or an index
                  array.

        Raises
        ------
        TypeError
            If the transform is a single transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> t = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]  # 3 translations
        >>> tf = Tf.from_translation(t)

        A single index returns a single transform:

        >>> tf[0].as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        A slice returns a stack of transforms:

        >>> tf[1:3].translation
        array([[1., 0., 0.],
               [2., 0., 0.]])

        An index array returns a stack of transforms:

        >>> tf[[0, 2]].translation
        array([[0., 0., 0.],
               [2., 0., 0.]])
        """
        if self.single:
            raise TypeError("Single transform is not subscriptable.")

        is_array = isinstance(indexer, type(self._matrix))
        xp = self._xp
        # Masking is only specified in the Array API when the array is the sole index
        # This special case handling is necessary to support boolean indexing and
        # integer array indexing with take (see
        # https://github.com/data-apis/array-api/pull/900#issuecomment-2674432480)
        # Ideally we would converge to [indexer, ...] indexing, but this is not
        # supported for now.
        if is_array and indexer.dtype == xp.bool:
            return RigidTransform(self._matrix[indexer], normalize=False)
        if is_array and xp.isdtype(indexer.dtype, "integral"):
            if self._matrix.shape[0] == 0:
                raise IndexError("cannot take from an empty array")
            return RigidTransform(
                xp.take(self._matrix, indexer, axis=0), normalize=False
            )
        return RigidTransform(self._matrix[indexer, ...], normalize=False)

    @xp_capabilities(
        jax_jit=False,
        skip_backends=[("dask.array", "cannot handle zero-length rigid transforms")],
    )
    def __setitem__(
        self,
        indexer: int | slice | EllipsisType | None | ArrayLike,
        value: RigidTransform,
    ):
        """Set transform(s) at given index(es) in this object.

        Parameters
        ----------
        indexer : int or slice or array_like
            Specifies which transform(s) to replace. A single indexer must be
            specified, i.e. as if indexing a 1 dimensional array or list.

        value : `RigidTransform` instance
            The transform(s) to set.

        Raises
        ------
        TypeError
            If the transform is a single transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> t = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]  # 3 translations
        >>> tf = Tf.from_translation(t)

        Set a single transform:

        >>> tf[0] = Tf.from_translation([9, 9, 9])
        >>> tf.translation
        array([[9., 9., 9.],
               [1., 0., 0.],
               [2., 0., 0.]])
        """
        if self.single:
            raise TypeError("Single transform is not subscriptable.")

        if not isinstance(value, RigidTransform):
            raise TypeError("value must be a RigidTransform object")

        self._matrix = self._backend.setitem(self._matrix, indexer, value.as_matrix())

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def __mul__(self, other: RigidTransform) -> RigidTransform | NotImplementedType:
        """Compose this transform with the other.

        If ``p`` and ``q`` are two transforms, then the composition of '``q``
        followed by ``p``' is equivalent to ``p * q``. In terms of
        transformation matrices, the composition can be expressed as
        ``p.as_matrix() @ q.as_matrix()``.

        In terms of translations and rotations, the composition when applied to
        a vector ``v`` is equivalent to
        ``p.translation + p.rotation.apply(q.translation)
        + (p.rotation * q.rotation).apply(v)``.

        This function supports composition of multiple transforms at a time using
        broadcasting rules. The resulting shape for two `RigidTransform` instances
        ``p`` and ``q`` is `np.broadcast_shapes(p.shape, q.shape)`.

        Parameters
        ----------
        other : `RigidTransform` instance
            Transform(s) to be composed with this one. The shapes must be
            broadcastable.

        Returns
        -------
        `RigidTransform` instance
            The composed transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Compose two transforms:

        >>> tf1 = Tf.from_translation([1, 0, 0])
        >>> tf2 = Tf.from_translation([0, 1, 0])
        >>> tf = tf1 * tf2
        >>> tf.translation
        array([1., 1., 0.])
        >>> tf.single
        True

        When applied to a vector, the composition of two transforms is applied
        in right-to-left order.

        >>> t1, r1 = [1, 2, 3], R.from_euler('z', 60, degrees=True)
        >>> t2, r2 = [0, 1, 0], R.from_euler('x', 30, degrees=True)
        >>> tf1 = Tf.from_components(t1, r1)
        >>> tf2 = Tf.from_components(t2, r2)
        >>> tf = tf1 * tf2
        >>> tf.apply([1, 0, 0])
        array([0.6339746, 3.3660254, 3.       ])
        >>> tf1.apply(tf2.apply([1, 0, 0]))
        array([0.6339746, 3.3660254, 3.       ])

        When at least one of the transforms is not single, the result is a stack
        of transforms.

        >>> tf1 = Tf.from_translation([1, 0, 0])
        >>> tf2 = Tf.from_translation([[0, 2, 0], [0, 0, 3]])
        >>> tf = tf1 * tf2
        >>> tf.translation
        array([[1., 2., 0.],
               [1., 0., 3.]])
        >>> tf.single
        False
        >>> len(tf)
        2

        Broadcasting rules apply when composing multiple transforms at a time.

        >>> tf1 = Tf.from_translation(np.ones((5, 1, 3)))  # Shape (5, 1, 3)
        >>> tf2 = Tf.from_translation(np.ones((1, 4, 3)))  # Shape (1, 4, 3)
        >>> tf = tf1 * tf2  # Shape (5, 4, 3)
        >>> tf.translation.shape
        (5, 4, 3)
        """
        if not isinstance(other, RigidTransform):
            # If other is not a RigidTransform, we return NotImplemented to allow other
            # types to implement __rmul__
            return NotImplemented
        if not broadcastable(self._matrix.shape, other._matrix.shape):
            raise ValueError(
                f"Cannot broadcast {self._matrix.shape[:-2]} transforms in "
                f"first to {other._matrix.shape[:-2]} transforms in second object."
            )
        cython_compatible = self._matrix.ndim < 4 and other._matrix.ndim < 4
        backend = select_backend(self._xp, cython_compatible=cython_compatible)
        matrix = backend.compose_transforms(self._matrix, other._matrix)
        # Only necessary for cython. Array API broadcasting handles this by default
        if self._single and other._single:
            matrix = matrix[0, ...]
        return RigidTransform(matrix, normalize=True, copy=False)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def __pow__(self, n: float) -> RigidTransform:
        """Compose this transform with itself `n` times.

        A rigid transform `p` when raised to non-integer powers can be thought
        of as finding a fraction of the transformation. For example, a power of
        0.5 finds a "halfway" transform from the identity to `p`.

        This is implemented by applying screw linear interpolation (ScLERP)
        between `p` and the identity transform, where the angle of the rotation
        component is scaled by `n`, and the translation is proportionally
        adjusted along the screw axis.

        ``q = p ** n`` can also be expressed as
        ``q = RigidTransform.from_exp_coords(p.as_exp_coords() * n)``.

        If `n` is negative, then the transform is inverted before the power
        is applied. In other words, ``p ** -abs(n) == p.inv() ** abs(n)``.

        Parameters
        ----------
        n : float
            The number of times to compose the transform with itself.

        Returns
        -------
        `RigidTransform` instance
            If the input Rotation `p` contains `N` multiple rotations, then
            the output will contain `N` rotations where the `i` th rotation
            is equal to ``p[i] ** n``.

        Notes
        -----
        There are three notable cases: if ``n == 1`` then a copy of the original
        transform is returned, if ``n == 0`` then the identity transform is
        returned, and if ``n == -1`` then the inverse transform is returned.

        Note that fractional powers ``n`` which effectively take a root of
        rotation, do so using the shortest path smallest representation of that
        angle (the principal root). This means that powers of ``n`` and ``1/n``
        are not necessarily inverses of each other. For example, a 0.5 power of
        a +240 degree rotation will be calculated as the 0.5 power of a -120
        degree rotation, with the result being a rotation of -60 rather than
        +120 degrees.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> import numpy as np

        A power of 2 returns the transform composed with itself:

        >>> tf = Tf.from_translation([1, 2, 3])
        >>> (tf ** 2).translation
        array([2., 4., 6.])
        >>> (tf ** 2).as_matrix()
        array([[1., 0., 0., 2.],
               [0., 1., 0., 4.],
               [0., 0., 1., 6.],
               [0., 0., 0., 1.]])

        A negative power returns the inverse of the transform raised to the
        absolute value of `n`:

        >>> (tf ** -2).translation
        array([-2., -4., -6.])
        >>> np.allclose((tf ** -2).as_matrix(), (tf.inv() ** 2).as_matrix(),
        ...             atol=1e-12)
        True

        A power of 0 returns the identity transform:

        >>> (tf ** 0).as_matrix()
        array([[1., 0., 0., 0.],
               [0., 1., 0., 0.],
               [0., 0., 1., 0.],
               [0., 0., 0., 1.]])

        A power of 1 returns a copy of the original transform:

        >>> (tf ** 1).as_matrix()
        array([[1., 0., 0., 1.],
               [0., 1., 0., 2.],
               [0., 0., 1., 3.],
               [0., 0., 0., 1.]])

        A fractional power returns a transform with a scaled rotation and
        translated along the screw axis. Here we take the square root of the
        transform, which when squared recovers the original transform:

        >>> tf_half = (tf ** 0.5)
        >>> tf_half.translation
        array([0.5, 1., 1.5])
        >>> (tf_half ** 2).as_matrix()
        array([[1., 0., 0., 1.],
               [0., 1., 0., 2.],
               [0., 0., 1., 3.],
               [0., 0., 0., 1.]])
        """
        matrix = self._backend.pow(self._matrix, n)
        if self._single:
            matrix = matrix[0, ...]
        return RigidTransform._from_raw_matrix(matrix, self._xp, self._backend)

    @xp_capabilities(
        skip_backends=[("dask.array", "missing linalg.cross/det functions")]
    )
    def inv(self) -> RigidTransform:
        """Invert this transform.

        Composition of a transform with its inverse results in an identity
        transform.

        A rigid transform is a composition of a rotation and a translation,
        where the rotation is applied first, followed by the translation. So the
        inverse transform is equivalent to the inverse translation followed by
        the inverse rotation.

        Returns
        -------
        `RigidTransform` instance
            The inverse of this transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        A transform composed with its inverse results in an identity transform:

        >>> rng = np.random.default_rng(seed=123)
        >>> t = rng.random(3)
        >>> r = R.random(rng=rng)
        >>> tf = Tf.from_components(t, r)
        >>> tf.as_matrix()
        array([[-0.45431291,  0.67276178, -0.58394466,  0.68235186],
               [-0.23272031,  0.54310598,  0.80676958,  0.05382102],
               [ 0.85990758,  0.50242162, -0.09017473,  0.22035987],
               [ 0.        ,  0.        ,  0.        ,  1.        ]])

        >>> (tf.inv() * tf).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])

        The inverse rigid transform is the same as the inverse translation
        followed by the inverse rotation:

        >>> t, r = tf.as_components()
        >>> r_inv = r.inv()  # inverse rotation
        >>> t_inv = -t  # inverse translation
        >>> tf_r_inv = Tf.from_rotation(r_inv)
        >>> tf_t_inv = Tf.from_translation(t_inv)
        >>> np.allclose((tf_r_inv * tf_t_inv).as_matrix(),
        ...             tf.inv().as_matrix(),
        ...             atol=1e-12)
        True
        >>> (tf_r_inv * tf_t_inv * tf).as_matrix()
        array([[[1., 0., 0., 0.],
                [0., 1., 0., 0.],
                [0., 0., 1., 0.],
                [0., 0., 0., 1.]]])
        """
        matrix = self._backend.inv(self._matrix)
        if self._single:
            matrix = matrix[0, ...]
        return RigidTransform._from_raw_matrix(matrix, self._xp, self._backend)

    @xp_capabilities(
        skip_backends=[
            ("dask.array", "missing linalg.cross/det functions"),
            ("cupy", "missing .mT attribute in cupy<14.*"),
        ]
    )
    def apply(self, vector: ArrayLike, inverse: bool = False) -> Array:
        """Apply the transform to a vector.

        If the original frame transforms to the final frame by this transform,
        then its application to a vector can be seen in two ways:

            - As a projection of vector components expressed in the final frame
              to the original frame.
            - As the physical transformation of a vector being glued to the
              original frame as it transforms. In this case the vector
              components are expressed in the original frame before and after
              the transformation.

        In terms of rotation matrices and translation vectors, this application
        is the same as
        ``self.translation + self.rotation.as_matrix() @ vector``.

        Parameters
        ----------
        vector : array_like, shape (..., 3)
            Vector(s) to be transformed. Each vector[..., :] represents a 3D
            vector.
        inverse : bool, optional
            If True, the inverse of the transform is applied to the vector.

        Returns
        -------
        transformed_vector : numpy.ndarray, shape (..., 3)
            The transformed vector(s) with shape determined by broadcasting
            the transform and vector shapes together.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        Apply a single transform to a vector. Here the transform is just a
        translation, so the result is the vector added to the translation
        vector.

        >>> t = np.array([1, 2, 3])
        >>> tf = Tf.from_translation(t)
        >>> t + np.array([1, 0, 0])
        array([2, 2, 3])
        >>> tf.apply([1, 0, 0])
        array([2., 2., 3.])

        Apply a single transform to a stack of vectors:

        >>> tf.apply([[1, 0, 0], [0, 1, 0]])
        array([[2., 2., 3.],
               [1., 3., 3.]])

        Apply the inverse of a transform to a vector, so the result is the
        negative of the translation vector added to the vector.

        >>> -t + np.array([1, 0, 0])
        array([0, -2, -3])
        >>> tf.apply([1, 0, 0], inverse=True)
        array([0., -2., -3.])

        Broadcasting is supported when applying multiple transforms to an N-dimensional
        array of vectors.

        >>> tf = Tf.from_translation(np.ones((4, 5, 1, 3)))
        >>> vectors = np.zeros((2, 3))
        >>> tf.apply(vectors).shape
        (4, 5, 2, 3)

        For transforms which are not just pure translations, applying it to a
        vector is the same as applying the rotation component to the vector and
        then adding the translation component.

        >>> r = R.from_euler('z', 60, degrees=True)
        >>> tf = Tf.from_components(t, r)
        >>> t + r.apply([1, 0, 0])
        array([1.5,       2.8660254, 3.       ])
        >>> tf.apply([1, 0, 0])
        array([1.5,       2.8660254, 3.       ])

        When applying the inverse of a transform, the result is the negative of
        the translation vector added to the vector, and then rotated by the
        inverse rotation.

        >>> r.inv().apply(-t + np.array([1, 0, 0]))
        array([-1.73205081, -1.        , -3.        ])
        >>> tf.apply([1, 0, 0], inverse=True)
        array([-1.73205081, -1.        , -3.        ])
        """
        # We do not use the cached xp here to catch cases where matrix and vector have
        # different xp instances.
        xp = array_namespace(self._matrix, vector)
        vector = xp.asarray(
            vector, dtype=self._matrix.dtype, device=device(self._matrix)
        )
        cython_compatible = self._matrix.ndim < 4 and vector.ndim < 3
        backend = select_backend(xp, cython_compatible=cython_compatible)
        result = backend.apply(self._matrix, vector, inverse)
        if self._single and vector.ndim == 1:
            result = result[0, ...]
        return result

    @property
    def rotation(self) -> Rotation:
        """Return the rotation component of the transform.

        A transform is a composition of a rotation and a translation, such that
        when applied to a vector, the vector is first rotated and then
        translated. This property returns the rotation part of the transform.

        Returns
        -------
        rotation : `Rotation` instance
            A single rotation or a stack of rotations.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        The rotation component is extracted from the transform:

        >>> t = np.array([1, 0, 0])
        >>> r = R.random(3)
        >>> tf = Tf.from_components(t, r)
        >>> np.allclose(tf.rotation.as_matrix(), r.as_matrix())
        True
        """
        if self._single:
            return Rotation.from_matrix(self._matrix[0, :3, :3])
        return Rotation.from_matrix(self._matrix[..., :3, :3])

    @property
    def translation(self) -> Array:
        """Return the translation component of the transform.

        A transform is a composition of a rotation and a translation, such that
        when applied to a vector, the vector is first rotated and then
        translated. This property returns the translation part of the transform.

        Returns
        -------
        translation : numpy.ndarray, shape (..., 3)
            Translation vectors with the same leading dimensions as the transform.

        Examples
        --------
        >>> from scipy.spatial.transform import RigidTransform as Tf
        >>> from scipy.spatial.transform import Rotation as R
        >>> import numpy as np

        The translation component is extracted from the transform:

        >>> t = np.array([[1, 0, 0], [2, 0, 0], [3, 0, 0]])
        >>> r = R.random()
        >>> tf = Tf.from_components(t, r)
        >>> np.allclose(tf.translation, t)
        True
        """
        if self._single:
            return self._xp.asarray(self._matrix[0, :3, 3], copy=True)
        return self._xp.asarray(self._matrix[..., :3, 3], copy=True)

    @property
    def single(self) -> bool:
        """Whether this instance represents a single transform.

        Single transforms are not subscriptable, and do not have a length.

        Returns
        -------
        single : bool
            True if this instance represents a single transform, False
            otherwise.
        """
        return self._single or self._matrix.ndim == 2

    @property
    def shape(self) -> tuple[int, ...]:
        """The shape of the transform's leading dimensions."""
        if self._single or self._matrix.ndim == 2:
            return ()
        return self._matrix.shape[:-2]

    def __reduce__(self) -> tuple[Callable, tuple]:
        """Reduce the RigidTransform for pickling.

        We store modules inside RigidTransforms which cannot be pickled. To circumvent
        this, we pickle only the matrix and restore the cached modules from the matrix
        type in `from_matrix`.
        """
        matrix = self._matrix
        if self._single:
            matrix = matrix[0, ...]
        return (self.__class__.from_matrix, (matrix,))

    @xp_capabilities()
    def __iter__(self) -> Iterator[RigidTransform]:
        """Iterate over transforms."""
        if self._single or self._matrix.ndim == 2:
            raise TypeError("Single transform is not iterable.")
        # We return a generator that yields a new RigidTransform object for each
        # transform in the current object. We cannot rely on the default implementation
        # because jax will not raise an IndexError for out-of-bounds indices.
        for i in range(self._matrix.shape[0]):
            yield RigidTransform._from_raw_matrix(
                self._matrix[i, ...], self._xp, self._backend
            )

    @staticmethod
    def _from_raw_matrix(
        matrix: Array, xp: ModuleType, backend: ModuleType | None = None
    ) -> RigidTransform:
        """Create a RigidTransform skipping all sanitization steps.

        This method is is intended for internal, performant creation of RigidTransforms
        with matrices that are guaranteed to be valid.
        """
        tf = RigidTransform.__new__(RigidTransform)
        tf._single = matrix.ndim == 2 and is_numpy(xp)
        if tf._single:
            matrix = xpx.atleast_nd(matrix, ndim=3, xp=xp)
        tf._matrix = matrix
        tf._xp = xp
        if backend is None:
            backend = select_backend(xp, matrix.ndim < 4)
        tf._backend = backend
        return tf
