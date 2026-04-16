"""Array API backend for the `Rotation` class.

This module provides generic, functional implementations of the `Rotation` class methods
that work with any Array API-compatible backend.
"""

# Parts of the implementation are adapted from the cython backend and
# https://github.com/jax-ml/jax/blob/d695aa4c63ffcebefce52794427c46bad576680c/jax/_src/scipy/spatial/transform.py.
import re
import warnings
from types import EllipsisType

import numpy as np
from scipy._lib._array_api import (
    array_namespace,
    Array,
    ArrayLike,
    is_lazy_array,
    xp_vector_norm,
    xp_result_type,
    xp_promote,
    is_jax,
)
from scipy._lib._util import broadcastable
from scipy._lib.array_api_compat import device as xp_device
from scipy._lib.array_api_compat import is_array_api_obj
import scipy._lib.array_api_extra as xpx


def from_quat(
    quat: Array,
    normalize: bool = True,
    copy: bool = True,
    *,
    scalar_first: bool = False,
) -> Array:
    xp = array_namespace(quat)
    if scalar_first:
        quat = xp.roll(quat, -1, axis=-1)
    # Normalize will always copy, so we avoid the extra copy if we normalize
    if copy and not normalize:
        quat = xp.asarray(quat, copy=True)
    if normalize:
        quat = _normalize_quaternion(quat)
    return quat


def from_matrix(matrix: Array, assume_valid: bool = False) -> Array:
    xp = array_namespace(matrix)
    device = xp_device(matrix)

    if not assume_valid:
        mask = xp.linalg.det(matrix) <= 0
        lazy = is_lazy_array(mask)
        # Only non-lazy backends raise an error for non-positive determinants.
        if not lazy and xp.any(mask):
            ind = int(xp.nonzero(xpx.atleast_nd(mask, ndim=1, xp=xp))[0][0])
            raise ValueError(
                "Non-positive determinant (left-handed or null coordinate frame) in "
                f"rotation matrix {ind}: {matrix[ind, ...]}."
            )
        elif lazy:
            matrix = xp.where(mask[..., None, None], xp.nan, matrix)

        gramians = matrix @ xp.matrix_transpose(matrix)
        eye = xp.eye(3, dtype=matrix.dtype, device=device)
        is_orthogonal = xp.all(
            xpx.isclose(gramians, eye, atol=1e-12, xp=xp), axis=(-2, -1)
        )

        if lazy:
            # Lazy backends do not support non-concrete boolean indexing or any form of
            # computation without statically known shapes, so we always compute SVD and
            # use xp.where to select the result.
            U, _, Vt = xp.linalg.svd(matrix, full_matrices=False)
            matrix = xp.where(is_orthogonal[..., None, None], matrix, U @ Vt)
        elif not xp.all(is_orthogonal):
            # For eager frameworks, only compute SVD if needed.
            is_not_orthogonal = ~is_orthogonal
            U, _, Vt = xp.linalg.svd(matrix[is_not_orthogonal], full_matrices=False)
            matrix = xpx.at(matrix)[is_not_orthogonal].set(U @ Vt)

    return _from_matrix_orthogonal(matrix)


def _from_matrix_orthogonal(matrix: Array) -> Array:
    """Convert known orthogonal rotation matrix to quaternion"""
    xp = array_namespace(matrix)
    device = xp_device(matrix)

    matrix_trace = matrix[..., 0, 0] + matrix[..., 1, 1] + matrix[..., 2, 2]
    decision = xp.stack(
        [matrix[..., 0, 0], matrix[..., 1, 1], matrix[..., 2, 2], matrix_trace],
        axis=-1,
    )
    choice = xp.argmax(decision, axis=-1, keepdims=True)
    quat = xp.empty((*matrix.shape[:-2], 4), dtype=matrix.dtype, device=device)

    # The Array API does not support mixing integer indexing with ellipsis, so we cannot
    # follow the same pattern as the cython backend. Instead, we compute each case
    # explicitly and assemble the final result with `xp.where`.
    # TODO: Revisit this implementation if the array API supports mixed integer and
    # ellipsis indexing.

    # Case 0
    quat_0 = xp.stack(
        [
            1 - matrix_trace[...] + 2 * matrix[..., 0, 0],
            matrix[..., 1, 0] + matrix[..., 0, 1],
            matrix[..., 2, 0] + matrix[..., 0, 2],
            matrix[..., 2, 1] - matrix[..., 1, 2],
        ],
        axis=-1,
    )
    quat = xp.where(choice == 0, quat_0, quat)

    # Case 1
    quat_1 = xp.stack(
        [
            matrix[..., 1, 0] + matrix[..., 0, 1],
            1 - matrix_trace[...] + 2 * matrix[..., 1, 1],
            matrix[..., 2, 1] + matrix[..., 1, 2],
            matrix[..., 0, 2] - matrix[..., 2, 0],
        ],
        axis=-1,
    )
    quat = xp.where(choice == 1, quat_1, quat)

    # Case 2
    quat_2 = xp.stack(
        [
            matrix[..., 2, 0] + matrix[..., 0, 2],
            matrix[..., 2, 1] + matrix[..., 1, 2],
            1 - matrix_trace[...] + 2 * matrix[..., 2, 2],
            matrix[..., 1, 0] - matrix[..., 0, 1],
        ],
        axis=-1,
    )
    quat = xp.where(choice == 2, quat_2, quat)

    # Case 3
    quat_3 = xp.stack(
        [
            matrix[..., 2, 1] - matrix[..., 1, 2],
            matrix[..., 0, 2] - matrix[..., 2, 0],
            matrix[..., 1, 0] - matrix[..., 0, 1],
            1 + matrix_trace[...],
        ],
        axis=-1,
    )
    quat = xp.where(choice == 3, quat_3, quat)

    return _normalize_quaternion(quat)


def from_rotvec(rotvec: Array, degrees: bool = False) -> Array:
    xp = array_namespace(rotvec)
    if rotvec.shape[-1] != 3:
        raise ValueError(
            f"Expected `rot_vec` to have shape (..., 3), got {rotvec.shape}"
        )
    rotvec = _deg2rad(rotvec) if degrees else rotvec

    angle = xp_vector_norm(rotvec, axis=-1, keepdims=True, xp=xp)
    small_angle = angle <= 1e-3
    angle2 = angle**2
    small_scale = 0.5 - angle2 / 48 + angle2**2 / 3840
    # We need to handle the case where angle is 0 to avoid division by zero. We use the
    # value of the Taylor series approximation, but non-branching operations require
    # that we still divide by the angle. Since we do not use the result where the angle
    # is close to 0, this is safe.
    div_angle = angle + xp.asarray(small_angle, dtype=angle.dtype)
    large_scale = xp.sin(angle / 2) / div_angle
    scale = xp.where(small_angle, small_scale, large_scale)
    quat = xp.concat([rotvec * scale, xp.cos(angle / 2)], axis=-1)
    return quat


def from_mrp(mrp: Array) -> Array:
    xp = array_namespace(mrp)
    if mrp.shape[-1] != 3:
        raise ValueError(f"Expected `mrp` to have shape (..., 3), got {mrp.shape}")
    mrp2_plus_1 = xp.linalg.vecdot(mrp, mrp, axis=-1)[..., None] + 1
    q_no_norm = xp.concat([2 * mrp[..., :3], (2 - mrp2_plus_1)], axis=-1)
    quat = q_no_norm / mrp2_plus_1
    return quat


def from_euler(seq: str, angles: Array, degrees: bool = False) -> Array:
    xp = array_namespace(angles)
    num_axes = len(seq)
    if num_axes < 1 or num_axes > 3:
        raise ValueError(
            "Expected axis specification to be a non-empty "
            f"string of up to 3 characters, got {seq}"
        )

    intrinsic = re.match(r"^[XYZ]{1,3}$", seq) is not None
    extrinsic = re.match(r"^[xyz]{1,3}$", seq) is not None
    if not (intrinsic or extrinsic):
        raise ValueError(
            "Expected axes from `seq` to be from ['x', 'y', "
            f"'z'] or ['X', 'Y', 'Z'], got {seq}"
        )

    if any(seq[i] == seq[i + 1] for i in range(num_axes - 1)):
        raise ValueError(f"Expected consecutive axes to be different, got {seq}")

    if degrees:
        angles = _deg2rad(angles)

    angles = xpx.atleast_nd(angles, ndim=1, xp=xp)

    if angles.shape[-1] != num_axes:
        raise ValueError(
            "Expected last dimension of `angles` to match number of sequence axes "
            f"specified, got {angles.shape[-1]}."
        )
    axes = [_elementary_basis_index(x) for x in seq.lower()]
    q = _elementary_quat_compose(axes, angles, intrinsic)
    return q


def from_davenport(
    axes: Array, order: str, angles: Array | float, degrees: bool = False
) -> Array:
    xp = array_namespace(axes)
    device = xp_device(axes)
    if order in ["e", "extrinsic"]:  # Must be static, cannot be jitted
        extrinsic = True
    elif order in ["i", "intrinsic"]:
        extrinsic = False
    else:
        raise ValueError(
            "order should be 'e'/'extrinsic' for extrinsic sequences or 'i'/"
            f"'intrinsic' for intrinsic sequences, got {order}"
        )

    if axes.shape[-1] != 3:
        raise ValueError("Axes must be vectors of length 3.")

    axes = xpx.atleast_nd(axes, ndim=2, xp=xp)
    angles = xpx.atleast_nd(angles, ndim=1, xp=xp)
    num_axes = axes.shape[-2]
    if num_axes < 1 or num_axes > 3:
        raise ValueError(f"Expected up to 3 axes, got {num_axes}")

    axes = axes / xp_vector_norm(axes, axis=-1, keepdims=True, xp=xp)

    # Check if axes are orthogonal. Shape checks also work for lazy backends.
    axes_not_orthogonal = xp.zeros(axes.shape[:-2], dtype=xp.bool, device=device)
    if num_axes > 1:
        # Cannot be True yet, so we do not need to use xp.logical_or
        axes_not_orthogonal = axes_not_orthogonal | (
            xp.abs(xp.vecdot(axes[..., 0, :], axes[..., 1, :])) > 1e-7
        )
    if num_axes > 2:
        axes_not_orthogonal = axes_not_orthogonal | (
            xp.abs(xp.vecdot(axes[..., 1, :], axes[..., 2, :])) > 1e-7
        )
    if not is_lazy_array(axes_not_orthogonal) and xp.any(axes_not_orthogonal):
        raise ValueError("Consecutive axes must be orthogonal.")
    else:
        axes = xp.where(axes_not_orthogonal[..., None, None], xp.nan, axes)

    if degrees:
        angles = _deg2rad(angles)

    if (
        not broadcastable(axes.shape[:-1], angles.shape)
        or axes.shape[-2] != angles.shape[-1]
    ):
        raise ValueError(
            f"Expected `angles` to match number of axes, got {angles.shape} angles "
            f"and {axes.shape} axes."
        )

    q_shape = angles.shape[:-1] + (4,)
    q = xp.zeros(q_shape, dtype=angles.dtype, device=xp_device(angles))
    q = xpx.at(q)[..., 3].set(1)

    for i in range(num_axes):
        qi = from_rotvec(angles[..., i, None] * axes[..., i, :])
        q = compose_quat(qi, q) if extrinsic else compose_quat(q, qi)
    return q


def as_quat(
    quat: Array, canonical: bool = False, *, scalar_first: bool = False
) -> Array:
    xp = array_namespace(quat)
    if canonical:
        quat = _quat_canonical(quat)
    if scalar_first:
        quat = xp.roll(quat, 1, axis=-1)
    return quat


def as_matrix(quat: Array) -> Array:
    xp = array_namespace(quat)
    x = quat[..., 0]
    y = quat[..., 1]
    z = quat[..., 2]
    w = quat[..., 3]

    x2 = x * x
    y2 = y * y
    z2 = z * z
    w2 = w * w

    xy = x * y
    zw = z * w
    xz = x * z
    yw = y * w
    yz = y * z
    xw = x * w

    matrix_elements = [
        x2 - y2 - z2 + w2,
        2 * (xy - zw),
        2 * (xz + yw),
        2 * (xy + zw),
        -x2 + y2 - z2 + w2,
        2 * (yz - xw),
        2 * (xz - yw),
        2 * (yz + xw),
        -x2 - y2 + z2 + w2,
    ]
    matrix = xp.reshape(xp.stack(matrix_elements, axis=-1), (*quat.shape[:-1], 3, 3))
    return matrix


def as_rotvec(quat: Array, degrees: bool = False) -> Array:
    xp = array_namespace(quat)
    quat = _quat_canonical(quat)
    ax_norm = xp_vector_norm(quat[..., :3], axis=-1, keepdims=True, xp=xp)
    angle = 2 * xp.atan2(ax_norm, quat[..., 3][..., None])
    small_angle = angle <= 1e-3
    angle2 = angle**2
    small_scale = 2 + angle2 / 12 + 7 * angle2**2 / 2880
    # We need to handle the case where sin(angle/2) is 0 to avoid division by zero. We
    # use the value of the Taylor series approximation, but non-branching operations
    # require that we still divide by the sin. Since we do not use the result where the
    # angle is close to 0, adding one to the sin where we discard the result is safe.
    div_sin = xp.sin(angle / 2.0) + xp.asarray(small_angle, dtype=angle.dtype)
    large_scale = angle / div_sin
    scale = xp.where(small_angle, small_scale, large_scale)
    if degrees:
        scale = _rad2deg(scale)
    rotvec = scale * quat[..., :3]
    return rotvec


def as_mrp(quat: Array) -> Array:
    xp = array_namespace(quat)
    one = xp.asarray(1.0, device=xp_device(quat), dtype=quat.dtype)
    sign = xp.where(quat[..., 3, None] < 0, -1, one)
    denominator = 1.0 + sign * quat[..., 3, None]
    return sign * quat[..., :3] / denominator


def as_euler(
    quat: Array, seq: str, degrees: bool = False, *, suppress_warnings: bool = False
) -> Array:
    xp = array_namespace(quat)

    # Sanitize the sequence
    if len(seq) != 3:
        raise ValueError(f"Expected 3 axes, got {seq}.")
    intrinsic = re.match(r"^[XYZ]{1,3}$", seq) is not None
    extrinsic = re.match(r"^[xyz]{1,3}$", seq) is not None
    if not (intrinsic or extrinsic):
        raise ValueError(
            "Expected axes from `seq` to be from ['x', 'y', 'z'] or ['X', 'Y', 'Z'], "
            f"got {seq}"
        )
    if any(seq[i] == seq[i + 1] for i in range(2)):
        raise ValueError(f"Expected consecutive axes to be different, got {seq}")

    device = xp_device(quat)
    axes = [_elementary_basis_index(x) for x in seq.lower()]
    axes = axes if extrinsic else axes[::-1]
    i, j, k = axes
    symmetric = i == k
    k = 3 - i - j if symmetric else k

    mask = xp.asarray(symmetric, device=device)
    sign = xp.asarray((i - j) * (j - k) * (k - i) // 2, dtype=quat.dtype, device=device)
    # Permute quaternion elements
    a = xp.where(mask, quat[..., 3], quat[..., 3] - quat[..., j])
    b = xp.where(mask, quat[..., i], quat[..., i] + quat[..., k] * sign)
    c = xp.where(mask, quat[..., j], quat[..., j] + quat[..., 3])
    d = xp.where(mask, quat[..., k] * sign, quat[..., k] * sign - quat[..., i])

    angles = _get_angles(
        extrinsic, symmetric, sign, xp.pi / 2, a, b, c, d, suppress_warnings
    )
    return _rad2deg(angles) if degrees else angles


def as_davenport(
    quat: Array,
    axes: ArrayLike,
    order: str,
    degrees: bool = False,
    *,
    suppress_warnings: bool = False,
) -> Array:
    xp = array_namespace(quat)

    # Check argument validity
    if order in ["e", "extrinsic"]:
        extrinsic = True
    elif order in ["i", "intrinsic"]:
        extrinsic = False
    else:
        raise ValueError(
            "order should be 'e'/'extrinsic' for extrinsic sequences or 'i'/'intrinsic'"
            f" for intrinsic sequences, got {order}"
        )
    if axes.shape[-2] != 3:
        raise ValueError(f"Expected 3 axes, got {axes.shape}.")
    if axes.shape[-1] != 3:
        raise ValueError("Axes must be vectors of length 3.")
    if not broadcastable(axes.shape[:-2], quat.shape[:-1]):
        raise ValueError(
            f"Expected `axes` to match number of rotations, got {axes.shape} axes "
            f"and {quat.shape} rotations."
        )

    # normalize axes
    axes = axes / xp_vector_norm(axes, axis=-1, keepdims=True, xp=xp)
    vdot_ax0_ax1 = xp.vecdot(axes[..., 0, :], axes[..., 1, :])
    vdot_ax1_ax2 = xp.vecdot(axes[..., 1, :], axes[..., 2, :])
    is_invalid = (vdot_ax0_ax1 >= 1e-7) | (vdot_ax1_ax2 >= 1e-7)
    if is_lazy_array(is_invalid):
        axes = xp.where(is_invalid[..., None, None], xp.nan, axes)
    elif xp.any(is_invalid):
        raise ValueError("Consecutive axes must be orthogonal.")

    angles = _compute_davenport_from_quat(
        quat,
        axes[..., 0, :],
        axes[..., 1, :],
        axes[..., 2, :],
        extrinsic,
        suppress_warnings,
    )
    if degrees:
        angles = _rad2deg(angles)
    return angles


def inv(quat: Array) -> Array:
    return xpx.at(quat)[..., :3].multiply(-1, copy=True)


def magnitude(quat: Array) -> Array:
    xp = array_namespace(quat)
    sin_q = xp_vector_norm(quat[..., :3], axis=-1, xp=xp)
    cos_q = xp.abs(quat[..., 3])
    angles = 2 * xp.atan2(sin_q, cos_q)
    return angles


def approx_equal(
    quat: Array, other: Array, atol: float | None = None, degrees: bool = False
) -> Array:
    if atol is None:
        if degrees:
            warnings.warn(
                "atol must be set to use the degrees flag, defaulting to 1e-8 radians.",
                stacklevel=2,
            )
        atol = 1e-8
    elif degrees:
        atol = _deg2rad(atol)

    if not broadcastable(quat.shape, other.shape):
        raise ValueError(
            f"Expected broadcastable shapes in both rotations, got {quat.shape[:-1]} "
            f"rotations in first and {other.shape[:-1]} rotations in second object."
        )

    quat_result = compose_quat(other, inv(quat))
    angles = magnitude(quat_result)
    return angles < atol


def mean(
    quat: Array,
    weights: ArrayLike | None = None,
    axis: None | int | tuple[int, ...] = None,
) -> Array:
    xp = array_namespace(quat)
    device = xp_device(quat)
    dtype = xp_result_type(quat, force_floating=True, xp=xp)
    if quat.shape[0] == 0:
        raise ValueError("Mean of an empty rotation set is undefined.")
    # Axis logic: For None, we reduce over all axes. For int, we only reduce over that
    # axis. For tuple, we reduce over all specified axes.
    all_axes = tuple(range(quat.ndim - 1))
    if axis is None:
        axis = all_axes
    elif isinstance(axis, int):
        axis = (axis,)
    if not isinstance(axis, tuple):
        raise ValueError("`axis` must be None, int, or tuple of ints.")
    # Ensure all axes are within bounds
    if axis != () and (min(axis) < -(quat.ndim - 1) or max(axis) > (quat.ndim - 2)):
        raise ValueError(
            f"axis {axis} is out of bounds for rotation with shape {quat.shape[:-1]}."
        )
    # Ensure all axes are positive and unique
    axis = tuple(sorted(set(x % (quat.ndim - 1) for x in axis)))

    lazy = is_lazy_array(quat)
    # Branching code is okay for checks that include meta info such as shapes and types
    quat_expand = quat[..., None, :]
    if weights is None:
        K = xp.matrix_transpose(quat_expand) @ quat_expand
    else:
        weights = xp.asarray(weights, dtype=dtype, device=device)
        neg_weights = weights < 0
        if not lazy and xp.any(neg_weights):
            raise ValueError("`weights` must be non-negative.")
        elif lazy:
            # We cannot check for negative weights because jit code needs to be
            # non-branching. We return NaN instead
            weights = xp.where(neg_weights, xp.nan, weights)

        if not broadcastable(quat.shape[:-1], weights.shape):
            raise ValueError(
                "Expected `weights` to be broadcastable to rotation shape, got shape "
                f"{weights.shape} for {quat.shape[:-1]} rotations."
            )

        # Make sure we can transpose quat
        weighted_quat = weights[..., None, None] * quat_expand
        K = xp.matrix_transpose(weighted_quat) @ quat_expand

    # Move reduction axes to the end
    keep_axes = tuple(i for i in all_axes if i not in axis)
    axes_order = keep_axes + axis
    K_reordered = xp.moveaxis(K, axes_order, all_axes)
    # Reshape to flatten reduction axes
    new_shape = K_reordered.shape[: len(keep_axes)] + (-1, 4, 4)
    K = xp.mean(xp.reshape(K_reordered, new_shape), axis=-3)
    _, v = xp.linalg.eigh(K)
    return v[..., -1]


def reduce(
    quat: Array,
    left: Array | None = None,
    right: Array | None = None,
) -> tuple[Array, Array | None, Array | None]:
    if left is None and right is None:
        return quat, None, None
    # DECISION: We cannot have variable number of return arguments for jit compiled
    # functions. We therefore always return the indices, and filter out later.
    # TOOD: Properly support broadcasting.
    xp = array_namespace(quat)
    quat = xpx.atleast_nd(quat, ndim=2, xp=xp)
    if left is None:
        left = xp.ones_like(quat)
    if right is None:
        right = xp.ones_like(quat)

    # We want to calculate the real components of q = l * p * r. It can
    # be shown that:
    #     qs = ls * ps * rs - ls * dot(pv, rv) - ps * dot(lv, rv)
    #          - rs * dot(lv, pv) - dot(cross(lv, pv), rv)
    # where ls and lv denote the scalar and vector components of l.

    p = quat
    ps, pv = _split_rotation(p, xp)
    ls, lv = _split_rotation(left, xp)
    rs, rv = _split_rotation(right, xp)

    # Compute each term without einsum (not accessible in the Array API)
    # First term: np.einsum("i,j,k", ls, ps, rs)
    term1 = ls[..., :, None, None] * ps[..., None, :, None] * rs[..., None, None, :]
    # Second term: np.einsum('i,jx,kx', ls, pv, rv)
    prv = xp.sum(pv[..., :, None, :] * rv[..., None, :, :], axis=-1)
    term2 = ls[..., :, None, None] * prv[..., None, :, :]
    # Third term: np.einsum('ix,j,kx', lv, ps, rv)
    lrv = xp.sum(lv[..., :, None, :] * rv[..., None, :, :], axis=-1)
    term3 = ps[..., None, :, None] * lrv[..., :, None, :]
    # Fourth term: np.einsum('ix,jx,k', lv, pv, rs)
    lpv = xp.sum(lv[..., :, None, :] * pv[..., None, :, :], axis=-1)
    term4 = rs[..., None, None, :] * lpv[..., :, :, None]
    # Fifth term: np.einsum('xyz,ix,jy,kz', e, lv, pv, rv). We want to avoid expanding
    # the einsum into a 6D tensor to avoid excessive memory usage. Instead, we compute
    # the cross product between lv and pv and then compute the dot product with rv.
    # First compute cross products between lv and pv
    lv_expanded = lv[..., :, None, :]
    pv_expanded = pv[..., None, :, :]
    cross_lp = xp.linalg.cross(lv_expanded, pv_expanded)
    # Then compute dot product with rv
    term5 = xp.sum(cross_lp[..., :, :, None, :] * rv[..., None, None, :, :], axis=-1)
    # Combine all terms with proper shape alignment
    qs = xp.abs(term1 - term2 - term3 - term4 - term5)
    qs = xp.reshape(xp.moveaxis(qs, 1, 0), (qs.shape[1], -1))

    # Find best indices from scalar components
    max_ind = xp.argmax(xp.reshape(qs, (qs.shape[0], -1)), axis=1)
    left_best = max_ind // rv.shape[0]
    right_best = max_ind % rv.shape[0]
    # Array API limitation: Integer index arrays are only allowed with integer indices
    # TODO: Can we somehow avoid this?
    all_idx = xp.reshape(xp.arange(left.shape[-1]), (1, -1))
    left_idx = xp.reshape(left_best, (-1, 1))
    left = left[left_idx, all_idx]
    right_idx = xp.reshape(right_best, (-1, 1))
    right = right[right_idx, all_idx]

    # Reduce the rotation using the best indices
    reduced = compose_quat(left, compose_quat(p, right))

    if left is None:
        left_best = None
    if right is None:
        right_best = None
    return reduced, left_best, right_best


def apply(quat: Array, points: Array, inverse: bool = False) -> Array:
    xp = array_namespace(quat)
    mat = as_matrix(quat)
    # We do not have access to einsum. To avoid broadcasting issues, we add a singleton
    # dimension to the points array and remove it after the operation.
    points = points[..., None]
    if not broadcastable(mat.shape, points.shape):
        raise ValueError(
            f"Cannot broadcast {quat.shape[:-1]} rotations to {points.shape[:-1]} "
            "vectors."
        )
    if inverse:
        # TODO: Replace with .mT once numpy 2.0 is the minimum supported version
        return (xp.matrix_transpose(mat) @ points)[..., 0]
    return (mat @ points)[..., 0]


def setitem(
    quat: Array, value: Array, indexer: int | slice | EllipsisType | None
) -> Array:
    return xpx.at(quat)[indexer, ...].set(value)


def align_vectors(
    a: Array, b: Array, weights: Array | None = None, return_sensitivity: bool = False
) -> tuple[Array, Array, Array]:
    xp = array_namespace(a)
    # Check input vectors
    dtype = xp_result_type(a, b, force_floating=True, xp=xp)
    a_original = xp.asarray(a, dtype=dtype)
    b_original = xp.asarray(b, dtype=dtype)
    a = xpx.atleast_nd(a_original, ndim=2, xp=xp)
    b = xpx.atleast_nd(b_original, ndim=2, xp=xp)
    if a.shape[-1] != 3:
        raise ValueError(
            f"Expected input `a` to have shape (3,) or (N, 3), got {a_original.shape}"
        )
    if b.shape[-1] != 3:
        raise ValueError(
            f"Expected input `b` to have shape (3,) or (N, 3), got {b_original.shape}"
        )
    if a.shape != b.shape:
        raise ValueError(
            "Expected inputs `a` and `b` to have same shapes"
            f", got {a_original.shape} and {b_original.shape} respectively."
        )
    if a.ndim > 2 or b.ndim > 2:  # This function does not support broadcasting
        raise ValueError(
            "Expected inputs `a` and `b` to have shape (3,) or (N, 3), got "
            f"{a_original.shape} and {b_original.shape} respectively."
        )
    N = a.shape[0]

    # Check weights
    if weights is None:
        weights = xp.ones(N, device=xp_device(a), dtype=a.dtype)
    else:
        weights = xp.asarray(weights, device=xp_device(a), dtype=a.dtype)
        if weights.ndim != 1:
            raise ValueError(
                f"Expected `weights` to be 1 dimensional, got shape {weights.shape}."
            )
        if N > 1 and (weights.shape[0] != N):
            raise ValueError(
                "Expected `weights` to have number of values equal to number of input "
                f"vectors, got {weights.shape[0]} values and {N} vectors."
            )
        # We can only check for negative weights in eager execution models. Lazy
        # backends will return NaNs instead.
        negative_weights = weights < 0
        if not is_lazy_array(negative_weights) and xp.any(negative_weights):
            raise ValueError("`weights` may not contain negative values")
        weights = xp.where(negative_weights, xp.nan, weights)

    # For the special case of a single vector pair, we use the infinite
    # weight code path
    weight_is_inf = xp.asarray([True]) if N == 1 else weights == xp.inf
    n_inf = xp.sum(xp.astype(weight_is_inf, a.dtype))
    # We can only error out on multiple infinite weights or sensitivity return with
    # infinite weights in eager execution models. Lazy backends will return NaNs.
    if not is_lazy_array(n_inf):
        if n_inf > 1:
            raise ValueError("Only one infinite weight is allowed")
        if n_inf == 1 and return_sensitivity:
            raise ValueError(
                "Cannot return sensitivity matrix with an "
                "infinite weight or one vector pair"
            )

    weights = xp.where(n_inf > 1, xp.nan, weights)

    inf_branch = xp.any(weight_is_inf, axis=-1)
    # DECISION: We cannot compute both branches for all frameworks. There are two main
    # reasons:
    # 1. Computing both for eager execution models is expensive.
    # 2. Some operations will fail when running the unused branch because of numerical
    # and algorithmical issues. Numpy e.g. will raise an exception when trying to
    # compute the svd of a matrix with infinite weights. To prevent this, we only
    # compute the branch that is needed. Lazy backends however require us to take the
    # full compute graph. Therefore, we use xp.where for lazy backends and a branching
    # version for eager frameworks.
    #
    # Note that we could also solve this by exploiting the externals of xpx.apply_where.
    # However, we'd have to rely on the implementation details of apply_where, which is
    # something we should avoid.
    # See https://github.com/scipy/scipy/pull/22777#discussion_r2028868364
    if is_lazy_array(inf_branch):
        q_opt, rssd, sensitivity = _align_vectors(a, b, weights)
        q_opt_inf, rssd_inf, sensitivity_inf = _align_vectors_fixed(a, b, weights)
        q_opt = xp.where(inf_branch, q_opt_inf, q_opt)
        rssd = xp.where(inf_branch, rssd_inf, rssd)
        sensitivity = xp.where(inf_branch, sensitivity_inf, sensitivity)
    else:
        if xp.any(inf_branch):
            q_opt, rssd, sensitivity = _align_vectors_fixed(a, b, weights)
        else:
            q_opt, rssd, sensitivity = _align_vectors(a, b, weights)
    return q_opt, rssd, sensitivity


def _align_vectors(a: Array, b: Array, weights: Array) -> tuple[Array, Array, Array]:
    xp = array_namespace(a)
    device = xp_device(a)
    B = (weights[:, None] * a).mT @ b
    u, s, vh = xp.linalg.svd(B)

    # Correct improper rotation if necessary (as in Kabsch algorithm)
    neg_det = xp.linalg.det(u @ vh) < 0
    s = xpx.at(s)[..., -1].set(xp.where(neg_det, -s[..., -1], s[..., -1]))
    u = xpx.at(u)[..., :, -1].set(xp.where(neg_det, -u[..., :, -1], u[..., :, -1]))

    C = u @ vh

    # DECISION: We cannot branch on the condition because jit code needs to be
    # non-branching. Hence, we omit the check for uniqueness
    # (s[1] + s[2] < 1e-16 * s[0])
    ssd = xp.sum(weights * xp.sum(b**2 + a**2, axis=-1), axis=-1) - 2 * xp.sum(
        s, axis=-1
    )
    rssd = xp.sqrt(xp.maximum(ssd, xp.zeros(1, device=device)))[..., 0]

    # TODO: We currently need to always compute the sensitivity matrix because lazy code
    # needs to be non-branching. We should check if compilers can optimize the where
    # statement (e.g. in jax) and check if we can have an eager version that only
    # evaluates the branch that is needed.
    # See xpx.apply_where, issue: https://github.com/data-apis/array-api-extra/pull/141
    zeta = (s[..., 0] + s[..., 1]) * (s[..., 1] + s[..., 2]) * (s[..., 2] + s[..., 0])
    kappa = s[..., 0] * s[..., 1] + s[..., 1] * s[..., 2] + s[..., 2] * s[..., 0]
    eye = xp.eye(3, dtype=a.dtype, device=device)
    sensitivity = xp.mean(weights) / zeta * (kappa * eye + B @ B.mT)
    q_opt = _from_matrix_orthogonal(C)
    return q_opt, rssd, sensitivity


def _align_vectors_fixed(
    a: Array, b: Array, weights: Array
) -> tuple[Array, Array, Array]:
    xp = array_namespace(a)
    device = xp_device(a)
    N = a.shape[0]
    weight_is_inf = xp.asarray([True], device=device) if N == 1 else weights == xp.inf
    # We cannot use boolean masks for indexing because of jax. For the same reason, we
    # also cannot use dynamic slices. As a workaround, we roll the array so that the
    # infinitely weighted vector is at index 0. We then use static slices to get the
    # primary and secondary vectors.
    #
    # Note that argmax fulfils a secondary purpose here:
    # When we trace this function with jax, this function might get executed even if
    # weight_is_inf does not have a single valid entry. Argmax will still give us a
    # valid index which allows us to proceed with the function (the results are
    # discarded anyways), whereas boolean indexing would result in invalid, zero-shaped
    # arrays.

    inf_idx = xp.argmax(xp.astype(weight_is_inf, xp.uint8))
    # xp.argmax returns an Array, but the specification of xp.roll does not support
    # Arrays as shifts. For lazy execution models we cannot convert to int because this
    # raises a concretization error. However, jax does accept Arrays as shifts which
    # allows jit compiling the function. Hence we do not convert to int for jax.
    # See https://github.com/data-apis/array-api/issues/914#issuecomment-2918959918
    if not is_jax(xp):
        inf_idx = int(inf_idx)
    a_sorted = xp.roll(a, shift=-inf_idx, axis=0)
    b_sorted = xp.roll(b, shift=-inf_idx, axis=0)
    weights_sorted = xp.roll(weights, shift=-inf_idx, axis=0)

    a_pri = a_sorted[0, ...][None, ...]  # Effectively [[0], ...]
    b_pri = b_sorted[0, ...][None, ...]
    a_pri_norm = xp_vector_norm(a_pri, axis=-1, keepdims=True, xp=xp)
    b_pri_norm = xp_vector_norm(b_pri, axis=-1, keepdims=True, xp=xp)
    if not is_lazy_array(a_pri_norm):
        if xp.any(a_pri_norm == 0) or xp.any(b_pri_norm == 0):
            raise ValueError("Cannot align zero length primary vectors")

    # We cannot error out on zero length vectors. We set the norm to NaN to avoid
    # division by zero and mark the result as invalid.
    a_pri_norm = xpx.at(a_pri_norm)[a_pri_norm == 0].set(xp.nan)
    b_pri_norm = xpx.at(b_pri_norm)[b_pri_norm == 0].set(xp.nan)

    a_pri = a_pri / a_pri_norm
    b_pri = b_pri / b_pri_norm

    # We first find the minimum angle rotation between the primary
    # vectors.
    cross = xp.linalg.cross(b_pri[..., 0, :], a_pri[..., 0, :])
    cross_norm = xp_vector_norm(cross, axis=-1, xp=xp)
    theta = xp.atan2(cross_norm, xp.sum(a_pri[..., 0, :] * b_pri[..., 0, :], axis=-1))
    tolerance = 1e-3  # tolerance for small angle approximation (rad)
    q_flip = xp.asarray([0.0, 0.0, 0.0, 1.0], dtype=a_pri.dtype, device=device)

    # Near pi radians, the Taylor series approximation of x/sin(x) diverges, so for
    # numerical stability we flip pi and then rotate back by the small angle pi - theta
    flip = xp.pi - theta < tolerance
    # For antiparallel vectors, cross = [0, 0, 0] so we need to manually set an
    # arbitrary orthogonal axis of rotation
    i = xp.argmin(xp.abs(a_pri[..., 0, :]))
    # TODO: Array API does not support fancy indexing __setitem__. The code is
    # equivalent to doing the following:
    # r_components = xp.zeros(3)
    # r_components = xpx.at(r_components)[i - 1].set(a_pri[0, i - 2])
    # r_components = xpx.at(r_components)[i - 2].set(-a_pri[0, i - 1])
    r_components = xp.stack(
        [
            xp.where(i == 1, a_pri[0, 2], xp.where(i == 2, -a_pri[0, 1], 0)),
            xp.where(i == 2, a_pri[0, 0], xp.where(i == 0, -a_pri[0, 2], 0)),
            xp.where(i == 0, a_pri[0, 1], xp.where(i == 1, -a_pri[0, 0], 0)),
        ]
    )
    r = xp.where(cross_norm == 0, r_components, cross)

    q_flip = xp.where(flip, from_rotvec(r / xp_vector_norm(r, xp=xp) * xp.pi), q_flip)
    theta = xp.where(flip, xp.pi - theta, theta)
    cross = xp.where(flip, -cross, cross)

    # Small angle Taylor series approximation for numerical stability
    theta2 = theta * theta
    small_scale = xp.abs(theta) < tolerance
    r_small_scale = cross * (1 + theta2 / 6 + theta2 * theta2 * 7 / 360)
    # We need to handle the case where theta is 0 to avoid division by zero. We use the
    # value of the Taylor series approximation, but non-branching operations require
    # that we still divide by the angle. Since we do not use the result where the angle
    # is close to 0, this is safe.
    theta = theta + xp.asarray(small_scale, dtype=theta.dtype)
    r_large_scale = cross * theta / xp.sin(theta)
    r = xp.where(small_scale, r_small_scale, r_large_scale)
    q_pri = compose_quat(from_rotvec(r), q_flip)

    # Case 1): No secondary vectors, q_opt is q_pri -> Immediately done
    # Case 2): Secondary vectors exist
    # We cannot use boolean masks here because of jax
    a_sec = a_sorted[1:, ...]
    b_sec = b_sorted[1:, ...]
    weights_sec = weights_sorted[1:]

    # We apply the first rotation to the b vectors to align the
    # primary vectors, resulting in vectors c. The secondary
    # vectors must now be rotated about that primary vector to best
    # align c to a.
    c_sec = apply(q_pri, b_sec)

    # Calculate vector components for the angle calculation. The
    # angle phi to rotate a single 2D vector C to align to 2D
    # vector A in the xy plane can be found with the equation
    # phi = atan2(cross(C, A), dot(C, A))
    #     = atan2(|C|*|A|*sin(phi), |C|*|A|*cos(phi))
    # The below equations perform the same operation, but with the
    # 3D rotation restricted to the 2D plane normal to a_pri, where
    # the secondary vectors are projected into that plane. We then
    # find the composite angle from the weighted sum of the
    # axial components in that plane, minimizing the 2D alignment
    # problem.
    sin_term = xp.linalg.vecdot(xp.linalg.cross(c_sec, a_sec), a_pri)
    cos_term = xp.linalg.vecdot(c_sec, a_sec) - (
        xp.linalg.vecdot(c_sec, a_pri) * xp.linalg.vecdot(a_sec, a_pri)
    )
    phi = xp.atan2(xp.sum(weights_sec * sin_term), xp.sum(weights_sec * cos_term))
    q_sec = from_rotvec(phi * a_pri[0, ...])

    # Compose these to get the optimal rotation
    q_opt = q_pri if N == 1 else compose_quat(q_sec, q_pri)

    # Calculate the root sum squared distance. We force the error to
    # be zero for the infinite weight vectors since they will align
    # exactly.
    mask = ((N > 1) | (weights[0] == xp.inf)) & weight_is_inf
    # Skip non-infinite weight single vectors pairs, we used the
    # infinite weight code path but don't want to zero that weight
    weights_inf_zero = xpx.at(weights, mask).set(0, copy=True, xp=xp)
    a_est = apply(q_opt, b)
    rssd = xp.sqrt(xp.sum(weights_inf_zero @ (a - a_est) ** 2))

    mask = xp.any(xp.isnan(weights), axis=-1)
    q_opt = xp.where(mask, xp.nan, q_opt)
    return q_opt, rssd, xp.asarray(xp.nan, device=device)


def pow(quat: Array, n: float | Array) -> Array:
    xp = array_namespace(quat)
    device = xp_device(quat)
    # If n is an array, we sanitize it to a scalar and promote quat and n to
    # the same dtype.
    if is_array_api_obj(n):
        if n.shape == (1,):
            n = n[0]
        elif n.ndim != 0:
            raise ValueError("Array exponent must be a scalar")
        quat, n = xp_promote(quat, n, force_floating=True, xp=xp)

    # If n is a lazy array, we cannot take fast paths for special cases.
    if is_lazy_array(n):
        result = from_rotvec(n * as_rotvec(quat))  # general scaling of rotation angle
        # Special cases 0 -> identity, -1 -> inv, 1 -> copy
        identity = xp.zeros((*quat.shape[:-1], 4), dtype=quat.dtype, device=device)
        identity = xpx.at(identity)[..., 3].set(1)
        result = xp.where(n == 0, identity, result)
        result = xp.where(n == -1, inv(quat), result)
        result = xp.where(n == 1, quat, result)
        return result
    if n == 0:
        identity = xp.zeros((*quat.shape[:-1], 4), dtype=quat.dtype, device=device)
        return xpx.at(identity)[..., 3].set(1)
    if n == -1:
        return inv(quat)
    if n == 1:
        return quat
    return from_rotvec(n * as_rotvec(quat))


def _normalize_quaternion(quat: Array) -> Array:
    xp = array_namespace(quat)
    quat_norm = xp_vector_norm(quat, axis=-1, keepdims=True, xp=xp)
    zero_norm = quat_norm == 0
    if is_lazy_array(quat_norm):
        quat = xp.where(zero_norm, xp.nan, quat)
    elif xp.any(zero_norm):
        raise ValueError("Found zero norm quaternions in `quat`.")
    return quat / quat_norm


def _quat_canonical(quat: Array) -> Array:
    xp = array_namespace(quat)
    mask = quat[..., 3] < 0
    zero_w = quat[..., 3] == 0
    mask = mask | (zero_w & (quat[..., 0] < 0))
    zero_wx = zero_w & (quat[..., 0] == 0)
    mask = mask | (zero_wx & (quat[..., 1] < 0))
    zero_wxy = zero_wx & (quat[..., 1] == 0)
    mask = mask | (zero_wxy & (quat[..., 2] < 0))
    return xp.where(mask[..., None], -quat, quat)


def _elementary_basis_index(axis: str) -> int:
    if axis == "x":
        return 0
    elif axis == "y":
        return 1
    elif axis == "z":
        return 2
    raise ValueError(f"Expected axis to be from ['x', 'y', 'z'], got {axis}")


def _compute_davenport_from_quat(
    quat: Array,
    n1: Array,
    n2: Array,
    n3: Array,
    extrinsic: bool,
    suppress_warnings: bool,
) -> Array:
    # The algorithm assumes extrinsic frame transformations. The algorithm
    # in the paper is formulated for rotation quaternions, which are stored
    # directly by Rotation.
    # Adapt the algorithm for our case by reversing both axis sequence and
    # angles for intrinsic rotations when needed
    xp = array_namespace(quat)
    n1, n3 = (n1, n3) if extrinsic else (n3, n1)

    n_cross = xp.linalg.cross(n1, n2)
    lamb = xp.atan2(xp.vecdot(n3, n_cross), xp.vecdot(n3, n1))

    # alternative set of angles compatible with as_euler implementation
    mask = lamb < 0
    n2 = xp.where(mask[..., None], -n2, n2)
    lamb = xp.where(mask, -lamb, lamb)
    n_cross = xp.where(mask[..., None], -n_cross, n_cross)
    correct_set = mask

    quat_lamb = xp.concat(
        (xp.sin(lamb / 2)[..., None] * n2, xp.cos(lamb / 2)[..., None]), axis=-1
    )

    q_trans = compose_quat(quat_lamb, quat)
    a = q_trans[..., 3]
    b = xp.linalg.vecdot(q_trans[..., :3], n1)
    c = xp.linalg.vecdot(q_trans[..., :3], n2)
    d = xp.linalg.vecdot(q_trans[..., :3], n_cross)

    angles = _get_angles(extrinsic, False, 1, lamb, a, b, c, d, suppress_warnings)
    angles = xpx.at(angles)[..., 1].set(
        xp.where(correct_set, -angles[..., 1], angles[..., 1])
    )

    return angles


def _elementary_quat_compose(axes: list[int], angles: Array, intrinsic: bool) -> Array:
    xp = array_namespace(angles)
    device = xp_device(angles)
    quat = _make_elementary_quat(axes[0], angles[..., 0], device=device, xp=xp)
    for i in range(1, len(axes)):
        ax_quat = _make_elementary_quat(axes[i], angles[..., i], device=device, xp=xp)
        quat = compose_quat(quat, ax_quat) if intrinsic else compose_quat(ax_quat, quat)
    return quat


def _make_elementary_quat(axis: int, angle: Array, device, xp) -> Array:
    quat = xp.zeros((*angle.shape, 4), dtype=angle.dtype, device=device)
    quat = xpx.at(quat)[..., 3].set(xp.cos(angle / 2.0))
    quat = xpx.at(quat)[..., axis].set(xp.sin(angle / 2.0))
    return quat


def _get_angles(
    extrinsic: bool,
    symmetric: bool,
    sign: int,
    lamb: float,
    a: Array,
    b: Array,
    c: Array,
    d: Array,
    suppress_warnings: bool,
) -> Array:
    xp = array_namespace(a)
    device = xp_device(a)
    eps = 1e-7
    half_sum = xp.atan2(b, a)
    half_diff = xp.atan2(d, c)
    # We zero-initialize to automatically cover singular cases where the second angle is
    # not defined uniquely.
    angles = xp.zeros((*a.shape, 3), dtype=a.dtype, device=device)

    angles = xpx.at(angles)[..., 1].set(2 * xp.atan2(xp.hypot(c, d), xp.hypot(a, b)))

    angle_first = 0 if extrinsic else 2
    angle_third = 2 if extrinsic else 0

    # Check if the second angle is close to 0 or pi, causing a singularity.
    # - Case 0: Second angle is neither close to 0 nor pi.
    # - Case 1: Second angle is close to 0.
    # - Case 2: Second angle is close to pi.
    case1 = xp.abs(angles[..., 1]) <= eps
    case2 = xp.abs(angles[..., 1] - xp.pi) <= eps
    case0 = ~(case1 | case2)
    if not suppress_warnings and not is_lazy_array(case0) and xp.any(~case0):
        warnings.warn(
            "Gimbal lock detected. Setting third angle to zero "
            "since it is not possible to uniquely determine "
            "all angles.",
            stacklevel=3,
        )

    # This writes case1 into a0 where True and case2 everywhere else. This is sound
    # because we later overwrite any values without singularity with the regular value
    # of case0. The second angle is covered by default since we zero-initialized the
    # second dimension.
    a0 = xp.where(case1, 2 * half_sum, 2 * half_diff * (-1 if extrinsic else 1))
    angles = xpx.at(angles)[..., 0].set(a0)

    # We overwrite the values of angles without singularities (case0)
    a1 = xp.where(case0, half_sum - half_diff, angles[..., angle_first])
    angles = xpx.at(angles)[..., angle_first].set(a1)

    # Same as above but for the third angle. We overwrite the non-singular case0 values
    a3 = xp.where(case0, half_sum + half_diff, angles[..., angle_third])
    if not symmetric:
        a3 = a3 * sign
        angles = xpx.at(angles)[..., 1].set(angles[..., 1] - lamb)
    angles = xpx.at(angles)[..., angle_third].set(a3)

    angles = (angles + xp.pi) % (2 * xp.pi) - xp.pi
    return angles


def compose_quat(p: Array, q: Array) -> Array:
    xp = array_namespace(p)
    cross = xp.linalg.cross(p[..., :3], q[..., :3])
    qx = p[..., 3] * q[..., 0] + q[..., 3] * p[..., 0] + cross[..., 0]
    qy = p[..., 3] * q[..., 1] + q[..., 3] * p[..., 1] + cross[..., 1]
    qz = p[..., 3] * q[..., 2] + q[..., 3] * p[..., 2] + cross[..., 2]
    qw = (
        p[..., 3] * q[..., 3]
        - p[..., 0] * q[..., 0]
        - p[..., 1] * q[..., 1]
        - p[..., 2] * q[..., 2]
    )
    quat = xp.stack([qx, qy, qz, qw], axis=-1)
    return quat


def _split_rotation(q: Array, xp) -> tuple[Array, Array]:
    q = xpx.atleast_nd(q, ndim=2, xp=xp)
    return q[..., -1], q[..., :-1]


def _deg2rad(angles: Array) -> Array:
    return angles * (np.pi / 180.0)


def _rad2deg(angles: Array) -> Array:
    return angles * (180.0 / np.pi)
