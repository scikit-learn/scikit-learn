import pickle
from itertools import product

import pytest

import numpy as np
from scipy.spatial.transform import Rotation, RigidTransform
from scipy.spatial.transform._rigid_transform import normalize_dual_quaternion
from scipy._lib._array_api import (
    is_lazy_array,
    xp_vector_norm,
    is_numpy,
    xp_assert_close,
    make_xp_test_case,
    xp_assert_equal,
    xp_promote
)
import scipy._lib.array_api_extra as xpx

lazy_xp_modules = [RigidTransform]


def rotation_to_xp(r: Rotation, xp):
    dtype = xpx.default_dtype(xp)
    return Rotation.from_quat(xp.asarray(r.as_quat(), dtype=dtype))


def rigid_transform_to_xp(r: RigidTransform, xp):
    dtype = xpx.default_dtype(xp)
    return RigidTransform.from_matrix(xp.asarray(r.as_matrix(), dtype=dtype))


@make_xp_test_case(RigidTransform.as_matrix)
def test_repr(xp):
    actual = repr(RigidTransform.from_matrix(xp.eye(4)))
    expected = """\
RigidTransform.from_matrix(array([[1., 0., 0., 0.],
                                  [0., 1., 0., 0.],
                                  [0., 0., 1., 0.],
                                  [0., 0., 0., 1.]]))"""
    if is_numpy(xp):
        assert actual == expected
    else:
        assert actual.startswith("RigidTransform.from_matrix(")

    tf = RigidTransform.from_matrix(xp.asarray(RigidTransform.identity(2).as_matrix()))
    actual = repr(tf)
    expected = """\
RigidTransform.from_matrix(array([[[1., 0., 0., 0.],
                                   [0., 1., 0., 0.],
                                   [0., 0., 1., 0.],
                                   [0., 0., 0., 1.]],
                           
                                  [[1., 0., 0., 0.],
                                   [0., 1., 0., 0.],
                                   [0., 0., 1., 0.],
                                   [0., 0., 0., 1.]]]))"""
    if is_numpy(xp):
        assert actual == expected
    else:
        assert actual.startswith("RigidTransform.from_matrix(")


@make_xp_test_case(RigidTransform.from_rotation)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_from_rotation(xp, ndim: int):
    atol = 1e-12
    rng = np.random.default_rng(0)
    shape = (ndim,) * (ndim - 1) + (4,)
    r = rotation_to_xp(Rotation.from_quat(rng.normal(size=shape)), xp=xp)
    tf = RigidTransform.from_rotation(r)
    xp_assert_close(tf.as_matrix()[..., :3, :3], r.as_matrix(), atol=atol)
    xp_assert_close(tf.as_matrix()[..., :3, 3], xp.zeros(shape[:-1] + (3,)), atol=atol)
    xp_assert_close(tf.as_matrix()[..., 3, :3], xp.zeros(shape[:-1] + (3,)), atol=atol)
    xp_assert_close(tf.as_matrix()[..., 3, 3], xp.ones(shape[:-1]), atol=atol)
    assert tf.single == (ndim == 1)


@make_xp_test_case(RigidTransform.from_translation)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_from_translation(xp, ndim: int):
    shape = (ndim,) * (ndim - 1)
    t = xp.reshape(xp.arange(ndim ** (ndim-1) * 3), shape + (3,))
    tf = RigidTransform.from_translation(t)

    expected = xp.tile(xp.eye(4), shape + (1, 1))
    t_float = xp_promote(t, force_floating=True, xp=xp)
    expected = xpx.at(expected)[..., :3, 3].set(t_float)
    xp_assert_close(tf.as_matrix(), expected)
    assert tf.single == (ndim == 1)


def test_from_translation_array_like():
    # Test single translation
    t = [1, 2, 3]
    tf = RigidTransform.from_translation(t)
    tf_expected = RigidTransform.from_translation(np.array(t))
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix())
    assert tf.single

    # Test multiple translations
    t = [[1, 2, 3], [4, 5, 6]]
    tf = RigidTransform.from_translation(t)
    tf_expected = RigidTransform.from_translation(np.array(t))
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix())
    assert not tf.single


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.as_matrix)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_from_matrix(xp, ndim: int):
    atol = 1e-12
    shape = (ndim,) * (ndim - 1) + (4, 4)
    dtype = xpx.default_dtype(xp)

    matrix = xp.tile(xp.eye(4), shape[:-2] + (1, 1))
    t = xp.reshape(xp.arange(ndim ** (ndim-1) * 3, dtype=dtype), shape[:-2] + (3,))
    matrix = xpx.at(matrix)[..., :3, 3].set(t)

    tf = RigidTransform.from_matrix(matrix)
    xp_assert_close(tf.as_matrix(), matrix, atol=atol)
    assert tf.single == (ndim == 1)

    # Test non-1 determinant
    matrix = xp.tile(xp.eye(4), shape[:-2] + (1, 1))
    matrix = xpx.at(matrix)[..., :3, :3].set(xp.eye(3) * 2.0)
    tf = RigidTransform.from_matrix(matrix)
    expected = xp.tile(xp.eye(4), shape[:-2] + (1, 1))
    xp_assert_close(tf.as_matrix(), expected, atol=atol)

    # Test non-orthogonal rotation matrix
    matrix = xp.tile(xp.eye(4), shape[:-2] + (1, 1))
    # matrix is equivalent to [[1, 1, 0, 0],
    #                          [0, 1, 0, 0],
    #                          [0, 0, 1, 0],
    #                          [0, 0, 0, 1]]
    matrix = xpx.at(matrix)[..., 0, 1].set(1.0)
    tf = RigidTransform.from_matrix(matrix)
    expected = xp.tile(xp.eye(4), shape[:-2] + (1, 1))
    expected = xpx.at(expected)[..., 0, 0].set(0.894427)
    expected = xpx.at(expected)[..., 0, 1].set(0.447214)
    expected = xpx.at(expected)[..., 1, 0].set(-0.447214)
    expected = xpx.at(expected)[..., 1, 1].set(0.894427)
    xp_assert_close(tf.as_matrix(), expected, atol=1e-6)

    # Test invalid matrix
    invalid = xp.tile(xp.eye(4), shape[:-2] + (1, 1))
    invalid = xpx.at(invalid)[..., 3, 3].set(2)  # Invalid last row
    if is_lazy_array(invalid):
        tf = RigidTransform.from_matrix(invalid)
        assert xp.all(xp.isnan(tf.as_matrix()))
    else:
        with pytest.raises(ValueError):
            RigidTransform.from_matrix(invalid)


def test_from_matrix_array_like():
    # Test single transform matrix
    matrix = [[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]]
    expected = np.eye(4)
    tf = RigidTransform.from_matrix(matrix)
    xp_assert_close(tf.as_matrix(), expected)
    assert tf.single

    # Test multiple transform matrices
    matrices = [matrix, matrix]
    tf = RigidTransform.from_matrix(matrices)
    for i in range(len(matrices)):
        xp_assert_close(tf.as_matrix()[i, ...], expected)
    assert not tf.single


@make_xp_test_case(RigidTransform.from_components)
@pytest.mark.parametrize("r_ndim", range(1, 3))
@pytest.mark.parametrize("t_ndim", range(1, 3))
def test_from_components(xp, r_ndim: int, t_ndim: int):
    atol = 1e-12
    dims = (6, 5, 4, 3)  # Common shape
    q_shape = dims[:r_ndim - 1][::-1] + (4,)
    t_shape = dims[:t_ndim - 1][::-1] + (3,)
    tf_shape = np.broadcast_shapes(q_shape[:-1], t_shape[:-1]) + (4, 4)
    rng = np.random.default_rng(0)

    t = xp.reshape(xp.arange(np.prod(t_shape[:-1]) * 3), t_shape)
    r = rotation_to_xp(Rotation.from_quat(rng.random(size=q_shape)), xp=xp)
    tf = RigidTransform.from_components(t, r)

    expected = xp.zeros(tf_shape)
    expected = xpx.at(expected)[..., :3, :3].set(r.as_matrix())
    t_float = xp_promote(t, force_floating=True, xp=xp)
    expected = xpx.at(expected)[..., :3, 3].set(t_float)
    expected = xpx.at(expected)[..., 3, 3].set(1)
    xp_assert_close(tf.as_matrix(), expected, atol=atol)
    assert tf.single == (r_ndim == 1 and t_ndim == 1)


def test_from_components_array_like():
    rng = np.random.default_rng(123)
    # Test single rotation and translation
    t = [1, 2, 3]
    r = Rotation.random(rng=rng)
    tf = RigidTransform.from_components(t, r)
    tf_expected = RigidTransform.from_components(np.array(t), r)
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix(), atol=1e-12)
    assert tf.single

    # Test multiple rotations and translations
    t = [[1, 2, 3], [4, 5, 6]]
    r = Rotation.random(len(t), rng=rng)
    tf = RigidTransform.from_components(t, r)
    tf_expected = RigidTransform.from_components(np.array(t), r)
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix(), atol=1e-12)
    assert not tf.single


@make_xp_test_case(RigidTransform.as_components)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_as_components(xp, ndim):
    dtype = xpx.default_dtype(xp)
    atol = 1e-12 if dtype == xp.float64 else 1e-6
    shape = (ndim,) * (ndim - 1)
    rng = np.random.default_rng(123)
    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    r = rotation_to_xp(Rotation.from_quat(rng.random(shape + (4,))), xp=xp)
    tf = RigidTransform.from_components(t, r)
    new_t, new_r = tf.as_components()
    assert xp.all(new_r.approx_equal(r, atol=atol))
    xp_assert_close(new_t, t, atol=atol)


@make_xp_test_case(RigidTransform.from_exp_coords)
@pytest.mark.parametrize("dim", range(1, 4))
def test_from_exp_coords(xp, dim: int):
    shape = (dim,) * (dim - 1)
    # example from 3.3 of
    # https://hades.mech.northwestern.edu/images/2/25/MR-v2.pdf
    dtype = xpx.default_dtype(xp)
    angle1 = np.deg2rad(30.0)
    mat = xp.asarray([
        [np.cos(angle1), -np.sin(angle1), 0.0, 1.0],
        [np.sin(angle1), np.cos(angle1), 0.0, 2.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ], dtype=dtype)
    mat = xp.tile(mat, shape + (1, 1))
    tf1 = RigidTransform.from_matrix(mat)
    angle2 = np.deg2rad(60.0)
    mat = xp.asarray([
        [np.cos(angle2), -np.sin(angle2), 0.0, 2.0],
        [np.sin(angle2), np.cos(angle2), 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0]
    ], dtype=dtype)
    mat = xp.tile(mat, shape + (1, 1))
    tf2 = RigidTransform.from_matrix(mat)
    expected = tf2 * tf1.inv()
    deg2rag = xp.asarray(np.deg2rad(30.0))
    exp_coords = deg2rag * xp.asarray([0.0, 0.0, 1.0, 3.37, -3.37, 0.0])
    exp_coords = xp.tile(exp_coords, shape + (1,))
    actual = RigidTransform.from_exp_coords(exp_coords)
    xp_assert_close(actual.as_matrix(), expected.as_matrix(), atol=1e-2)

    # test cases generated by comparison to pytransform3d
    exp_coords = xp.asarray([
        [-2.01041204, -0.52983629, 0.65773501,
         0.10386614, 0.05855009, 0.54959179],
        [-0.22537438, -0.24132627, -2.4747121,
         -0.09158594,  1.88075832, -0.03197204]
    ])
    exp_coords = xp.tile(exp_coords, shape + (1,))
    expected_matrix = xp.asarray([
        [[0.76406621, 0.10504613, -0.63652819, -0.10209961],
         [0.59956454, -0.47987325, 0.64050295, 0.40158789],
         [-0.2381705, -0.87102639, -0.42963687, 0.19637636],
         [0., 0., 0., 1.]],
        [[-0.78446989, 0.61157488, 0.10287448, 1.33330055],
         [-0.58017785, -0.78232107, 0.22664378, 0.52660831],
         [0.21909052, 0.11810973, 0.96852952, -0.02968529],
         [0., 0., 0., 1.]]
    ])
    expected_matrix = xp.tile(expected_matrix, shape + (1, 1))
    xp_assert_close(
        RigidTransform.from_exp_coords(exp_coords).as_matrix(),
        expected_matrix, atol=1e-8)

    # identity
    expected_matrix = xp.tile(xp.eye(4), shape + (1, 1))
    exp_coords = xp.zeros(shape + (6,), dtype=dtype)
    xp_assert_close(
        RigidTransform.from_exp_coords(exp_coords).as_matrix(),
        expected_matrix, atol=1e-12)

    # only translation
    expected_matrix = xp.asarray([
        [[1.0, 0.0, 0.0, 3.0],
         [0.0, 1.0, 0.0, -5.4],
         [0.0, 0.0, 1.0, 100.2],
         [0.0, 0.0, 0.0, 1.0]],
        [[1.0, 0.0, 0.0, -3.0],
         [0.0, 1.0, 0.0, 13.3],
         [0.0, 0.0, 1.0, 1.3],
         [0.0, 0.0, 0.0, 1.0]]
    ])
    expected_matrix = xp.tile(expected_matrix, shape + (1, 1, 1))
    exp_coords = xp.asarray([
        [0.0, 0.0, 0.0, 3.0, -5.4, 100.2],
        [0.0, 0.0, 0.0, -3.0, 13.3, 1.3],
    ])
    exp_coords = xp.tile(exp_coords, shape + (1, 1))
    actual = RigidTransform.from_exp_coords(exp_coords)
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=1e-12)

    # only rotation
    angles = xp.asarray([[34, -12, 0.5], [-102, -55, 30]])
    angles = xp.tile(angles, shape + (1, 1))
    rot = Rotation.from_euler('zyx', angles, degrees=True)
    rotvec = rot.as_rotvec()
    expected_matrix = xp.tile(xp.eye(4), shape + (2, 1, 1))
    expected_matrix = xpx.at(expected_matrix)[..., :3, :3].set(rot.as_matrix())
    exp_coords = xp.concat((rotvec, xp.zeros_like(rotvec)), axis=-1)
    actual = RigidTransform.from_exp_coords(exp_coords)
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=1e-12)


def test_from_exp_coords_array_like():
    rng = np.random.default_rng(123)
    # Test single transform
    t = np.array([1, 2, 3])
    r = Rotation.random(rng=rng)
    tf_expected = RigidTransform.from_components(t, r)
    exp_coords = tf_expected.as_exp_coords().tolist()
    assert isinstance(exp_coords, list)
    tf = RigidTransform.from_exp_coords(exp_coords)
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix(), atol=1e-12)

    # Test multiple transforms
    t = [[1, 2, 3], [4, 5, 6]]
    r = Rotation.random(len(t), rng=rng)
    tf_expected = RigidTransform.from_components(t, r)
    exp_coords = tf_expected.as_exp_coords().tolist()
    assert isinstance(exp_coords, list)
    tf = RigidTransform.from_exp_coords(exp_coords)
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix(), atol=1e-12)


@make_xp_test_case(RigidTransform.as_exp_coords)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_as_exp_coords(xp, ndim: int):
    shape = (ndim,) * (ndim - 1)
    # identity
    expected = xp.zeros(shape + (6,))
    actual = RigidTransform.from_exp_coords(expected).as_exp_coords()
    xp_assert_close(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # pure rotation
    rot_vec = xp.asarray(rng.normal(scale=0.1, size=shape + (1000, 3)))
    tf = RigidTransform.from_rotation(Rotation.from_rotvec(rot_vec))
    exp_coords = tf.as_exp_coords()
    xp_assert_close(exp_coords[..., :3], rot_vec, atol=1e-12)
    expected = xp.zeros_like(rot_vec)
    xp_assert_close(exp_coords[..., 3:], expected, atol=1e-16)

    # pure translation
    translation = xp.asarray(rng.normal(scale=100.0, size=shape + (1000, 3)))
    tf = RigidTransform.from_translation(translation)
    exp_coords = tf.as_exp_coords()
    xp_assert_close(exp_coords[..., :3], expected, atol=1e-16)
    xp_assert_close(exp_coords[..., 3:], translation, atol=1e-15)


@make_xp_test_case(RigidTransform.from_dual_quat)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_from_dual_quat(xp, ndim: int):
    dtype = xpx.default_dtype(xp)
    atol = 1e-12 if dtype == xp.float64 else 1e-7
    shape = (ndim,) * (ndim - 1)

    # identity
    dq = xp.asarray([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=dtype)
    dq = xp.tile(dq, shape + (1,))
    expected = xp.tile(xp.eye(4), shape + (1, 1))
    xp_assert_close(RigidTransform.from_dual_quat(dq).as_matrix(), expected, atol=atol)
    dq = xp.asarray([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], dtype=dtype)
    dq = xp.tile(dq, shape + (1,))
    xp_assert_close(RigidTransform.from_dual_quat(dq, scalar_first=True).as_matrix(),
                    expected, atol=atol)

    # only translation
    dq = xp.asarray([0, 0, 0, 1, 0.25, 0.15, -0.7, 0], dtype=dtype)
    dq = xp.tile(dq, shape + (1,))
    actual = RigidTransform.from_dual_quat(dq)
    expected_matrix = xp.asarray([
        [1, 0, 0, 0.5],
        [0, 1, 0, 0.3],
        [0, 0, 1, -1.4],
        [0, 0, 0, 1]
    ])
    expected_matrix = xp.tile(expected_matrix, shape + (1, 1))
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=atol)
    dq = xp.asarray([1, 0, 0, 0, 0, 0.25, 0.15, -0.7], dtype=dtype)
    dq = xp.tile(dq, shape + (1,))
    actual = RigidTransform.from_dual_quat(dq, scalar_first=True)
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=atol)

    # only rotation
    angles = xp.asarray([65, -13, 90], dtype=dtype)
    angles = xp.tile(angles, shape + (1,))
    actual_rot = Rotation.from_euler("xyz", angles, degrees=True)
    qrot = actual_rot.as_quat()
    dq = xp.concat((qrot, xp.zeros_like(qrot)), axis=-1)
    actual = RigidTransform.from_dual_quat(dq)
    expected_matrix = xp.tile(xp.eye(4), shape + (1, 1))
    expected_matrix = xpx.at(expected_matrix)[..., :3, :3].set(actual_rot.as_matrix())
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=atol)

    qrot = actual_rot.as_quat(scalar_first=True)
    dq = xp.concat((qrot, xp.zeros_like(qrot)), axis=-1)
    actual = RigidTransform.from_dual_quat(dq, scalar_first=True)
    expected_matrix = xp.tile(xp.eye(4), shape + (1, 1))
    expected_matrix = xpx.at(expected_matrix)[..., :3, :3].set(actual_rot.as_matrix())
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=atol)

    # rotation and translation
    # rtol is set to 1e-7 because xp_assert_close deviates from
    # np.testing.assert_allclose in that it does not automatically default to 1e-7 for
    # floating point inputs.
    # See https://numpy.org/doc/2.2/reference/generated/numpy.testing.assert_allclose.html
    dq = xp.asarray(
        [[0.0617101, -0.06483886, 0.31432811, 0.94508498,
          0.04985168, -0.26119618, 0.1691491, -0.07743254],
         [0.19507259, 0.49404931, -0.06091285, 0.8450749,
          0.65049656, -0.30782513, 0.16566752, 0.04174544]])
    dq = xp.tile(dq, shape + (1, 1))
    actual = RigidTransform.from_dual_quat(dq)
    expected_matrix = xp.asarray(
        [[[0.79398752, -0.60213598, -0.08376202, 0.24605262],
          [0.58613113, 0.79477941, -0.15740392, -0.4932833],
          [0.16135089, 0.07588122, 0.98397557, 0.34262676],
          [0., 0., 0., 1.]],
         [[0.50440981, 0.2957028, 0.81125249, 1.20934468],
          [0.08979911, 0.91647262, -0.3898898, -0.70540077],
          [-0.8587822, 0.26951399, 0.43572393, -0.47776265],
          [0., 0., 0., 1.]]])
    expected_matrix = xp.tile(expected_matrix, shape + (1, 1, 1))
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=atol, rtol=1e-7)

    dq = xp.asarray(
        [[0.94508498, 0.0617101, -0.06483886, 0.31432811,
          -0.07743254, 0.04985168, -0.26119618, 0.1691491],
         [0.8450749, 0.19507259, 0.49404931, -0.06091285,
          0.04174544, 0.65049656, -0.30782513, 0.16566752]])
    dq = xp.tile(dq, shape + (1, 1))
    actual = RigidTransform.from_dual_quat(dq, scalar_first=True)
    xp_assert_close(actual.as_matrix(), expected_matrix, atol=atol, rtol=1e-7)

    # unnormalized dual quaternions

    # invalid real quaternion with norm 0
    dq = xp.zeros(shape + (8,))
    actual = RigidTransform.from_dual_quat(dq)
    expected = xp.tile(xp.eye(4), shape + (1, 1))
    xp_assert_close(actual.as_matrix(), expected, atol=atol)

    # real quaternion with norm != 1
    unnormalized_dual_quat = xp.asarray(
        [-0.2547655, 1.23506123, 0.20230088, 0.24247194,  # norm 1.3
         0.38559628, 0.08184063, 0.1755943, -0.1582222]  # orthogonal
    )
    xp_assert_close(xp_vector_norm(unnormalized_dual_quat[:4]), xp.asarray(1.3)[()],
                    atol=atol)
    xp_assert_close(xp.vecdot(unnormalized_dual_quat[:4],
                              unnormalized_dual_quat[4:])[()],
                    xp.asarray(0.0)[()], atol=1e-8)

    dq = xp.tile(unnormalized_dual_quat, shape + (1,))
    dual_quat = RigidTransform.from_dual_quat(dq).as_dual_quat()

    expected_ones = xp.ones(shape) if shape != () else xp.asarray(1.0)[()]
    expected_zeros = xp.zeros(shape) if shape != () else xp.asarray(0.0)[()]
    xp_assert_close(xp_vector_norm(dual_quat[..., :4], axis=-1), expected_ones,
                    atol=1e-12)
    vecdot = xp.vecdot(dual_quat[..., :4], dual_quat[..., 4:])
    vecdot = vecdot[()] if vecdot.shape == () else vecdot
    xp_assert_close(vecdot, expected_zeros, atol=atol)

    # real and dual quaternion are not orthogonal
    unnormalized_dual_quat = xp.asarray(
        [0.20824458, 0.75098079, 0.54542913, -0.30849493,  # unit norm
         -0.16051025, 0.10742978, 0.21277201, 0.20596935]  # not orthogonal
    )
    xp_assert_close(xp_vector_norm(unnormalized_dual_quat[:4]), xp.asarray(1.0)[()],
                    atol=atol)
    assert xp.vecdot(unnormalized_dual_quat[:4], unnormalized_dual_quat[4:]) != 0.0
    dq = xp.tile(unnormalized_dual_quat, shape + (1,))
    dual_quat = RigidTransform.from_dual_quat(dq).as_dual_quat()

    xp_assert_close(xp_vector_norm(dual_quat[..., :4], axis=-1), expected_ones,
                    atol=1e-12)
    vecdot = xp.vecdot(dual_quat[..., :4], dual_quat[..., 4:])
    vecdot = vecdot[()] if vecdot.shape == () else vecdot
    xp_assert_close(vecdot, expected_zeros, atol=atol)

    # invalid real quaternion with norm 0, non-orthogonal dual quaternion
    unnormalized_dual_quat = xp.asarray(
        [0.0, 0.0, 0.0, 0.0, -0.16051025, 0.10742978, 0.21277201, 0.20596935])
    assert xp.vecdot(xp.asarray([0.0, 0, 0, 1]), unnormalized_dual_quat[4:]) != 0.0
    dq = xp.tile(unnormalized_dual_quat, shape + (1,))
    dual_quat = RigidTransform.from_dual_quat(dq).as_dual_quat()

    xp_assert_close(xp_vector_norm(dual_quat[..., :4], axis=-1), expected_ones,
                    atol=1e-12)
    vecdot = xp.vecdot(dual_quat[..., :4], dual_quat[..., 4:])
    vecdot = vecdot[()] if vecdot.shape == () else vecdot
    xp_assert_close(vecdot, expected_zeros, atol=atol)

    # compensation for precision loss in real quaternion
    rng = np.random.default_rng(1000)
    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    q = xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype)
    r = Rotation.from_quat(q)
    random_dual_quats = RigidTransform.from_components(t, r).as_dual_quat()

    # ensure that random quaternions are not normalized
    random_dual_quats = xpx.at(random_dual_quats)[..., :4].add(0.01)
    assert not xp.any(xpx.isclose(xp_vector_norm(random_dual_quats[..., :4], axis=-1),
                                  1.0, atol=0.0001))
    dual_quat_norm = RigidTransform.from_dual_quat(
        random_dual_quats).as_dual_quat()
    xp_assert_close(xp_vector_norm(dual_quat_norm[..., :4], axis=-1), expected_ones,
                    atol=atol)

    # compensation for precision loss in dual quaternion, results in violation
    # of orthogonality constraint
    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    q = xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype)
    r = Rotation.from_quat(q)
    random_dual_quats = RigidTransform.from_components(t, r).as_dual_quat()

    # ensure that random quaternions are not normalized
    random_dual_quats = xpx.at(random_dual_quats)[..., 4:].add(0.1)
    q_norm = xp.vecdot(random_dual_quats[..., :4], random_dual_quats[..., 4:])
    assert not xp.any(xpx.isclose(q_norm, 0.0, atol=0.0001))
    dual_quat_norm = RigidTransform.from_dual_quat(
        random_dual_quats).as_dual_quat()
    vecdot = xp.vecdot(dual_quat[..., :4], dual_quat[..., 4:])
    vecdot = vecdot[()] if vecdot.shape == () else vecdot
    xp_assert_close(vecdot, expected_zeros, atol=atol)
    xp_assert_close(random_dual_quats[..., :4], dual_quat_norm[..., :4], atol=atol)


def test_from_dual_quat_array_like():
    rng = np.random.default_rng(123)
    # Test single transform
    t = np.array([1, 2, 3])
    r = Rotation.random(rng=rng)
    tf_expected = RigidTransform.from_components(t, r)
    dual_quat = tf_expected.as_dual_quat().tolist()
    assert isinstance(dual_quat, list)
    tf = RigidTransform.from_dual_quat(dual_quat)
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix(), atol=1e-12)

    # Test multiple transforms
    t = [[1, 2, 3], [4, 5, 6]]
    r = Rotation.random(len(t), rng=rng)
    tf_expected = RigidTransform.from_components(t, r)
    dual_quat = tf_expected.as_dual_quat().tolist()
    assert isinstance(dual_quat, list)
    tf = RigidTransform.from_dual_quat(dual_quat)
    xp_assert_close(tf.as_matrix(), tf_expected.as_matrix(), atol=1e-12)


@make_xp_test_case(RigidTransform.as_dual_quat)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_as_dual_quat(xp, ndim: int):
    dtype = xpx.default_dtype(xp)
    shape = (ndim,) * (ndim - 1)
    # identity
    expected = xp.asarray([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0], dtype=dtype)
    actual = rigid_transform_to_xp(RigidTransform.identity(), xp).as_dual_quat()
    xp_assert_close(actual, expected, atol=1e-12)

    expected = xp.asarray([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    tf = rigid_transform_to_xp(RigidTransform.identity(), xp)
    actual = tf.as_dual_quat(scalar_first=True)
    xp_assert_close(actual, expected, atol=1e-12)

    rng = np.random.default_rng(10)

    # only rotation
    for _ in range(10):
        q = xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype)
        real_part = Rotation.from_quat(q).as_quat()
        dual_part = xp.zeros_like(real_part)
        expected = xp.concat((real_part, dual_part), axis=-1)
        actual = RigidTransform.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        actual = actual * xp.sign(actual[..., 0, None])
        expected = expected * xp.sign(expected[..., 0, None])
        xp_assert_close(actual, expected, atol=1e-12)

    # only translation
    for _ in range(10):
        tf = 0.5 * xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
        expected = xp.zeros(shape + (8,), dtype=dtype)
        expected = xpx.at(expected)[..., 3].set(1.0)
        expected = xpx.at(expected)[..., 4:7].set(tf)
        actual = RigidTransform.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        actual = actual * xp.sign(actual[..., 0, None])
        expected = expected * xp.sign(expected[..., 0, None])
        xp_assert_close(actual, expected, atol=1e-12)

    # rotation and translation
    for _ in range(10):
        t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
        r = Rotation.from_quat(xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype))
        expected = RigidTransform.from_components(t, r).as_dual_quat()
        actual = RigidTransform.from_dual_quat(expected).as_dual_quat()
        # because of double cover:
        actual = actual * xp.sign(actual[..., 0, None])
        expected = expected * xp.sign(expected[..., 0, None])
        xp_assert_close(actual, expected, atol=1e-12)


@make_xp_test_case(RigidTransform.from_components, RigidTransform.as_components,
                   RigidTransform.from_exp_coords, RigidTransform.as_exp_coords,
                   RigidTransform.from_matrix, RigidTransform.as_matrix,
                   RigidTransform.from_dual_quat, RigidTransform.as_dual_quat)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_from_as_internal_consistency(xp, ndim: int):
    dtype = xpx.default_dtype(xp)
    atol = 1e-12
    n = 10
    rng = np.random.default_rng(10)
    shape = (n,) + (ndim,) * (ndim - 1)

    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    r = Rotation.from_quat(xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype))
    tf0 = RigidTransform.from_components(t, r)

    tf1 = RigidTransform.from_components(*tf0.as_components())
    xp_assert_close(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_components(tf0.translation, tf0.rotation)
    xp_assert_close(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_exp_coords(tf0.as_exp_coords())
    xp_assert_close(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_matrix(tf0.as_matrix())
    xp_assert_close(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    tf1 = RigidTransform.from_dual_quat(tf0.as_dual_quat())
    xp_assert_close(tf0.as_matrix(), tf1.as_matrix(), atol=atol)

    # exp_coords small rotation
    t = xp.asarray(rng.normal(scale=1000.0, size=shape + (3,)), dtype=dtype)
    rotvec = xp.asarray(rng.normal(scale=1e-10, size=shape + (3,)), dtype=dtype)
    r = Rotation.from_rotvec(rotvec)
    tf0 = RigidTransform.from_components(t, r)
    tf1 = RigidTransform.from_exp_coords(tf0.as_exp_coords())
    xp_assert_close(tf0.as_matrix(), tf1.as_matrix(), atol=atol)


def test_identity():
    # We do not use xp here because identity always returns numpy arrays
    atol = 1e-12
    # Test single identity
    tf = RigidTransform.identity()
    xp_assert_close(tf.as_matrix(), np.eye(4), atol=atol)
    # Test multiple identities
    tf = RigidTransform.identity(5)
    xp_assert_close(tf.as_matrix(), np.array([np.eye(4)] * 5), atol=atol)
    # Test shape
    tf = RigidTransform.identity(shape=3)
    expected = np.tile(np.eye(4), (3, 1, 1))
    xp_assert_close(tf.as_matrix(), expected, atol=atol)
    tf = RigidTransform.identity(shape=(2, 3))
    expected = np.tile(np.eye(4), (2, 3, 1, 1))
    xp_assert_close(tf.as_matrix(), expected, atol=atol)
    # Test errors
    with pytest.raises(ValueError, match="Only one of `num` and `shape` can be."):
        RigidTransform.identity(10, shape=(2, 3))
    with pytest.raises(TypeError, match="takes from 0 to 1 positional arguments"):
        RigidTransform.identity(None, (-1, 3))  # Shape is kwarg only
    with pytest.raises(ValueError, match="`shape` must be an int or a tuple of ints"):
        RigidTransform.identity(shape="invalid")


@make_xp_test_case(RigidTransform.apply)
def test_apply(xp):
    atol = 1e-12
    # Broadcast shape: (6, 5, 4, 2) ( + (3,) for vectors, + (4,) for rotations)
    vector_shapes = [(), (1,), (2,), (1, 2), (5, 1, 2)]
    tf_shapes = [(), (1,), (2,), (1, 2), (4, 2), (1, 4, 2), (5, 4, 2), (6, 5, 4, 2)]
    rng = np.random.default_rng(123)

    for tf_shape, v_shape in product(tf_shapes, vector_shapes):
        # Random rotation and translation
        t = xp.asarray(rng.normal(size=tf_shape + (3,)))
        q = xp.asarray(rng.normal(size=tf_shape + (4,)))
        r = Rotation.from_quat(q)
        tf = RigidTransform.from_components(t, r)

        vecs = xp.asarray(rng.normal(size=v_shape + (3,)))
        expected = t + r.apply(vecs)
        res = tf.apply(vecs)
        assert res.shape == np.broadcast_shapes(tf_shape, v_shape) + (3,)
        xp_assert_close(res, expected, atol=atol)


def test_apply_array_like():
    rng = np.random.default_rng(123)
    # Single vector
    t = np.array([1, 2, 3])
    r = Rotation.random(rng=rng)
    tf = RigidTransform.from_components(t, r)
    vec = [1, 0, 0]
    expected = t + r.apply(vec)
    xp_assert_close(tf.apply(vec), expected, atol=1e-12)

    # Multiple vectors
    t = np.array([[1, 2, 3], [4, 5, 6]])
    r = Rotation.random(len(t), rng=rng)
    tf = RigidTransform.from_components(t, r)
    vec = [[1, 0, 0], [0, 1, 0]]
    expected = t + r.apply(vec)
    xp_assert_close(tf.apply(vec), expected, atol=1e-12)


@make_xp_test_case(RigidTransform.apply)
def test_inverse_apply(xp):
    atol = 1e-12
    # Broadcast shape: (6, 5, 4, 2) ( + (3,) for vectors, + (4,) for rotations)
    vector_shapes = [(), (1,), (2,), (1, 2), (5, 1, 2)]
    tf_shapes = [(), (1,), (2,), (1, 2), (4, 2), (1, 4, 2), (5, 4, 2), (6, 5, 4, 2)]
    rng = np.random.default_rng(123)

    for tf_shape, v_shape in product(tf_shapes, vector_shapes):
        # Random rotation and translation
        t = xp.asarray(rng.normal(size=tf_shape + (3,)))
        q = xp.asarray(rng.normal(size=tf_shape + (4,)))
        r = Rotation.from_quat(q)
        tf = RigidTransform.from_components(t, r)

        vecs = xp.asarray(rng.normal(size=v_shape + (3,)))
        expected = tf.inv().apply(vecs)
        res = tf.apply(vecs, inverse=True)
        assert res.shape == np.broadcast_shapes(tf_shape, v_shape) + (3,)
        xp_assert_close(res, expected, atol=atol)


@make_xp_test_case(RigidTransform.apply)
def test_rotation_alone(xp):
    atol = 1e-12

    r = Rotation.from_euler('z', xp.asarray(90), degrees=True)
    tf = RigidTransform.from_rotation(r)
    vec = xp.asarray([1, 0, 0])
    expected = r.apply(vec)
    xp_assert_close(tf.apply(vec), expected, atol=atol)


@make_xp_test_case(RigidTransform.apply)
def test_translation_alone(xp):
    atol = 1e-12
    t = xp.asarray([1.0, 2, 3])
    tf = RigidTransform.from_translation(t)
    vec = xp.asarray([5.0, 6, 7])
    expected = t + vec
    xp_assert_close(tf.apply(vec), expected, atol=atol)


@make_xp_test_case(RigidTransform.apply, RigidTransform.__mul__)
def test_composition(xp):
    atol = 1e-12
    tf_shapes = [(), (1,), (2,), (1, 2), (4, 2), (5, 4, 2)]
    dtype = xpx.default_dtype(xp)
    rng = np.random.default_rng(123)

    for tf_shape1, tf_shape2 in product(tf_shapes, repeat=2):
        # Random rotation and translation
        t1 = xp.asarray(rng.normal(size=tf_shape1 + (3,)))
        q1 = xp.asarray(rng.normal(size=tf_shape1 + (4,)))
        r1 = Rotation.from_quat(q1)
        tf1 = RigidTransform.from_components(t1, r1)

        t2 = xp.asarray(rng.normal(size=tf_shape2 + (3,)))
        q2 = xp.asarray(rng.normal(size=tf_shape2 + (4,)))
        r2 = Rotation.from_quat(q2)
        tf2 = RigidTransform.from_components(t2, r2)

        composed = tf2 * tf1
        vec = xp.asarray(rng.normal(size=(3,)), dtype=dtype)
        expected = tf2.apply(tf1.apply(vec))
        res = composed.apply(vec)
        assert res.shape == np.broadcast_shapes(tf_shape1, tf_shape2) + (3,)
        xp_assert_close(res, expected, atol=atol)

        expected = t2 + r2.apply(t1 + r1.apply(vec))
        xp_assert_close(composed.apply(vec), expected, atol=atol)
        assert composed.single == (tf1.single and tf2.single)


@make_xp_test_case(RigidTransform.__pow__, RigidTransform.__mul__)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_pow(xp, ndim: int):
    dtype = xpx.default_dtype(xp)
    atol = 1e-12 if dtype == xp.float64 else 1e-6
    num = 10
    rng = np.random.default_rng(100)
    shape = (num,) + (ndim,) * (ndim - 1)
    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    q = xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype)
    r = Rotation.from_quat(q)
    p = RigidTransform.from_components(t, r)
    p_inv = p.inv()

    # Test the short-cuts and other integers
    for n in [-5, -2, -1, 0, 1, 2, 5]:
        q = p**n
        # Regression test for gh-24436
        assert isinstance(q._matrix, type(p._matrix))
        r = RigidTransform.from_matrix(xp.tile(xp.eye(4), shape + (1, 1)))
        for _ in range(abs(n)):
            if n > 0:
                r = r * p
            else:
                r = r * p_inv
        xp_assert_close(q.as_matrix(), r.as_matrix(), atol=atol)

    # Test shape preservation of single
    single_tf = RigidTransform.identity()
    assert (single_tf**n).as_matrix().shape == (4, 4)

    # Test fractional powers
    q = p**0.5
    xp_assert_close((q * q).as_matrix(), p.as_matrix(), atol=atol)
    q = p**-0.5
    xp_assert_close((q * q).as_matrix(), p.inv().as_matrix(), atol=atol)
    q = p** 1.5
    xp_assert_close((q * q).as_matrix(), (p**3).as_matrix(), atol=atol)
    q = p** -1.5
    xp_assert_close((q * q).as_matrix(), (p**-3).as_matrix(), atol=atol)

    # pow function
    tf = pow(RigidTransform.from_matrix(xp.eye(4)), 2)
    xp_assert_close(tf.as_matrix(), xp.eye(4), atol=atol)


@make_xp_test_case(RigidTransform.__pow__)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_pow_equivalence_with_rotation(xp, ndim: int):
    atol = 1e-12
    num = 10
    rng = np.random.default_rng(100)
    dtype = xpx.default_dtype(xp)
    shape = (num,) + (ndim,) * (ndim - 1)

    r = Rotation.from_quat(xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype))
    p = RigidTransform.from_rotation(r)
    for n in [-5, -2, -1.5, -1, -0.5, 0.0, 0.5, 1, 1.5, 2, 5]:
        xp_assert_close((p**n).rotation.as_matrix(), (r**n).as_matrix(), atol=atol)


@make_xp_test_case(RigidTransform.inv, RigidTransform.__mul__)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_inverse(xp, ndim: int):
    dtype = xpx.default_dtype(xp)
    atol = 1e-12 if dtype == xp.float64 else 1e-6
    rng = np.random.default_rng(100)
    shape = (ndim,) * (ndim - 1)

    # Test inverse transform
    r = Rotation.from_quat(xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype))
    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    tf = RigidTransform.from_components(t, r)

    # Test that tf * tf.inv() equals identity
    tf_inv = tf.inv()
    composed = tf * tf_inv
    expected = xp.tile(xp.eye(4), shape + (1, 1))
    xp_assert_close(composed.as_matrix(), expected, atol=atol)


@make_xp_test_case(RigidTransform.as_matrix)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_properties(xp, ndim: int):
    atol = 1e-12 if xpx.default_dtype(xp) == xp.float64 else 1e-6
    shape = (ndim,) * (ndim - 1)
    dtype = xpx.default_dtype(xp)
    rng = np.random.default_rng(100)

    # Test rotation and translation properties
    r = Rotation.from_quat(xp.asarray(rng.normal(size=shape + (4,)), dtype=dtype))
    t = xp.asarray(rng.normal(size=shape + (3,)), dtype=dtype)
    tf = RigidTransform.from_components(t, r)

    xp_assert_close(tf.rotation.as_matrix(), r.as_matrix(), atol=atol)
    assert xp.all(tf.rotation.approx_equal(r, atol=atol))
    xp_assert_close(tf.translation, t, atol=atol)
    # Test that we don't return views that would modify the original array
    xpx.at(tf.translation)[..., 0].set(0.0)
    xp_assert_close(tf.translation, t, atol=atol)
    assert tf.single == (shape == ())


@make_xp_test_case(RigidTransform.__getitem__)
def test_indexing(xp):
    atol = 1e-12

    # Test indexing for multiple transforms
    r = Rotation.from_euler('zyx', xp.asarray([[90, 0, 0], [0, 90, 0]]), degrees=True)
    t = xp.asarray([[1.0, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)

    # Test single index
    xp_assert_close(tf[0].as_matrix()[:3, :3], r[0].as_matrix(), atol=atol)
    xp_assert_close(tf[0].as_matrix()[:3, 3], t[0, ...], atol=atol)

    # Test slice
    tf_slice = tf[0:2]
    xp_assert_close(tf_slice.as_matrix()[:, :3, :3], r[0:2].as_matrix(), atol=atol)
    xp_assert_close(tf_slice.as_matrix()[:, :3, 3], t[0:2, ...], atol=atol)

    # Test boolean indexing
    tf_masked = tf[xp.asarray([True, True])]
    xp_assert_close(tf_masked.as_matrix()[:, :3, :3], r.as_matrix(), atol=atol)
    xp_assert_close(tf_masked.as_matrix()[:, :3, 3], t, atol=atol)

    tf_masked = tf[xp.asarray([False, True])]
    xp_assert_close(tf_masked.as_matrix()[:, :3, :3],
                    r[xp.asarray([False, True])].as_matrix(), atol=atol)
    xp_assert_close(tf_masked.as_matrix()[:, :3, 3], t[xp.asarray([False, True])],
                    atol=atol)

    tf_masked = tf[xp.asarray([False, False])]
    assert len(tf_masked) == 0

    # Test integer array indexing
    idx = xp.asarray([0, 1])
    xp_assert_close(tf[idx].as_matrix()[:, :3, :3], r[idx].as_matrix(), atol=atol)
    xp_assert_close(tf[idx].as_matrix()[:, :3, 3], t, atol=atol)


def test_indexing_array_like():
    atol = 1e-12

    r = Rotation.from_euler('zyx', np.array([[90, 0, 0], [0, 90, 0]]), degrees=True)
    t = np.array([[1.0, 2, 3], [4, 5, 6]])
    tf = RigidTransform.from_components(t, r)

    tf_masked = tf[[False, True]]
    xp_assert_close(tf_masked.as_matrix()[:, :3, :3], r[[False, True]].as_matrix(),
                    atol=atol)
    xp_assert_close(tf_masked.as_matrix()[:, :3, 3], t[[False, True]], atol=atol)
    tf_masked = tf[[False, False]]
    assert len(tf_masked) == 0


@make_xp_test_case(RigidTransform.concatenate)
def test_concatenate(xp):
    atol = 1e-12

    # Test concatenation of transforms
    t1 = xp.asarray([1, 0, 0])
    r1 = Rotation.from_euler('z', xp.asarray(90), degrees=True)
    tf1 = RigidTransform.from_components(t1, r1)

    t2 = xp.asarray([0, 1, 0])
    r2 = Rotation.from_euler('x', xp.asarray(90), degrees=True)
    tf2 = RigidTransform.from_components(t2, r2)

    # Concatenate single transforms
    concatenated1 = RigidTransform.concatenate([tf1, tf2])
    xp_assert_close(concatenated1[0].as_matrix(), tf1.as_matrix(), atol=atol)
    xp_assert_close(concatenated1[1].as_matrix(), tf2.as_matrix(), atol=atol)

    # Concatenate multiple transforms
    concatenated2 = RigidTransform.concatenate([tf1, concatenated1])
    xp_assert_close(concatenated2[0].as_matrix(), tf1.as_matrix(), atol=atol)
    xp_assert_close(concatenated2[1].as_matrix(), tf1.as_matrix(), atol=atol)
    xp_assert_close(concatenated2[2].as_matrix(), tf2.as_matrix(), atol=atol)

    # Test ND concatenation
    tf3 = RigidTransform.from_translation(xp.reshape(xp.arange(18), (3, 2, 3)))
    tf4 = RigidTransform.from_translation(xp.reshape(xp.arange(18) + 18, (3, 2, 3)))
    concatenated3 = RigidTransform.concatenate([tf3, tf4])
    xp_assert_close(concatenated3.as_matrix()[:3, ...], tf3.as_matrix(), atol=atol)
    xp_assert_close(concatenated3.as_matrix()[3:, ...], tf4.as_matrix(), atol=atol)


@make_xp_test_case(RigidTransform.from_matrix)
def test_input_validation(xp):
    # Test invalid matrix shapes
    inputs = [xp.eye(3), xp.zeros((4, 3)), []]
    for input in inputs:
        with pytest.raises(ValueError, match="Expected `matrix` to have shape"):
            RigidTransform.from_matrix(input)

    # Test invalid last row
    for ndim in range(3):
        shape = (ndim,) * (ndim - 1)
        matrix = xp.zeros(shape + (4, 4))
        matrix = xpx.at(matrix)[...].set(xp.eye(4))
        matrix = xpx.at(matrix)[..., 3, :].set(xp.asarray([1.0, 0, 0, 1]))
        if is_lazy_array(matrix):
            matrix = RigidTransform.from_matrix(matrix).as_matrix()
            assert xp.all(xp.isnan(matrix[..., 3, :]))
        else:
            with pytest.raises(ValueError, match="last row of transformation matrix"):
                RigidTransform.from_matrix(matrix)

    # Test left handed rotation matrix
    matrix = xp.eye(4)
    matrix = xpx.at(matrix)[0, 0].set(-1)
    if is_lazy_array(matrix):
        matrix = RigidTransform.from_matrix(matrix).as_matrix()
        assert xp.all(xp.isnan(matrix[..., :3, :3]))
    else:
        with pytest.raises(ValueError, match="Non-positive determinant"):
            RigidTransform(matrix, normalize=True)

    # Test non-Rotation input
    with pytest.raises(TypeError,
                       match="Expected `rotation` to be a `Rotation` instance"):
        RigidTransform.from_rotation(xp.eye(3))


@make_xp_test_case(RigidTransform.mean)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_mean(xp, ndim: int):
    atol = 1e-12
    rng = np.random.default_rng(123)

    dtype = xpx.default_dtype(xp)
    t = xp.asarray(rng.normal(size=(ndim,) * (ndim - 1) + (3,)), dtype=dtype)
    q = xp.asarray(rng.normal(size=(ndim,) * (ndim - 1) + (4,)), dtype=dtype)
    r = Rotation.from_quat(q)
    tf = RigidTransform.from_components(t, r)

    # Unweighted mean
    axis = tuple(range(t.ndim - 1))
    t_mean = xp.mean(t, axis=axis)
    r_mean = r.mean()
    tf_mean = tf.mean()
    assert tf_mean.shape == ()
    xp_assert_close(tf_mean.as_matrix(),
                    RigidTransform.from_components(t_mean, r_mean).as_matrix(),
                    atol=atol)

    # Weighted mean
    if ndim == 1:
        weights = None
        t_mean = t
    else:
        weights = xp.asarray(rng.random(size=(ndim,) * (ndim - 1)), dtype=dtype)
        norm = xp.sum(weights[..., None], axis=axis)
        wsum = xp.sum(t * weights[..., None], axis=axis)
        t_mean = wsum/norm
    r_mean = r.mean(weights=weights)
    tf_mean = tf.mean(weights=weights)
    assert tf_mean.shape == ()
    xp_assert_close(tf_mean.as_matrix(),
                    RigidTransform.from_components(t_mean, r_mean).as_matrix(),
                    atol=atol)


@make_xp_test_case(
    RigidTransform.from_rotation, RigidTransform.mean, Rotation.magnitude
)
@pytest.mark.parametrize("ndim", range(1, 5))
def test_mean_axis(xp, ndim: int):
    axes = xp.tile(xp.concat((-xp.eye(3), xp.eye(3))), (3,) * (ndim - 1) + (1, 1))
    theta = xp.pi / 4
    r = Rotation.from_rotvec(theta * axes)
    tf = RigidTransform.from_rotation(r)

    # Test mean over last axis
    desired = xp.full(axes.shape[:-2], 0.0)
    if ndim == 1:
        desired = desired[()]
    atol = 1e-6 if xpx.default_dtype(xp) is xp.float32 else 1e-10
    xp_assert_close(tf.mean(axis=-1).rotation.magnitude(), desired, atol=atol)

    # Test tuple axes
    desired = xp.full(axes.shape[1:-2], 0.0)
    if ndim < 3:
        desired = desired[()]
    xp_assert_close(tf.mean(axis=(0, -1)).rotation.magnitude(), desired, atol=atol)

    # Empty axis tuple should return RigidTransform unchanged
    tf_mean = tf.mean(axis=())
    xp_assert_close(tf_mean.as_matrix(), tf.as_matrix(), atol=atol)


@make_xp_test_case(RigidTransform.mean, Rotation.magnitude)
def test_mean_compare_axis(xp):
    # Create a random set of transforms and compare the mean over an axis with
    # the mean without axis of the sliced transform
    atol = 1e-10 if xpx.default_dtype(xp) == xp.float64 else 1e-6
    rng = np.random.default_rng(0)
    q = xp.asarray(rng.normal(size=(4, 5, 6, 4)), dtype=xpx.default_dtype(xp))
    r = Rotation.from_quat(q)
    t = xp.asarray(rng.normal(size=(4, 5, 6, 3)), dtype=xpx.default_dtype(xp))
    tf = RigidTransform.from_components(t, r)

    mean_0 = tf.mean(axis=0)
    for i in range(q.shape[1]):
        for j in range(q.shape[2]):
            r_slice = Rotation.from_quat(q[:, i, j, ...])
            t_slice = t[:, i, j, ...]
            mean_slice_tf = RigidTransform.from_components(t_slice, r_slice).mean()
            xp_assert_close(
                (mean_0[i][j].rotation * mean_slice_tf.rotation.inv()).magnitude(),
                xp.asarray(0.0)[()], atol=atol,
            )
            xp_assert_close(
                mean_0[i][j].translation, mean_slice_tf.translation, atol=atol,
            )
    mean_1_2 = tf.mean(axis=(1, 2))
    for i in range(q.shape[0]):
        r_slice = Rotation.from_quat(q[i, ...])
        t_slice = t[i, ...]
        mean_slice_tf = RigidTransform.from_components(t_slice, r_slice).mean()
        xp_assert_close(
            (mean_1_2[i].rotation * mean_slice_tf.rotation.inv()).magnitude(),
            xp.asarray(0.0)[()], atol=atol,
        )
        xp_assert_close(
            mean_1_2[i].translation, mean_slice_tf.translation, atol=atol,
        )


@make_xp_test_case(RigidTransform.mean)
def test_mean_invalid_weights(xp):
    tf = RigidTransform.from_matrix(xp.tile(xp.eye(4), (4, 1, 1)))
    if is_lazy_array(tf.as_matrix()):
        m = tf.mean(weights=-xp.ones(4))
        assert xp.all(xp.isnan(m.as_matrix()))
    else:
        with pytest.raises(ValueError, match="non-negative"):
            tf.mean(weights=-xp.ones(4))

    # Test weight shape mismatch
    tf = RigidTransform.from_matrix(xp.eye(4))
    with pytest.raises(ValueError, match="Expected `weights` to"):
        tf.mean(weights=xp.ones((2,)))
    tf = RigidTransform.from_matrix(xp.tile(xp.eye(4), (3, 2, 1, 1, 1)))
    with pytest.raises(ValueError, match="Expected `weights` to"):
        tf.mean(weights=xp.ones((2, 1)))


@make_xp_test_case(RigidTransform.from_translation)
def test_translation_validation(xp):
    # Test invalid translation shapes
    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransform.from_translation(xp.asarray([1, 2]))

    with pytest.raises(ValueError, match="Expected `translation` to have shape"):
        RigidTransform.from_translation(xp.zeros((2, 2)))


@make_xp_test_case(RigidTransform.apply)
def test_vector_validation(xp):
    tf = rigid_transform_to_xp(RigidTransform.identity(2), xp=xp)

    # Test invalid vector shapes
    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(xp.asarray([1, 2]))

    with pytest.raises(ValueError, match="Expected vector to have shape"):
        tf.apply(xp.zeros((2, 2)))

    with pytest.raises(ValueError, match="operands could not be broadcast"):
        tf.apply(xp.zeros((1, 4, 3)))


@make_xp_test_case(RigidTransform.__getitem__)
def test_indexing_validation(xp):
    tf = RigidTransform.from_matrix(xp.eye(4))

    # Test indexing on single transform
    with pytest.raises(TypeError, match="Single transform is not subscriptable"):
        tf[0]

    with pytest.raises(TypeError, match="Single transform is not subscriptable"):
        tf[0:1]

    # Test length on single transform
    with pytest.raises(TypeError, match="Single transform has no len"):
        len(tf)


@make_xp_test_case(RigidTransform.__mul__)
def test_composition_validation(xp):
    tf2 = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6]]))
    tf3 = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))

    # Test incompatible shapes
    with pytest.raises(ValueError, match="Cannot broadcast"):
        tf2 * tf3

    tf4 = RigidTransform.from_matrix(xp.tile(xp.eye(4), (1, 4, 1, 1)))
    # Test invalid broadcasting shape
    with pytest.raises(ValueError, match="Cannot broadcast"):
        tf2 * tf4


@make_xp_test_case(RigidTransform.concatenate)
def test_concatenate_validation(xp):
    tf = RigidTransform.from_matrix(xp.eye(4))

    # Test invalid inputs
    with pytest.raises(TypeError,
                       match="input must contain RigidTransform objects"):
        RigidTransform.concatenate([tf, xp.eye(4)])

    # Test incompatible shapes
    tf2 = RigidTransform.from_translation(xp.ones((1, 1, 1, 3)))
    # Frameworks have a highly heterogeneous way of reporting errors for this case
    with pytest.raises((ValueError, TypeError, RuntimeError)):
        RigidTransform.concatenate([tf, tf2])


@make_xp_test_case(RigidTransform.__setitem__)
def test_setitem(xp):
    tf = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
    single = RigidTransform.from_translation(xp.asarray([1, 1, 1]))
    double = RigidTransform.from_translation(xp.asarray([[2, 2, 2], [3, 3, 3]]))
    triple = RigidTransform.from_translation(xp.asarray([[3, 3, 3],
                                                         [4, 4, 4],
                                                         [5, 5, 5]]))

    # Test indexing with integer index
    tf[0] = single
    xp_assert_close(tf.translation, xp.asarray([[1.0, 1, 1], [4, 5, 6], [7, 8, 9]]))

    # Test indexing with slice
    tf = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
    tf[:2] = double
    xp_assert_close(tf.translation, xp.asarray([[2.0, 2, 2], [3, 3, 3], [7, 8, 9]]))

    # Test indexing with ellipsis
    tf = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
    tf[...] = triple
    xp_assert_close(tf.translation, xp.asarray([[3.0, 3, 3], [4, 4, 4], [5, 5, 5]]))

    # Test indexing with boolean array
    tf = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
    mask = xp.asarray([True, False, True])
    tf[mask] = double
    xp_assert_close(tf.translation, xp.asarray([[2.0, 2, 2], [4, 5, 6], [3, 3, 3]]))


@make_xp_test_case(RigidTransform.__setitem__)
@pytest.mark.skip_xp_backends("array_api_strict",
                              reason="doesn't support fancy indexing __setitem__")
def test_setitem_fancy_indexing(xp):
    double = RigidTransform.from_translation(xp.asarray([[2, 2, 2], [3, 3, 3]]))
    tf = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]))
    idx = xp.asarray([0, 2])
    tf[idx] = double
    xp_assert_close(tf.translation, xp.asarray([[2.0, 2, 2], [4, 5, 6], [3, 3, 3]]))


@make_xp_test_case(RigidTransform.__setitem__)
def test_setitem_validation(xp):
    tf = RigidTransform.from_translation(xp.asarray([[1, 2, 3], [4, 5, 6]]))
    single = RigidTransform.from_matrix(xp.eye(4))

    # Test setting item on single transform
    with pytest.raises(TypeError, match="Single transform is not subscriptable"):
        single[0] = tf

    # Test invalid value type
    with pytest.raises(TypeError, match="value must be a RigidTransform"):
        tf[0] = xp.eye(4)


@pytest.mark.skip_xp_backends("jax.numpy",
                              reason="JAX does not support memory sharing")
@make_xp_test_case(RigidTransform.as_matrix)
def test_copy_flag(xp):
    # Test that copy=True creates new memory
    matrix = xp.eye(4)
    tf = RigidTransform(matrix, normalize=False, copy=True)
    matrix[0, 0] = 2
    assert tf.as_matrix()[0, 0] == 1

    # Test that copy=False shares memory
    matrix = xp.eye(4)
    tf = RigidTransform(matrix, normalize=False, copy=False)
    matrix[0, 0] = 2
    assert tf.as_matrix()[0, 0] == 2


@make_xp_test_case(normalize_dual_quaternion)
@pytest.mark.parametrize("ndim", range(1, 4))
def test_normalize_dual_quaternion(xp, ndim: int):
    dtype = xpx.default_dtype(xp)
    atol = 1e-12 if dtype == xp.float64 else 1e-6
    rng = np.random.default_rng(100)
    shape = (ndim,) * (ndim - 1)

    dual_quat = normalize_dual_quaternion(xp.zeros((1, 8)))
    xp_assert_close(xp_vector_norm(dual_quat[0, :4], axis=-1), xp.asarray(1.0)[()],
                    atol=1e-12)
    xp_assert_close(xp.vecdot(dual_quat[0, :4], dual_quat[0, 4:])[()],
                    xp.asarray(0.0)[()], atol=1e-12)

    dual_quat = xp.asarray(rng.normal(size=shape + (8,)), dtype=dtype)
    dual_quat = normalize_dual_quaternion(dual_quat)
    expected = xp.ones(shape) if shape != () else xp.asarray(1.0)[()]
    xp_assert_close(xp_vector_norm(dual_quat[..., :4], axis=-1), expected, atol=atol)
    expected = xp.zeros(shape) if shape != () else xp.asarray(0.0)[()]
    vecdot = xp.vecdot(dual_quat[..., :4], dual_quat[..., 4:])
    vecdot = vecdot[()] if vecdot.shape == () else vecdot
    xp_assert_close(vecdot, expected, atol=atol)


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.from_rotation,
                   RigidTransform.from_translation, RigidTransform.from_components,
                   RigidTransform.from_exp_coords, RigidTransform.from_dual_quat)
def test_empty_transform_construction(xp):
    tf = RigidTransform.from_matrix(xp.empty((0, 4, 4)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_rotation(Rotation.from_quat(xp.zeros((0, 4))))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_translation(xp.empty((0, 3)))
    assert len(tf) == 0
    assert not tf.single

    empty_rot = Rotation.from_quat(xp.zeros((0, 4)))
    tf = RigidTransform.from_components(xp.empty((0, 3)), empty_rot)
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_exp_coords(xp.empty((0, 6)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.from_dual_quat(xp.empty((0, 8)))
    assert len(tf) == 0
    assert not tf.single

    tf = RigidTransform.identity(0)
    assert len(tf) == 0
    assert not tf.single


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.as_components,
                   RigidTransform.as_exp_coords, RigidTransform.as_dual_quat)
def test_empty_transform_representation(xp):
    tf = RigidTransform.from_matrix(xp.empty((0, 4, 4)))

    assert len(tf.rotation) == 0
    assert tf.translation.shape == (0, 3)

    t, r = tf.as_components()
    assert t.shape == (0, 3)
    assert len(r) == 0

    assert tf.as_matrix().shape == (0, 4, 4)
    assert tf.as_exp_coords().shape == (0, 6)
    assert tf.as_dual_quat().shape == (0, 8)


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.apply)
def test_empty_transform_application(xp):
    tf = RigidTransform.from_matrix(xp.empty((0, 4, 4)))

    assert tf.apply(xp.zeros((3,))).shape == (0, 3)
    assert tf.apply(xp.empty((0, 3))).shape == (0, 3)

    with pytest.raises(ValueError, match="operands could not be broadcast together"):
        tf.apply(xp.zeros((2, 3)))


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.__mul__)
def test_empty_transform_composition(xp):
    tf_empty = RigidTransform.from_matrix(xp.empty((0, 4, 4)))
    tf_single = RigidTransform.from_matrix(xp.eye(4))
    tf_many = rigid_transform_to_xp(RigidTransform.identity(3), xp=xp)

    assert len(tf_empty * tf_empty) == 0
    assert len(tf_empty * tf_single) == 0
    assert len(tf_single * tf_empty) == 0

    with pytest.raises(ValueError, match="Cannot broadcast"):
        tf_many * tf_empty

    with pytest.raises(ValueError, match="Cannot broadcast"):
        tf_empty * tf_many


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.concatenate)
def test_empty_transform_concatenation(xp):
    tf_empty = RigidTransform.from_matrix(xp.empty((0, 4, 4)))
    tf_single = RigidTransform.from_matrix(xp.eye(4))
    tf_many = rigid_transform_to_xp(RigidTransform.identity(2), xp=xp)

    assert len(RigidTransform.concatenate([tf_empty, tf_empty])) == 0
    assert len(RigidTransform.concatenate([tf_empty, tf_single])) == 1
    assert len(RigidTransform.concatenate([tf_single, tf_empty])) == 1
    assert len(RigidTransform.concatenate([tf_empty, tf_many])) == 2
    assert len(RigidTransform.concatenate([tf_many, tf_empty])) == 2
    assert len(RigidTransform.concatenate([tf_many, tf_empty, tf_single])) == 3


@make_xp_test_case(RigidTransform.from_matrix, RigidTransform.inv,
                   RigidTransform.__pow__)
def test_empty_transform_inv_and_pow(xp):
    tf = RigidTransform.from_matrix(xp.empty((0, 4, 4)))
    assert len(tf.inv()) == 0
    assert len(tf ** 0) == 0
    assert len(tf ** 1) == 0
    assert len(tf ** -1) == 0
    assert len(tf ** 0.5) == 0


@make_xp_test_case(RigidTransform.__getitem__)
def test_empty_transform_indexing(xp):
    tf_many = rigid_transform_to_xp(RigidTransform.identity(3), xp=xp)
    tf_zero = tf_many[xp.asarray([], dtype=xp.int32)]
    assert len(tf_zero) == 0

    # Array API does not specify out-of-bounds indexing. Only check for numpy.
    if is_numpy(xp):
        assert len(tf_zero[:5]) == 0  # Slices can go out of bounds.

    with pytest.raises(IndexError):
        tf_zero[0]

    with pytest.raises(IndexError):
        tf_zero[xp.asarray([0, 2])]

    with pytest.raises(IndexError):
        tf_zero[xp.asarray([False, True])]


@make_xp_test_case(RigidTransform.from_matrix)
@pytest.mark.skip_xp_backends("array_api_strict",
                              reason="array API doesn't support pickling")
def test_pickling(xp):
    # Note: Array API makes no provision for arrays to be pickleable, so
    # it's OK to skip this test for the backends that don't support it
    mat = xp.eye(4)
    mat = xpx.at(mat)[0, 3].set(2.0)
    tf = RigidTransform.from_matrix(mat)
    pkl = pickle.dumps(tf)
    unpickled = pickle.loads(pkl)
    xp_assert_close(tf.as_matrix(), unpickled.as_matrix(), atol=1e-15)


@make_xp_test_case(RigidTransform.as_matrix, RigidTransform.__iter__)
def test_rigid_transform_iter(xp):
    r = rigid_transform_to_xp(RigidTransform.identity(3), xp)
    for i, r_i in enumerate(r):
        assert isinstance(r_i, RigidTransform)
        xp_assert_equal(r_i.as_matrix(), r[i].as_matrix())
        if i > len(r):
            raise RuntimeError("Iteration exceeded length of transforms")


@make_xp_test_case()
@pytest.mark.parametrize("dim", range(1, 5))
def test_shape_property(xp, dim: int):
    shape = (dim,) * (dim - 1)
    tf = RigidTransform.from_translation(xp.zeros(shape + (3,)))
    assert tf.shape == shape


def test_non_writeable():
    mat = np.eye(4)
    mat.flags.writeable = False
    RigidTransform.from_matrix(mat)  # Regression test against gh-24378
