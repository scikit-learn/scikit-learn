from __future__ import division, print_function, absolute_import

import pytest

import numpy as np
from numpy.testing import assert_equal, assert_array_almost_equal
from numpy.testing import assert_allclose
from scipy.spatial.transform import Rotation, Slerp
from scipy.stats import special_ortho_group
from itertools import permutations


def test_generic_quat_matrix():
    x = np.array([[3, 4, 0, 0], [5, 12, 0, 0]])
    r = Rotation.from_quat(x)
    expected_quat = x / np.array([[5], [13]])
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_from_single_1d_quaternion():
    x = np.array([3, 4, 0, 0])
    r = Rotation.from_quat(x)
    expected_quat = x / 5
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_from_single_2d_quaternion():
    x = np.array([[3, 4, 0, 0]])
    r = Rotation.from_quat(x)
    expected_quat = x / 5
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_from_square_quat_matrix():
    # Ensure proper norm array broadcasting
    x = np.array([
        [3, 0, 0, 4],
        [5, 0, 12, 0],
        [0, 0, 0, 1],
        [0, 0, 0, -1]
        ])
    r = Rotation.from_quat(x)
    expected_quat = x / np.array([[5], [13], [1], [1]])
    assert_array_almost_equal(r.as_quat(), expected_quat)


def test_malformed_1d_from_quat():
    with pytest.raises(ValueError):
        Rotation.from_quat(np.array([1, 2, 3]))


def test_malformed_2d_from_quat():
    with pytest.raises(ValueError):
        Rotation.from_quat(np.array([
            [1, 2, 3, 4, 5],
            [4, 5, 6, 7, 8]
            ]))


def test_zero_norms_from_quat():
    x = np.array([
            [3, 4, 0, 0],
            [0, 0, 0, 0],
            [5, 0, 12, 0]
            ])
    with pytest.raises(ValueError):
        Rotation.from_quat(x)


def test_as_dcm_single_1d_quaternion():
    quat = [0, 0, 0, 1]
    mat = Rotation.from_quat(quat).as_dcm()
    # mat.shape == (3,3) due to 1d input
    assert_array_almost_equal(mat, np.eye(3))


def test_as_dcm_single_2d_quaternion():
    quat = [[0, 0, 1, 1]]
    mat = Rotation.from_quat(quat).as_dcm()
    assert_equal(mat.shape, (1, 3, 3))
    expected_mat = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
        ])
    assert_array_almost_equal(mat[0], expected_mat)


def test_as_dcm_from_square_input():
    quats = [
            [0, 0, 1, 1],
            [0, 1, 0, 1],
            [0, 0, 0, 1],
            [0, 0, 0, -1]
            ]
    mat = Rotation.from_quat(quats).as_dcm()
    assert_equal(mat.shape, (4, 3, 3))

    expected0 = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
        ])
    assert_array_almost_equal(mat[0], expected0)

    expected1 = np.array([
        [0, 0, 1],
        [0, 1, 0],
        [-1, 0, 0]
        ])
    assert_array_almost_equal(mat[1], expected1)

    assert_array_almost_equal(mat[2], np.eye(3))
    assert_array_almost_equal(mat[3], np.eye(3))


def test_as_dcm_from_generic_input():
    quats = [
            [0, 0, 1, 1],
            [0, 1, 0, 1],
            [1, 2, 3, 4]
            ]
    mat = Rotation.from_quat(quats).as_dcm()
    assert_equal(mat.shape, (3, 3, 3))

    expected0 = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
        ])
    assert_array_almost_equal(mat[0], expected0)

    expected1 = np.array([
        [0, 0, 1],
        [0, 1, 0],
        [-1, 0, 0]
        ])
    assert_array_almost_equal(mat[1], expected1)

    expected2 = np.array([
        [0.4, -2, 2.2],
        [2.8, 1, 0.4],
        [-1, 2, 2]
        ]) / 3
    assert_array_almost_equal(mat[2], expected2)


def test_from_single_2d_dcm():
    dcm = [
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0]
            ]
    expected_quat = [0.5, 0.5, 0.5, 0.5]
    assert_array_almost_equal(
            Rotation.from_dcm(dcm).as_quat(),
            expected_quat)


def test_from_single_3d_dcm():
    dcm = np.array([
        [0, 0, 1],
        [1, 0, 0],
        [0, 1, 0]
        ]).reshape((1, 3, 3))
    expected_quat = np.array([0.5, 0.5, 0.5, 0.5]).reshape((1, 4))
    assert_array_almost_equal(
            Rotation.from_dcm(dcm).as_quat(),
            expected_quat)


def test_from_dcm_calculation():
    expected_quat = np.array([1, 1, 6, 1]) / np.sqrt(39)
    dcm = np.array([
            [-0.8974359, -0.2564103, 0.3589744],
            [0.3589744, -0.8974359, 0.2564103],
            [0.2564103, 0.3589744, 0.8974359]
            ])
    assert_array_almost_equal(
            Rotation.from_dcm(dcm).as_quat(),
            expected_quat)
    assert_array_almost_equal(
            Rotation.from_dcm(dcm.reshape((1, 3, 3))).as_quat(),
            expected_quat.reshape((1, 4)))


def test_dcm_calculation_pipeline():
    dcm = special_ortho_group.rvs(3, size=10, random_state=0)
    assert_array_almost_equal(Rotation.from_dcm(dcm).as_dcm(), dcm)


def test_from_dcm_ortho_output():
    np.random.seed(0)
    dcm = np.random.random((100, 3, 3))
    ortho_dcm = Rotation.from_dcm(dcm).as_dcm()

    mult_result = np.einsum('...ij,...jk->...ik', ortho_dcm,
                            ortho_dcm.transpose((0, 2, 1)))

    eye3d = np.zeros((100, 3, 3))
    for i in range(3):
        eye3d[:, i, i] = 1.0

    assert_array_almost_equal(mult_result, eye3d)


def test_from_1d_single_rotvec():
    rotvec = [1, 0, 0]
    expected_quat = np.array([0.4794255, 0, 0, 0.8775826])
    result = Rotation.from_rotvec(rotvec)
    assert_array_almost_equal(result.as_quat(), expected_quat)


def test_from_2d_single_rotvec():
    rotvec = [[1, 0, 0]]
    expected_quat = np.array([[0.4794255, 0, 0, 0.8775826]])
    result = Rotation.from_rotvec(rotvec)
    assert_array_almost_equal(result.as_quat(), expected_quat)


def test_from_generic_rotvec():
    rotvec = [
            [1, 2, 2],
            [1, -1, 0.5],
            [0, 0, 0]
            ]
    expected_quat = np.array([
        [0.3324983, 0.6649967, 0.6649967, 0.0707372],
        [0.4544258, -0.4544258, 0.2272129, 0.7316889],
        [0, 0, 0, 1]
        ])
    assert_array_almost_equal(
            Rotation.from_rotvec(rotvec).as_quat(),
            expected_quat)


def test_from_rotvec_small_angle():
    rotvec = np.array([
        [5e-4 / np.sqrt(3), -5e-4 / np.sqrt(3), 5e-4 / np.sqrt(3)],
        [0.2, 0.3, 0.4],
        [0, 0, 0]
        ])

    quat = Rotation.from_rotvec(rotvec).as_quat()
    # cos(theta/2) ~~ 1 for small theta
    assert_allclose(quat[0, 3], 1)
    # sin(theta/2) / theta ~~ 0.5 for small theta
    assert_allclose(quat[0, :3], rotvec[0] * 0.5)

    assert_allclose(quat[1, 3], 0.9639685)
    assert_allclose(
            quat[1, :3],
            np.array([
                0.09879603932153465,
                0.14819405898230198,
                0.19759207864306931
                ]))

    assert_equal(quat[2], np.array([0, 0, 0, 1]))


def test_malformed_1d_from_rotvec():
    with pytest.raises(ValueError, match='Expected `rot_vec` to have shape'):
        Rotation.from_rotvec([1, 2])


def test_malformed_2d_from_rotvec():
    with pytest.raises(ValueError, match='Expected `rot_vec` to have shape'):
        Rotation.from_rotvec([
            [1, 2, 3, 4],
            [5, 6, 7, 8]
            ])


def test_as_generic_rotvec():
    quat = np.array([
            [1, 2, -1, 0.5],
            [1, -1, 1, 0.0003],
            [0, 0, 0, 1]
            ])
    quat /= np.linalg.norm(quat, axis=1)[:, None]

    rotvec = Rotation.from_quat(quat).as_rotvec()
    angle = np.linalg.norm(rotvec, axis=1)

    assert_allclose(quat[:, 3], np.cos(angle/2))
    assert_allclose(np.cross(rotvec, quat[:, :3]), np.zeros((3, 3)))


def test_as_rotvec_single_1d_input():
    quat = np.array([1, 2, -3, 2])
    expected_rotvec = np.array([0.5772381, 1.1544763, -1.7317144])

    actual_rotvec = Rotation.from_quat(quat).as_rotvec()

    assert_equal(actual_rotvec.shape, (3,))
    assert_allclose(actual_rotvec, expected_rotvec)


def test_as_rotvec_single_2d_input():
    quat = np.array([[1, 2, -3, 2]])
    expected_rotvec = np.array([[0.5772381, 1.1544763, -1.7317144]])

    actual_rotvec = Rotation.from_quat(quat).as_rotvec()

    assert_equal(actual_rotvec.shape, (1, 3))
    assert_allclose(actual_rotvec, expected_rotvec)


def test_rotvec_calc_pipeline():
    # Include small angles
    rotvec = np.array([
        [0, 0, 0],
        [1, -1, 2],
        [-3e-4, 3.5e-4, 7.5e-5]
        ])
    assert_allclose(Rotation.from_rotvec(rotvec).as_rotvec(), rotvec)


def test_from_euler_single_rotation():
    quat = Rotation.from_euler('z', 90, degrees=True).as_quat()
    expected_quat = np.array([0, 0, 1, 1]) / np.sqrt(2)
    assert_allclose(quat, expected_quat)


def test_single_intrinsic_extrinsic_rotation():
    extrinsic = Rotation.from_euler('z', 90, degrees=True).as_dcm()
    intrinsic = Rotation.from_euler('Z', 90, degrees=True).as_dcm()
    assert_allclose(extrinsic, intrinsic)


def test_from_euler_rotation_order():
    # Intrinsic rotation is same as extrinsic with order reversed
    np.random.seed(0)
    a = np.random.randint(low=0, high=180, size=(6, 3))
    b = a[:, ::-1]
    x = Rotation.from_euler('xyz', a, degrees=True).as_quat()
    y = Rotation.from_euler('ZYX', b, degrees=True).as_quat()
    assert_allclose(x, y)


def test_from_euler_elementary_extrinsic_rotation():
    # Simple test to check if extrinsic rotations are implemented correctly
    dcm = Rotation.from_euler('zx', [90, 90], degrees=True).as_dcm()
    expected_dcm = np.array([
        [0, -1, 0],
        [0, 0, -1],
        [1, 0, 0]
    ])
    assert_array_almost_equal(dcm, expected_dcm)


def test_from_euler_intrinsic_rotation_312():
    angles = [
        [30, 60, 45],
        [30, 60, 30],
        [45, 30, 60]
        ]
    dcm = Rotation.from_euler('ZXY', angles, degrees=True).as_dcm()

    assert_array_almost_equal(dcm[0], np.array([
        [0.3061862, -0.2500000, 0.9185587],
        [0.8838835, 0.4330127, -0.1767767],
        [-0.3535534, 0.8660254, 0.3535534]
    ]))

    assert_array_almost_equal(dcm[1], np.array([
        [0.5334936, -0.2500000, 0.8080127],
        [0.8080127, 0.4330127, -0.3995191],
        [-0.2500000, 0.8660254, 0.4330127]
    ]))

    assert_array_almost_equal(dcm[2], np.array([
        [0.0473672, -0.6123725, 0.7891491],
        [0.6597396, 0.6123725, 0.4355958],
        [-0.7500000, 0.5000000, 0.4330127]
    ]))


def test_from_euler_intrinsic_rotation_313():
    angles = [
        [30, 60, 45],
        [30, 60, 30],
        [45, 30, 60]
        ]
    dcm = Rotation.from_euler('ZXZ', angles, degrees=True).as_dcm()

    assert_array_almost_equal(dcm[0], np.array([
        [0.43559574, -0.78914913, 0.4330127],
        [0.65973961, -0.04736717, -0.750000],
        [0.61237244, 0.61237244, 0.500000]
    ]))

    assert_array_almost_equal(dcm[1], np.array([
        [0.6250000, -0.64951905, 0.4330127],
        [0.64951905, 0.1250000, -0.750000],
        [0.4330127, 0.750000, 0.500000]
    ]))

    assert_array_almost_equal(dcm[2], np.array([
        [-0.1767767, -0.91855865, 0.35355339],
        [0.88388348, -0.30618622, -0.35355339],
        [0.4330127, 0.25000000, 0.8660254]
    ]))


def test_from_euler_extrinsic_rotation_312():
    angles = [
        [30, 60, 45],
        [30, 60, 30],
        [45, 30, 60]
        ]
    dcm = Rotation.from_euler('zxy', angles, degrees=True).as_dcm()

    assert_array_almost_equal(dcm[0], np.array([
        [0.91855865, 0.1767767, 0.35355339],
        [0.25000000, 0.4330127, -0.8660254],
        [-0.30618622, 0.88388348, 0.35355339]
    ]))

    assert_array_almost_equal(dcm[1], np.array([
        [0.96650635, -0.0580127, 0.2500000],
        [0.25000000, 0.4330127, -0.8660254],
        [-0.0580127, 0.89951905, 0.4330127]
    ]))

    assert_array_almost_equal(dcm[2], np.array([
        [0.65973961, -0.04736717, 0.7500000],
        [0.61237244, 0.61237244, -0.5000000],
        [-0.43559574, 0.78914913, 0.4330127]
    ]))


def test_from_euler_extrinsic_rotation_313():
    angles = [
        [30, 60, 45],
        [30, 60, 30],
        [45, 30, 60]
        ]
    dcm = Rotation.from_euler('zxz', angles, degrees=True).as_dcm()

    assert_array_almost_equal(dcm[0], np.array([
        [0.43559574, -0.65973961, 0.61237244],
        [0.78914913, -0.04736717, -0.61237244],
        [0.4330127, 0.75000000, 0.500000]
    ]))

    assert_array_almost_equal(dcm[1], np.array([
        [0.62500000, -0.64951905, 0.4330127],
        [0.64951905, 0.12500000, -0.750000],
        [0.4330127, 0.75000000, 0.500000]
    ]))

    assert_array_almost_equal(dcm[2], np.array([
        [-0.1767767, -0.88388348, 0.4330127],
        [0.91855865, -0.30618622, -0.250000],
        [0.35355339, 0.35355339, 0.8660254]
    ]))


def test_as_euler_asymmetric_axes():
    np.random.seed(0)
    n = 10
    angles = np.empty((n, 3))
    angles[:, 0] = np.random.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles[:, 1] = np.random.uniform(low=-np.pi / 2, high=np.pi / 2, size=(n,))
    angles[:, 2] = np.random.uniform(low=-np.pi, high=np.pi, size=(n,))

    for seq_tuple in permutations('xyz'):
        # Extrinsic rotations
        seq = ''.join(seq_tuple)
        assert_allclose(angles, Rotation.from_euler(seq, angles).as_euler(seq))
        # Intrinsic rotations
        seq = seq.upper()
        assert_allclose(angles, Rotation.from_euler(seq, angles).as_euler(seq))


def test_as_euler_symmetric_axes():
    np.random.seed(0)
    n = 10
    angles = np.empty((n, 3))
    angles[:, 0] = np.random.uniform(low=-np.pi, high=np.pi, size=(n,))
    angles[:, 1] = np.random.uniform(low=0, high=np.pi, size=(n,))
    angles[:, 2] = np.random.uniform(low=-np.pi, high=np.pi, size=(n,))

    for axis1 in ['x', 'y', 'z']:
        for axis2 in ['x', 'y', 'z']:
            if axis1 == axis2:
                continue
            # Extrinsic rotations
            seq = axis1 + axis2 + axis1
            assert_allclose(
                angles, Rotation.from_euler(seq, angles).as_euler(seq))
            # Intrinsic rotations
            seq = seq.upper()
            assert_allclose(
                angles, Rotation.from_euler(seq, angles).as_euler(seq))


def test_as_euler_degenerate_asymmetric_axes():
    # Since we cannot check for angle equality, we check for dcm equality
    angles = np.array([
        [45, 90, 35],
        [35, -90, 20],
        [35, 90, 25],
        [25, -90, 15]
        ])

    with pytest.warns(UserWarning, match="Gimbal lock"):
        for seq_tuple in permutations('xyz'):
            # Extrinsic rotations
            seq = ''.join(seq_tuple)
            rotation = Rotation.from_euler(seq, angles, degrees=True)
            dcm_expected = rotation.as_dcm()

            angle_estimates = rotation.as_euler(seq, degrees=True)
            dcm_estimated = Rotation.from_euler(
                seq, angle_estimates, degrees=True
                ).as_dcm()

            assert_array_almost_equal(dcm_expected, dcm_estimated)

            # Intrinsic rotations
            seq = seq.upper()
            rotation = Rotation.from_euler(seq, angles, degrees=True)
            dcm_expected = rotation.as_dcm()

            angle_estimates = rotation.as_euler(seq, degrees=True)
            dcm_estimated = Rotation.from_euler(
                seq, angle_estimates, degrees=True
                ).as_dcm()

            assert_array_almost_equal(dcm_expected, dcm_estimated)


def test_as_euler_degenerate_symmetric_axes():
    # Since we cannot check for angle equality, we check for dcm equality
    angles = np.array([
        [15, 0, 60],
        [35, 0, 75],
        [60, 180, 35],
        [15, -180, 25],
        ])

    with pytest.warns(UserWarning, match="Gimbal lock"):
        for axis1 in ['x', 'y', 'z']:
            for axis2 in ['x', 'y', 'z']:
                if axis1 == axis2:
                    continue

                # Extrinsic rotations
                seq = axis1 + axis2 + axis1
                rotation = Rotation.from_euler(seq, angles, degrees=True)
                dcm_expected = rotation.as_dcm()

                angle_estimates = rotation.as_euler(seq, degrees=True)
                dcm_estimated = Rotation.from_euler(
                    seq, angle_estimates, degrees=True
                    ).as_dcm()

                assert_array_almost_equal(dcm_expected, dcm_estimated)

                # Intrinsic rotations
                seq = seq.upper()
                rotation = Rotation.from_euler(seq, angles, degrees=True)
                dcm_expected = rotation.as_dcm()

                angle_estimates = rotation.as_euler(seq, degrees=True)
                dcm_estimated = Rotation.from_euler(
                    seq, angle_estimates, degrees=True
                    ).as_dcm()

                assert_array_almost_equal(dcm_expected, dcm_estimated)


def test_inv():
    np.random.seed(0)
    n = 10
    p = Rotation.from_quat(np.random.normal(size=(n, 4)))
    q = p.inv()

    p_dcm = p.as_dcm()
    q_dcm = q.as_dcm()
    result1 = np.einsum('...ij,...jk->...ik', p_dcm, q_dcm)
    result2 = np.einsum('...ij,...jk->...ik', q_dcm, p_dcm)

    eye3d = np.empty((n, 3, 3))
    eye3d[:] = np.eye(3)

    assert_array_almost_equal(result1, eye3d)
    assert_array_almost_equal(result2, eye3d)


def test_inv_single_rotation():
    np.random.seed(0)
    p = Rotation.from_quat(np.random.normal(size=(4,)))
    q = p.inv()

    p_dcm = p.as_dcm()
    q_dcm = q.as_dcm()
    res1 = np.dot(p_dcm, q_dcm)
    res2 = np.dot(q_dcm, p_dcm)

    eye = np.eye(3)

    assert_array_almost_equal(res1, eye)
    assert_array_almost_equal(res2, eye)

    x = Rotation.from_quat(np.random.normal(size=(1, 4)))
    y = x.inv()

    x_dcm = x.as_dcm()
    y_dcm = y.as_dcm()
    result1 = np.einsum('...ij,...jk->...ik', x_dcm, y_dcm)
    result2 = np.einsum('...ij,...jk->...ik', y_dcm, x_dcm)

    eye3d = np.empty((1, 3, 3))
    eye3d[:] = np.eye(3)

    assert_array_almost_equal(result1, eye3d)
    assert_array_almost_equal(result2, eye3d)


def test_apply_single_rotation_single_point():
    dcm = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
    ])
    r_1d = Rotation.from_dcm(dcm)
    r_2d = Rotation.from_dcm(np.expand_dims(dcm, axis=0))

    v_1d = np.array([1, 2, 3])
    v_2d = np.expand_dims(v_1d, axis=0)
    v1d_rotated = np.array([-2, 1, 3])
    v2d_rotated = np.expand_dims(v1d_rotated, axis=0)

    assert_allclose(r_1d.apply(v_1d), v1d_rotated)
    assert_allclose(r_1d.apply(v_2d), v2d_rotated)
    assert_allclose(r_2d.apply(v_1d), v2d_rotated)
    assert_allclose(r_2d.apply(v_2d), v2d_rotated)

    v1d_inverse = np.array([2, -1, 3])
    v2d_inverse = np.expand_dims(v1d_inverse, axis=0)

    assert_allclose(r_1d.apply(v_1d, inverse=True), v1d_inverse)
    assert_allclose(r_1d.apply(v_2d, inverse=True), v2d_inverse)
    assert_allclose(r_2d.apply(v_1d, inverse=True), v2d_inverse)
    assert_allclose(r_2d.apply(v_2d, inverse=True), v2d_inverse)


def test_apply_single_rotation_multiple_points():
    dcm = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
    ])
    r1 = Rotation.from_dcm(dcm)
    r2 = Rotation.from_dcm(np.expand_dims(dcm, axis=0))

    v = np.array([[1, 2, 3], [4, 5, 6]])
    v_rotated = np.array([[-2, 1, 3], [-5, 4, 6]])

    assert_allclose(r1.apply(v), v_rotated)
    assert_allclose(r2.apply(v), v_rotated)

    v_inverse = np.array([[2, -1, 3], [5, -4, 6]])

    assert_allclose(r1.apply(v, inverse=True), v_inverse)
    assert_allclose(r2.apply(v, inverse=True), v_inverse)


def test_apply_multiple_rotations_single_point():
    dcm = np.empty((2, 3, 3))
    dcm[0] = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
    ])
    dcm[1] = np.array([
        [1, 0, 0],
        [0, 0, -1],
        [0, 1, 0]
    ])
    r = Rotation.from_dcm(dcm)

    v1 = np.array([1, 2, 3])
    v2 = np.expand_dims(v1, axis=0)

    v_rotated = np.array([[-2, 1, 3], [1, -3, 2]])

    assert_allclose(r.apply(v1), v_rotated)
    assert_allclose(r.apply(v2), v_rotated)

    v_inverse = np.array([[2, -1, 3], [1, 3, -2]])

    assert_allclose(r.apply(v1, inverse=True), v_inverse)
    assert_allclose(r.apply(v2, inverse=True), v_inverse)


def test_apply_multiple_rotations_multiple_points():
    dcm = np.empty((2, 3, 3))
    dcm[0] = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
    ])
    dcm[1] = np.array([
        [1, 0, 0],
        [0, 0, -1],
        [0, 1, 0]
    ])
    r = Rotation.from_dcm(dcm)

    v = np.array([[1, 2, 3], [4, 5, 6]])
    v_rotated = np.array([[-2, 1, 3], [4, -6, 5]])
    assert_allclose(r.apply(v), v_rotated)

    v_inverse = np.array([[2, -1, 3], [4, 6, -5]])
    assert_allclose(r.apply(v, inverse=True), v_inverse)


def test_getitem():
    dcm = np.empty((2, 3, 3))
    dcm[0] = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
    ])
    dcm[1] = np.array([
        [1, 0, 0],
        [0, 0, -1],
        [0, 1, 0]
    ])
    r = Rotation.from_dcm(dcm)

    assert_allclose(r[0].as_dcm(), dcm[0])
    assert_allclose(r[1].as_dcm(), dcm[1])
    assert_allclose(r[:-1].as_dcm(), np.expand_dims(dcm[0], axis=0))


def test_n_rotations():
    dcm = np.empty((2, 3, 3))
    dcm[0] = np.array([
        [0, -1, 0],
        [1, 0, 0],
        [0, 0, 1]
    ])
    dcm[1] = np.array([
        [1, 0, 0],
        [0, 0, -1],
        [0, 1, 0]
    ])
    r = Rotation.from_dcm(dcm)

    assert_equal(len(r), 2)
    assert_equal(len(r[0]), 1)
    assert_equal(len(r[1]), 1)
    assert_equal(len(r[:-1]), 1)


def test_quat_ownership():
    # Ensure that users cannot accidentally corrupt object
    quat = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0]
    ])
    r = Rotation.from_quat(quat, normalized=True)
    s = r[0:2]

    r._quat[0] = np.array([0, -1, 0, 0])
    assert_allclose(s._quat[0], np.array([1, 0, 0, 0]))


def test_match_vectors_no_rotation():
    x = np.array([[1, 2, 3], [4, 5, 6]])
    y = x.copy()

    r, p = Rotation.match_vectors(x, y)
    assert_array_almost_equal(r.as_dcm(), np.eye(3))


def test_match_vectors_no_noise():
    np.random.seed(0)
    c = Rotation.from_quat(np.random.normal(size=4))
    b = np.random.normal(size=(5, 3))
    a = c.apply(b)

    est, cov = Rotation.match_vectors(a, b)
    assert_allclose(c.as_quat(), est.as_quat())


def test_match_vectors_noise():
    np.random.seed(0)
    n_vectors = 100
    rot = Rotation.from_euler('xyz', np.random.normal(size=3))
    vectors = np.random.normal(size=(n_vectors, 3))
    result = rot.apply(vectors)

    # The paper adds noise as indepedently distributed angular errors
    sigma = np.deg2rad(1)
    tolerance = 1.5 * sigma
    noise = Rotation.from_rotvec(
        np.random.normal(
            size=(n_vectors, 3),
            scale=sigma
        )
    )

    # Attitude errors must preserve norm. Hence apply individual random
    # rotations to each vector.
    noisy_result = noise.apply(result)

    est, cov = Rotation.match_vectors(noisy_result, vectors)

    # Use rotation compositions to find out closeness
    error_vector = (rot * est.inv()).as_rotvec()
    assert_allclose(error_vector[0], 0, atol=tolerance)
    assert_allclose(error_vector[1], 0, atol=tolerance)
    assert_allclose(error_vector[2], 0, atol=tolerance)

    # Check error bounds using covariance matrix
    cov *= sigma
    assert_allclose(cov[0, 0], 0, atol=tolerance)
    assert_allclose(cov[1, 1], 0, atol=tolerance)
    assert_allclose(cov[2, 2], 0, atol=tolerance)


def test_random_rotation_shape():
    assert_equal(Rotation.random().as_quat().shape, (4,))
    assert_equal(Rotation.random(None).as_quat().shape, (4,))

    assert_equal(Rotation.random(1).as_quat().shape, (1, 4))
    assert_equal(Rotation.random(5).as_quat().shape, (5, 4))


def test_slerp():
    np.random.seed(0)

    key_rots = Rotation.from_quat(np.random.uniform(size=(5, 4)))
    key_quats = key_rots.as_quat()

    key_times = [0, 1, 2, 3, 4]
    interpolator = Slerp(key_times, key_rots)

    times = [0, 0.5, 0.25, 1, 1.5, 2, 2.75, 3, 3.25, 3.60, 4]
    interp_rots = interpolator(times)
    interp_quats = interp_rots.as_quat()

    # Dot products are affected by sign of quaternions
    interp_quats[interp_quats[:, -1] < 0] *= -1
    # Checking for quaternion equality, perform same operation
    key_quats[key_quats[:, -1] < 0] *= -1

    # Equality at keyframes, including both endpoints
    assert_allclose(interp_quats[0], key_quats[0])
    assert_allclose(interp_quats[3], key_quats[1])
    assert_allclose(interp_quats[5], key_quats[2])
    assert_allclose(interp_quats[7], key_quats[3])
    assert_allclose(interp_quats[10], key_quats[4])

    # Constant angular velocity between keyframes. Check by equating
    # cos(theta) between quaternion pairs with equal time difference.
    cos_theta1 = np.sum(interp_quats[0] * interp_quats[2])
    cos_theta2 = np.sum(interp_quats[2] * interp_quats[1])
    assert_allclose(cos_theta1, cos_theta2)

    cos_theta4 = np.sum(interp_quats[3] * interp_quats[4])
    cos_theta5 = np.sum(interp_quats[4] * interp_quats[5])
    assert_allclose(cos_theta4, cos_theta5)

    # theta1: 0 -> 0.25, theta3 : 0.5 -> 1
    # Use double angle formula for double the time difference
    cos_theta3 = np.sum(interp_quats[1] * interp_quats[3])
    assert_allclose(cos_theta3, 2 * (cos_theta1**2) - 1)

    # Miscellaneous checks
    assert_equal(len(interp_rots), len(times))


def test_slerp_single_rot():
    with pytest.raises(ValueError, match="at least 2 rotations"):
        r = Rotation.from_quat([1, 2, 3, 4])
        Slerp([1], r)


def test_slerp_time_dim_mismatch():
    with pytest.raises(ValueError,
                       match="times to be specified in a 1 dimensional array"):
        np.random.seed(0)
        r = Rotation.from_quat(np.random.uniform(size=(2, 4)))
        t = np.array([[1],
                      [2]])
        Slerp(t, r)


def test_slerp_num_rotations_mismatch():
    with pytest.raises(ValueError, match="number of rotations to be equal to "
                                         "number of timestamps"):
        np.random.seed(0)
        r = Rotation.from_quat(np.random.uniform(size=(5, 4)))
        t = np.arange(7)
        Slerp(t, r)


def test_slerp_equal_times():
    with pytest.raises(ValueError, match="strictly increasing order"):
        np.random.seed(0)
        r = Rotation.from_quat(np.random.uniform(size=(5, 4)))
        t = [0, 1, 2, 2, 4]
        Slerp(t, r)


def test_slerp_decreasing_times():
    with pytest.raises(ValueError, match="strictly increasing order"):
        np.random.seed(0)
        r = Rotation.from_quat(np.random.uniform(size=(5, 4)))
        t = [0, 1, 3, 2, 4]
        Slerp(t, r)


def test_slerp_call_time_dim_mismatch():
    np.random.seed(0)
    r = Rotation.from_quat(np.random.uniform(size=(5, 4)))
    t = np.arange(5)
    s = Slerp(t, r)

    with pytest.raises(ValueError,
                       match="times to be specified in a 1 dimensional array"):
        interp_times = np.array([[3.5],
                                 [4.2]])
        s(interp_times)


def test_slerp_call_time_out_of_range():
    np.random.seed(0)
    r = Rotation.from_quat(np.random.uniform(size=(5, 4)))
    t = np.arange(5) + 1
    s = Slerp(t, r)

    with pytest.raises(ValueError, match="times must be within the range"):
        s([0, 1, 2])
    with pytest.raises(ValueError, match="times must be within the range"):
        s([1, 2, 6])
