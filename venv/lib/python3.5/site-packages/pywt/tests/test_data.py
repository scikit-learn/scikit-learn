import numpy as np
from numpy.testing import assert_allclose

import pywt.data


def test_data_aero():
    aero = pywt.data.aero()

    ref = np.array([[178, 178, 179],
                    [170, 173, 171],
                    [185, 174, 171]])

    assert_allclose(aero[:3, :3], ref)


def test_data_ascent():
    ascent = pywt.data.ascent()

    ref = np.array([[83, 83, 83],
                    [82, 82, 83],
                    [80, 81, 83]])

    assert_allclose(ascent[:3, :3], ref)


def test_data_camera():
    ascent = pywt.data.camera()

    ref = np.array([[156, 157, 160],
                    [156, 157, 159],
                    [158, 157, 156]])

    assert_allclose(ascent[:3, :3], ref)


def test_data_ecg():
    ecg = pywt.data.ecg()

    ref = np.array([-86, -87, -87])

    assert_allclose(ecg[:3], ref)
