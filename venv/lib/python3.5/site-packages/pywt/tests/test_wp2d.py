#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (run_module_suite, assert_allclose, assert_,
                           assert_raises)

import pywt


def test_traversing_tree_2d():
    x = np.array([[1, 2, 3, 4, 5, 6, 7, 8]] * 8, dtype=np.float64)
    wp = pywt.WaveletPacket2D(data=x, wavelet='db1', mode='symmetric')

    assert_(np.all(wp.data == x))
    assert_(wp.path == '')
    assert_(wp.level == 0)
    assert_(wp.maxlevel == 3)

    assert_allclose(wp['a'].data, np.array([[3., 7., 11., 15.]] * 4),
                    rtol=1e-12)
    assert_allclose(wp['h'].data, np.zeros((4, 4)), rtol=1e-12, atol=1e-14)
    assert_allclose(wp['v'].data, -np.ones((4, 4)), rtol=1e-12, atol=1e-14)
    assert_allclose(wp['d'].data, np.zeros((4, 4)), rtol=1e-12, atol=1e-14)

    assert_allclose(wp['aa'].data, np.array([[10., 26.]] * 2), rtol=1e-12)

    assert_(wp['a']['a'].data is wp['aa'].data)
    assert_allclose(wp['aaa'].data, np.array([[36.]]), rtol=1e-12)

    assert_raises(IndexError, lambda: wp['aaaa'])
    assert_raises(ValueError, lambda: wp['f'])


def test_accessing_node_atributes_2d():
    x = np.array([[1, 2, 3, 4, 5, 6, 7, 8]] * 8, dtype=np.float64)
    wp = pywt.WaveletPacket2D(data=x, wavelet='db1', mode='symmetric')

    assert_allclose(wp['av'].data, np.zeros((2, 2)) - 4, rtol=1e-12)
    assert_(wp['av'].path == 'av')
    assert_(wp['av'].node_name == 'v')
    assert_(wp['av'].parent.path == 'a')

    assert_allclose(wp['av'].parent.data, np.array([[3., 7., 11., 15.]] * 4),
                    rtol=1e-12)
    assert_(wp['av'].level == 2)
    assert_(wp['av'].maxlevel == 3)
    assert_(wp['av'].mode == 'symmetric')


def test_collecting_nodes_2d():
    x = np.array([[1, 2, 3, 4, 5, 6, 7, 8]] * 8, dtype=np.float64)
    wp = pywt.WaveletPacket2D(data=x, wavelet='db1', mode='symmetric')

    assert_(len(wp.get_level(0)) == 1)
    assert_(wp.get_level(0)[0].path == '')

    # First level
    assert_(len(wp.get_level(1)) == 4)
    assert_([node.path for node in wp.get_level(1)] == ['a', 'h', 'v', 'd'])

    # Second level
    assert_(len(wp.get_level(2)) == 16)
    paths = [node.path for node in wp.get_level(2)]
    expected_paths = ['aa', 'ah', 'av', 'ad', 'ha', 'hh', 'hv', 'hd', 'va',
                      'vh', 'vv', 'vd', 'da', 'dh', 'dv', 'dd']
    assert_(paths == expected_paths)

    # Third level.
    assert_(len(wp.get_level(3)) == 64)
    paths = [node.path for node in wp.get_level(3)]
    expected_paths = ['aaa', 'aah', 'aav', 'aad', 'aha', 'ahh', 'ahv', 'ahd',
                      'ava', 'avh', 'avv', 'avd', 'ada', 'adh', 'adv', 'add',
                      'haa', 'hah', 'hav', 'had', 'hha', 'hhh', 'hhv', 'hhd',
                      'hva', 'hvh', 'hvv', 'hvd', 'hda', 'hdh', 'hdv', 'hdd',
                      'vaa', 'vah', 'vav', 'vad', 'vha', 'vhh', 'vhv', 'vhd',
                      'vva', 'vvh', 'vvv', 'vvd', 'vda', 'vdh', 'vdv', 'vdd',
                      'daa', 'dah', 'dav', 'dad', 'dha', 'dhh', 'dhv', 'dhd',
                      'dva', 'dvh', 'dvv', 'dvd', 'dda', 'ddh', 'ddv', 'ddd']

    assert_(paths == expected_paths)


def test_data_reconstruction_2d():
    x = np.array([[1, 2, 3, 4, 5, 6, 7, 8]] * 8, dtype=np.float64)
    wp = pywt.WaveletPacket2D(data=x, wavelet='db1', mode='symmetric')

    new_wp = pywt.WaveletPacket2D(data=None, wavelet='db1', mode='symmetric')
    new_wp['vh'] = wp['vh'].data
    new_wp['vv'] = wp['vh'].data
    new_wp['vd'] = np.zeros((2, 2), dtype=np.float64)
    new_wp['a'] = [[3.0, 7.0, 11.0, 15.0]] * 4
    new_wp['d'] = np.zeros((4, 4), dtype=np.float64)
    new_wp['h'] = wp['h']       # all zeros

    assert_allclose(new_wp.reconstruct(update=False),
                    np.array([[1.5, 1.5, 3.5, 3.5, 5.5, 5.5, 7.5, 7.5]] * 8),
                    rtol=1e-12)
    assert_allclose(wp['va'].data, np.zeros((2, 2)) - 2, rtol=1e-12)

    new_wp['va'] = wp['va'].data
    assert_allclose(new_wp.reconstruct(update=False), x, rtol=1e-12)


def test_data_reconstruction_delete_nodes_2d():
    x = np.array([[1, 2, 3, 4, 5, 6, 7, 8]] * 8, dtype=np.float64)
    wp = pywt.WaveletPacket2D(data=x, wavelet='db1', mode='symmetric')

    new_wp = pywt.WaveletPacket2D(data=None, wavelet='db1', mode='symmetric')
    new_wp['vh'] = wp['vh'].data
    new_wp['vv'] = wp['vh'].data
    new_wp['vd'] = np.zeros((2, 2), dtype=np.float64)
    new_wp['a'] = [[3.0, 7.0, 11.0, 15.0]] * 4
    new_wp['d'] = np.zeros((4, 4), dtype=np.float64)
    new_wp['h'] = wp['h']       # all zeros

    assert_allclose(new_wp.reconstruct(update=False),
                    np.array([[1.5, 1.5, 3.5, 3.5, 5.5, 5.5, 7.5, 7.5]] * 8),
                    rtol=1e-12)

    new_wp['va'] = wp['va'].data
    assert_allclose(new_wp.reconstruct(update=False), x, rtol=1e-12)

    del(new_wp['va'])
    new_wp['va'] = wp['va'].data
    assert_(new_wp.data is None)

    assert_allclose(new_wp.reconstruct(update=True), x, rtol=1e-12)
    assert_allclose(new_wp.data, x, rtol=1e-12)

    # TODO: decompose=True


def test_lazy_evaluation_2D():
    # Note: internal implementation detail not to be relied on.  Testing for
    # now for backwards compatibility, but this test may be broken in needed.
    x = np.array([[1, 2, 3, 4, 5, 6, 7, 8]] * 8)
    wp = pywt.WaveletPacket2D(data=x, wavelet='db1', mode='symmetric')

    assert_(wp.a is None)
    assert_allclose(wp['a'].data, np.array([[3., 7., 11., 15.]] * 4),
                    rtol=1e-12)
    assert_allclose(wp.a.data, np.array([[3., 7., 11., 15.]] * 4), rtol=1e-12)
    assert_allclose(wp.d.data, np.zeros((4, 4)), rtol=1e-12, atol=1e-12)


if __name__ == '__main__':
    run_module_suite()
