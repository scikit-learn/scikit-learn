#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (run_module_suite, assert_allclose, assert_,
                           assert_raises)

import pywt


def test_wavelet_packet_structure():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    assert_(wp.data == [1, 2, 3, 4, 5, 6, 7, 8])
    assert_(wp.path == '')
    assert_(wp.level == 0)
    assert_(wp['ad'].maxlevel == 3)


def test_traversing_wp_tree():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    assert_(wp.maxlevel == 3)

    # First level
    assert_allclose(wp['a'].data, np.array([2.12132034356, 4.949747468306,
                                           7.778174593052, 10.606601717798]),
                    rtol=1e-12)

    # Second level
    assert_allclose(wp['aa'].data, np.array([5., 13.]), rtol=1e-12)

    # Third level
    assert_allclose(wp['aaa'].data, np.array([12.727922061358]), rtol=1e-12)


def test_acess_path():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    assert_(wp['a'].path == 'a')
    assert_(wp['aa'].path == 'aa')
    assert_(wp['aaa'].path == 'aaa')

    # Maximum level reached:
    assert_raises(IndexError, lambda: wp['aaaa'].path)

    # Wrong path
    assert_raises(ValueError, lambda: wp['ac'].path)


def test_access_node_atributes():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    assert_allclose(wp['ad'].data, np.array([-2., -2.]), rtol=1e-12)
    assert_(wp['ad'].path == 'ad')
    assert_(wp['ad'].node_name == 'd')
    assert_(wp['ad'].parent.path == 'a')
    assert_(wp['ad'].level == 2)
    assert_(wp['ad'].maxlevel == 3)
    assert_(wp['ad'].mode == 'symmetric')


def test_collecting_nodes():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    # All nodes in natural order
    assert_([node.path for node in wp.get_level(3, 'natural')] ==
            ['aaa', 'aad', 'ada', 'add', 'daa', 'dad', 'dda', 'ddd'])

    # and in frequency order.
    assert_([node.path for node in wp.get_level(3, 'freq')] ==
            ['aaa', 'aad', 'add', 'ada', 'dda', 'ddd', 'dad', 'daa'])


def test_reconstructing_data():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')

    # Create another Wavelet Packet and feed it with some data.
    new_wp = pywt.WaveletPacket(data=None, wavelet='db1', mode='symmetric')
    new_wp['aa'] = wp['aa'].data
    new_wp['ad'] = [-2., -2.]

    # For convenience, :attr:`Node.data` gets automatically extracted
    # from the :class:`Node` object:
    new_wp['d'] = wp['d']

    # Reconstruct data from aa, ad, and d packets.
    assert_allclose(new_wp.reconstruct(update=False), x, rtol=1e-12)

    # The node's :attr:`~Node.data` will not be updated
    assert_(new_wp.data is None)

    # When `update` is True:
    assert_allclose(new_wp.reconstruct(update=True), x, rtol=1e-12)
    assert_allclose(new_wp.data, np.arange(1, 9), rtol=1e-12)

    assert_([n.path for n in new_wp.get_leaf_nodes(False)] ==
            ['aa', 'ad', 'd'])
    assert_([n.path for n in new_wp.get_leaf_nodes(True)] ==
            ['aaa', 'aad', 'ada', 'add', 'daa', 'dad', 'dda', 'ddd'])


def test_removing_nodes():
    x = [1, 2, 3, 4, 5, 6, 7, 8]
    wp = pywt.WaveletPacket(data=x, wavelet='db1', mode='symmetric')
    wp.get_level(2)

    dataleafs = [n.data for n in wp.get_leaf_nodes(False)]
    expected = np.array([[5., 13.], [-2, -2], [-1, -1], [0, 0]])

    for i in range(4):
        assert_allclose(dataleafs[i], expected[i, :], atol=1e-12)

    node = wp['ad']
    del(wp['ad'])
    dataleafs = [n.data for n in wp.get_leaf_nodes(False)]
    expected = np.array([[5., 13.], [-1, -1], [0, 0]])

    for i in range(3):
        assert_allclose(dataleafs[i], expected[i, :], atol=1e-12)

    wp.reconstruct()
    # The reconstruction is:
    assert_allclose(wp.reconstruct(),
                    np.array([2., 3., 2., 3., 6., 7., 6., 7.]), rtol=1e-12)

    # Restore the data
    wp['ad'].data = node.data

    dataleafs = [n.data for n in wp.get_leaf_nodes(False)]
    expected = np.array([[5., 13.], [-2, -2], [-1, -1], [0, 0]])

    for i in range(4):
        assert_allclose(dataleafs[i], expected[i, :], atol=1e-12)

    assert_allclose(wp.reconstruct(), np.arange(1, 9), rtol=1e-12)


if __name__ == '__main__':
    run_module_suite()
