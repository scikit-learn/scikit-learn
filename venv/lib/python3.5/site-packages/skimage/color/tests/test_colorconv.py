#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for color conversion functions.

Authors
-------
- the rgb2hsv test was written by Nicolas Pinto, 2009
- other tests written by Ralf Gommers, 2009

:license: modified BSD
"""

from __future__ import division
import os.path

import numpy as np
from skimage._shared.testing import assert_equal, assert_almost_equal
from skimage._shared.testing import assert_array_almost_equal
from skimage._shared.testing import TestCase

from skimage import img_as_float, img_as_ubyte
from skimage.io import imread
from skimage.color import (rgb2hsv, hsv2rgb,
                           rgb2xyz, xyz2rgb,
                           rgb2hed, hed2rgb,
                           separate_stains,
                           combine_stains,
                           rgb2rgbcie, rgbcie2rgb,
                           convert_colorspace,
                           rgb2grey, gray2rgb,
                           xyz2lab, lab2xyz,
                           lab2rgb, rgb2lab,
                           xyz2luv, luv2xyz,
                           luv2rgb, rgb2luv,
                           lab2lch, lch2lab,
                           rgb2yuv, yuv2rgb,
                           rgb2yiq, yiq2rgb,
                           rgb2ypbpr, ypbpr2rgb,
                           rgb2ycbcr, ycbcr2rgb,
                           rgb2ydbdr, ydbdr2rgb,
                           rgba2rgb,
                           guess_spatial_dimensions)

from skimage import data_dir
from skimage._shared._warnings import expected_warnings
from skimage._shared import testing
import colorsys


def test_guess_spatial_dimensions():
    im1 = np.zeros((5, 5))
    im2 = np.zeros((5, 5, 5))
    im3 = np.zeros((5, 5, 3))
    im4 = np.zeros((5, 5, 5, 3))
    im5 = np.zeros((5,))
    assert_equal(guess_spatial_dimensions(im1), 2)
    assert_equal(guess_spatial_dimensions(im2), 3)
    assert_equal(guess_spatial_dimensions(im3), None)
    assert_equal(guess_spatial_dimensions(im4), 3)
    with testing.raises(ValueError):
        guess_spatial_dimensions(im5)


class TestColorconv(TestCase):

    img_rgb = imread(os.path.join(data_dir, 'color.png'))
    img_grayscale = imread(os.path.join(data_dir, 'camera.png'))
    img_rgba = np.array([[[0, 0.5, 1, 0],
                          [0, 0.5, 1, 1],
                          [0, 0.5, 1, 0.5]]]).astype(np.float)

    colbars = np.array([[1, 1, 0, 0, 1, 1, 0, 0],
                        [1, 1, 1, 1, 0, 0, 0, 0],
                        [1, 0, 1, 0, 1, 0, 1, 0]]).astype(np.float)
    colbars_array = np.swapaxes(colbars.reshape(3, 4, 2), 0, 2)
    colbars_point75 = colbars * 0.75
    colbars_point75_array = np.swapaxes(colbars_point75.reshape(3, 4, 2), 0, 2)

    xyz_array = np.array([[[0.4124, 0.21260, 0.01930]],    # red
                          [[0, 0, 0]],    # black
                          [[.9505, 1., 1.089]],    # white
                          [[.1805, .0722, .9505]],    # blue
                          [[.07719, .15438, .02573]],    # green
                          ])
    lab_array = np.array([[[53.233, 80.109, 67.220]],    # red
                          [[0., 0., 0.]],    # black
                          [[100.0, 0.005, -0.010]],    # white
                          [[32.303, 79.197, -107.864]],    # blue
                          [[46.229, -51.7, 49.898]],    # green
                          ])

    luv_array = np.array([[[53.233, 175.053, 37.751]],   # red
                          [[0., 0., 0.]],   # black
                          [[100., 0.001, -0.017]],   # white
                          [[32.303, -9.400, -130.358]],   # blue
                          [[46.228, -43.774, 56.589]],   # green
                          ])

    # RGBA to RGB
    def test_rgba2rgb_conversion(self):
        rgba = self.img_rgba
        rgb = rgba2rgb(rgba)
        expected = np.array([[[1, 1, 1],
                              [0, 0.5, 1],
                              [0.5, 0.75, 1]]]).astype(np.float)
        self.assertEqual(rgb.shape, expected.shape)
        assert_almost_equal(rgb, expected)

    def test_rgba2rgb_error_grayscale(self):
        self.assertRaises(ValueError, rgba2rgb, self.img_grayscale)

    def test_rgba2rgb_error_rgb(self):
        self.assertRaises(ValueError, rgba2rgb, self.img_rgb)

    # RGB to HSV
    def test_rgb2hsv_conversion(self):
        rgb = img_as_float(self.img_rgb)[::16, ::16]
        hsv = rgb2hsv(rgb).reshape(-1, 3)
        # ground truth from colorsys
        gt = np.array([colorsys.rgb_to_hsv(pt[0], pt[1], pt[2])
                       for pt in rgb.reshape(-1, 3)]
                      )
        assert_almost_equal(hsv, gt)

    def test_rgb2hsv_error_grayscale(self):
        self.assertRaises(ValueError, rgb2hsv, self.img_grayscale)

    def test_rgb2hsv_error_one_element(self):
        self.assertRaises(ValueError, rgb2hsv, self.img_rgb[0, 0])

    # HSV to RGB
    def test_hsv2rgb_conversion(self):
        rgb = self.img_rgb.astype("float32")[::16, ::16]
        # create HSV image with colorsys
        hsv = np.array([colorsys.rgb_to_hsv(pt[0], pt[1], pt[2])
                        for pt in rgb.reshape(-1, 3)]).reshape(rgb.shape)
        # convert back to RGB and compare with original.
        # relative precision for RGB -> HSV roundtrip is about 1e-6
        assert_almost_equal(rgb, hsv2rgb(hsv), decimal=4)

    def test_hsv2rgb_error_grayscale(self):
        self.assertRaises(ValueError, hsv2rgb, self.img_grayscale)

    def test_hsv2rgb_error_one_element(self):
        self.assertRaises(ValueError, hsv2rgb, self.img_rgb[0, 0])

    # RGB to XYZ
    def test_rgb2xyz_conversion(self):
        gt = np.array([[[0.950456, 1.      , 1.088754],
                        [0.538003, 0.787329, 1.06942 ],
                        [0.592876, 0.28484 , 0.969561],
                        [0.180423, 0.072169, 0.950227]],
                       [[0.770033, 0.927831, 0.138527],
                        [0.35758 , 0.71516 , 0.119193],
                        [0.412453, 0.212671, 0.019334],
                        [0.      , 0.      , 0.      ]]])
        assert_almost_equal(rgb2xyz(self.colbars_array), gt)

    # stop repeating the "raises" checks for all other functions that are
    # implemented with color._convert()
    def test_rgb2xyz_error_grayscale(self):
        self.assertRaises(ValueError, rgb2xyz, self.img_grayscale)

    def test_rgb2xyz_error_one_element(self):
        self.assertRaises(ValueError, rgb2xyz, self.img_rgb[0, 0])

    # XYZ to RGB
    def test_xyz2rgb_conversion(self):
        assert_almost_equal(xyz2rgb(rgb2xyz(self.colbars_array)),
                            self.colbars_array)

    # RGB<->XYZ roundtrip on another image
    def test_xyz_rgb_roundtrip(self):
        img_rgb = img_as_float(self.img_rgb)
        assert_array_almost_equal(xyz2rgb(rgb2xyz(img_rgb)), img_rgb)

    # RGB<->HED roundtrip with ubyte image
    def test_hed_rgb_roundtrip(self):
        img_rgb = img_as_ubyte(self.img_rgb)
        with expected_warnings(['precision loss']):
            new = img_as_ubyte(hed2rgb(rgb2hed(img_rgb)))
        assert_equal(new, img_rgb)

    # RGB<->HED roundtrip with float image
    def test_hed_rgb_float_roundtrip(self):
        img_rgb = img_as_float(self.img_rgb)
        assert_array_almost_equal(hed2rgb(rgb2hed(img_rgb)), img_rgb)

    # RGB<->HDX roundtrip with ubyte image
    def test_hdx_rgb_roundtrip(self):
        from skimage.color.colorconv import hdx_from_rgb, rgb_from_hdx
        img_rgb = self.img_rgb
        conv = combine_stains(separate_stains(img_rgb, hdx_from_rgb),
                              rgb_from_hdx)
        assert_equal(img_as_ubyte(conv), img_rgb)

    # RGB<->HDX roundtrip with ubyte image
    def test_hdx_rgb_roundtrip(self):
        from skimage.color.colorconv import hdx_from_rgb, rgb_from_hdx
        img_rgb = img_as_float(self.img_rgb)
        conv = combine_stains(separate_stains(img_rgb, hdx_from_rgb),
                              rgb_from_hdx)
        assert_array_almost_equal(conv, img_rgb)

    # RGB to RGB CIE
    def test_rgb2rgbcie_conversion(self):
        gt = np.array([[[ 0.1488856 ,  0.18288098,  0.19277574],
                        [ 0.01163224,  0.16649536,  0.18948516],
                        [ 0.12259182,  0.03308008,  0.17298223],
                        [-0.01466154,  0.01669446,  0.16969164]],
                       [[ 0.16354714,  0.16618652,  0.0230841 ],
                        [ 0.02629378,  0.1498009 ,  0.01979351],
                        [ 0.13725336,  0.01638562,  0.00329059],
                        [ 0.        ,  0.        ,  0.        ]]])
        assert_almost_equal(rgb2rgbcie(self.colbars_array), gt)

    # RGB CIE to RGB
    def test_rgbcie2rgb_conversion(self):
        # only roundtrip test, we checked rgb2rgbcie above already
        assert_almost_equal(rgbcie2rgb(rgb2rgbcie(self.colbars_array)),
                            self.colbars_array)

    def test_convert_colorspace(self):
        colspaces = ['HSV', 'RGB CIE', 'XYZ', 'YCbCr', 'YPbPr', 'YDbDr']
        colfuncs_from = [
            hsv2rgb, rgbcie2rgb, xyz2rgb,
            ycbcr2rgb, ypbpr2rgb, ydbdr2rgb
        ]
        colfuncs_to = [
            rgb2hsv, rgb2rgbcie, rgb2xyz,
            rgb2ycbcr, rgb2ypbpr, rgb2ydbdr
        ]

        assert_almost_equal(
            convert_colorspace(self.colbars_array, 'RGB', 'RGB'),
            self.colbars_array)

        for i, space in enumerate(colspaces):
            gt = colfuncs_from[i](self.colbars_array)
            assert_almost_equal(
                convert_colorspace(self.colbars_array, space, 'RGB'), gt)
            gt = colfuncs_to[i](self.colbars_array)
            assert_almost_equal(
                convert_colorspace(self.colbars_array, 'RGB', space), gt)

        self.assertRaises(ValueError, convert_colorspace,
                          self.colbars_array, 'nokey', 'XYZ')
        self.assertRaises(ValueError, convert_colorspace,
                          self.colbars_array, 'RGB', 'nokey')

    def test_rgb2grey(self):
        x = np.array([1, 1, 1]).reshape((1, 1, 3)).astype(np.float)
        g = rgb2grey(x)
        assert_array_almost_equal(g, 1)

        assert_equal(g.shape, (1, 1))

    def test_rgb2grey_contiguous(self):
        x = np.random.rand(10, 10, 3)
        assert rgb2grey(x).flags["C_CONTIGUOUS"]
        assert rgb2grey(x[:5, :5]).flags["C_CONTIGUOUS"]

    def test_rgb2grey_alpha(self):
        x = np.random.rand(10, 10, 4)
        assert rgb2grey(x).ndim == 2

    def test_rgb2grey_on_grey(self):
        rgb2grey(np.random.rand(5, 5))

    # test matrices for xyz2lab and lab2xyz generated using
    # http://www.easyrgb.com/index.php?X=CALC
    # Note: easyrgb website displays xyz*100
    def test_xyz2lab(self):
        assert_array_almost_equal(xyz2lab(self.xyz_array),
                                  self.lab_array, decimal=3)

        # Test the conversion with the rest of the illuminants.
        for I in ["d50", "d55", "d65", "d75"]:
            for obs in ["2", "10"]:
                fname = "lab_array_{0}_{1}.npy".format(I, obs)
                lab_array_I_obs = np.load(
                    os.path.join(os.path.dirname(__file__), 'data', fname))
                assert_array_almost_equal(lab_array_I_obs,
                                          xyz2lab(self.xyz_array, I, obs),
                                          decimal=2)
        for I in ["a", "e"]:
            fname = "lab_array_{0}_2.npy".format(I)
            lab_array_I_obs = np.load(
                os.path.join(os.path.dirname(__file__), 'data', fname))
            assert_array_almost_equal(lab_array_I_obs,
                                      xyz2lab(self.xyz_array, I, "2"),
                                      decimal=2)

    def test_lab2xyz(self):
        assert_array_almost_equal(lab2xyz(self.lab_array),
                                  self.xyz_array, decimal=3)

        # Test the conversion with the rest of the illuminants.
        for I in ["d50", "d55", "d65", "d75"]:
            for obs in ["2", "10"]:
                fname = "lab_array_{0}_{1}.npy".format(I, obs)
                lab_array_I_obs = np.load(
                    os.path.join(os.path.dirname(__file__), 'data', fname))
                assert_array_almost_equal(lab2xyz(lab_array_I_obs, I, obs),
                                          self.xyz_array, decimal=3)
        for I in ["a", "e"]:
            fname = "lab_array_{0}_2.npy".format(I, obs)
            lab_array_I_obs = np.load(
                os.path.join(os.path.dirname(__file__), 'data', fname))
            assert_array_almost_equal(lab2xyz(lab_array_I_obs, I, "2"),
                                      self.xyz_array, decimal=3)

        # And we include a call to test the exception handling in the code.
        try:
            xs = lab2xyz(lab_array_I_obs, "NaI", "2")   # Not an illuminant
        except ValueError:
            pass

        try:
            xs = lab2xyz(lab_array_I_obs, "d50", "42")   # Not a degree
        except ValueError:
            pass

    def test_rgb2lab_brucelindbloom(self):
        """
        Test the RGB->Lab conversion by comparing to the calculator on the
        authoritative Bruce Lindbloom
        [website](http://brucelindbloom.com/index.html?ColorCalculator.html).
        """
        # Obtained with D65 white point, sRGB model and gamma
        gt_for_colbars = np.array([
            [100,0,0],
            [97.1393, -21.5537, 94.4780],
            [91.1132, -48.0875, -14.1312],
            [87.7347, -86.1827, 83.1793],
            [60.3242, 98.2343, -60.8249],
            [53.2408, 80.0925, 67.2032],
            [32.2970, 79.1875, -107.8602],
            [0,0,0]]).T
        gt_array = np.swapaxes(gt_for_colbars.reshape(3, 4, 2), 0, 2)
        assert_array_almost_equal(rgb2lab(self.colbars_array), gt_array, decimal=2)

    def test_lab_rgb_roundtrip(self):
        img_rgb = img_as_float(self.img_rgb)
        assert_array_almost_equal(lab2rgb(rgb2lab(img_rgb)), img_rgb)

    # test matrices for xyz2luv and luv2xyz generated using
    # http://www.easyrgb.com/index.php?X=CALC
    # Note: easyrgb website displays xyz*100
    def test_xyz2luv(self):
        assert_array_almost_equal(xyz2luv(self.xyz_array),
                                  self.luv_array, decimal=3)

        # Test the conversion with the rest of the illuminants.
        for I in ["d50", "d55", "d65", "d75"]:
            for obs in ["2", "10"]:
                fname = "luv_array_{0}_{1}.npy".format(I, obs)
                luv_array_I_obs = np.load(
                    os.path.join(os.path.dirname(__file__), 'data', fname))
                assert_array_almost_equal(luv_array_I_obs,
                                          xyz2luv(self.xyz_array, I, obs),
                                          decimal=2)
        for I in ["a", "e"]:
            fname = "luv_array_{0}_2.npy".format(I)
            luv_array_I_obs = np.load(
                os.path.join(os.path.dirname(__file__), 'data', fname))
            assert_array_almost_equal(luv_array_I_obs,
                                      xyz2luv(self.xyz_array, I, "2"),
                                      decimal=2)

    def test_luv2xyz(self):
        assert_array_almost_equal(luv2xyz(self.luv_array),
                                  self.xyz_array, decimal=3)

        # Test the conversion with the rest of the illuminants.
        for I in ["d50", "d55", "d65", "d75"]:
            for obs in ["2", "10"]:
                fname = "luv_array_{0}_{1}.npy".format(I, obs)
                luv_array_I_obs = np.load(
                    os.path.join(os.path.dirname(__file__), 'data', fname))
                assert_array_almost_equal(luv2xyz(luv_array_I_obs, I, obs),
                                          self.xyz_array, decimal=3)
        for I in ["a", "e"]:
            fname = "luv_array_{0}_2.npy".format(I, obs)
            luv_array_I_obs = np.load(
                os.path.join(os.path.dirname(__file__), 'data', fname))
            assert_array_almost_equal(luv2xyz(luv_array_I_obs, I, "2"),
                                      self.xyz_array, decimal=3)

    def test_rgb2luv_brucelindbloom(self):
        """
        Test the RGB->Lab conversion by comparing to the calculator on the
        authoritative Bruce Lindbloom
        [website](http://brucelindbloom.com/index.html?ColorCalculator.html).
        """
        # Obtained with D65 white point, sRGB model and gamma
        gt_for_colbars = np.array([
            [100, 0, 0],
            [97.1393, 7.7056, 106.7866],
            [91.1132, -70.4773, -15.2042],
            [87.7347, -83.0776, 107.3985],
            [60.3242, 84.0714, -108.6834],
            [53.2408, 175.0151, 37.7564],
            [32.2970, -9.4054, -130.3423],
            [0, 0, 0]]).T
        gt_array = np.swapaxes(gt_for_colbars.reshape(3, 4, 2), 0, 2)
        assert_array_almost_equal(rgb2luv(self.colbars_array),
                                  gt_array, decimal=2)

    def test_luv_rgb_roundtrip(self):
        img_rgb = img_as_float(self.img_rgb)
        assert_array_almost_equal(luv2rgb(rgb2luv(img_rgb)), img_rgb)

    def test_lab_rgb_outlier(self):
        lab_array = np.ones((3, 1, 3))
        lab_array[0] = [50, -12, 85]
        lab_array[1] = [50, 12, -85]
        lab_array[2] = [90, -4, -47]
        rgb_array = np.array([[[0.501, 0.481, 0]],
                              [[0, 0.482, 1.]],
                              [[0.578, 0.914, 1.]],
                              ])
        assert_almost_equal(lab2rgb(lab_array), rgb_array, decimal=3)

    def test_lab_full_gamut(self):
        a, b = np.meshgrid(np.arange(-100, 100), np.arange(-100, 100))
        L = np.ones(a.shape)
        lab = np.dstack((L, a, b))
        for value in [0, 10, 20]:
            lab[:, :, 0] = value
            with expected_warnings(['Color data out of range']):
                lab2xyz(lab)

    def test_lab_lch_roundtrip(self):
        rgb = img_as_float(self.img_rgb)
        lab = rgb2lab(rgb)
        lab2 = lch2lab(lab2lch(lab))
        assert_array_almost_equal(lab2, lab)

    def test_rgb_lch_roundtrip(self):
        rgb = img_as_float(self.img_rgb)
        lab = rgb2lab(rgb)
        lch = lab2lch(lab)
        lab2 = lch2lab(lch)
        rgb2 = lab2rgb(lab2)
        assert_array_almost_equal(rgb, rgb2)

    def test_lab_lch_0d(self):
        lab0 = self._get_lab0()
        lch0 = lab2lch(lab0)
        lch2 = lab2lch(lab0[None, None, :])
        assert_array_almost_equal(lch0, lch2[0, 0, :])

    def test_lab_lch_1d(self):
        lab0 = self._get_lab0()
        lch0 = lab2lch(lab0)
        lch1 = lab2lch(lab0[None, :])
        assert_array_almost_equal(lch0, lch1[0, :])

    def test_lab_lch_3d(self):
        lab0 = self._get_lab0()
        lch0 = lab2lch(lab0)
        lch3 = lab2lch(lab0[None, None, None, :])
        assert_array_almost_equal(lch0, lch3[0, 0, 0, :])

    def _get_lab0(self):
        rgb = img_as_float(self.img_rgb[:1, :1, :])
        return rgb2lab(rgb)[0, 0, :]

    def test_yuv(self):
        rgb = np.array([[[1.0, 1.0, 1.0]]])
        assert_array_almost_equal(rgb2yuv(rgb), np.array([[[1, 0, 0]]]))
        assert_array_almost_equal(rgb2yiq(rgb), np.array([[[1, 0, 0]]]))
        assert_array_almost_equal(rgb2ypbpr(rgb), np.array([[[1, 0, 0]]]))
        assert_array_almost_equal(rgb2ycbcr(rgb), np.array([[[235, 128, 128]]]))
        assert_array_almost_equal(rgb2ydbdr(rgb), np.array([[[1, 0, 0]]]))
        rgb = np.array([[[0.0, 1.0, 0.0]]])
        assert_array_almost_equal(rgb2yuv(rgb), np.array([[[0.587, -0.28886916, -0.51496512]]]))
        assert_array_almost_equal(rgb2yiq(rgb), np.array([[[0.587, -0.27455667, -0.52273617]]]))
        assert_array_almost_equal(rgb2ypbpr(rgb), np.array([[[0.587, -0.331264, -0.418688]]]))
        assert_array_almost_equal(rgb2ycbcr(rgb), np.array([[[144.553,   53.797,   34.214]]]))
        assert_array_almost_equal(rgb2ydbdr(rgb), np.array([[[0.587, -0.883, 1.116]]]))

    def test_yuv_roundtrip(self):
        img_rgb = img_as_float(self.img_rgb)[::16, ::16]
        assert_array_almost_equal(yuv2rgb(rgb2yuv(img_rgb)), img_rgb)
        assert_array_almost_equal(yiq2rgb(rgb2yiq(img_rgb)), img_rgb)
        assert_array_almost_equal(ypbpr2rgb(rgb2ypbpr(img_rgb)), img_rgb)
        assert_array_almost_equal(ycbcr2rgb(rgb2ycbcr(img_rgb)), img_rgb)
        assert_array_almost_equal(ydbdr2rgb(rgb2ydbdr(img_rgb)), img_rgb)

    def test_rgb2yiq_conversion(self):
        rgb = img_as_float(self.img_rgb)[::16, ::16]
        yiq = rgb2yiq(rgb).reshape(-1, 3)
        gt = np.array([colorsys.rgb_to_yiq(pt[0], pt[1], pt[2])
                       for pt in rgb.reshape(-1, 3)]
                      )
        assert_almost_equal(yiq, gt, decimal=2)


def test_gray2rgb():
    x = np.array([0, 0.5, 1])
    w = gray2rgb(x)
    expected_output = np.array([[ 0, 0, 0 ],
                                [ 0.5, 0.5, 0.5, ],
                                [ 1, 1, 1 ]])

    assert_equal(w, expected_output)

    x = x.reshape((3, 1))
    y = gray2rgb(x)

    assert_equal(y.shape, (3, 1, 3))
    assert_equal(y.dtype, x.dtype)
    assert_equal(y[..., 0], x)
    assert_equal(y[0, 0, :], [0, 0, 0])

    x = np.array([[0, 128, 255]], dtype=np.uint8)
    z = gray2rgb(x)

    assert_equal(z.shape, (1, 3, 3))
    assert_equal(z[..., 0], x)
    assert_equal(z[0, 1, :], [128, 128, 128])


def test_gray2rgb_rgb():
    x = np.random.rand(5, 5, 4)
    y = gray2rgb(x)
    assert_equal(x, y)


def test_gray2rgb_alpha():
    x = np.random.random((5, 5, 4))
    assert_equal(gray2rgb(x, alpha=None).shape, (5, 5, 4))
    assert_equal(gray2rgb(x, alpha=False).shape, (5, 5, 3))
    assert_equal(gray2rgb(x, alpha=True).shape, (5, 5, 4))

    x = np.random.random((5, 5, 3))
    assert_equal(gray2rgb(x, alpha=None).shape, (5, 5, 3))
    assert_equal(gray2rgb(x, alpha=False).shape, (5, 5, 3))
    assert_equal(gray2rgb(x, alpha=True).shape, (5, 5, 4))

    assert_equal(gray2rgb(np.array([[1, 2], [3, 4.]]),
                          alpha=True)[0, 0, 3], 1)
    assert_equal(gray2rgb(np.array([[1, 2], [3, 4]], dtype=np.uint8),
                          alpha=True)[0, 0, 3], 255)
