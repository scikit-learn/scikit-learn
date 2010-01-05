#! /usr/bin/env python
# Last Change: Mon Jul 02 05:00 PM 2007 J

# For now, just test that all mode/dim execute correctly

import sys
from numpy.testing import *

import numpy as N

set_package_path()
from pyem import GM
from pyem.densities import multiple_gauss_den
restore_path()

class test_BasicFunc(NumpyTestCase):
    """Check that basic functionalities work."""
    def test_conf_ellip(self):
        """Only test whether the call succeed. To check wether the result is
        OK, you have to plot the results."""
        d = 3
        k = 3
        w, mu, va = GM.gen_param(d, k)
        gm = GM.fromvalues(w, mu, va)
        gm.conf_ellipses()

    def test_1d_bogus(self):
        """Check that functions which do not make sense for 1d fail nicely."""
        d = 1
        k = 2
        w, mu, va = GM.gen_param(d, k)
        gm = GM.fromvalues(w, mu, va)
        try:
            gm.conf_ellipses()
            raise AssertionError("This should not work !")
        except ValueError, e:
            print "Ok, conf_ellipses failed as expected (with msg: " + str(e) + ")"

        try:
            gm.density_on_grid()
            raise AssertionError("This should not work !")
        except ValueError, e:
            print "Ok, density_grid failed as expected (with msg: " + str(e) + ")"

    def test_get_va(self):
        """Test _get_va for diag and full mode."""
        d = 3
        k = 2
        ld = 2
        dim = [0, 2]
        w, mu, va = GM.gen_param(d, k, 'full')
        va = N.arange(d*d*k).reshape(d*k, d)
        gm = GM.fromvalues(w, mu, va)

        tva = N.empty(ld * ld * k)
        for i in range(k * ld * ld):
            tva[i] = dim[i%ld] + (i % 4)/ ld  * dim[1] * d + d*d * (i / (ld*ld))
        tva = tva.reshape(ld * k, ld)
        sva = gm._get_va(dim)
        assert N.all(sva == tva)

    def test_2d_diag_pdf(self):
        d = 2
        w = N.array([0.4, 0.6])
        mu = N.array([[0., 2], [-1, -2]])
        va = N.array([[1, 0.5], [0.5, 1]])
        x = N.random.randn(100, 2)
        gm = GM.fromvalues(w, mu, va)
        y1 = N.sum(multiple_gauss_den(x, mu, va) * w, 1)
        y2 = gm.pdf(x)
        assert_array_almost_equal(y1, y2)

    def test_2d_diag_logpdf(self):
        d = 2
        w = N.array([0.4, 0.6])
        mu = N.array([[0., 2], [-1, -2]])
        va = N.array([[1, 0.5], [0.5, 1]])
        x = N.random.randn(100, 2)
        gm = GM.fromvalues(w, mu, va)
        y1 = N.sum(multiple_gauss_den(x, mu, va) * w, 1)
        y2 = gm.pdf(x, log = True)
        assert_array_almost_equal(N.log(y1), y2)

if __name__ == "__main__":
    NumpyTest().run()
