# Created by Pearu Peterson, June 2003
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import (assert_equal, assert_almost_equal, assert_array_equal,
        assert_array_almost_equal, assert_allclose)
from scipy._lib._numpy_compat import suppress_warnings
from pytest import raises as assert_raises

from numpy import array, diff, linspace, meshgrid, ones, pi, shape
from scipy.interpolate.fitpack import bisplrep, bisplev
from scipy.interpolate.fitpack2 import (UnivariateSpline,
        LSQUnivariateSpline, InterpolatedUnivariateSpline,
        LSQBivariateSpline, SmoothBivariateSpline, RectBivariateSpline,
        LSQSphereBivariateSpline, SmoothSphereBivariateSpline,
        RectSphereBivariateSpline)


class TestUnivariateSpline(object):
    def test_linear_constant(self):
        x = [1,2,3]
        y = [3,3,3]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[3,3,3])

    def test_preserve_shape(self):
        x = [1, 2, 3]
        y = [0, 2, 4]
        lut = UnivariateSpline(x, y, k=1)
        arg = 2
        assert_equal(shape(arg), shape(lut(arg)))
        assert_equal(shape(arg), shape(lut(arg, nu=1)))
        arg = [1.5, 2, 2.5]
        assert_equal(shape(arg), shape(lut(arg)))
        assert_equal(shape(arg), shape(lut(arg, nu=1)))

    def test_linear_1d(self):
        x = [1,2,3]
        y = [0,2,4]
        lut = UnivariateSpline(x,y,k=1)
        assert_array_almost_equal(lut.get_knots(),[1,3])
        assert_array_almost_equal(lut.get_coeffs(),[0,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2]),[0,1,2])

    def test_subclassing(self):
        # See #731

        class ZeroSpline(UnivariateSpline):
            def __call__(self, x):
                return 0*array(x)

        sp = ZeroSpline([1,2,3,4,5], [3,2,3,2,3], k=2)
        assert_array_equal(sp([1.5, 2.5]), [0., 0.])

    def test_empty_input(self):
        # Test whether empty input returns an empty output. Ticket 1014
        x = [1,3,5,7,9]
        y = [0,4,9,12,21]
        spl = UnivariateSpline(x, y, k=3)
        assert_array_equal(spl([]), array([]))

    def test_resize_regression(self):
        """Regression test for #1375."""
        x = [-1., -0.65016502, -0.58856235, -0.26903553, -0.17370892,
             -0.10011001, 0., 0.10011001, 0.17370892, 0.26903553, 0.58856235,
             0.65016502, 1.]
        y = [1.,0.62928599, 0.5797223, 0.39965815, 0.36322694, 0.3508061,
             0.35214793, 0.3508061, 0.36322694, 0.39965815, 0.5797223,
             0.62928599, 1.]
        w = [1.00000000e+12, 6.88875973e+02, 4.89314737e+02, 4.26864807e+02,
             6.07746770e+02, 4.51341444e+02, 3.17480210e+02, 4.51341444e+02,
             6.07746770e+02, 4.26864807e+02, 4.89314737e+02, 6.88875973e+02,
             1.00000000e+12]
        spl = UnivariateSpline(x=x, y=y, w=w, s=None)
        desired = array([0.35100374, 0.51715855, 0.87789547, 0.98719344])
        assert_allclose(spl([0.1, 0.5, 0.9, 0.99]), desired, atol=5e-4)

    def test_out_of_range_regression(self):
        # Test different extrapolation modes. See ticket 3557
        x = np.arange(5, dtype=float)
        y = x**3

        xp = linspace(-8, 13, 100)
        xp_zeros = xp.copy()
        xp_zeros[np.logical_or(xp_zeros < 0., xp_zeros > 4.)] = 0
        xp_clip = xp.copy()
        xp_clip[xp_clip < x[0]] = x[0]
        xp_clip[xp_clip > x[-1]] = x[-1]

        for cls in [UnivariateSpline, InterpolatedUnivariateSpline]:
            spl = cls(x=x, y=y)
            for ext in [0, 'extrapolate']:
                assert_allclose(spl(xp, ext=ext), xp**3, atol=1e-16)
                assert_allclose(cls(x, y, ext=ext)(xp), xp**3, atol=1e-16)
            for ext in [1, 'zeros']:
                assert_allclose(spl(xp, ext=ext), xp_zeros**3, atol=1e-16)
                assert_allclose(cls(x, y, ext=ext)(xp), xp_zeros**3, atol=1e-16)
            for ext in [2, 'raise']:
                assert_raises(ValueError, spl, xp, **dict(ext=ext))
            for ext in [3, 'const']:
                assert_allclose(spl(xp, ext=ext), xp_clip**3, atol=1e-16)
                assert_allclose(cls(x, y, ext=ext)(xp), xp_clip**3, atol=1e-16)

        # also test LSQUnivariateSpline [which needs explicit knots]
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        spl = LSQUnivariateSpline(x, y, t)
        assert_allclose(spl(xp, ext=0), xp**3, atol=1e-16)
        assert_allclose(spl(xp, ext=1), xp_zeros**3, atol=1e-16)
        assert_raises(ValueError, spl, xp, **dict(ext=2))
        assert_allclose(spl(xp, ext=3), xp_clip**3, atol=1e-16)

        # also make sure that unknown values for `ext` are caught early
        for ext in [-1, 'unknown']:
            spl = UnivariateSpline(x, y)
            assert_raises(ValueError, spl, xp, **dict(ext=ext))
            assert_raises(ValueError, UnivariateSpline,
                    **dict(x=x, y=y, ext=ext))

    def test_lsq_fpchec(self):
        xs = np.arange(100) * 1.
        ys = np.arange(100) * 1.
        knots = np.linspace(0, 99, 10)
        bbox = (-1, 101)
        assert_raises(ValueError, LSQUnivariateSpline, xs, ys, knots,
                      bbox=bbox)

    def test_derivative_and_antiderivative(self):
        # Thin wrappers to splder/splantider, so light smoke test only.
        x = np.linspace(0, 1, 70)**3
        y = np.cos(x)

        spl = UnivariateSpline(x, y, s=0)
        spl2 = spl.antiderivative(2).derivative(2)
        assert_allclose(spl(0.3), spl2(0.3))

        spl2 = spl.antiderivative(1)
        assert_allclose(spl2(0.6) - spl2(0.2),
                        spl.integral(0.2, 0.6))

    def test_nan(self):
        # bail out early if the input data contains nans
        x = np.arange(10, dtype=float)
        y = x**3
        w = np.ones_like(x)
        # also test LSQUnivariateSpline [which needs explicit knots]
        spl = UnivariateSpline(x, y, check_finite=True)
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        y_end = y[-1]
        for z in [np.nan, np.inf, -np.inf]:
            y[-1] = z
            assert_raises(ValueError, UnivariateSpline,
                    **dict(x=x, y=y, check_finite=True))
            assert_raises(ValueError, InterpolatedUnivariateSpline,
                    **dict(x=x, y=y, check_finite=True))
            assert_raises(ValueError, LSQUnivariateSpline,
                    **dict(x=x, y=y, t=t, check_finite=True))
            y[-1] = y_end  # check valid y but invalid w
            w[-1] = z
            assert_raises(ValueError, UnivariateSpline,
                    **dict(x=x, y=y, w=w, check_finite=True))
            assert_raises(ValueError, InterpolatedUnivariateSpline,
                    **dict(x=x, y=y, w=w, check_finite=True))
            assert_raises(ValueError, LSQUnivariateSpline,
                    **dict(x=x, y=y, t=t, w=w, check_finite=True))

    def test_increasing_x(self):
        xx = np.arange(10, dtype=float)
        yy = xx**3
        x = np.arange(10, dtype=float)
        x[1] = x[0]
        y = x**3
        w = np.ones_like(x)
        # also test LSQUnivariateSpline [which needs explicit knots]
        spl = UnivariateSpline(xx, yy, check_finite=True)
        t = spl.get_knots()[3:4]  # interior knots w/ default k=3
        assert_raises(ValueError, UnivariateSpline,
                **dict(x=x, y=y, check_finite=True))
        assert_raises(ValueError, InterpolatedUnivariateSpline,
                **dict(x=x, y=y, check_finite=True))
        assert_raises(ValueError, LSQUnivariateSpline,
                **dict(x=x, y=y, t=t, w=w, check_finite=True))


class TestLSQBivariateSpline(object):
    # NOTE: The systems in this test class are rank-deficient
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with suppress_warnings() as sup:
            r = sup.record(UserWarning, "\nThe coefficients of the spline")
            lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)
            assert_equal(len(r), 1)

        assert_almost_equal(lut(2,2), 3.)

    def test_bilinearity(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,7,8,3,4,7,1,3,4]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with suppress_warnings() as sup:
            # This seems to fail (ier=1, see ticket 1642).
            sup.filter(UserWarning, "\nThe coefficients of the spline")
            lut = LSQBivariateSpline(x,y,z,tx,ty,kx=1,ky=1)

        tx, ty = lut.get_knots()
        for xa, xb in zip(tx[:-1], tx[1:]):
            for ya, yb in zip(ty[:-1], ty[1:]):
                for t in [0.1, 0.5, 0.9]:
                    for s in [0.3, 0.4, 0.7]:
                        xp = xa*(1-t) + xb*t
                        yp = ya*(1-s) + yb*s
                        zp = (+ lut(xa, ya)*(1-t)*(1-s)
                              + lut(xb, ya)*t*(1-s)
                              + lut(xa, yb)*(1-t)*s
                              + lut(xb, yb)*t*s)
                        assert_almost_equal(lut(xp,yp), zp)

    def test_integral(self):
        x = [1,1,1,2,2,2,8,8,8]
        y = [1,2,3,1,2,3,1,2,3]
        z = array([0,7,8,3,4,7,1,3,4])

        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with suppress_warnings() as sup:
            r = sup.record(UserWarning, "\nThe coefficients of the spline")
            lut = LSQBivariateSpline(x, y, z, tx, ty, kx=1, ky=1)
            assert_equal(len(r), 1)
        tx, ty = lut.get_knots()
        tz = lut(tx, ty)
        trpz = .25*(diff(tx)[:,None]*diff(ty)[None,:]
                    * (tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()

        assert_almost_equal(lut.integral(tx[0], tx[-1], ty[0], ty[-1]),
                            trpz)

    def test_empty_input(self):
        # Test whether empty inputs returns an empty output. Ticket 1014
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        s = 0.1
        tx = [1+s,3-s]
        ty = [1+s,3-s]
        with suppress_warnings() as sup:
            r = sup.record(UserWarning, "\nThe coefficients of the spline")
            lut = LSQBivariateSpline(x, y, z, tx, ty, kx=1, ky=1)
            assert_equal(len(r), 1)

        assert_array_equal(lut([], []), np.zeros((0,0)))
        assert_array_equal(lut([], [], grid=False), np.zeros((0,)))


class TestSmoothBivariateSpline(object):
    def test_linear_constant(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [3,3,3,3,3,3,3,3,3]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[3,3,3,3])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[3,3],[3,3],[3,3]])

    def test_linear_1d(self):
        x = [1,1,1,2,2,2,3,3,3]
        y = [1,2,3,1,2,3,1,2,3]
        z = [0,0,0,2,2,2,4,4,4]
        lut = SmoothBivariateSpline(x,y,z,kx=1,ky=1)
        assert_array_almost_equal(lut.get_knots(),([1,1,3,3],[1,1,3,3]))
        assert_array_almost_equal(lut.get_coeffs(),[0,0,4,4])
        assert_almost_equal(lut.get_residual(),0.0)
        assert_array_almost_equal(lut([1,1.5,2],[1,1.5]),[[0,0],[1,1],[2,2]])

    def test_integral(self):
        x = [1,1,1,2,2,2,4,4,4]
        y = [1,2,3,1,2,3,1,2,3]
        z = array([0,7,8,3,4,7,1,3,4])

        with suppress_warnings() as sup:
            # This seems to fail (ier=1, see ticket 1642).
            sup.filter(UserWarning, "\nThe required storage space")
            lut = SmoothBivariateSpline(x, y, z, kx=1, ky=1, s=0)

        tx = [1,2,4]
        ty = [1,2,3]

        tz = lut(tx, ty)
        trpz = .25*(diff(tx)[:,None]*diff(ty)[None,:]
                    * (tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()
        assert_almost_equal(lut.integral(tx[0], tx[-1], ty[0], ty[-1]), trpz)

        lut2 = SmoothBivariateSpline(x, y, z, kx=2, ky=2, s=0)
        assert_almost_equal(lut2.integral(tx[0], tx[-1], ty[0], ty[-1]), trpz,
                            decimal=0)  # the quadratures give 23.75 and 23.85

        tz = lut(tx[:-1], ty[:-1])
        trpz = .25*(diff(tx[:-1])[:,None]*diff(ty[:-1])[None,:]
                    * (tz[:-1,:-1]+tz[1:,:-1]+tz[:-1,1:]+tz[1:,1:])).sum()
        assert_almost_equal(lut.integral(tx[0], tx[-2], ty[0], ty[-2]), trpz)

    def test_rerun_lwrk2_too_small(self):
        # in this setting, lwrk2 is too small in the default run. Here we
        # check for equality with the bisplrep/bisplev output because there,
        # an automatic re-run of the spline representation is done if ier>10.
        x = np.linspace(-2, 2, 80)
        y = np.linspace(-2, 2, 80)
        z = x + y
        xi = np.linspace(-1, 1, 100)
        yi = np.linspace(-2, 2, 100)
        tck = bisplrep(x, y, z)
        res1 = bisplev(xi, yi, tck)
        interp_ = SmoothBivariateSpline(x, y, z)
        res2 = interp_(xi, yi)
        assert_almost_equal(res1, res2)


class TestLSQSphereBivariateSpline(object):
    def setup_method(self):
        # define the input data and coordinates
        ntheta, nphi = 70, 90
        theta = linspace(0.5/(ntheta - 1), 1 - 0.5/(ntheta - 1), ntheta) * pi
        phi = linspace(0.5/(nphi - 1), 1 - 0.5/(nphi - 1), nphi) * 2. * pi
        data = ones((theta.shape[0], phi.shape[0]))
        # define knots and extract data values at the knots
        knotst = theta[::5]
        knotsp = phi[::5]
        knotdata = data[::5, ::5]
        # calculate spline coefficients
        lats, lons = meshgrid(theta, phi)
        lut_lsq = LSQSphereBivariateSpline(lats.ravel(), lons.ravel(),
                                           data.T.ravel(), knotst, knotsp)
        self.lut_lsq = lut_lsq
        self.data = knotdata
        self.new_lons, self.new_lats = knotsp, knotst

    def test_linear_constant(self):
        assert_almost_equal(self.lut_lsq.get_residual(), 0.0)
        assert_array_almost_equal(self.lut_lsq(self.new_lats, self.new_lons),
                                  self.data)

    def test_empty_input(self):
        assert_array_almost_equal(self.lut_lsq([], []), np.zeros((0,0)))
        assert_array_almost_equal(self.lut_lsq([], [], grid=False), np.zeros((0,)))


class TestSmoothSphereBivariateSpline(object):
    def setup_method(self):
        theta = array([.25*pi, .25*pi, .25*pi, .5*pi, .5*pi, .5*pi, .75*pi,
                       .75*pi, .75*pi])
        phi = array([.5 * pi, pi, 1.5 * pi, .5 * pi, pi, 1.5 * pi, .5 * pi, pi,
                     1.5 * pi])
        r = array([3, 3, 3, 3, 3, 3, 3, 3, 3])
        self.lut = SmoothSphereBivariateSpline(theta, phi, r, s=1E10)

    def test_linear_constant(self):
        assert_almost_equal(self.lut.get_residual(), 0.)
        assert_array_almost_equal(self.lut([1, 1.5, 2],[1, 1.5]),
                                  [[3, 3], [3, 3], [3, 3]])

    def test_empty_input(self):
        assert_array_almost_equal(self.lut([], []), np.zeros((0,0)))
        assert_array_almost_equal(self.lut([], [], grid=False), np.zeros((0,)))


class TestRectBivariateSpline(object):
    def test_defaults(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y),z)

    def test_evaluate(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)

        xi = [1, 2.3, 5.3, 0.5, 3.3, 1.2, 3]
        yi = [1, 3.3, 1.2, 4.0, 5.0, 1.0, 3]
        zi = lut.ev(xi, yi)
        zi2 = array([lut(xp, yp)[0,0] for xp, yp in zip(xi, yi)])

        assert_almost_equal(zi, zi2)

    def test_derivatives_grid(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        dx = array([[0,0,-20,0,0],[0,0,13,0,0],[0,0,4,0,0],
            [0,0,-11,0,0],[0,0,4,0,0]])/6.
        dy = array([[4,-1,0,1,-4],[4,-1,0,1,-4],[0,1.5,0,-1.5,0],
            [2,.25,0,-.25,-2],[4,-1,0,1,-4]])
        dxdy = array([[40,-25,0,25,-40],[-26,16.25,0,-16.25,26],
            [-8,5,0,-5,8],[22,-13.75,0,13.75,-22],[-8,5,0,-5,8]])/6.
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y,dx=1),dx)
        assert_array_almost_equal(lut(x,y,dy=1),dy)
        assert_array_almost_equal(lut(x,y,dx=1,dy=1),dxdy)

    def test_derivatives(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        dx = array([0,0,2./3,0,0])
        dy = array([4,-1,0,-.25,-4])
        dxdy = array([160,65,0,55,32])/24.
        lut = RectBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y,dx=1,grid=False),dx)
        assert_array_almost_equal(lut(x,y,dy=1,grid=False),dy)
        assert_array_almost_equal(lut(x,y,dx=1,dy=1,grid=False),dxdy)

    def test_broadcast(self):
        x = array([1,2,3,4,5])
        y = array([1,2,3,4,5])
        z = array([[1,2,1,2,1],[1,2,1,2,1],[1,2,3,2,1],[1,2,2,2,1],[1,2,1,2,1]])
        lut = RectBivariateSpline(x,y,z)
        assert_allclose(lut(x, y), lut(x[:,None], y[None,:], grid=False))


class TestRectSphereBivariateSpline(object):
    def test_defaults(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])
        lut = RectSphereBivariateSpline(x,y,z)
        assert_array_almost_equal(lut(x,y),z)

    def test_evaluate(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])
        lut = RectSphereBivariateSpline(x,y,z)
        yi = [0.2, 1, 2.3, 2.35, 3.0, 3.99, 5.25]
        xi = [1.5, 0.4, 1.1, 0.45, 0.2345, 1., 0.0001]
        zi = lut.ev(xi, yi)
        zi2 = array([lut(xp, yp)[0,0] for xp, yp in zip(xi, yi)])
        assert_almost_equal(zi, zi2)

    def test_derivatives_grid(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])

        lut = RectSphereBivariateSpline(x,y,z)

        y = linspace(0.02, 2*pi-0.02, 7)
        x = linspace(0.02, pi-0.02, 7)

        assert_allclose(lut(x, y, dtheta=1), _numdiff_2d(lut, x, y, dx=1),
                        rtol=1e-4, atol=1e-4)
        assert_allclose(lut(x, y, dphi=1), _numdiff_2d(lut, x, y, dy=1),
                        rtol=1e-4, atol=1e-4)
        assert_allclose(lut(x, y, dtheta=1, dphi=1), _numdiff_2d(lut, x, y, dx=1, dy=1, eps=1e-6),
                        rtol=1e-3, atol=1e-3)

    def test_derivatives(self):
        y = linspace(0.01, 2*pi-0.01, 7)
        x = linspace(0.01, pi-0.01, 7)
        z = array([[1,2,1,2,1,2,1],[1,2,1,2,1,2,1],[1,2,3,2,1,2,1],
                   [1,2,2,2,1,2,1],[1,2,1,2,1,2,1],[1,2,2,2,1,2,1],
                   [1,2,1,2,1,2,1]])

        lut = RectSphereBivariateSpline(x,y,z)

        y = linspace(0.02, 2*pi-0.02, 7)
        x = linspace(0.02, pi-0.02, 7)

        assert_equal(lut(x, y, dtheta=1, grid=False).shape, x.shape)
        assert_allclose(lut(x, y, dtheta=1, grid=False),
                        _numdiff_2d(lambda x,y: lut(x,y,grid=False), x, y, dx=1),
                        rtol=1e-4, atol=1e-4)
        assert_allclose(lut(x, y, dphi=1, grid=False),
                        _numdiff_2d(lambda x,y: lut(x,y,grid=False), x, y, dy=1),
                        rtol=1e-4, atol=1e-4)
        assert_allclose(lut(x, y, dtheta=1, dphi=1, grid=False),
                        _numdiff_2d(lambda x,y: lut(x,y,grid=False), x, y, dx=1, dy=1, eps=1e-6),
                        rtol=1e-3, atol=1e-3)


def _numdiff_2d(func, x, y, dx=0, dy=0, eps=1e-8):
    if dx == 0 and dy == 0:
        return func(x, y)
    elif dx == 1 and dy == 0:
        return (func(x + eps, y) - func(x - eps, y)) / (2*eps)
    elif dx == 0 and dy == 1:
        return (func(x, y + eps) - func(x, y - eps)) / (2*eps)
    elif dx == 1 and dy == 1:
        return (func(x + eps, y + eps) - func(x - eps, y + eps)
                - func(x + eps, y - eps) + func(x - eps, y - eps)) / (2*eps)**2
    else:
        raise ValueError("invalid derivative order")
