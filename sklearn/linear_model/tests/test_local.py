from unittest import TestCase
from sklearn.linear_model._local import LocalRegression
import numpy as np


class TestLocalRegression(TestCase):

    def setUp(self):
        self.n = 10
        self.x = np.linspace(-10, 10, self.n)
        self.X = self.x[:, np.newaxis]
        np.random.seed(1)
        self.e = np.random.normal(0, 0.01, self.n)

    def test_nw_analytic_vs_qr_vs_leastsq_1d(self):
        self.y = self.x**2 + self.e
        self.nw = LocalRegression(degree=0, warm_start=False)
        self.nw.fit(self.X, self.y)
        X_eval = np.linspace(-10, 10, 2 * self.n)[:, np.newaxis]
        y_pred_analytic = self.nw.predict(X_eval, method='analytic')
        y_pred_qr = self.nw.predict(X_eval, method='qr')
        y_pred_leastsq = self.nw.predict(X_eval, method='leastsq')
        np.testing.assert_almost_equal(y_pred_analytic, y_pred_qr, decimal=6)
        np.testing.assert_almost_equal(y_pred_qr, y_pred_leastsq, decimal=6)

    def test_nw_analytic_vs_qr_vs_leastsq_2d(self):
        grid = np.linspace(-10, 10, self.n)
        x1, x2 = np.meshgrid(grid, grid)
        self.X = np.vstack((x1.flatten(), x2.flatten())).T
        self.y = self.X[:, 0] * self.X[:, 1] + np.random.normal(0, 0.01, self.X.shape[0])
        self.nw = LocalRegression(degree=0, warm_start=False)
        self.nw.fit(self.X, self.y)
        grid_eval = np.linspace(-10, 10, 2 * self.n)[:, np.newaxis]
        x1e, x2e = np.meshgrid(grid_eval, grid_eval)
        X_eval = np.vstack((x1e.flatten(), x2e.flatten())).T
        y_pred_analytic = self.nw.predict(X_eval, method='analytic')
        y_pred_qr = self.nw.predict(X_eval, method='qr')
        y_pred_leastsq = self.nw.predict(X_eval, method='leastsq')
        np.testing.assert_almost_equal(y_pred_qr, y_pred_leastsq, decimal=5)
        np.testing.assert_almost_equal(y_pred_analytic, y_pred_qr, decimal=5)

    def test_ll_analytic_vs_cost(self):
        self.y = self.x + self.e
        self.ll = LocalRegression(degree=1, warm_start=False)
        self.ll.fit(self.X, self.y)
        X_eval = np.linspace(-10, 10, 2 * self.n)[:, np.newaxis]
        y_pred_analytic = self.ll.predict(X_eval, method='analytic')
        y_pred_qr = self.ll.predict(X_eval, method='qr')
        y_pred_leastsq = self.ll.predict(X_eval, method='leastsq')
        np.testing.assert_almost_equal(y_pred_analytic, y_pred_qr)
        np.testing.assert_almost_equal(y_pred_qr, y_pred_leastsq)

    def test_nw_analytic_loc_const_const_1d(self):
        self.y = np.ones_like(self.x) + self.e
        self.nw = LocalRegression(degree=0)
        self.nw.fit(self.X, self.y)
        np.testing.assert_almost_equal(self.nw.predict(self.X), self.y)

    def test_ll_analytic_loc_lin_lin_1d(self):
        self.y = self.x + self.e
        self.ll = LocalRegression(degree=1)
        self.ll.fit(self.X, self.y)
        np.testing.assert_almost_equal(self.ll.predict(self.X), self.y, decimal=1)

    def test_2d_loc_const_const(self):
        n = 5
        c = 7.
        x = np.linspace(-5, 5, n)
        y = np.linspace(-5, 5, n)
        z = np.array([[c for xx in x] for yy in y])
        x, y = np.meshgrid(x, y)
        X = np.array(list(zip(x.flatten(), y.flatten())))
        np.random.seed(1)
        e = np.random.normal(0, 0.01, (n, n))
        z = z + e
        locreg3d = LocalRegression(degree=0).fit(X, z.flatten())
        z_pred = locreg3d.predict(X)
        self.assertTrue(np.all(np.abs(z_pred / c - 1) <= 0.01))

    def test_2d_loc_lin_lin(self):
        n = 5
        x = np.linspace(1, 10, n)
        y = np.linspace(1, 10, n)
        z = np.array([[xx + yy for xx in x] for yy in y])
        res = z.flatten()
        x, y = np.meshgrid(x, y)
        X = np.array(list(zip(x.flatten(), y.flatten())))
        np.random.seed(1)
        e = np.random.normal(0, 0.01, (n, n))
        z = z + e
        locreg3d = LocalRegression(degree=1).fit(X, z.flatten())
        z_pred = locreg3d.predict(X)
        self.assertTrue(np.all(np.abs(z_pred / res - 1) <= 0.01))

    def test_2d_least_sq_vs_qr(self):
        n = 5
        x1 = np.linspace(1, 10, n)
        x2 = np.linspace(1, 10, n)
        y = np.array([[xx1 ** 2 + xx2 ** 2 for xx1 in x1] for xx2 in x2])
        x1, x2 = np.meshgrid(x1, x2)
        X = np.array(list(zip(x1.flatten(), x2.flatten())))
        np.random.seed(1)
        e = np.random.normal(0, 0.01, (n, n))
        y = y + e
        locreg3d = LocalRegression(degree=2).fit(X, y.flatten())
        y_pred_leastsq = locreg3d.predict(X, method='leastsq')
        y_pred_qr = locreg3d.predict(X, method='qr')
        np.testing.assert_array_almost_equal(y_pred_leastsq, y_pred_qr)

    def test_fit_partial(self):
        self.y = self.x + self.e
        self.ll = LocalRegression(degree=1, warm_start=True)
        self.ll.fit(self.X, self.y)
        X_eval = np.linspace(-10, 10, 2 * self.n)[:, np.newaxis]
        _ = self.ll.predict(X_eval, method='analytic')
        self.y = self.x ** 2 + self.e
        h = self.ll.bandwidth
        self.ll.fit_partial(self.y, h)
        y_pred_partial = self.ll.predict(X_eval, method='analytic')
        self.ll2 = LocalRegression(degree=1, warm_start=False)
        self.ll2.fit(self.X, self.y, h)
        y_pred = self.ll2.predict(X_eval, method='analytic')
        np.testing.assert_almost_equal(y_pred, y_pred_partial)
        # test against new instance to validate correct state
        self.ll3 = LocalRegression(degree=1, warm_start=False)
        self.ll3.fit(self.X, self.y, h)
        y_pred_new = self.ll3.predict(X_eval, method='analytic')
        np.testing.assert_array_almost_equal(y_pred_new, y_pred)
