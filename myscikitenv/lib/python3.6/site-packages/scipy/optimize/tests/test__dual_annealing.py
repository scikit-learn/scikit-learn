# Dual annealing unit tests implementation.
# Copyright (c) 2018 Sylvain Gubian <sylvain.gubian@pmi.com>,
# Yang Xiang <yang.xiang@pmi.com>
# Author: Sylvain Gubian, PMP S.A.
"""
Unit tests for the dual annealing global optimizer
"""
from scipy.optimize import dual_annealing
from scipy.optimize._dual_annealing import VisitingDistribution
from scipy.optimize._dual_annealing import ObjectiveFunWrapper
from scipy.optimize._dual_annealing import EnergyState
from scipy.optimize._dual_annealing import LocalSearchWrapper
from scipy.optimize import rosen, rosen_der
import numpy as np
from numpy.testing import (assert_equal, TestCase, assert_allclose,
                           assert_array_less)
from pytest import raises as assert_raises
from scipy._lib._util import check_random_state


class TestDualAnnealing(TestCase):

    def setUp(self):
        # A function that returns always infinity for initialization tests
        self.weirdfunc = lambda x: np.inf
        # 2-D bounds for testing function
        self.ld_bounds = [(-5.12, 5.12)] * 2
        # 4-D bounds for testing function
        self.hd_bounds = self.ld_bounds * 4
        # Number of values to be generated for testing visit function
        self.nbtestvalues = 5000
        self.high_temperature = 5230
        self.low_temperature = 0.1
        self.qv = 2.62
        self.seed = 1234
        self.rs = check_random_state(self.seed)
        self.nb_fun_call = 0
        self.ngev = 0

    def tearDown(self):
        pass

    def callback(self, x, f, context):
        # For testing callback mechanism. Should stop for e <= 1 as
        # the callback function returns True
        if f <= 1.0:
            return True

    def func(self, x, args=()):
        # Using Rastrigin function for performing tests
        if args:
            shift = args
        else:
            shift = 0
        y = np.sum((x - shift) ** 2 - 10 * np.cos(2 * np.pi * (
            x - shift))) + 10 * np.size(x) + shift
        self.nb_fun_call += 1
        return y

    def rosen_der_wrapper(self, x, args=()):
        self.ngev += 1
        return rosen_der(x, *args)

    def test_visiting_stepping(self):
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        dim = lower.size
        vd = VisitingDistribution(lower, upper, self.qv, self.rs)
        values = np.zeros(dim)
        x_step_low = vd.visiting(values, 0, self.high_temperature)
        # Make sure that only the first component is changed
        assert_equal(np.not_equal(x_step_low, 0), True)
        values = np.zeros(dim)
        x_step_high = vd.visiting(values, dim, self.high_temperature)
        # Make sure that component other than at dim has changed
        assert_equal(np.not_equal(x_step_high[0], 0), True)

    def test_visiting_dist_high_temperature(self):
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        vd = VisitingDistribution(lower, upper, self.qv, self.rs)
        # values = np.zeros(self.nbtestvalues)
        # for i in np.arange(self.nbtestvalues):
        #     values[i] = vd.visit_fn(self.high_temperature)
        values = vd.visit_fn(self.high_temperature, self.nbtestvalues)

        # Visiting distribution is a distorted version of Cauchy-Lorentz
        # distribution, and as no 1st and higher moments (no mean defined,
        # no variance defined).
        # Check that big tails values are generated
        assert_array_less(np.min(values), 1e-10)
        assert_array_less(1e+10, np.max(values))

    def test_reset(self):
        owf = ObjectiveFunWrapper(self.weirdfunc)
        lu = list(zip(*self.ld_bounds))
        lower = np.array(lu[0])
        upper = np.array(lu[1])
        es = EnergyState(lower, upper)
        assert_raises(ValueError, es.reset, owf, check_random_state(None))

    def test_low_dim(self):
        ret = dual_annealing(
            self.func, self.ld_bounds, seed=self.seed)
        assert_allclose(ret.fun, 0., atol=1e-12)
        assert ret.success

    def test_high_dim(self):
        ret = dual_annealing(self.func, self.hd_bounds)
        assert_allclose(ret.fun, 0., atol=1e-12)
        assert ret.success

    def test_low_dim_no_ls(self):
        ret = dual_annealing(self.func, self.ld_bounds,
                             no_local_search=True)
        assert_allclose(ret.fun, 0., atol=1e-4)

    def test_high_dim_no_ls(self):
        ret = dual_annealing(self.func, self.hd_bounds,
                             no_local_search=True)
        assert_allclose(ret.fun, 0., atol=1e-4)

    def test_nb_fun_call(self):
        ret = dual_annealing(self.func, self.ld_bounds)
        assert_equal(self.nb_fun_call, ret.nfev)

    def test_nb_fun_call_no_ls(self):
        ret = dual_annealing(self.func, self.ld_bounds,
                no_local_search=True)
        assert_equal(self.nb_fun_call, ret.nfev)

    def test_max_reinit(self):
        assert_raises(ValueError, dual_annealing, self.weirdfunc,
                self.ld_bounds)

    def test_reproduce(self):
        seed = 1234
        res1 = dual_annealing(self.func, self.ld_bounds, seed=seed)
        res2 = dual_annealing(self.func, self.ld_bounds, seed=seed)
        res3 = dual_annealing(self.func, self.ld_bounds, seed=seed)
        # If we have reproducible results, x components found has to
        # be exactly the same, which is not the case with no seeding
        assert_equal(res1.x, res2.x)
        assert_equal(res1.x, res3.x)

    def test_bounds_integrity(self):
        wrong_bounds = [(-5.12, 5.12), (1, 0), (5.12, 5.12)]
        assert_raises(ValueError, dual_annealing, self.func,
                wrong_bounds)

    def test_bound_validity(self):
        invalid_bounds = [(-5, 5), (-np.inf, 0), (-5, 5)]
        assert_raises(ValueError, dual_annealing, self.func,
                invalid_bounds)
        invalid_bounds = [(-5, 5), (0, np.inf), (-5, 5)]
        assert_raises(ValueError, dual_annealing, self.func,
                invalid_bounds)
        invalid_bounds = [(-5, 5), (0, np.nan), (-5, 5)]
        assert_raises(ValueError, dual_annealing, self.func,
                invalid_bounds)

    def test_max_fun_ls(self):
        ret = dual_annealing(self.func, self.ld_bounds, maxfun=100)

        ls_max_iter = min(max(
            len(self.ld_bounds) * LocalSearchWrapper.LS_MAXITER_RATIO,
            LocalSearchWrapper.LS_MAXITER_MIN),
            LocalSearchWrapper.LS_MAXITER_MAX)
        assert ret.nfev <= 100 + ls_max_iter
        assert not ret.success

    def test_max_fun_no_ls(self):
        ret = dual_annealing(self.func, self.ld_bounds,
                             no_local_search=True, maxfun=500)
        assert ret.nfev <= 500
        assert not ret.success

    def test_maxiter(self):
        ret = dual_annealing(self.func, self.ld_bounds, maxiter=700)
        assert ret.nit <= 700

    # Testing that args are passed correctly for dual_annealing
    def test_fun_args_ls(self):
        ret = dual_annealing(self.func, self.ld_bounds,
                             args=((3.14159, )))
        assert_allclose(ret.fun, 3.14159, atol=1e-6)

    # Testing that args are passed correctly for pure simulated annealing
    def test_fun_args_no_ls(self):
        ret = dual_annealing(self.func, self.ld_bounds,
                             args=((3.14159, )), no_local_search=True)
        assert_allclose(ret.fun, 3.14159, atol=1e-4)

    def test_callback_stop(self):
        # Testing that callback make the algorithm stop for
        # fun value <= 1.0 (see callback method)
        ret = dual_annealing(self.func, self.ld_bounds,
                             callback=self.callback)
        assert ret.fun <= 1.0
        assert 'stop early' in ret.message[0]
        assert not ret.success

    def test_neldermed_ls_minimizer(self):
        minimizer_opts = {
            'method': 'Nelder-Mead',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-6)

    def test_powell_ls_minimizer(self):
        minimizer_opts = {
            'method': 'Powell',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-8)

    def test_cg_ls_minimizer(self):
        minimizer_opts = {
            'method': 'CG',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-8)

    def test_bfgs_ls_minimizer(self):
        minimizer_opts = {
            'method': 'BFGS',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-8)

    def test_tnc_ls_minimizer(self):
        minimizer_opts = {
            'method': 'TNC',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-8)

    def test_colyba_ls_minimizer(self):
        minimizer_opts = {
            'method': 'COBYLA',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-5)

    def test_slsqp_ls_minimizer(self):
        minimizer_opts = {
            'method': 'SLSQP',
        }
        ret = dual_annealing(self.func, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert_allclose(ret.fun, 0., atol=1e-7)

    def test_wrong_restart_temp(self):
        assert_raises(ValueError, dual_annealing, self.func,
                      self.ld_bounds, restart_temp_ratio=1)
        assert_raises(ValueError, dual_annealing, self.func,
                      self.ld_bounds, restart_temp_ratio=0)

    def test_gradient_gnev(self):
        minimizer_opts = {
            'jac': self.rosen_der_wrapper,
        }
        ret = dual_annealing(rosen, self.ld_bounds,
                             local_search_options=minimizer_opts)
        assert ret.njev == self.ngev
