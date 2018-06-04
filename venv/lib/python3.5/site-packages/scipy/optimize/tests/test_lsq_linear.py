import numpy as np
from numpy.linalg import lstsq
from numpy.testing import assert_allclose, assert_equal, assert_
from pytest import raises as assert_raises

from scipy.sparse import rand
from scipy.sparse.linalg import aslinearoperator
from scipy.optimize import lsq_linear


A = np.array([
    [0.171, -0.057],
    [-0.049, -0.248],
    [-0.166, 0.054],
])
b = np.array([0.074, 1.014, -0.383])


class BaseMixin(object):
    def setup_method(self):
        self.rnd = np.random.RandomState(0)

    def test_dense_no_bounds(self):
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, method=self.method, lsq_solver=lsq_solver)
            assert_allclose(res.x, lstsq(A, b, rcond=-1)[0])

    def test_dense_bounds(self):
        # Solutions for comparison are taken from MATLAB.
        lb = np.array([-1, -10])
        ub = np.array([1, 0])
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (lb, ub), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, lstsq(A, b, rcond=-1)[0])

        lb = np.array([0.0, -np.inf])
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (lb, np.inf), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, np.array([0.0, -4.084174437334673]),
                            atol=1e-6)

        lb = np.array([-1, 0])
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (lb, np.inf), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, np.array([0.448427311733504, 0]),
                            atol=1e-15)

        ub = np.array([np.inf, -5])
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (-np.inf, ub), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, np.array([-0.105560998682388, -5]))

        ub = np.array([-1, np.inf])
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (-np.inf, ub), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, np.array([-1, -4.181102129483254]))

        lb = np.array([0, -4])
        ub = np.array([1, 0])
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (lb, ub), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, np.array([0.005236663400791, -4]))

    def test_dense_rank_deficient(self):
        A = np.array([[-0.307, -0.184]])
        b = np.array([0.773])
        lb = [-0.1, -0.1]
        ub = [0.1, 0.1]
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (lb, ub), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.x, [-0.1, -0.1])

        A = np.array([
            [0.334, 0.668],
            [-0.516, -1.032],
            [0.192, 0.384],
        ])
        b = np.array([-1.436, 0.135, 0.909])
        lb = [0, -1]
        ub = [1, -0.5]
        for lsq_solver in self.lsq_solvers:
            res = lsq_linear(A, b, (lb, ub), method=self.method,
                             lsq_solver=lsq_solver)
            assert_allclose(res.optimality, 0, atol=1e-11)

    def test_full_result(self):
        lb = np.array([0, -4])
        ub = np.array([1, 0])
        res = lsq_linear(A, b, (lb, ub), method=self.method)

        assert_allclose(res.x, [0.005236663400791, -4])

        r = A.dot(res.x) - b
        assert_allclose(res.cost, 0.5 * np.dot(r, r))
        assert_allclose(res.fun, r)

        assert_allclose(res.optimality, 0.0, atol=1e-12)
        assert_equal(res.active_mask, [0, -1])
        assert_(res.nit < 15)
        assert_(res.status == 1 or res.status == 3)
        assert_(isinstance(res.message, str))
        assert_(res.success)


class SparseMixin(object):
    def test_sparse_and_LinearOperator(self):
        m = 5000
        n = 1000
        A = rand(m, n, random_state=0)
        b = self.rnd.randn(m)
        res = lsq_linear(A, b)
        assert_allclose(res.optimality, 0, atol=1e-6)

        A = aslinearoperator(A)
        res = lsq_linear(A, b)
        assert_allclose(res.optimality, 0, atol=1e-6)

    def test_sparse_bounds(self):
        m = 5000
        n = 1000
        A = rand(m, n, random_state=0)
        b = self.rnd.randn(m)
        lb = self.rnd.randn(n)
        ub = lb + 1
        res = lsq_linear(A, b, (lb, ub))
        assert_allclose(res.optimality, 0.0, atol=1e-8)

        res = lsq_linear(A, b, (lb, ub), lsmr_tol=1e-13)
        assert_allclose(res.optimality, 0.0, atol=1e-8)

        res = lsq_linear(A, b, (lb, ub), lsmr_tol='auto')
        assert_allclose(res.optimality, 0.0, atol=1e-8)


class TestTRF(BaseMixin, SparseMixin):
    method = 'trf'
    lsq_solvers = ['exact', 'lsmr']


class TestBVLS(BaseMixin):
    method = 'bvls'
    lsq_solvers = ['exact']

