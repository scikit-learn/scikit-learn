from __future__ import division, print_function, absolute_import
import pytest
import numpy as np
from numpy.testing import TestCase, assert_array_equal
import scipy.sparse as sps
from scipy.optimize._constraints import (
    Bounds, LinearConstraint, NonlinearConstraint, PreparedConstraint,
    new_bounds_to_old, old_bound_to_new, strict_bounds)


class TestStrictBounds(TestCase):
    def test_scalarvalue_unique_enforce_feasibility(self):
        m = 3
        lb = 2
        ub = 4
        enforce_feasibility = False
        strict_lb, strict_ub = strict_bounds(lb, ub,
                                             enforce_feasibility,
                                             m)
        assert_array_equal(strict_lb, [-np.inf, -np.inf, -np.inf])
        assert_array_equal(strict_ub, [np.inf, np.inf, np.inf])

        enforce_feasibility = True
        strict_lb, strict_ub = strict_bounds(lb, ub,
                                             enforce_feasibility,
                                             m)
        assert_array_equal(strict_lb, [2, 2, 2])
        assert_array_equal(strict_ub, [4, 4, 4])

    def test_vectorvalue_unique_enforce_feasibility(self):
        m = 3
        lb = [1, 2, 3]
        ub = [4, 5, 6]
        enforce_feasibility = False
        strict_lb, strict_ub = strict_bounds(lb, ub,
                                              enforce_feasibility,
                                              m)
        assert_array_equal(strict_lb, [-np.inf, -np.inf, -np.inf])
        assert_array_equal(strict_ub, [np.inf, np.inf, np.inf])

        enforce_feasibility = True
        strict_lb, strict_ub = strict_bounds(lb, ub,
                                              enforce_feasibility,
                                              m)
        assert_array_equal(strict_lb, [1, 2, 3])
        assert_array_equal(strict_ub, [4, 5, 6])

    def test_scalarvalue_vector_enforce_feasibility(self):
        m = 3
        lb = 2
        ub = 4
        enforce_feasibility = [False, True, False]
        strict_lb, strict_ub = strict_bounds(lb, ub,
                                             enforce_feasibility,
                                             m)
        assert_array_equal(strict_lb, [-np.inf, 2, -np.inf])
        assert_array_equal(strict_ub, [np.inf, 4, np.inf])

    def test_vectorvalue_vector_enforce_feasibility(self):
        m = 3
        lb = [1, 2, 3]
        ub = [4, 6, np.inf]
        enforce_feasibility = [True, False, True]
        strict_lb, strict_ub = strict_bounds(lb, ub,
                                             enforce_feasibility,
                                             m)
        assert_array_equal(strict_lb, [1, -np.inf, 3])
        assert_array_equal(strict_ub, [4, np.inf, np.inf])


def test_prepare_constraint_infeasible_x0():
    lb = np.array([0, 20, 30])
    ub = np.array([0.5, np.inf, 70])
    x0 = np.array([1, 2, 3])
    enforce_feasibility = np.array([False, True, True], dtype=bool)
    bounds = Bounds(lb, ub, enforce_feasibility)
    pytest.raises(ValueError, PreparedConstraint, bounds, x0)

    x0 = np.array([1, 2, 3, 4])
    A = np.array([[1, 2, 3, 4], [5, 0, 0, 6], [7, 0, 8, 0]])
    enforce_feasibility = np.array([True, True, True], dtype=bool)
    linear = LinearConstraint(A, -np.inf, 0, enforce_feasibility)
    pytest.raises(ValueError, PreparedConstraint, linear, x0)

    def fun(x):
        return A.dot(x)

    def jac(x):
        return A

    def hess(x, v):
        return sps.csr_matrix((4, 4))

    nonlinear = NonlinearConstraint(fun, -np.inf, 0, jac, hess,
                                    enforce_feasibility)
    pytest.raises(ValueError, PreparedConstraint, nonlinear, x0)


def test_new_bounds_to_old():
    lb = np.array([-np.inf, 2, 3])
    ub = np.array([3, np.inf, 10])

    bounds = [(None, 3), (2, None), (3, 10)]
    assert_array_equal(new_bounds_to_old(lb, ub, 3), bounds)

    bounds_single_lb = [(-1, 3), (-1, None), (-1, 10)]
    assert_array_equal(new_bounds_to_old(-1, ub, 3), bounds_single_lb)

    bounds_no_lb = [(None, 3), (None, None), (None, 10)]
    assert_array_equal(new_bounds_to_old(-np.inf, ub, 3), bounds_no_lb)

    bounds_single_ub = [(None, 20), (2, 20), (3, 20)]
    assert_array_equal(new_bounds_to_old(lb, 20, 3), bounds_single_ub)

    bounds_no_ub = [(None, None), (2, None), (3, None)]
    assert_array_equal(new_bounds_to_old(lb, np.inf, 3), bounds_no_ub)

    bounds_single_both = [(1, 2), (1, 2), (1, 2)]
    assert_array_equal(new_bounds_to_old(1, 2, 3), bounds_single_both)

    bounds_no_both = [(None, None), (None, None), (None, None)]
    assert_array_equal(new_bounds_to_old(-np.inf, np.inf, 3), bounds_no_both)


def test_old_bounds_to_new():
    bounds = ([1, 2], (None, 3), (-1, None))
    lb_true = np.array([1, -np.inf, -1])
    ub_true = np.array([2, 3, np.inf])

    lb, ub = old_bound_to_new(bounds)
    assert_array_equal(lb, lb_true)
    assert_array_equal(ub, ub_true)
