"""
Unit test for Linear Programming via Simplex Algorithm.
"""
from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.testing import assert_, assert_allclose
from pytest import raises as assert_raises
from scipy.optimize._linprog_ip import _clean_inputs
from copy import deepcopy


def test_aliasing():
    c = 1
    A_ub = [[1]]
    b_ub = [1]
    A_eq = [[1]]
    b_eq = [1]
    bounds = (-np.inf, np.inf)

    c_copy = deepcopy(c)
    A_ub_copy = deepcopy(A_ub)
    b_ub_copy = deepcopy(b_ub)
    A_eq_copy = deepcopy(A_eq)
    b_eq_copy = deepcopy(b_eq)
    bounds_copy = deepcopy(bounds)

    _clean_inputs(c, A_ub, b_ub, A_eq, b_eq, bounds)

    assert_(c == c_copy, "c modified by _clean_inputs")
    assert_(A_ub == A_ub_copy, "A_ub modified by _clean_inputs")
    assert_(b_ub == b_ub_copy, "b_ub modified by _clean_inputs")
    assert_(A_eq == A_eq_copy, "A_eq modified by _clean_inputs")
    assert_(b_eq == b_eq_copy, "b_eq modified by _clean_inputs")
    assert_(bounds == bounds_copy, "bounds modified by _clean_inputs")


def test_aliasing2():
    c = np.array([1, 1])
    A_ub = np.array([[1, 1], [2, 2]])
    b_ub = np.array([[1], [1]])
    A_eq = np.array([[1, 1]])
    b_eq = np.array([1])
    bounds = [(-np.inf, np.inf), (None, 1)]

    c_copy = c.copy()
    A_ub_copy = A_ub.copy()
    b_ub_copy = b_ub.copy()
    A_eq_copy = A_eq.copy()
    b_eq_copy = b_eq.copy()
    bounds_copy = deepcopy(bounds)

    _clean_inputs(c, A_ub, b_ub, A_eq, b_eq, bounds)

    assert_allclose(c, c_copy, err_msg="c modified by _clean_inputs")
    assert_allclose(A_ub, A_ub_copy, err_msg="A_ub modified by _clean_inputs")
    assert_allclose(b_ub, b_ub_copy, err_msg="b_ub modified by _clean_inputs")
    assert_allclose(A_eq, A_eq_copy, err_msg="A_eq modified by _clean_inputs")
    assert_allclose(b_eq, b_eq_copy, err_msg="b_eq modified by _clean_inputs")
    assert_(bounds == bounds_copy, "bounds modified by _clean_inputs")


def test_missing_inputs():
    c = [1, 2]
    A_ub = np.array([[1, 1], [2, 2]])
    b_ub = np.array([1, 1])
    A_eq = np.array([[1, 1], [2, 2]])
    b_eq = np.array([1, 1])

    assert_raises(TypeError, _clean_inputs)
    assert_raises(TypeError, _clean_inputs, c=None)
    assert_raises(ValueError, _clean_inputs, c=c, A_ub=A_ub)
    assert_raises(ValueError, _clean_inputs, c=c, A_ub=A_ub, b_ub=None)
    assert_raises(ValueError, _clean_inputs, c=c, b_ub=b_ub)
    assert_raises(ValueError, _clean_inputs, c=c, A_ub=None, b_ub=b_ub)
    assert_raises(ValueError, _clean_inputs, c=c, A_eq=A_eq)
    assert_raises(ValueError, _clean_inputs, c=c, A_eq=A_eq, b_eq=None)
    assert_raises(ValueError, _clean_inputs, c=c, b_eq=b_eq)
    assert_raises(ValueError, _clean_inputs, c=c, A_eq=None, b_eq=b_eq)


def test_too_many_dimensions():
    cb = [1, 2, 3, 4]
    A = np.random.rand(4, 4)
    bad2D = [[1, 2], [3, 4]]
    bad3D = np.random.rand(4, 4, 4)
    assert_raises(ValueError, _clean_inputs, c=bad2D, A_ub=A, b_ub=cb)
    assert_raises(ValueError, _clean_inputs, c=cb, A_ub=bad3D, b_ub=cb)
    assert_raises(ValueError, _clean_inputs, c=cb, A_ub=A, b_ub=bad2D)
    assert_raises(ValueError, _clean_inputs, c=cb, A_eq=bad3D, b_eq=cb)
    assert_raises(ValueError, _clean_inputs, c=cb, A_eq=A, b_eq=bad2D)


def test_too_few_dimensions():
    bad = np.random.rand(4, 4).ravel()
    cb = np.random.rand(4)
    assert_raises(ValueError, _clean_inputs, c=cb, A_ub=bad, b_ub=cb)
    assert_raises(ValueError, _clean_inputs, c=cb, A_eq=bad, b_eq=cb)


def test_inconsistent_dimensions():
    m = 2
    n = 4
    c = [1, 2, 3, 4]

    Agood = np.random.rand(m, n)
    Abad = np.random.rand(m, n + 1)
    bgood = np.random.rand(m)
    bbad = np.random.rand(m + 1)
    boundsbad = [(0, 1)] * (n + 1)
    assert_raises(ValueError, _clean_inputs, c=c, A_ub=Abad, b_ub=bgood)
    assert_raises(ValueError, _clean_inputs, c=c, A_ub=Agood, b_ub=bbad)
    assert_raises(ValueError, _clean_inputs, c=c, A_eq=Abad, b_eq=bgood)
    assert_raises(ValueError, _clean_inputs, c=c, A_eq=Agood, b_eq=bbad)
    assert_raises(ValueError, _clean_inputs, c=c, bounds=boundsbad)


def test_type_errors():
    bad = "hello"
    c = [1, 2]
    A_ub = np.array([[1, 1], [2, 2]])
    b_ub = np.array([1, 1])
    A_eq = np.array([[1, 1], [2, 2]])
    b_eq = np.array([1, 1])
    bounds = [(0, 1)]
    assert_raises(
        TypeError,
        _clean_inputs,
        c=bad,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=bad,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=bad,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=bad,
        b_eq=b_eq,
        bounds=bounds)

    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bad)
    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds="hi")
    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=["hi"])
    assert_raises(
        TypeError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=[
            ("hi")])
    assert_raises(TypeError, _clean_inputs, c=c, A_ub=A_ub,
                  b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(1, "")])
    assert_raises(TypeError, _clean_inputs, c=c, A_ub=A_ub,
                  b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(1, 2), (1, "")])


def test_non_finite_errors():
    c = [1, 2]
    A_ub = np.array([[1, 1], [2, 2]])
    b_ub = np.array([1, 1])
    A_eq = np.array([[1, 1], [2, 2]])
    b_eq = np.array([1, 1])
    bounds = [(0, 1)]
    assert_raises(
        ValueError, _clean_inputs, c=[0, None], A_ub=A_ub, b_ub=b_ub,
        A_eq=A_eq, b_eq=b_eq, bounds=bounds)
    assert_raises(
        ValueError, _clean_inputs, c=[np.inf, 0], A_ub=A_ub, b_ub=b_ub,
        A_eq=A_eq, b_eq=b_eq, bounds=bounds)
    assert_raises(
        ValueError, _clean_inputs, c=[0, -np.inf], A_ub=A_ub, b_ub=b_ub,
        A_eq=A_eq, b_eq=b_eq, bounds=bounds)
    assert_raises(
        ValueError, _clean_inputs, c=[np.nan, 0], A_ub=A_ub, b_ub=b_ub,
        A_eq=A_eq, b_eq=b_eq, bounds=bounds)

    assert_raises(ValueError, _clean_inputs, c=c, A_ub=[[1, 2], [None, 1]],
                  b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds)
    assert_raises(
        ValueError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=[
            np.inf,
            1],
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_raises(ValueError, _clean_inputs, c=c, A_ub=A_ub, b_ub=b_ub, A_eq=[
                  [1, 2], [1, -np.inf]], b_eq=b_eq, bounds=bounds)
    assert_raises(
        ValueError,
        _clean_inputs,
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=[
            1,
            np.nan],
        bounds=bounds)


def test__clean_inputs1():
    c = [1, 2]
    A_ub = [[1, 1], [2, 2]]
    b_ub = [1, 1]
    A_eq = [[1, 1], [2, 2]]
    b_eq = [1, 1]
    bounds = None
    outputs = _clean_inputs(
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_allclose(outputs[0], np.array(c))
    assert_allclose(outputs[1], np.array(A_ub))
    assert_allclose(outputs[2], np.array(b_ub))
    assert_allclose(outputs[3], np.array(A_eq))
    assert_allclose(outputs[4], np.array(b_eq))
    assert_(outputs[5] == [(0, None)] * 2, "")

    assert_(outputs[0].shape == (2,), "")
    assert_(outputs[1].shape == (2, 2), "")
    assert_(outputs[2].shape == (2,), "")
    assert_(outputs[3].shape == (2, 2), "")
    assert_(outputs[4].shape == (2,), "")


def test__clean_inputs2():
    c = 1
    A_ub = [[1]]
    b_ub = 1
    A_eq = [[1]]
    b_eq = 1
    bounds = (0, 1)
    outputs = _clean_inputs(
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_allclose(outputs[0], np.array(c))
    assert_allclose(outputs[1], np.array(A_ub))
    assert_allclose(outputs[2], np.array(b_ub))
    assert_allclose(outputs[3], np.array(A_eq))
    assert_allclose(outputs[4], np.array(b_eq))
    assert_(outputs[5] == [(0, 1)], "")

    assert_(outputs[0].shape == (1,), "")
    assert_(outputs[1].shape == (1, 1), "")
    assert_(outputs[2].shape == (1,), "")
    assert_(outputs[3].shape == (1, 1), "")
    assert_(outputs[4].shape == (1,), "")


def test__clean_inputs3():
    c = [[1, 2]]
    A_ub = np.random.rand(2, 2)
    b_ub = [[1], [2]]
    A_eq = np.random.rand(2, 2)
    b_eq = [[1], [2]]
    bounds = [(0, 1)]
    outputs = _clean_inputs(
        c=c,
        A_ub=A_ub,
        b_ub=b_ub,
        A_eq=A_eq,
        b_eq=b_eq,
        bounds=bounds)
    assert_allclose(outputs[0], np.array([1, 2]))
    assert_allclose(outputs[2], np.array([1, 2]))
    assert_allclose(outputs[4], np.array([1, 2]))
    assert_(outputs[5] == [(0, 1)] * 2, "")

    assert_(outputs[0].shape == (2,), "")
    assert_(outputs[2].shape == (2,), "")
    assert_(outputs[4].shape == (2,), "")


def test_bad_bounds():
    c = [1, 2]
    assert_raises(ValueError, _clean_inputs, c=c, bounds=(1, -2))
    assert_raises(ValueError, _clean_inputs, c=c, bounds=[(1, -2)])
    assert_raises(ValueError, _clean_inputs, c=c, bounds=[(1, -2), (1, 2)])

    assert_raises(ValueError, _clean_inputs, c=c, bounds=(1, 2, 2))
    assert_raises(ValueError, _clean_inputs, c=c, bounds=[(1, 2, 2)])
    assert_raises(ValueError, _clean_inputs, c=c, bounds=[(1, 2), (1, 2, 2)])
    assert_raises(ValueError, _clean_inputs, c=c,
                  bounds=[(1, 2), (1, 2), (1, 2)])


def test_good_bounds():
    c = [1, 2]
    outputs = _clean_inputs(c=c, bounds=None)
    assert_(outputs[5] == [(0, None)] * 2, "")

    outputs = _clean_inputs(c=c, bounds=(1, 2))
    assert_(outputs[5] == [(1, 2)] * 2, "")

    outputs = _clean_inputs(c=c, bounds=[(1, 2)])
    assert_(outputs[5] == [(1, 2)] * 2, "")

    outputs = _clean_inputs(c=c, bounds=[(1, np.inf)])
    assert_(outputs[5] == [(1, None)] * 2, "")

    outputs = _clean_inputs(c=c, bounds=[(-np.inf, 1)])
    assert_(outputs[5] == [(None, 1)] * 2, "")

    outputs = _clean_inputs(c=c, bounds=[(-np.inf, np.inf), (-np.inf, np.inf)])
    assert_(outputs[5] == [(None, None)] * 2, "")
