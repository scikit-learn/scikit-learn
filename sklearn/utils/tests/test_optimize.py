import warnings

import numpy as np
import pytest
from scipy.optimize import fmin_ncg

from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._bunch import Bunch
from sklearn.utils._testing import assert_allclose
from sklearn.utils.optimize import _check_optimize_result, _newton_cg


def test_newton_cg(global_random_seed):
    # Test that newton_cg gives same result as scipy's fmin_ncg

    rng = np.random.RandomState(global_random_seed)
    A = rng.normal(size=(10, 10))
    x0 = np.ones(10)

    def func(x):
        Ax = A.dot(x)
        return 0.5 * (Ax).dot(Ax)

    def grad(x):
        return A.T.dot(A.dot(x))

    def hess(x, p):
        return p.dot(A.T.dot(A.dot(x.all())))

    def grad_hess(x):
        return grad(x), lambda x: A.T.dot(A.dot(x))

    # func is a definite positive quadratic form, so the minimum is at x = 0
    # hence the use of absolute tolerance.
    assert np.all(np.abs(_newton_cg(grad_hess, func, grad, x0, tol=1e-10)[0]) <= 1e-7)
    assert_allclose(
        _newton_cg(grad_hess, func, grad, x0, tol=1e-7)[0],
        fmin_ncg(f=func, x0=x0, fprime=grad, fhess_p=hess),
        atol=1e-5,
    )


@pytest.mark.parametrize("verbose", [0, 1, 2])
def test_newton_cg_verbosity(capsys, verbose):
    """Test the std output of verbose newton_cg solver."""
    A = np.eye(2)
    b = np.array([1, 2], dtype=float)

    _newton_cg(
        grad_hess=lambda x: (A @ x - b, lambda z: A @ z),
        func=lambda x: 0.5 * x @ A @ x - b @ x,
        grad=lambda x: A @ x - b,
        x0=np.zeros(A.shape[0]),
        verbose=verbose,
    )  # returns array([1., 2])
    captured = capsys.readouterr()

    if verbose == 0:
        assert captured.out == ""
    else:
        msg = [
            "Newton-CG iter = 1",
            "Check Convergence",
            "max |gradient|",
            "Solver did converge at loss = ",
        ]
        for m in msg:
            assert m in captured.out

    if verbose >= 2:
        msg = [
            "Inner CG solver iteration 1 stopped with",
            "sum(|residuals|) <= tol",
            "Line Search",
            "try line search wolfe1",
            "wolfe1 line search was successful",
        ]
        for m in msg:
            assert m in captured.out

    if verbose >= 2:
        # Set up a badly scaled singular Hessian with a completely wrong starting
        # position. This should trigger 2nd line search check
        A = np.array([[1.0, 2], [2, 4]]) * 1e30  # collinear columns
        b = np.array([1.0, 2.0])
        # Note that scipy.optimize._linesearch LineSearchWarning inherits from
        # RuntimeWarning, but we do not want to import from non public APIs.
        with pytest.warns(RuntimeWarning):
            _newton_cg(
                grad_hess=lambda x: (A @ x - b, lambda z: A @ z),
                func=lambda x: 0.5 * x @ A @ x - b @ x,
                grad=lambda x: A @ x - b,
                x0=np.array([-2.0, 1]),  # null space of hessian
                verbose=verbose,
            )
        captured = capsys.readouterr()
        msg = [
            "wolfe1 line search was not successful",
            "check loss |improvement| <= eps * |loss_old|:",
            "check sum(|gradient|) < sum(|gradient_old|):",
            "last resort: try line search wolfe2",
        ]
        for m in msg:
            assert m in captured.out

        # Set up a badly conditioned Hessian that leads to tiny curvature.
        # X.T @ X have singular values array([1.00000400e+01, 1.00008192e-11])
        A = np.array([[1.0, 2], [1, 2 + 1e-15]])
        b = np.array([-2.0, 1])
        with pytest.warns(ConvergenceWarning):
            _newton_cg(
                grad_hess=lambda x: (A @ x - b, lambda z: A @ z),
                func=lambda x: 0.5 * x @ A @ x - b @ x,
                grad=lambda x: A @ x - b,
                x0=b,
                verbose=verbose,
                maxiter=2,
            )
        captured = capsys.readouterr()
        msg = [
            "tiny_|p| = eps * ||p||^2",
        ]
        for m in msg:
            assert m in captured.out

        # Test for a case with negative Hessian.
        # We do not trigger "Inner CG solver iteration {i} stopped with negative
        # curvature", but that is very hard to trigger.
        A = np.eye(2)
        b = np.array([-2.0, 1])
        with pytest.warns(RuntimeWarning):
            _newton_cg(
                # Note the wrong sign in the hessian product.
                grad_hess=lambda x: (A @ x - b, lambda z: -A @ z),
                func=lambda x: 0.5 * x @ A @ x - b @ x,
                grad=lambda x: A @ x - b,
                x0=np.array([1.0, 1.0]),
                verbose=verbose,
                maxiter=3,
            )
        captured = capsys.readouterr()
        msg = [
            "Inner CG solver iteration 0 fell back to steepest descent",
        ]
        for m in msg:
            assert m in captured.out

        A = np.diag([1e-3, 1, 1e3])
        b = np.array([-2.0, 1, 2.0])
        with pytest.warns(ConvergenceWarning):
            _newton_cg(
                grad_hess=lambda x: (A @ x - b, lambda z: A @ z),
                func=lambda x: 0.5 * x @ A @ x - b @ x,
                grad=lambda x: A @ x - b,
                x0=np.ones_like(b),
                verbose=verbose,
                maxiter=2,
                maxinner=1,
            )
        captured = capsys.readouterr()
        msg = [
            "Inner CG solver stopped reaching maxiter=1",
        ]
        for m in msg:
            assert m in captured.out


def test_check_optimize():
    # Mock some lbfgs output using a Bunch instance:
    result = Bunch()

    # First case: no warnings
    result.nit = 1
    result.status = 0
    result.message = "OK"

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _check_optimize_result("lbfgs", result)

    # Second case: warning about implicit `max_iter`: do not recommend the user
    # to increase `max_iter` this is not a user settable parameter.
    result.status = 1
    result.message = "STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT"
    with pytest.warns(ConvergenceWarning) as record:
        _check_optimize_result("lbfgs", result)

    assert len(record) == 1
    warn_msg = record[0].message.args[0]
    assert "lbfgs failed to converge after 1 iteration(s)" in warn_msg
    assert result.message in warn_msg
    assert "Increase the number of iterations" not in warn_msg
    assert "scale the data" in warn_msg

    # Third case: warning about explicit `max_iter`: recommend user to increase
    # `max_iter`.
    with pytest.warns(ConvergenceWarning) as record:
        _check_optimize_result("lbfgs", result, max_iter=1)

    assert len(record) == 1
    warn_msg = record[0].message.args[0]
    assert "lbfgs failed to converge after 1 iteration(s)" in warn_msg
    assert result.message in warn_msg
    assert "Increase the number of iterations" in warn_msg
    assert "scale the data" in warn_msg

    # Fourth case: other convergence problem before reaching `max_iter`: do not
    # recommend increasing `max_iter`.
    result.nit = 2
    result.status = 2
    result.message = "ABNORMAL"
    with pytest.warns(ConvergenceWarning) as record:
        _check_optimize_result("lbfgs", result, max_iter=10)

    assert len(record) == 1
    warn_msg = record[0].message.args[0]
    assert "lbfgs failed to converge after 2 iteration(s)" in warn_msg
    assert result.message in warn_msg
    assert "Increase the number of iterations" not in warn_msg
    assert "scale the data" in warn_msg
