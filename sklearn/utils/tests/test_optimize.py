import numpy as np
import pytest
from scipy.optimize import fmin_ncg

from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._testing import assert_array_almost_equal
from sklearn.utils.optimize import _newton_cg


def test_newton_cg():
    # Test that newton_cg gives same result as scipy's fmin_ncg

    rng = np.random.RandomState(0)
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

    assert_array_almost_equal(
        _newton_cg(grad_hess, func, grad, x0, tol=1e-10)[0],
        fmin_ncg(f=func, x0=x0, fprime=grad, fhess_p=hess),
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
