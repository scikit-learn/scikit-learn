import numpy as np
import scipy

from ...utils.optimize import _check_optimize_result


class LbfgsSolver:
    def __init__(
        self,
        linear_loss,
        coef,
        max_iter,
        verbose,
        tol,
        l2_reg_strength,
        n_threads,
    ):
        self.linear_loss = linear_loss
        self.coef = coef
        self.max_iter = max_iter
        self.verbose = verbose
        self.tol = tol
        self.l2_reg_strength = l2_reg_strength
        self.n_threads = n_threads

    def solve(self, X, y, sample_weight):
        opt_res = scipy.optimize.minimize(
            self.linear_loss.loss_gradient,
            self.coef,
            method="L-BFGS-B",
            jac=True,
            options={
                "maxiter": self.max_iter,
                "maxls": 50,  # default is 20
                "iprint": self.verbose - 1,
                "gtol": self.tol,
                # The constant 64 was found empirically to pass the test suite.
                # The point is that ftol is very small, but a bit larger than
                # machine precision for float64, which is the dtype used by lbfgs.
                "ftol": 64 * np.finfo(float).eps,
            },
            args=(X, y, sample_weight, self.l2_reg_strength, self.n_threads),
        )
        self.iteration = _check_optimize_result("lbfgs", opt_res)
        return opt_res.x
