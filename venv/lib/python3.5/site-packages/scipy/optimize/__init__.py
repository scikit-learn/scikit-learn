"""
=====================================================
Optimization and root finding (:mod:`scipy.optimize`)
=====================================================

.. currentmodule:: scipy.optimize

Optimization
============

Local Optimization
------------------

.. autosummary::
   :toctree: generated/

   minimize - Interface for minimizers of multivariate functions
   minimize_scalar - Interface for minimizers of univariate functions
   OptimizeResult - The optimization result returned by some optimizers
   OptimizeWarning - The optimization encountered problems

The `minimize` function supports the following methods:

.. toctree::

   optimize.minimize-neldermead
   optimize.minimize-powell
   optimize.minimize-cg
   optimize.minimize-bfgs
   optimize.minimize-newtoncg
   optimize.minimize-lbfgsb
   optimize.minimize-tnc
   optimize.minimize-cobyla
   optimize.minimize-slsqp
   optimize.minimize-trustconstr
   optimize.minimize-dogleg
   optimize.minimize-trustncg
   optimize.minimize-trustkrylov
   optimize.minimize-trustexact

Constraints are passed to `minimize` function as a single object or
as a list of objects from the following classes:

.. autosummary::
   :toctree: generated/

   NonlinearConstraint - Class defining general nonlinear constraints.
   LinearConstraint - Class defining general linear constraints.

Simple bound constraints are handled separately and there is a special class
for them:

.. autosummary::
   :toctree: generated/

   Bounds - Bound constraints.

Quasi-Newton strategies implementing `HessianUpdateStrategy`
interface can be used to approximate the Hessian in `minimize`
function (available only for the 'trust-constr' method). Available
quasi-Newton methods implementing this interface are:

.. autosummary::
   :toctree: generated/

   BFGS - Broyden-Fletcher-Goldfarb-Shanno (BFGS) Hessian update strategy.
   SR1 - Symmetric-rank-1 Hessian update strategy.

The `minimize_scalar` function supports the following methods:

.. toctree::

   optimize.minimize_scalar-brent
   optimize.minimize_scalar-bounded
   optimize.minimize_scalar-golden

The specific optimization method interfaces below in this subsection are
not recommended for use in new scripts; all of these methods are accessible
via a newer, more consistent interface provided by the functions above.

General-purpose multivariate methods:

.. autosummary::
   :toctree: generated/

   fmin - Nelder-Mead Simplex algorithm
   fmin_powell - Powell's (modified) level set method
   fmin_cg - Non-linear (Polak-Ribiere) conjugate gradient algorithm
   fmin_bfgs - Quasi-Newton method (Broydon-Fletcher-Goldfarb-Shanno)
   fmin_ncg - Line-search Newton Conjugate Gradient

Constrained multivariate methods:

.. autosummary::
   :toctree: generated/

   fmin_l_bfgs_b - Zhu, Byrd, and Nocedal's constrained optimizer
   fmin_tnc - Truncated Newton code
   fmin_cobyla - Constrained optimization by linear approximation
   fmin_slsqp - Minimization using sequential least-squares programming
   differential_evolution - stochastic minimization using differential evolution

Univariate (scalar) minimization methods:

.. autosummary::
   :toctree: generated/

   fminbound - Bounded minimization of a scalar function
   brent - 1-D function minimization using Brent method
   golden - 1-D function minimization using Golden Section method

Equation (Local) Minimizers
---------------------------

.. autosummary::
   :toctree: generated/

   leastsq - Minimize the sum of squares of M equations in N unknowns
   least_squares - Feature-rich least-squares minimization.
   nnls - Linear least-squares problem with non-negativity constraint
   lsq_linear - Linear least-squares problem with bound constraints

Global Optimization
-------------------

.. autosummary::
   :toctree: generated/

   basinhopping - Basinhopping stochastic optimizer
   brute - Brute force searching optimizer
   differential_evolution - stochastic minimization using differential evolution

Rosenbrock function
-------------------

.. autosummary::
   :toctree: generated/

   rosen - The Rosenbrock function.
   rosen_der - The derivative of the Rosenbrock function.
   rosen_hess - The Hessian matrix of the Rosenbrock function.
   rosen_hess_prod - Product of the Rosenbrock Hessian with a vector.

Fitting
=======

.. autosummary::
   :toctree: generated/

   curve_fit -- Fit curve to a set of points

Root finding
============

Scalar functions
----------------
.. autosummary::
   :toctree: generated/

   brentq - quadratic interpolation Brent method
   brenth - Brent method, modified by Harris with hyperbolic extrapolation
   ridder - Ridder's method
   bisect - Bisection method
   newton - Secant method or Newton's method

Fixed point finding:

.. autosummary::
   :toctree: generated/

   fixed_point - Single-variable fixed-point solver

Multidimensional
----------------

General nonlinear solvers:

.. autosummary::
   :toctree: generated/

   root - Unified interface for nonlinear solvers of multivariate functions
   fsolve - Non-linear multi-variable equation solver
   broyden1 - Broyden's first method
   broyden2 - Broyden's second method

The `root` function supports the following methods:

.. toctree::

   optimize.root-hybr
   optimize.root-lm
   optimize.root-broyden1
   optimize.root-broyden2
   optimize.root-anderson
   optimize.root-linearmixing
   optimize.root-diagbroyden
   optimize.root-excitingmixing
   optimize.root-krylov
   optimize.root-dfsane

Large-scale nonlinear solvers:

.. autosummary::
   :toctree: generated/

   newton_krylov
   anderson

Simple iterations:

.. autosummary::
   :toctree: generated/

   excitingmixing
   linearmixing
   diagbroyden

:mod:`Additional information on the nonlinear solvers <scipy.optimize.nonlin>`

Linear Programming
==================

General linear programming solver:

.. autosummary::
   :toctree: generated/

   linprog -- Unified interface for minimizers of linear programming problems

The `linprog` function supports the following methods:

.. toctree::

   optimize.linprog-simplex
   optimize.linprog-interior-point

The simplex method supports callback functions, such as:

.. autosummary::
   :toctree: generated/

   linprog_verbose_callback -- Sample callback function for linprog (simplex)

Assignment problems:

.. autosummary::
   :toctree: generated/

   linear_sum_assignment -- Solves the linear-sum assignment problem

Utilities
=========

.. autosummary::
   :toctree: generated/

   approx_fprime - Approximate the gradient of a scalar function
   bracket - Bracket a minimum, given two starting points
   check_grad - Check the supplied derivative using finite differences
   line_search - Return a step that satisfies the strong Wolfe conditions

   show_options - Show specific options optimization solvers
   LbfgsInvHessProduct - Linear operator for L-BFGS approximate inverse Hessian
   HessianUpdateStrategy - Interface for implementing Hessian update strategies

"""

from __future__ import division, print_function, absolute_import

from .optimize import *
from ._minimize import *
from ._root import *
from .minpack import *
from .zeros import *
from .lbfgsb import fmin_l_bfgs_b, LbfgsInvHessProduct
from .tnc import fmin_tnc
from .cobyla import fmin_cobyla
from .nonlin import *
from .slsqp import fmin_slsqp
from .nnls import nnls
from ._basinhopping import basinhopping
from ._linprog import linprog, linprog_verbose_callback
from ._hungarian import linear_sum_assignment
from ._differentialevolution import differential_evolution
from ._lsq import least_squares, lsq_linear
from ._constraints import (NonlinearConstraint,
                           LinearConstraint,
                           Bounds)
from ._hessian_update_strategy import HessianUpdateStrategy, BFGS, SR1

__all__ = [s for s in dir() if not s.startswith('_')]

from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
