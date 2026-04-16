"""
=================================================
Orthogonal distance regression (:mod:`scipy.odr`)
=================================================

.. currentmodule:: scipy.odr

.. deprecated:: 1.17.0
    `scipy.odr` is deprecated and will be removed in SciPy 1.19.0. Please use
    `pypi.org/project/odrpack/ <https://pypi.org/project/odrpack/>`_
    instead.

    The following example shows a brief comparison of the APIs::

        import numpy as np
        import scipy.odr
        import odrpack

        # Classic "Pearson data" that motivates ODR.
        # Errors are in both variables, and if you don't account for this,
        # doing a linear fit of X vs. Y or Y vs. X will give you quite
        # different results.
        p_x = np.array([0., .9, 1.8, 2.6, 3.3, 4.4, 5.2, 6.1, 6.5, 7.4])
        p_y = np.array([5.9, 5.4, 4.4, 4.6, 3.5, 3.7, 2.8, 2.8, 2.4, 1.5])
        p_sx = np.array([.03, .03, .04, .035, .07, .11, .13, .22, .74, 1.])
        p_sy = np.array([1., .74, .5, .35, .22, .22, .12, .12, .1, .04])

        # Old-style
        # The RealData class takes care of details like turning
        # standard-deviation error bars into weights.
        p_dat = scipy.odr.RealData(p_x, p_y, sx=p_sx, sy=p_sy)
        # Note, parameters come before `x` in scipy.odr
        p_mod = scipy.odr.Model(lambda beta, x: beta[0] + beta[1]*x)
        p_odr = scipy.odr.ODR(p_dat, p_mod, beta0=[1., 1.])
        old_out = p_odr.run()

        # New-style
        # Parameters come after data, in the new API.
        # We must convert the error bars into weights ourselves.
        new_out = odrpack.odr_fit(lambda x, beta: beta[0] + beta[1] * x,
            p_x, p_y, beta0=np.array([1.0, 1.0]),
            weight_x=p_sx**-2, weight_y=p_sy**-2)

        assert np.isclose(old_out.beta, new_out.beta).all()


Package Content
===============

.. autosummary::
   :toctree: generated/

   Data          -- The data to fit.
   RealData      -- Data with weights as actual std. dev.s and/or covariances.
   Model         -- Stores information about the function to be fit.
   ODR           -- Gathers all info & manages the main fitting routine.
   Output        -- Result from the fit.
   odr           -- Low-level function for ODR.

   OdrWarning    -- Warning about potential problems when running ODR.
   OdrError      -- Error exception.
   OdrStop       -- Stop exception.

   polynomial    -- Factory function for a general polynomial model.
   exponential   -- Exponential model
   multilinear   -- Arbitrary-dimensional linear model
   unilinear     -- Univariate linear model
   quadratic     -- Quadratic model

Usage information
=================

Introduction
------------

Why Orthogonal Distance Regression (ODR)? Sometimes one has
measurement errors in the explanatory (a.k.a., "independent")
variable(s), not just the response (a.k.a., "dependent") variable(s).
Ordinary Least Squares (OLS) fitting procedures treat the data for
explanatory variables as fixed, i.e., not subject to error of any kind.
Furthermore, OLS procedures require that the response variables be an
explicit function of the explanatory variables; sometimes making the
equation explicit is impractical and/or introduces errors.  ODR can
handle both of these cases with ease, and can even reduce to the OLS
case if that is sufficient for the problem.

ODRPACK is a FORTRAN-77 library for performing ODR with possibly
non-linear fitting functions. It uses a modified trust-region
Levenberg-Marquardt-type algorithm [1]_ to estimate the function
parameters.  The fitting functions are provided by Python functions
operating on NumPy arrays. The required derivatives may be provided
by Python functions as well, or may be estimated numerically. ODRPACK
can do explicit or implicit ODR fits, or it can do OLS. Input and
output variables may be multidimensional. Weights can be provided to
account for different variances of the observations, and even
covariances between dimensions of the variables.

The `scipy.odr` package offers an object-oriented interface to
ODRPACK, in addition to the low-level `odr` function.

Additional background information about ODRPACK can be found in the
`ODRPACK User's Guide
<https://docs.scipy.org/doc/external/odrpack_guide.pdf>`_, reading
which is recommended.

Basic usage
-----------

1. Define the function you want to fit against.::

       def f(B, x):
           '''Linear function y = m*x + b'''
           # B is a vector of the parameters.
           # x is an array of the current x values.
           # x is in the same format as the x passed to Data or RealData.
           #
           # Return an array in the same format as y passed to Data or RealData.
           return B[0]*x + B[1]

2. Create a Model.::

       linear = Model(f)

3. Create a Data or RealData instance.::

       mydata = Data(x, y, wd=1./power(sx,2), we=1./power(sy,2))

   or, when the actual covariances are known::

       mydata = RealData(x, y, sx=sx, sy=sy)

4. Instantiate ODR with your data, model and initial parameter estimate.::

       myodr = ODR(mydata, linear, beta0=[1., 2.])

5. Run the fit.::

       myoutput = myodr.run()

6. Examine output.::

       myoutput.pprint()


References
----------
.. [1] P. T. Boggs and J. E. Rogers, "Orthogonal Distance Regression,"
   in "Statistical analysis of measurement error models and
   applications: proceedings of the AMS-IMS-SIAM joint summer research
   conference held June 10-16, 1989," Contemporary Mathematics,
   vol. 112, pg. 186, 1990.

"""
# version: 0.7
# author: Robert Kern <robert.kern@gmail.com>
# date: 2006-09-21

from ._odrpack import *
from ._models import *
from . import _add_newdocs

# Deprecated namespaces, to be removed in v2.0.0
from . import models, odrpack

__all__ = [s for s in dir()
           if not (s.startswith('_') or s in ('odr_stop', 'odr_error'))]

import warnings
msg = ("`scipy.odr` is deprecated as of version 1.17.0 and will be removed in "
        "SciPy 1.19.0. Please use `https://pypi.org/project/odrpack/` instead.")
warnings.warn(msg, DeprecationWarning, stacklevel=2)
del warnings


from scipy._lib._testutils import PytestTester
test = PytestTester(__name__)
del PytestTester
