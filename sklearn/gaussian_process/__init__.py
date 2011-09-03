#!/usr/bin/python
# -*- coding: utf-8 -*-

# Author: Vincent Dubourg <vincent.dubourg@gmail.com>
#         (mostly translation, see implementation details)
# License: BSD style

"""
A module that implements scalar Gaussian Process based prediction (also
known as Kriging).

Contains
--------
GaussianProcess: The main class of the module that implements the Gaussian
                 Process prediction theory.
regression_models: A submodule that contains the built-in regression models.
correlation_models: A submodule that contains the built-in correlation models.

Implementation details
----------------------
The presentation implementation is based on a translation of the DACE
Matlab toolbox, see reference [1].

References
----------
[1] H.B. Nielsen, S.N. Lophaven, H. B. Nielsen and J. Sondergaard (2002).
    DACE - A MATLAB Kriging Toolbox.
    http://www2.imm.dtu.dk/~hbn/dace/dace.pdf

[2] W.J. Welch, R.J. Buck, J. Sacks, H.P. Wynn, T.J. Mitchell, and M.D.
    Morris (1992). Screening, predicting, and computer experiments.
    Technometrics, 34(1) 15--25.
"""

from .gaussian_process import GaussianProcess
from . import correlation_models
from . import regression_models
