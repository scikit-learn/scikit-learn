#! /usr/bin/env python
# Last Change: Thu Oct 12 11:00 PM 2006 J

from info import __doc__

from numpy.testing import NumpyTest
test = NumpyTest().test

from gauss_mix import GmParamError, GM
from gmm_em import GmmParamError, GMM, EM
