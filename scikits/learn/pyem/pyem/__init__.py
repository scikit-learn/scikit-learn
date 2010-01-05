#! /usr/bin/env python
# Last Change: Thu Sep 28 01:00 PM 2006 J

from info import __doc__

from numpy.testing import NumpyTest
test = NumpyTest().test

version = '0.5.2'

from gauss_mix import GmParamError, GM
from gmm_em import GmmParamError, GMM, EM
