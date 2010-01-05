#! /usr/bin/env python
# Last Change: Tue Oct 03 06:00 PM 2006 J

from info import __doc__

from numpy.testing import NumpyTest
test = NumpyTest().test

version = '0.5.3'

from gauss_mix import GmParamError, GM
from gmm_em import GmmParamError, GMM, EM
