#! /usr/bin/env python
# Last Change: Sat Jun 09 10:00 PM 2007 J

from info import __doc__

from gauss_mix import GmParamError, GM
from gmm_em import GmmParamError, GMM, EM
#from online_em import OnGMM as _OnGMM
#import examples as _examples

__all__ = filter(lambda s:not s.startswith('_'), dir())

from numpy.testing import NumpyTest
test = NumpyTest().test

