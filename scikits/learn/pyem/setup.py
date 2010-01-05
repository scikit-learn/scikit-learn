#! /usr/bin/env python
# Last Change: Thu Jul 13 01:00 PM 2006 J
""" pyem is a small python package to estimate Gaussian Mixtures Models
from data, using Expectation Maximization"""
from distutils.core import setup, Extension
from pyem import version as pyem_version

# distutils does not update MANIFEST correctly, removes it
import os
if os.path.exists('MANIFEST'): os.remove('MANIFEST')

from numpy.distutils.misc_util import get_numpy_include_dirs
NUMPYINC    = get_numpy_include_dirs()[0]
print NUMPYINC

# General variables:
#   - DISTNAME: name of the distributed package
#   - VERSION: the version reference is in pyem/__init__.py file
#   - other upper cased variables are the same than the corresponding 
#   keywords in setup call
DISTNAME    = 'pyem' 
VERSION     = pyem_version
DESCRIPTION ='A python module for Expectation Maximization learning of mixtures pdf',
AUTHOR      ='David Cournapeau',
AUTHOR_EMAIL='david@ar.media.kyoto-u.ac.jp',
URL         ='http://ar.media.kyoto-u.ac.jp/members/david',

def setup_pyrex():
    setup(name=DISTNAME,
        version=VERSION,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        packages=['pyem'], 
        ext_modules=[Extension('pyem/c_gmm', ['pyem/src/c_gmm.pyx'], 
            include_dirs=[NUMPYINC])],
        data_files=['c_numpy.pxd', 'c_python.pxd'],
        cmdclass = {'build_ext': build_ext},
    )

def setup_nopyrex():
    setup(name=DISTNAME,
        version=VERSION,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        packages=['pyem'], 
        #ext_modules=[Extension('_hello', ['hellomodule.c'])],
    )

try:
    from Pyrex.Distutils import build_ext
    setup_pyrex()
except:
    print "Pyrex not found, C extension won't be available"
    setup_nopyrex()

