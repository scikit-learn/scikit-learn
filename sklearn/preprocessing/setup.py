#-----------------------------------------------------------------------------
# Copyright (c) 2012, Enthought, Inc.
# All rights reserved.  See LICENSE.txt for details.
#
# Author: Kurt W. Smith
# Date: 26 March 2012
#-----------------------------------------------------------------------------

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
exts = [Extension("_mdlp", ["_mdlp.pyx"], include_dirs=[numpy.get_include()])]

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = exts,
)

#setup(
#    ext_modules=cythonize("_mdlp.pyx"),
#    include_dirs=[numpy.get_include()]
#)
