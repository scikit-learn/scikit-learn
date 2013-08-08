"""
Machine learning module for Python
==================================

sklearn is a Python module integrating classical machine
learning algorithms in the tightly-knit world of scientific Python
packages (numpy, scipy, matplotlib).

It aims to provide simple and efficient solutions to learning problems
that are accessible to everybody and reusable in various contexts:
machine-learning as a versatile tool for science and engineering.

See http://scikit-learn.org for complete documentation.
"""
import sys
__version__ = '0.14.1'

try:
    # This variable is injected in the __builtins__ by the build
    # process. It used to enable importing subpackages of sklearn when
    # the binaries are not built
    __SKLEARN_SETUP__
except NameError:
    __SKLEARN_SETUP__ = False

if __SKLEARN_SETUP__:
    sys.stderr.write('Partial import of sklearn during the build process.\n')
    # We are not importing the rest of the scikit during the build
    # process, as it may not be compiled yet
else:
    from . import __check_build
    from .base import clone

    def test(*args, **kwargs):
        import warnings
        # Not using a DeprecationWarning, as they are turned off by
        # default
        warnings.warn("""sklearn.test() is no longer supported to run the
scikit-learn test suite.

After installation, you can launch the test suite from outside the
source directory (you will need to have nosetests installed)::

   $ nosetests --exe sklearn

See the web page http://scikit-learn.org/stable/install.html#testing
for more information.

This function, `sklearn.test()` does not do anything. It does not run
the tests and will be removed in release 0.16.
""", stacklevel=2)

    # The following line is useful so that nosetests doesn't consider
    # "test" as a test function
    test.__test__ = False

    __all__ = ['cross_validation', 'cluster', 'covariance',
               'datasets', 'decomposition', 'feature_extraction',
               'feature_selection', 'semi_supervised',
               'gaussian_process', 'grid_search', 'hmm', 'lda', 'linear_model',
               'metrics', 'mixture', 'naive_bayes', 'neighbors', 'pipeline',
               'preprocessing', 'qda', 'svm', 'clone',
               'cross_decomposition',
               'isotonic', 'pls']


def setup_module(module):
    """Fixture for the tests to assure globally controllable seeding of RNGs
    """

    import os
    import numpy as np
    import random

    # It could have been provided in the environment
    _random_seed = os.environ.get('SKLEARN_SEED', None)
    if _random_seed is None:
        _random_seed = np.random.uniform() * (2 ** 31 - 1)
    _random_seed = int(_random_seed)
    print("I: Seeding RNGs with %r" % _random_seed)
    np.random.seed(_random_seed)
    random.seed(_random_seed)
