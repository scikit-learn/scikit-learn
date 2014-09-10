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
import re
import warnings
__version__ = '0.16-git'

# Make sure that DeprecationWarning within this package always gets printed
warnings.filterwarnings('always', category=DeprecationWarning,
                        module='^{0}\.'.format(re.escape(__name__)))

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
    __check_build  # avoid flakes unused variable error

    __all__ = ['cluster', 'covariance', 'cross_decomposition',
               'cross_validation', 'datasets', 'decomposition', 'dummy',
               'ensemble', 'externals', 'feature_extraction',
               'feature_selection', 'gaussian_process', 'grid_search', 'hmm',
               'isotonic', 'kernel_approximation', 'lda', 'learning_curve',
               'linear_model', 'manifold', 'metrics', 'mixture', 'multiclass',
               'naive_bayes', 'neighbors', 'neural_network', 'pipeline',
               'preprocessing', 'qda', 'random_projection', 'semi_supervised',
               'svm', 'tree',
               # Non-modules:
               'clone']


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
