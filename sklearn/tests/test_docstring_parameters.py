# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Raghav RV <rvraghav93@gmail.com>
# License: BSD 3 clause

from __future__ import print_function

import inspect
import sys
import warnings
import importlib

from pkgutil import walk_packages
from inspect import getsource

import sklearn
from sklearn.base import signature
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import check_docstring_parameters
from sklearn.utils.testing import _get_func_name
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.deprecation import _is_deprecated

PUBLIC_MODULES = set(['sklearn.' + modname
                      for _, modname, _ in walk_packages(sklearn.__path__)
                      if not modname.startswith('_') and
                      '.tests.' not in modname])

# TODO Uncomment all modules and fix doc inconsistencies everywhere
# The list of modules that are not tested for now
PUBLIC_MODULES -= set([
    'sklearn.ensemble',
    'sklearn.feature_selection',
    'sklearn.kernel_approximation',
    'sklearn.model_selection',
    'sklearn.multioutput',
    'sklearn.random_projection',
    'sklearn.setup',
    'sklearn.svm',
    'sklearn.utils',
    # Deprecated modules
    'sklearn.cross_validation',
    'sklearn.grid_search',
    'sklearn.learning_curve',
])

# functions to ignore args / docstring of
_DOCSTRING_IGNORES = [
    'sklearn.utils.deprecation.load_mlcomp',
    'sklearn.pipeline.make_pipeline',
    'sklearn.pipeline.make_union',
    'sklearn.utils.extmath.safe_sparse_dot',
]

# Methods where y param should be ignored if y=None by default
_METHODS_IGNORE_NONE_Y = [
    'fit',
    'score',
    'fit_predict',
    'fit_transform',
    'partial_fit',
    'predict'
]


def test_docstring_parameters():
    raise SkipTest('Not testing docstring parameter consistency')

    # Test module docstring formatting

    # Skip test if numpydoc is not found or if python version is < 3.5
    try:
        import numpydoc  # noqa
        assert sys.version_info >= (3, 5)
    except (ImportError, AssertionError):
        raise SkipTest("numpydoc is required to test the docstrings, "
                       "as well as python version >= 3.5")

    from numpydoc import docscrape

    incorrect = []
    for name in PUBLIC_MODULES:
        with warnings.catch_warnings(record=True):
            module = importlib.import_module(name)
        classes = inspect.getmembers(module, inspect.isclass)
        for cname, cls in classes:
            this_incorrect = []
            if cname in _DOCSTRING_IGNORES:
                continue
            if cname.startswith('_'):
                continue
            with warnings.catch_warnings(record=True) as w:
                cdoc = docscrape.ClassDoc(cls)
            if len(w):
                raise RuntimeError('Error for __init__ of %s in %s:\n%s'
                                   % (cls, name, w[0]))

            cls_init = getattr(cls, '__init__', None)

            if _is_deprecated(cls_init):
                continue

            elif cls_init is not None:
                this_incorrect += check_docstring_parameters(
                    cls.__init__, cdoc, class_name=cname)
            for method_name in cdoc.methods:
                method = getattr(cls, method_name)
                if _is_deprecated(method):
                    continue
                param_ignore = None
                # Now skip docstring test for y when y is None
                # by default for API reason
                if method_name in _METHODS_IGNORE_NONE_Y:
                    sig = signature(method)
                    if ('y' in sig.parameters and
                            sig.parameters['y'].default is None):
                        param_ignore = ['y']  # ignore y for fit and score
                result = check_docstring_parameters(
                    method, ignore=param_ignore, class_name=cname)
                this_incorrect += result

            incorrect += this_incorrect

        functions = inspect.getmembers(module, inspect.isfunction)
        for fname, func in functions:
            # Don't test private methods / functions
            if fname.startswith('_'):
                continue
            name_ = _get_func_name(func)
            if (not any(d in name_ for d in _DOCSTRING_IGNORES) and
                    not _is_deprecated(func)):
                incorrect += check_docstring_parameters(func)
    msg = '\n' + '\n'.join(sorted(list(set(incorrect))))
    if len(incorrect) > 0:
        raise AssertionError(msg)


@ignore_warnings(category=DeprecationWarning)
def test_tabs():
    # Test that there are no tabs in our source files
    for importer, modname, ispkg in walk_packages(sklearn.__path__,
                                                  prefix='sklearn.'):
        # because we don't import
        mod = importlib.import_module(modname)
        try:
            source = getsource(mod)
        except IOError:  # user probably should have run "make clean"
            continue
        assert '\t' not in source, ('"%s" has tabs, please remove them ',
                                    'or add it to theignore list'
                                    % modname)
