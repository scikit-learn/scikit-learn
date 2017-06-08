from __future__ import print_function

from nose.tools import assert_true
import inspect
import warnings
import importlib

from pkgutil import walk_packages
from inspect import getsource

import sklearn
from sklearn.base import signature
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import check_parameters_match, get_func_name

_doc_special_members = ('__contains__', '__getitem__', '__iter__', '__len__',
                        '__call__', '__add__', '__sub__', '__mul__', '__div__',
                        '__neg__', '__hash__')

public_modules = [
    # the list of modules users need to access for all functionality
    # 'sklearn',
    'sklearn.base',
    # 'sklearn.cluster',
    # 'sklearn.covariance',
    # 'sklearn.cross_decomposition',
    'sklearn.datasets',
    # 'sklearn.decomposition',
    # 'sklearn.ensemble',
    # 'sklearn.feature_extraction',
    # 'sklearn.feature_selection',
    # 'sklearn.gaussian_process',
    # 'sklearn.linear_model',
    # 'sklearn.manifold',
    # 'sklearn.metrics',
    # 'sklearn.mixture',
    # 'sklearn.model_selection',
    # 'sklearn.neighbors',
    # 'sklearn.neural_network',
    # 'sklearn.preprocessing',
    # 'sklearn.semi_supervised',
    # 'sklearn.tree',
    # 'sklearn.utils',
]


# functions to ignore args / docstring of
_docstring_ignores = [
    'sklearn.utils.deprecation.load_mlcomp',
]

_tab_ignores = [
]


def test_docstring_parameters():
    """Test module docstring formatting."""
    try:
        import numpydoc  # noqa
    except ImportError:
        raise SkipTest(
            "numpydoc is required to test the docstrings")

    from numpydoc import docscrape

    incorrect = []
    for name in public_modules:
        with warnings.catch_warnings(record=True):  # traits warnings
            module = __import__(name, globals())
        for submod in name.split('.')[1:]:
            module = getattr(module, submod)
        classes = inspect.getmembers(module, inspect.isclass)
        for cname, cls in classes:
            if cname.startswith('_') and cname not in _doc_special_members:
                continue
            with warnings.catch_warnings(record=True) as w:
                cdoc = docscrape.ClassDoc(cls)
            if len(w):
                raise RuntimeError('Error for __init__ of %s in %s:\n%s'
                                   % (cls, name, w[0]))
            if hasattr(cls, '__init__'):
                incorrect += check_parameters_match(cls.__init__, cdoc)
            for method_name in cdoc.methods:
                method = getattr(cls, method_name)
                param_ignore = None
                # Now skip docstring test for y when y is None
                # by default for API reason
                if method_name in ['fit', 'score', 'fit_predict']:
                    sig = signature(method)
                    if ('y' in sig.parameters and
                            sig.parameters['y'].default is None):
                        param_ignore = ['y']  # ignore y for fit and score
                incorrect += check_parameters_match(method,
                                                    ignore=param_ignore)
            if hasattr(cls, '__call__'):
                incorrect += check_parameters_match(cls.__call__)
        functions = inspect.getmembers(module, inspect.isfunction)
        for fname, func in functions:
            if fname.startswith('_'):
                # Don't test private methods / functions
                continue
            name_ = get_func_name(func)
            if not any(d in name_ for d in _docstring_ignores) and \
                    'deprecation_wrapped' not in func.__code__.co_name:
                incorrect += check_parameters_match(func)
    msg = '\n' + '\n'.join(sorted(list(set(incorrect))))
    if len(incorrect) > 0:
        raise AssertionError(msg)


def test_tabs():
    """Test that there are no tabs in our source files"""
    ignore = _tab_ignores[:]

    for importer, modname, ispkg in walk_packages(sklearn.__path__,
                                                  prefix='sklearn.'):
        # because we don't import e.g. mne.tests w/mne
        if not ispkg and modname not in ignore:
            mod = importlib.import_module(modname)
            try:
                source = getsource(mod)
            except IOError:  # user probably should have run "make clean"
                continue
            assert_true('\t' not in source,
                        '"%s" has tabs, please remove them or add it to the'
                        'ignore list' % modname)
