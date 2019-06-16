# Authors: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Raghav RV <rvraghav93@gmail.com>
# License: BSD 3 clause

import inspect
import warnings
import importlib

from pkgutil import walk_packages
from inspect import getsource, isabstract, signature

import sklearn
from sklearn.utils import IS_PYPY
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import check_docstring_parameters
from sklearn.utils.testing import _get_func_name
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.deprecation import _is_deprecated

import pytest

PUBLIC_MODULES = set([pckg[1] for pckg in walk_packages(prefix='sklearn.',
                                                        path=sklearn.__path__)
                      if not ("._" in pckg[1] or ".tests." in pckg[1])])


# functions to ignore args / docstring of
_DOCSTRING_IGNORES = [
    'sklearn.utils.deprecation.load_mlcomp',
    'sklearn.pipeline.make_pipeline',
    'sklearn.pipeline.make_union',
    'sklearn.utils.extmath.safe_sparse_dot',
    'sklearn.utils._joblib'
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


# numpydoc 0.8.0's docscrape tool raises because of collections.abc under
# Python 3.7
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.skipif(IS_PYPY, reason='test segfaults on PyPy')
def test_docstring_parameters():
    # Test module docstring formatting

    # Skip test if numpydoc is not found
    try:
        import numpydoc  # noqa
    except ImportError:
        raise SkipTest("numpydoc is required to test the docstrings")

    from numpydoc import docscrape

    incorrect = []
    for name in PUBLIC_MODULES:
        if name == 'sklearn.utils.fixes':
            # We cannot always control these docstrings
            continue
        with warnings.catch_warnings(record=True):
            module = importlib.import_module(name)
        classes = inspect.getmembers(module, inspect.isclass)
        # Exclude imported classes
        classes = [cls for cls in classes if cls[1].__module__ == name]
        for cname, cls in classes:
            this_incorrect = []
            if cname in _DOCSTRING_IGNORES or cname.startswith('_'):
                continue
            if isabstract(cls):
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
        # Exclude imported functions
        functions = [fn for fn in functions if fn[1].__module__ == name]
        for fname, func in functions:
            # Don't test private methods / functions
            if fname.startswith('_'):
                continue
            if fname == "configuration" and name.endswith("setup"):
                continue
            name_ = _get_func_name(func)
            if (not any(d in name_ for d in _DOCSTRING_IGNORES) and
                    not _is_deprecated(func)):
                incorrect += check_docstring_parameters(func)
    msg = '\n' + '\n'.join(sorted(list(set(incorrect))))
    if len(incorrect) > 0:
        raise AssertionError("Docstring Error: " + msg)


@ignore_warnings(category=DeprecationWarning)
def test_tabs():
    # Test that there are no tabs in our source files
    for importer, modname, ispkg in walk_packages(sklearn.__path__,
                                                  prefix='sklearn.'):

        if IS_PYPY and ('_svmlight_format' in modname or
                        'feature_extraction._hashing' in modname):
            continue

        # because we don't import
        mod = importlib.import_module(modname)
        try:
            source = getsource(mod)
        except IOError:  # user probably should have run "make clean"
            continue
        assert '\t' not in source, ('"%s" has tabs, please remove them ',
                                    'or add it to theignore list'
                                    % modname)
