from __future__ import print_function

from nose.tools import assert_true
import inspect
import warnings
import importlib

from pkgutil import walk_packages
from inspect import getsource

import sklearn
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import assert_raise_message

_doc_special_members = ('__contains__', '__getitem__', '__iter__', '__len__',
                        '__call__', '__add__', '__sub__', '__mul__', '__div__',
                        '__neg__', '__hash__')

public_modules = [
    # the list of modules users need to access for all functionality
    # 'sklearn',
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


# helpers to get function arguments
if hasattr(inspect, 'signature'):  # py35
    def _get_args(function, varargs=False):
        params = inspect.signature(function).parameters
        args = [key for key, param in params.items()
                if param.kind not in (param.VAR_POSITIONAL, param.VAR_KEYWORD)]
        if varargs:
            varargs = [param.name for param in params.values()
                       if param.kind == param.VAR_POSITIONAL]
            if len(varargs) == 0:
                varargs = None
            return args, varargs
        else:
            return args
else:
    def _get_args(function, varargs=False):
        out = inspect.getargspec(function)  # args, varargs, keywords, defaults
        if varargs:
            return out[:2]
        else:
            return out[0]


def get_name(func):
    parts = []
    module = inspect.getmodule(func)
    if module:
        parts.append(module.__name__)
    if hasattr(func, 'im_class'):
        parts.append(func.im_class.__name__)
    parts.append(func.__name__)
    return '.'.join(parts)


# functions to ignore args / docstring of
_docstring_ignores = [
    'sklearn.utils.deprecation.load_mlcomp',
]

_tab_ignores = [
]


def check_parameters_match(func, doc=None):
    """Helper to check docstring, returns list of incorrect results"""
    from numpydoc import docscrape
    incorrect = []
    name_ = get_name(func)
    if (not name_.startswith('sklearn.') or
            name_.startswith('sklearn.externals')):
        return incorrect
    if inspect.isdatadescriptor(func):
        return incorrect
    args = _get_args(func)
    # drop self
    if len(args) > 0 and args[0] == 'self':
        args = args[1:]

    if doc is None:
        with warnings.catch_warnings(record=True) as w:
            try:
                doc = docscrape.FunctionDoc(func)
            except Exception as exp:
                incorrect += [name_ + ' parsing error: ' + str(exp)]
                return incorrect
        if len(w):
            raise RuntimeError('Error for %s:\n%s' % (name_, w[0]))
    # check set
    param_names = [name for name, _, _ in doc['Parameters']]
    # clean up some docscrape output:
    param_names = [name.split(':')[0].strip('` ') for name in param_names]
    param_names = [name for name in param_names if '*' not in name]
    if len(param_names) != len(args):
        bad = str(sorted(list(set(param_names) - set(args)) +
                         list(set(args) - set(param_names))))
        if not any(d in name_ for d in _docstring_ignores) and \
                'deprecation_wrapped' not in func.__code__.co_name:
            incorrect += [name_ + ' arg mismatch: ' + bad]
    else:
        for n1, n2 in zip(param_names, args):
            if n1 != n2:
                incorrect += [name_ + ' ' + n1 + ' != ' + n2]
    return incorrect


def f_ok(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : float
        Parameter b

    Returns
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_bad_sections(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : float
        Parameter b

    Results
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_bad_order(b, a):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a
    b : float
        Parameter b

    Returns
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


def f_missing(a, b):
    """Function f

    Parameters
    ----------
    a : int
        Parameter a

    Returns
    -------
    c : list
        Parameter c
    """
    c = a + b
    return c


class Klass(object):
    def f_missing(self, X, y):
        pass

    def f_bad_sections(self, X, y):
        """Function f

        Parameter
        ----------
        a : int
            Parameter a
        b : float
            Parameter b

        Results
        -------
        c : list
            Parameter c
        """
        pass


def test_check_parameters_match():
    check_parameters_match(f_ok)
    assert_raise_message(RuntimeError, 'Unknown section Results',
                         check_parameters_match, f_bad_sections)
    assert_raise_message(RuntimeError, 'Unknown section Parameter',
                         check_parameters_match, Klass.f_bad_sections)

    messages = ['a != b']
    messages += ["arg mismatch: ['b']"]
    messages += ["arg mismatch: ['X', 'y']"]
    for mess, f in zip(messages, [f_bad_order, f_missing, Klass.f_missing]):
        incorrect = check_parameters_match(f)
        assert_true(len(incorrect) >= 1)
        assert_true(mess in incorrect[0],
                    '"%s" not in "%s"' % (mess, incorrect[0]))


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
                incorrect += check_parameters_match(method)
            if hasattr(cls, '__call__'):
                incorrect += check_parameters_match(cls.__call__)
        functions = inspect.getmembers(module, inspect.isfunction)
        for fname, func in functions:
            if fname.startswith('_'):
                continue
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
