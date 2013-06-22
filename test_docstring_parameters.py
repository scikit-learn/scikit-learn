# TODO inspect for Cython (see sagenb.misc.sageinspect)
from __future__ import print_function

import os.path
import inspect
import imp
import sys
import sklearn
from pkgutil import walk_packages

docscrape_path = 'doc/sphinxext/numpy_ext/docscrape.py'
docscrape = imp.load_source('docscrape',
                            os.path.join(os.path.dirname(__file__),
                                         docscrape_path))


def get_name(func):
    parts = []
    module = inspect.getmodule(func)
    if module:
        parts.append(module.__name__)
    if hasattr(func, 'im_class'):
        parts.append(func.im_class.__name__)
    parts.append(func.__name__)
    return '.'.join(parts)


def check_parameters_match(func, doc=None):
    name = get_name(func)
    if not name.startswith('sklearn.'):
        return
    if inspect.isdatadescriptor(func):
        return
    try:
        args, varargs, varkw, defaults = inspect.getargspec(func)
    except TypeError:
        return
    if inspect.ismethod(func):
        # drop self
        args = args[1:]

    if doc is None:
        try:
            doc = docscrape.FunctionDoc(func)
        except ValueError:
            print('Failed to parse doc of %r' % func, file=sys.stderr)
            return
    # check set
    param_names = [name for name, _, _ in doc['Parameters']]
    # clean up some docscrape output:
    param_names = [name.split(':')[0].strip('` ') for name in param_names]
    param_names = [name for name in param_names if '*' not in name]

    try:
        args_set = set(args)
    except TypeError:
        # TODO: handle arg tuples
        return
    extra_params = set(param_names) - args_set
    if extra_params and not varkw:
        print(get_name(func), 'in doc', sorted(extra_params))

    if defaults:
        none_defaults = [arg for arg, default in zip(args[-len(defaults):],
                                                     defaults)
                         if default is None]
    else:
        none_defaults = []
    extra_args = args_set - set(param_names) - set(none_defaults)
    if param_names and extra_args:
        print(get_name(func), 'in argspec', sorted(extra_args))
    # check order?


def test_docstring_parameters():
    # unique classes and functions
    encountered = [] # classes we've already processed, avoids duplicate checks
    
    for importer, modname, ispkg in walk_packages(
            path=sklearn.__path__, prefix='sklearn.', onerror=lambda x: None):            
        module = __import__(modname, fromlist="dummy")               
        classes = inspect.getmembers(module, inspect.isclass)  
        for cname, cls in classes:
            if cls in encountered:
                continue
            else: encountered.append(cls)
            cdoc = docscrape.ClassDoc(cls)
            # TODO: deal with __new__
            if hasattr(cls, '__init__'):
                check_parameters_match(cls.__init__, cdoc)
            for method_name in cdoc.methods:                                                           
                method = getattr(cls, method_name)
                if method in encountered:
                    continue
                else: 
                    encountered.append(method)
                    check_parameters_match(method)
            if hasattr(cls, '__call__'):
                check_parameters_match(cls.__call__)

        functions = inspect.getmembers(module, inspect.isfunction)
        for fname, func in functions:
            if func in encountered:
                continue
            else: 
                encountered.append(func)
                check_parameters_match(func)              


if __name__ == '__main__':
    test_docstring_parameters()