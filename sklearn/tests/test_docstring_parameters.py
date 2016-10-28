from pkgutil import walk_packages
from importlib import import_module
from numpydoc import docscrape

import inspect
import warnings

import sklearn


def get_name(func):
    parts = []
    module = inspect.getmodule(func)
    if module:
        parts.append(module.__name__)
    if hasattr(func, 'im_class'):
        parts.append(func.im_class.__name__)
    parts.append(func.__name__)
    return '.'.join(parts)


def get_all_modules():
    modules = []
    for importer, modname, ispkg in walk_packages(sklearn.__path__, prefix='sklearn.'):
        if ispkg:
            modules.append(modname)
    return modules


def check_parameters_match(func, doc=None):
    incorrect = []
    name_ = get_name(func)

    # skip deprecated and data descriptors
    if not name_.startswith('sklearn.') or 'deprecation' in name_ or inspect.isdatadescriptor(func):
        return incorrect

    args = inspect.getargspec(func)[0]

    # drop self
    if len(args) > 0 and args[0] == 'self':
        args = args[1:]

    if doc is None:
        with warnings.catch_warnings(record=True) as w:
            doc = docscrape.FunctionDoc(func)
        if len(w):
            raise RuntimeError('Error for %s:%s' % (name, w[0]))

    param_names = [name for name, _, _ in doc['Parameters']]
    param_names = [name.split(':')[0].strip('` ') for name in param_names]
    param_names = [name for name in param_names if '*' not in name]

    if len(param_names) != len(args):
        mismatch = str(sorted(list(set(args) ^ set(param_names))))
        incorrect += ['For function "' + name_ + '" arg mismatch: ' + mismatch]
    else:
        args.sort()
        param_names.sort()
        for n1, n2 in zip(param_names, args):
            if n1 != n2:
                incorrect += [name_ + ' ' + n1 + ' != ' + n2]

    return incorrect


def test_docstring_parameters():
    """Test module docstring formatting"""
    public_modules = get_all_modules()
    incorrect = []

    for name in public_modules:
        if name.endswith('tests'):
            continue

        module = import_module(name, globals())

        # check for classes in the module
        classes = inspect.getmembers(module, inspect.isclass)
        for cname, cls in classes:
            if cname.startswith('_'):
                continue

            with warnings.catch_warnings(record=True) as w:
                cdoc = docscrape.ClassDoc(cls)
            if len(w):
                raise RuntimeError('Error for __init__ of %s in %s:\n%s'
                                    % (cls, name, w[0]))

            # check for __init__ method of class
            if hasattr(cls, '__init__'):
                incorrect += check_parameters_match(cls.__init__, cdoc)

            # check for all methods of class
            for method_name in cdoc.methods:
                method = getattr(cls, method_name)
                incorrect += check_parameters_match(method)

            # check for _call__ method of class
            if hasattr(cls, '__call__'):
                incorrect += check_parameters_match(cls.__call__)

        # check for functions in module
        functions = inspect.getmembers(module, inspect.isfunction)
        for fname, func in functions:
            if fname.startswith('_'):
                continue
            incorrect += check_parameters_match(func)

    msg = '\n' + '\n'.join(sorted(list(set(incorrect))))
    if len(incorrect) > 0:
        raise AssertionError(msg)
