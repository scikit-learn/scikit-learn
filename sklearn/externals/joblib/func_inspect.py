"""
My own variation on function-specific inspect-like features.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import itertools
import inspect
import warnings
import os


def get_func_code(func):
    """ Attempts to retrieve a reliable function code hash.

        The reason we don't use inspect.getsource is that it caches the
        source, whereas we want this to be modified on the fly when the
        function is modified.

        Returns
        -------
        func_code: string
            The function code
        source_file: string
            The path to the file in which the function is defined.
        first_line: int
            The first line of the code in the source file.

        Notes
        ------
        This function does a bit more magic than inspect, and is thus
        more robust.
    """
    source_file = None
    try:
        # Try to retrieve the source code.
        source_file = func.func_code.co_filename
        source_file_obj = file(source_file)
        first_line = func.func_code.co_firstlineno
        # All the lines after the function definition:
        source_lines = list(itertools.islice(source_file_obj, first_line - 1,
                                             None))
        return ''.join(inspect.getblock(source_lines)), source_file, first_line
    except:
        # If the source code fails, we use the hash. This is fragile and
        # might change from one session to another.
        if hasattr(func, 'func_code'):
            return str(func.func_code.__hash__()), source_file, -1
        else:
            # Weird objects like numpy ufunc don't have func_code
            # This is fragile, as quite often the id of the object is
            # in the repr, so it might not persist accross sessions,
            # however it will work for ufuncs.
            return repr(func), source_file, -1


def get_func_name(func, resolv_alias=True, win_characters=True):
    """ Return the function import path (as a list of module names), and
        a name for the function.

        Parameters
        ----------
        func: callable
            The func to inspect
        resolv_alias: boolean, optional
            If true, possible local alias are indicated.
        win_characters: boolean, optional
            If true, substitute special characters using urllib.quote
            This is useful in Windows, as it cannot encode some filenames
    """
    if hasattr(func, '__module__'):
        module = func.__module__
    else:
        try:
            module = inspect.getmodule(func)
        except TypeError:
            if hasattr(func, '__class__'):
                module = func.__class__.__module__
            else:
                module = 'unkown'
    if module is None:
        # Happens in doctests, eg
        module = ''
    if module == '__main__':
        try:
            filename = inspect.getsourcefile(func)
        except:
            filename = None
        if filename is not None:
            # mangling of full path to filename
            filename = filename.replace(os.sep, '-')
            filename = filename.replace(":", "-")
            if filename.endswith('.py'):
                filename = filename[:-3]
            module = module + '-' + filename
    module = module.split('.')
    if hasattr(func, 'func_name'):
        name = func.func_name
    elif hasattr(func, '__name__'):
        name = func.__name__
    else:
        name = 'unknown'
    # Hack to detect functions not defined at the module-level
    if resolv_alias:
        # TODO: Maybe add a warning here?
        if hasattr(func, 'func_globals') and name in func.func_globals:
            if not func.func_globals[name] is func:
                name = '%s-alias' % name
    if inspect.ismethod(func):
        # We need to add the name of the class
        if hasattr(func, 'im_class'):
            klass = func.im_class
            module.append(klass.__name__)
    if os.name == 'nt' and win_characters:
        # Stupid windows can't encode certain characters in filenames
        import urllib
        for char in ('<', '>', '!', ':'):
            name = name.replace(char, urllib.quote(char))
    return module, name


def filter_args(func, ignore_lst, *args, **kwargs):
    """ Filters the given args and kwargs using a list of arguments to
        ignore, and a function specification.

        Parameters
        ----------
        func: callable
            Function giving the argument specification
        ignore_lst: list of strings
            List of arguments to ignore (either a name of an argument
            in the function spec, or '*', or '**')
        *args: list
            Positional arguments passed to the function.
        **kwargs: dict
            Keyword arguments passed to the function

        Returns
        -------
        filtered_args: list
            List of filtered positional arguments.
        filtered_kwdargs: dict
            List of filtered Keyword arguments.
    """
    args = list(args)
    if isinstance(ignore_lst, basestring):
        # Catch a common mistake
        raise ValueError('ignore_lst must be a list of parameters to ignore '
            '%s (type %s) was given' % (ignore_lst, type(ignore_lst)))
    # Special case for functools.partial objects
    if (not inspect.ismethod(func) and not inspect.isfunction(func)):
        if ignore_lst:
            warnings.warn('Cannot inspect object %s, ignore list will '
                'not work.' % func, stacklevel=2)
        return {'*': args, '**': kwargs}
    arg_spec = inspect.getargspec(func)
    # We need to if/them to account for different versions of Python
    if hasattr(arg_spec, 'args'):
        arg_names = arg_spec.args
        arg_defaults = arg_spec.defaults
        arg_keywords = arg_spec.keywords
        arg_varargs = arg_spec.varargs
    else:
        arg_names, arg_varargs, arg_keywords, arg_defaults = arg_spec
    arg_defaults = arg_defaults or {}
    if inspect.ismethod(func):
        # First argument is 'self', it has been removed by Python
        # we need to add it back:
        args = [func.im_self, ] + args
    # XXX: Maybe I need an inspect.isbuiltin to detect C-level methods, such
    # as on ndarrays.

    _, name = get_func_name(func, resolv_alias=False)
    arg_dict = dict()
    arg_position = -1
    for arg_position, arg_name in enumerate(arg_names):
        if arg_position < len(args):
            # Positional argument or keyword argument given as positional
            arg_dict[arg_name] = args[arg_position]
        else:
            position = arg_position - len(arg_names)
            if arg_name in kwargs:
                arg_dict[arg_name] = kwargs.pop(arg_name)
            else:
                try:
                    arg_dict[arg_name] = arg_defaults[position]
                except (IndexError, KeyError):
                    # Missing argument
                    raise ValueError('Wrong number of arguments for %s%s:\n'
                                     '     %s(%s, %s) was called.'
                        % (name,
                           inspect.formatargspec(*inspect.getargspec(func)),
                           name,
                           repr(args)[1:-1],
                           ', '.join('%s=%s' % (k, v)
                                    for k, v in kwargs.iteritems())
                           )
                        )

    varkwargs = dict()
    for arg_name, arg_value in kwargs.iteritems():
        if arg_name in arg_dict:
            arg_dict[arg_name] = arg_value
        elif arg_keywords is not None:
            varkwargs[arg_name] = arg_value
        else:
            raise TypeError("Ignore list for %s() contains an unexpected "
                            "keyword argument '%s'" % (name, arg_name))

    if arg_keywords is not None:
        arg_dict['**'] = varkwargs
    if arg_varargs is not None:
        varargs = args[arg_position + 1:]
        arg_dict['*'] = varargs

    # Now remove the arguments to be ignored
    for item in ignore_lst:
        if item in arg_dict:
            arg_dict.pop(item)
        else:
            raise ValueError("Ignore list: argument '%s' is not defined for "
            "function %s%s" %
                            (item, name,
                             inspect.formatargspec(arg_names,
                                                   arg_varargs,
                                                   arg_keywords,
                                                   arg_defaults,
                                                   )))
    # XXX: Return a sorted list of pairs?
    return arg_dict
