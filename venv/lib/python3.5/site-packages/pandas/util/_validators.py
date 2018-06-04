"""
Module that contains many useful utilities
for validating data or function arguments
"""
import warnings

from pandas.core.dtypes.common import is_bool


def _check_arg_length(fname, args, max_fname_arg_count, compat_args):
    """
    Checks whether 'args' has length of at most 'compat_args'. Raises
    a TypeError if that is not the case, similar to in Python when a
    function is called with too many arguments.

    """
    if max_fname_arg_count < 0:
        raise ValueError("'max_fname_arg_count' must be non-negative")

    if len(args) > len(compat_args):
        max_arg_count = len(compat_args) + max_fname_arg_count
        actual_arg_count = len(args) + max_fname_arg_count
        argument = 'argument' if max_arg_count == 1 else 'arguments'

        raise TypeError(
            "{fname}() takes at most {max_arg} {argument} "
            "({given_arg} given)".format(
                fname=fname, max_arg=max_arg_count,
                argument=argument, given_arg=actual_arg_count))


def _check_for_default_values(fname, arg_val_dict, compat_args):
    """
    Check that the keys in `arg_val_dict` are mapped to their
    default values as specified in `compat_args`.

    Note that this function is to be called only when it has been
    checked that arg_val_dict.keys() is a subset of compat_args

    """
    for key in arg_val_dict:
        # try checking equality directly with '=' operator,
        # as comparison may have been overridden for the left
        # hand object
        try:
            v1 = arg_val_dict[key]
            v2 = compat_args[key]

            # check for None-ness otherwise we could end up
            # comparing a numpy array vs None
            if (v1 is not None and v2 is None) or \
               (v1 is None and v2 is not None):
                match = False
            else:
                match = (v1 == v2)

            if not is_bool(match):
                raise ValueError("'match' is not a boolean")

        # could not compare them directly, so try comparison
        # using the 'is' operator
        except:
            match = (arg_val_dict[key] is compat_args[key])

        if not match:
            raise ValueError(("the '{arg}' parameter is not "
                              "supported in the pandas "
                              "implementation of {fname}()".
                              format(fname=fname, arg=key)))


def validate_args(fname, args, max_fname_arg_count, compat_args):
    """
    Checks whether the length of the `*args` argument passed into a function
    has at most `len(compat_args)` arguments and whether or not all of these
    elements in `args` are set to their default values.

    fname: str
        The name of the function being passed the `*args` parameter

    args: tuple
        The `*args` parameter passed into a function

    max_fname_arg_count: int
        The maximum number of arguments that the function `fname`
        can accept, excluding those in `args`. Used for displaying
        appropriate error messages. Must be non-negative.

    compat_args: OrderedDict
        A ordered dictionary of keys and their associated default values.
        In order to accommodate buggy behaviour in some versions of `numpy`,
        where a signature displayed keyword arguments but then passed those
        arguments **positionally** internally when calling downstream
        implementations, an ordered dictionary ensures that the original
        order of the keyword arguments is enforced. Note that if there is
        only one key, a generic dict can be passed in as well.

    Raises
    ------
    TypeError if `args` contains more values than there are `compat_args`
    ValueError if `args` contains values that do not correspond to those
    of the default values specified in `compat_args`

    """
    _check_arg_length(fname, args, max_fname_arg_count, compat_args)

    # We do this so that we can provide a more informative
    # error message about the parameters that we are not
    # supporting in the pandas implementation of 'fname'
    kwargs = dict(zip(compat_args, args))
    _check_for_default_values(fname, kwargs, compat_args)


def _check_for_invalid_keys(fname, kwargs, compat_args):
    """
    Checks whether 'kwargs' contains any keys that are not
    in 'compat_args' and raises a TypeError if there is one.

    """
    # set(dict) --> set of the dictionary's keys
    diff = set(kwargs) - set(compat_args)

    if diff:
        bad_arg = list(diff)[0]
        raise TypeError(("{fname}() got an unexpected "
                         "keyword argument '{arg}'".
                         format(fname=fname, arg=bad_arg)))


def validate_kwargs(fname, kwargs, compat_args):
    """
    Checks whether parameters passed to the **kwargs argument in a
    function `fname` are valid parameters as specified in `*compat_args`
    and whether or not they are set to their default values.

    Parameters
    ----------
    fname: str
        The name of the function being passed the `**kwargs` parameter

    kwargs: dict
        The `**kwargs` parameter passed into `fname`

    compat_args: dict
        A dictionary of keys that `kwargs` is allowed to have and their
        associated default values

    Raises
    ------
    TypeError if `kwargs` contains keys not in `compat_args`
    ValueError if `kwargs` contains keys in `compat_args` that do not
    map to the default values specified in `compat_args`

    """
    kwds = kwargs.copy()
    _check_for_invalid_keys(fname, kwargs, compat_args)
    _check_for_default_values(fname, kwds, compat_args)


def validate_args_and_kwargs(fname, args, kwargs,
                             max_fname_arg_count,
                             compat_args):
    """
    Checks whether parameters passed to the *args and **kwargs argument in a
    function `fname` are valid parameters as specified in `*compat_args`
    and whether or not they are set to their default values.

    Parameters
    ----------
    fname: str
        The name of the function being passed the `**kwargs` parameter

    args: tuple
        The `*args` parameter passed into a function

    kwargs: dict
        The `**kwargs` parameter passed into `fname`

    max_fname_arg_count: int
        The minimum number of arguments that the function `fname`
        requires, excluding those in `args`. Used for displaying
        appropriate error messages. Must be non-negative.

    compat_args: OrderedDict
        A ordered dictionary of keys that `kwargs` is allowed to
        have and their associated default values. Note that if there
        is only one key, a generic dict can be passed in as well.

    Raises
    ------
    TypeError if `args` contains more values than there are
    `compat_args` OR `kwargs` contains keys not in `compat_args`
    ValueError if `args` contains values not at the default value (`None`)
    `kwargs` contains keys in `compat_args` that do not map to the default
    value as specified in `compat_args`

    See Also
    --------
    validate_args : purely args validation
    validate_kwargs : purely kwargs validation

    """
    # Check that the total number of arguments passed in (i.e.
    # args and kwargs) does not exceed the length of compat_args
    _check_arg_length(fname, args + tuple(kwargs.values()),
                      max_fname_arg_count, compat_args)

    # Check there is no overlap with the positional and keyword
    # arguments, similar to what is done in actual Python functions
    args_dict = dict(zip(compat_args, args))

    for key in args_dict:
        if key in kwargs:
            raise TypeError("{fname}() got multiple values for keyword "
                            "argument '{arg}'".format(fname=fname, arg=key))

    kwargs.update(args_dict)
    validate_kwargs(fname, kwargs, compat_args)


def validate_bool_kwarg(value, arg_name):
    """ Ensures that argument passed in arg_name is of type bool. """
    if not (is_bool(value) or value is None):
        raise ValueError('For argument "{arg}" expected type bool, received '
                         'type {typ}.'.format(arg=arg_name,
                                              typ=type(value).__name__))
    return value


def validate_axis_style_args(data, args, kwargs, arg_name, method_name):
    """Argument handler for mixed index, columns / axis functions

    In an attempt to handle both `.method(index, columns)`, and
    `.method(arg, axis=.)`, we have to do some bad things to argument
    parsing. This translates all arguments to `{index=., columns=.}` style.

    Parameters
    ----------
    data : DataFrame or Panel
    arg : tuple
        All positional arguments from the user
    kwargs : dict
        All keyword arguments from the user
    arg_name, method_name : str
        Used for better error messages

    Returns
    -------
    kwargs : dict
        A dictionary of keyword arguments. Doesn't modify ``kwargs``
        inplace, so update them with the return value here.

    Examples
    --------
    >>> df._validate_axis_style_args((str.upper,), {'columns': id},
    ...                              'mapper', 'rename')
    {'columns': <function id>, 'index': <method 'upper' of 'str' objects>}

    This emits a warning
    >>> df._validate_axis_style_args((str.upper, id), {},
    ...                              'mapper', 'rename')
    {'columns': <function id>, 'index': <method 'upper' of 'str' objects>}
    """
    # TODO(PY3): Change to keyword-only args and remove all this

    out = {}
    # Goal: fill 'out' with index/columns-style arguments
    # like out = {'index': foo, 'columns': bar}

    # Start by validating for consistency
    if 'axis' in kwargs and any(x in kwargs for x in data._AXIS_NUMBERS):
        msg = "Cannot specify both 'axis' and any of 'index' or 'columns'."
        raise TypeError(msg)

    # First fill with explicit values provided by the user...
    if arg_name in kwargs:
        if args:
            msg = ("{} got multiple values for argument "
                   "'{}'".format(method_name, arg_name))
            raise TypeError(msg)

        axis = data._get_axis_name(kwargs.get('axis', 0))
        out[axis] = kwargs[arg_name]

    # More user-provided arguments, now from kwargs
    for k, v in kwargs.items():
        try:
            ax = data._get_axis_name(k)
        except ValueError:
            pass
        else:
            out[ax] = v

    # All user-provided kwargs have been handled now.
    # Now we supplement with positional arguments, emitting warnings
    # when there's ambiguity and raising when there's conflicts

    if len(args) == 0:
        pass  # It's up to the function to decide if this is valid
    elif len(args) == 1:
        axis = data._get_axis_name(kwargs.get('axis', 0))
        out[axis] = args[0]
    elif len(args) == 2:
        if 'axis' in kwargs:
            # Unambiguously wrong
            msg = ("Cannot specify both 'axis' and any of 'index' "
                   "or 'columns'")
            raise TypeError(msg)

        msg = ("Interpreting call\n\t'.{method_name}(a, b)' as "
               "\n\t'.{method_name}(index=a, columns=b)'.\nUse named "
               "arguments to remove any ambiguity. In the future, using "
               "positional arguments for 'index' or 'columns' will raise "
               " a 'TypeError'.")
        warnings.warn(msg.format(method_name=method_name,), FutureWarning,
                      stacklevel=4)
        out[data._AXIS_NAMES[0]] = args[0]
        out[data._AXIS_NAMES[1]] = args[1]
    else:
        msg = "Cannot specify all of '{}', 'index', 'columns'."
        raise TypeError(msg.format(arg_name))
    return out


def validate_fillna_kwargs(value, method, validate_scalar_dict_value=True):
    """Validate the keyword arguments to 'fillna'.

    This checks that exactly one of 'value' and 'method' is specified.
    If 'method' is specified, this validates that it's a valid method.

    Parameters
    ----------
    value, method : object
        The 'value' and 'method' keyword arguments for 'fillna'.
    validate_scalar_dict_value : bool, default True
        Whether to validate that 'value' is a scalar or dict. Specifically,
        validate that it is not a list or tuple.

    Returns
    -------
    value, method : object
    """
    from pandas.core.missing import clean_fill_method

    if value is None and method is None:
        raise ValueError("Must specify a fill 'value' or 'method'.")
    elif value is None and method is not None:
        method = clean_fill_method(method)

    elif value is not None and method is None:
        if validate_scalar_dict_value and isinstance(value, (list, tuple)):
            raise TypeError('"value" parameter must be a scalar or dict, but '
                            'you passed a "{0}"'.format(type(value).__name__))

    elif value is not None and method is not None:
        raise ValueError("Cannot specify both 'value' and 'method'.")

    return value, method
