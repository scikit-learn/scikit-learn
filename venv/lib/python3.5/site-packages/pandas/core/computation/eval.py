#!/usr/bin/env python

"""Top level ``eval`` module.
"""

import warnings
import tokenize
from pandas.io.formats.printing import pprint_thing
from pandas.core.computation.scope import _ensure_scope
from pandas.compat import string_types
from pandas.core.computation.engines import _engines
from pandas.util._validators import validate_bool_kwarg


def _check_engine(engine):
    """Make sure a valid engine is passed.

    Parameters
    ----------
    engine : str

    Raises
    ------
    KeyError
      * If an invalid engine is passed
    ImportError
      * If numexpr was requested but doesn't exist

    Returns
    -------
    string engine

    """
    from pandas.core.computation.check import _NUMEXPR_INSTALLED

    if engine is None:
        if _NUMEXPR_INSTALLED:
            engine = 'numexpr'
        else:
            engine = 'python'

    if engine not in _engines:
        valid = list(_engines.keys())
        raise KeyError('Invalid engine {engine!r} passed, valid engines are'
                       ' {valid}'.format(engine=engine, valid=valid))

    # TODO: validate this in a more general way (thinking of future engines
    # that won't necessarily be import-able)
    # Could potentially be done on engine instantiation
    if engine == 'numexpr':
        if not _NUMEXPR_INSTALLED:
            raise ImportError("'numexpr' is not installed or an "
                              "unsupported version. Cannot use "
                              "engine='numexpr' for query/eval "
                              "if 'numexpr' is not installed")

    return engine


def _check_parser(parser):
    """Make sure a valid parser is passed.

    Parameters
    ----------
    parser : str

    Raises
    ------
    KeyError
      * If an invalid parser is passed
    """
    from pandas.core.computation.expr import _parsers

    if parser not in _parsers:
        raise KeyError('Invalid parser {parser!r} passed, valid parsers are'
                       ' {valid}'.format(parser=parser, valid=_parsers.keys()))


def _check_resolvers(resolvers):
    if resolvers is not None:
        for resolver in resolvers:
            if not hasattr(resolver, '__getitem__'):
                name = type(resolver).__name__
                raise TypeError('Resolver of type {name!r} does not implement '
                                'the __getitem__ method'.format(name=name))


def _check_expression(expr):
    """Make sure an expression is not an empty string

    Parameters
    ----------
    expr : object
        An object that can be converted to a string

    Raises
    ------
    ValueError
      * If expr is an empty string
    """
    if not expr:
        raise ValueError("expr cannot be an empty string")


def _convert_expression(expr):
    """Convert an object to an expression.

    Thus function converts an object to an expression (a unicode string) and
    checks to make sure it isn't empty after conversion. This is used to
    convert operators to their string representation for recursive calls to
    :func:`~pandas.eval`.

    Parameters
    ----------
    expr : object
        The object to be converted to a string.

    Returns
    -------
    s : unicode
        The string representation of an object.

    Raises
    ------
    ValueError
      * If the expression is empty.
    """
    s = pprint_thing(expr)
    _check_expression(s)
    return s


def _check_for_locals(expr, stack_level, parser):
    from pandas.core.computation.expr import tokenize_string

    at_top_of_stack = stack_level == 0
    not_pandas_parser = parser != 'pandas'

    if not_pandas_parser:
        msg = "The '@' prefix is only supported by the pandas parser"
    elif at_top_of_stack:
        msg = ("The '@' prefix is not allowed in "
               "top-level eval calls, \nplease refer to "
               "your variables by name without the '@' "
               "prefix")

    if at_top_of_stack or not_pandas_parser:
        for toknum, tokval in tokenize_string(expr):
            if toknum == tokenize.OP and tokval == '@':
                raise SyntaxError(msg)


def eval(expr, parser='pandas', engine=None, truediv=True,
         local_dict=None, global_dict=None, resolvers=(), level=0,
         target=None, inplace=False):
    """Evaluate a Python expression as a string using various backends.

    The following arithmetic operations are supported: ``+``, ``-``, ``*``,
    ``/``, ``**``, ``%``, ``//`` (python engine only) along with the following
    boolean operations: ``|`` (or), ``&`` (and), and ``~`` (not).
    Additionally, the ``'pandas'`` parser allows the use of :keyword:`and`,
    :keyword:`or`, and :keyword:`not` with the same semantics as the
    corresponding bitwise operators.  :class:`~pandas.Series` and
    :class:`~pandas.DataFrame` objects are supported and behave as they would
    with plain ol' Python evaluation.

    Parameters
    ----------
    expr : str or unicode
        The expression to evaluate. This string cannot contain any Python
        `statements
        <https://docs.python.org/3/reference/simple_stmts.html#simple-statements>`__,
        only Python `expressions
        <https://docs.python.org/3/reference/simple_stmts.html#expression-statements>`__.
    parser : string, default 'pandas', {'pandas', 'python'}
        The parser to use to construct the syntax tree from the expression. The
        default of ``'pandas'`` parses code slightly different than standard
        Python. Alternatively, you can parse an expression using the
        ``'python'`` parser to retain strict Python semantics.  See the
        :ref:`enhancing performance <enhancingperf.eval>` documentation for
        more details.
    engine : string or None, default 'numexpr', {'python', 'numexpr'}

        The engine used to evaluate the expression. Supported engines are

        - None         : tries to use ``numexpr``, falls back to ``python``
        - ``'numexpr'``: This default engine evaluates pandas objects using
                         numexpr for large speed ups in complex expressions
                         with large frames.
        - ``'python'``: Performs operations as if you had ``eval``'d in top
                        level python. This engine is generally not that useful.

        More backends may be available in the future.

    truediv : bool, optional
        Whether to use true division, like in Python >= 3
    local_dict : dict or None, optional
        A dictionary of local variables, taken from locals() by default.
    global_dict : dict or None, optional
        A dictionary of global variables, taken from globals() by default.
    resolvers : list of dict-like or None, optional
        A list of objects implementing the ``__getitem__`` special method that
        you can use to inject an additional collection of namespaces to use for
        variable lookup. For example, this is used in the
        :meth:`~pandas.DataFrame.query` method to inject the
        ``DataFrame.index`` and ``DataFrame.columns``
        variables that refer to their respective :class:`~pandas.DataFrame`
        instance attributes.
    level : int, optional
        The number of prior stack frames to traverse and add to the current
        scope. Most users will **not** need to change this parameter.
    target : object, optional, default None
        This is the target object for assignment. It is used when there is
        variable assignment in the expression. If so, then `target` must
        support item assignment with string keys, and if a copy is being
        returned, it must also support `.copy()`.
    inplace : bool, default False
        If `target` is provided, and the expression mutates `target`, whether
        to modify `target` inplace. Otherwise, return a copy of `target` with
        the mutation.

    Returns
    -------
    ndarray, numeric scalar, DataFrame, Series

    Raises
    ------
    ValueError
        There are many instances where such an error can be raised:

        - `target=None`, but the expression is multiline.
        - The expression is multiline, but not all them have item assignment.
          An example of such an arrangement is this:

          a = b + 1
          a + 2

          Here, there are expressions on different lines, making it multiline,
          but the last line has no variable assigned to the output of `a + 2`.
        - `inplace=True`, but the expression is missing item assignment.
        - Item assignment is provided, but the `target` does not support
          string item assignment.
        - Item assignment is provided and `inplace=False`, but the `target`
          does not support the `.copy()` method

    Notes
    -----
    The ``dtype`` of any objects involved in an arithmetic ``%`` operation are
    recursively cast to ``float64``.

    See the :ref:`enhancing performance <enhancingperf.eval>` documentation for
    more details.

    See Also
    --------
    pandas.DataFrame.query
    pandas.DataFrame.eval
    """
    from pandas.core.computation.expr import Expr

    inplace = validate_bool_kwarg(inplace, "inplace")

    if isinstance(expr, string_types):
        _check_expression(expr)
        exprs = [e.strip() for e in expr.splitlines() if e.strip() != '']
    else:
        exprs = [expr]
    multi_line = len(exprs) > 1

    if multi_line and target is None:
        raise ValueError("multi-line expressions are only valid in the "
                         "context of data, use DataFrame.eval")

    ret = None
    first_expr = True
    target_modified = False

    for expr in exprs:
        expr = _convert_expression(expr)
        engine = _check_engine(engine)
        _check_parser(parser)
        _check_resolvers(resolvers)
        _check_for_locals(expr, level, parser)

        # get our (possibly passed-in) scope
        env = _ensure_scope(level + 1, global_dict=global_dict,
                            local_dict=local_dict, resolvers=resolvers,
                            target=target)

        parsed_expr = Expr(expr, engine=engine, parser=parser, env=env,
                           truediv=truediv)

        # construct the engine and evaluate the parsed expression
        eng = _engines[engine]
        eng_inst = eng(parsed_expr)
        ret = eng_inst.evaluate()

        if parsed_expr.assigner is None:
            if multi_line:
                raise ValueError("Multi-line expressions are only valid"
                                 " if all expressions contain an assignment")
            elif inplace:
                raise ValueError("Cannot operate inplace "
                                 "if there is no assignment")

        # assign if needed
        assigner = parsed_expr.assigner
        if env.target is not None and assigner is not None:
            target_modified = True

            # if returning a copy, copy only on the first assignment
            if not inplace and first_expr:
                try:
                    target = env.target.copy()
                except AttributeError:
                    raise ValueError("Cannot return a copy of the target")
            else:
                target = env.target

            # TypeError is most commonly raised (e.g. int, list), but you
            # get IndexError if you try to do this assignment on np.ndarray.
            # we will ignore numpy warnings here; e.g. if trying
            # to use a non-numeric indexer
            try:
                with warnings.catch_warnings(record=True):
                    target[assigner] = ret
            except (TypeError, IndexError):
                raise ValueError("Cannot assign expression output to target")

            if not resolvers:
                resolvers = ({assigner: ret},)
            else:
                # existing resolver needs updated to handle
                # case of mutating existing column in copy
                for resolver in resolvers:
                    if assigner in resolver:
                        resolver[assigner] = ret
                        break
                else:
                    resolvers += ({assigner: ret},)

            ret = None
            first_expr = False

    # We want to exclude `inplace=None` as being False.
    if inplace is False:
        return target if target_modified else ret
