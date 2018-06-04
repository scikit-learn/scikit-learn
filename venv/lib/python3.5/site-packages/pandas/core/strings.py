import numpy as np

from pandas.compat import zip
from pandas.core.dtypes.generic import ABCSeries, ABCIndex
from pandas.core.dtypes.missing import isna, notna
from pandas.core.dtypes.common import (
    is_bool_dtype,
    is_categorical_dtype,
    is_object_dtype,
    is_string_like,
    is_list_like,
    is_scalar,
    is_integer,
    is_re)

import pandas.core.common as com
from pandas.core.algorithms import take_1d
import pandas.compat as compat
from pandas.core.base import NoNewAttributesMixin
from pandas.util._decorators import Appender
import re
import pandas._libs.lib as lib
import pandas._libs.ops as libops
import warnings
import textwrap
import codecs

_cpython_optimized_encoders = (
    "utf-8", "utf8", "latin-1", "latin1", "iso-8859-1", "mbcs", "ascii"
)
_cpython_optimized_decoders = _cpython_optimized_encoders + (
    "utf-16", "utf-32"
)

_shared_docs = dict()


def _get_array_list(arr, others):
    """
    Auxiliary function for :func:`str_cat`

    Parameters
    ----------
    arr : ndarray
        The left-most ndarray of the concatenation
    others : list, ndarray, Series
        The rest of the content to concatenate. If list of list-likes,
        all elements must be passable to ``np.asarray``.

    Returns
    -------
    list
        List of all necessary arrays
    """
    from pandas.core.series import Series

    if len(others) and isinstance(com._values_from_object(others)[0],
                                  (list, np.ndarray, Series)):
        arrays = [arr] + list(others)
    else:
        arrays = [arr, others]

    return [np.asarray(x, dtype=object) for x in arrays]


def str_cat(arr, others=None, sep=None, na_rep=None):
    """
    Auxiliary function for :meth:`str.cat`

    If `others` is specified, this function concatenates the Series/Index
    and elements of `others` element-wise.
    If `others` is not being passed then all values in the Series are
    concatenated in a single string with a given `sep`.

    Parameters
    ----------
    others : list-like, or list of list-likes, optional
        List-likes (or a list of them) of the same length as calling object.
        If None, returns str concatenating strings of the Series.
    sep : string or None, default None
        If None, concatenates without any separator.
    na_rep : string or None, default None
        If None, NA in the series are ignored.

    Returns
    -------
    concat
        ndarray containing concatenated results (if `others is not None`)
        or str (if `others is None`)
    """
    if sep is None:
        sep = ''

    if others is not None:
        arrays = _get_array_list(arr, others)

        n = _length_check(arrays)
        masks = np.array([isna(x) for x in arrays])
        cats = None

        if na_rep is None:
            na_mask = np.logical_or.reduce(masks, axis=0)

            result = np.empty(n, dtype=object)
            np.putmask(result, na_mask, np.nan)

            notmask = ~na_mask

            tuples = zip(*[x[notmask] for x in arrays])
            cats = [sep.join(tup) for tup in tuples]

            result[notmask] = cats
        else:
            for i, x in enumerate(arrays):
                x = np.where(masks[i], na_rep, x)
                if cats is None:
                    cats = x
                else:
                    cats = cats + sep + x

            result = cats

        return result
    else:
        arr = np.asarray(arr, dtype=object)
        mask = isna(arr)
        if na_rep is None and mask.any():
            if sep == '':
                na_rep = ''
            else:
                return sep.join(arr[notna(arr)])
        return sep.join(np.where(mask, na_rep, arr))


def _length_check(others):
    n = None
    for x in others:
        try:
            if n is None:
                n = len(x)
            elif len(x) != n:
                raise ValueError('All arrays must be same length')
        except TypeError:
            raise ValueError('Must pass arrays containing strings to str_cat')
    return n


def _na_map(f, arr, na_result=np.nan, dtype=object):
    # should really _check_ for NA
    return _map(f, arr, na_mask=True, na_value=na_result, dtype=dtype)


def _map(f, arr, na_mask=False, na_value=np.nan, dtype=object):
    if not len(arr):
        return np.ndarray(0, dtype=dtype)

    if isinstance(arr, ABCSeries):
        arr = arr.values
    if not isinstance(arr, np.ndarray):
        arr = np.asarray(arr, dtype=object)
    if na_mask:
        mask = isna(arr)
        try:
            convert = not all(mask)
            result = lib.map_infer_mask(arr, f, mask.view(np.uint8), convert)
        except (TypeError, AttributeError) as e:
            # Reraise the exception if callable `f` got wrong number of args.
            # The user may want to be warned by this, instead of getting NaN
            if compat.PY2:
                p_err = r'takes (no|(exactly|at (least|most)) ?\d+) arguments?'
            else:
                p_err = (r'((takes)|(missing)) (?(2)from \d+ to )?\d+ '
                         r'(?(3)required )positional arguments?')

            if len(e.args) >= 1 and re.search(p_err, e.args[0]):
                raise e

            def g(x):
                try:
                    return f(x)
                except (TypeError, AttributeError):
                    return na_value

            return _map(g, arr, dtype=dtype)
        if na_value is not np.nan:
            np.putmask(result, mask, na_value)
            if result.dtype == object:
                result = lib.maybe_convert_objects(result)
        return result
    else:
        return lib.map_infer(arr, f)


def str_count(arr, pat, flags=0):
    """
    Count occurrences of pattern in each string of the Series/Index.

    This function is used to count the number of times a particular regex
    pattern is repeated in each of the string elements of the
    :class:`~pandas.Series`.

    Parameters
    ----------
    pat : str
        Valid regular expression.
    flags : int, default 0, meaning no flags
        Flags for the `re` module. For a complete list, `see here
        <https://docs.python.org/3/howto/regex.html#compilation-flags>`_.
    **kwargs
        For compatability with other string methods. Not used.

    Returns
    -------
    counts : Series or Index
        Same type as the calling object containing the integer counts.

    Notes
    -----
    Some characters need to be escaped when passing in `pat`.
    eg. ``'$'`` has a special meaning in regex and must be escaped when
    finding this literal character.

    See Also
    --------
    re : Standard library module for regular expressions.
    str.count : Standard library version, without regular expression support.

    Examples
    --------
    >>> s = pd.Series(['A', 'B', 'Aaba', 'Baca', np.nan, 'CABA', 'cat'])
    >>> s.str.count('a')
    0    0.0
    1    0.0
    2    2.0
    3    2.0
    4    NaN
    5    0.0
    6    1.0
    dtype: float64

    Escape ``'$'`` to find the literal dollar sign.

    >>> s = pd.Series(['$', 'B', 'Aab$', '$$ca', 'C$B$', 'cat'])
    >>> s.str.count('\$')
    0    1
    1    0
    2    1
    3    2
    4    2
    5    0
    dtype: int64

    This is also available on Index

    >>> pd.Index(['A', 'A', 'Aaba', 'cat']).str.count('a')
    Int64Index([0, 0, 2, 1], dtype='int64')
    """
    regex = re.compile(pat, flags=flags)
    f = lambda x: len(regex.findall(x))
    return _na_map(f, arr, dtype=int)


def str_contains(arr, pat, case=True, flags=0, na=np.nan, regex=True):
    """
    Test if pattern or regex is contained within a string of a Series or Index.

    Return boolean Series or Index based on whether a given pattern or regex is
    contained within a string of a Series or Index.

    Parameters
    ----------
    pat : str
        Character sequence or regular expression.
    case : bool, default True
        If True, case sensitive.
    flags : int, default 0 (no flags)
        Flags to pass through to the re module, e.g. re.IGNORECASE.
    na : default NaN
        Fill value for missing values.
    regex : bool, default True
        If True, assumes the pat is a regular expression.

        If False, treats the pat as a literal string.

    Returns
    -------
    Series or Index of boolean values
        A Series or Index of boolean values indicating whether the
        given pattern is contained within the string of each element
        of the Series or Index.

    See Also
    --------
    match : analogous, but stricter, relying on re.match instead of re.search

    Examples
    --------

    Returning a Series of booleans using only a literal pattern.

    >>> s1 = pd.Series(['Mouse', 'dog', 'house and parrot', '23', np.NaN])
    >>> s1.str.contains('og', regex=False)
    0    False
    1     True
    2    False
    3    False
    4      NaN
    dtype: object

    Returning an Index of booleans using only a literal pattern.

    >>> ind = pd.Index(['Mouse', 'dog', 'house and parrot', '23.0', np.NaN])
    >>> ind.str.contains('23', regex=False)
    Index([False, False, False, True, nan], dtype='object')

    Specifying case sensitivity using `case`.

    >>> s1.str.contains('oG', case=True, regex=True)
    0    False
    1    False
    2    False
    3    False
    4      NaN
    dtype: object

    Specifying `na` to be `False` instead of `NaN` replaces NaN values
    with `False`. If Series or Index does not contain NaN values
    the resultant dtype will be `bool`, otherwise, an `object` dtype.

    >>> s1.str.contains('og', na=False, regex=True)
    0    False
    1     True
    2    False
    3    False
    4    False
    dtype: bool

    Returning 'house' and 'parrot' within same string.

    >>> s1.str.contains('house|parrot', regex=True)
    0    False
    1    False
    2     True
    3    False
    4      NaN
    dtype: object

    Ignoring case sensitivity using `flags` with regex.

    >>> import re
    >>> s1.str.contains('PARROT', flags=re.IGNORECASE, regex=True)
    0    False
    1    False
    2     True
    3    False
    4      NaN
    dtype: object

    Returning any digit using regular expression.

    >>> s1.str.contains('\d', regex=True)
    0    False
    1    False
    2    False
    3     True
    4      NaN
    dtype: object

    Ensure `pat` is a not a literal pattern when `regex` is set to True.
    Note in the following example one might expect only `s2[1]` and `s2[3]` to
    return `True`. However, '.0' as a regex matches any character
    followed by a 0.

    >>> s2 = pd.Series(['40','40.0','41','41.0','35'])
    >>> s2.str.contains('.0', regex=True)
    0     True
    1     True
    2    False
    3     True
    4    False
    dtype: bool
    """
    if regex:
        if not case:
            flags |= re.IGNORECASE

        regex = re.compile(pat, flags=flags)

        if regex.groups > 0:
            warnings.warn("This pattern has match groups. To actually get the"
                          " groups, use str.extract.", UserWarning,
                          stacklevel=3)

        f = lambda x: bool(regex.search(x))
    else:
        if case:
            f = lambda x: pat in x
        else:
            upper_pat = pat.upper()
            f = lambda x: upper_pat in x
            uppered = _na_map(lambda x: x.upper(), arr)
            return _na_map(f, uppered, na, dtype=bool)
    return _na_map(f, arr, na, dtype=bool)


def str_startswith(arr, pat, na=np.nan):
    """
    Test if the start of each string element matches a pattern.

    Equivalent to :meth:`str.startswith`.

    Parameters
    ----------
    pat : str
        Character sequence. Regular expressions are not accepted.
    na : object, default NaN
        Object shown if element tested is not a string.

    Returns
    -------
    Series or Index of bool
        A Series of booleans indicating whether the given pattern matches
        the start of each string element.

    See Also
    --------
    str.startswith : Python standard library string method.
    Series.str.endswith : Same as startswith, but tests the end of string.
    Series.str.contains : Tests if string element contains a pattern.

    Examples
    --------
    >>> s = pd.Series(['bat', 'Bear', 'cat', np.nan])
    >>> s
    0     bat
    1    Bear
    2     cat
    3     NaN
    dtype: object

    >>> s.str.startswith('b')
    0     True
    1    False
    2    False
    3      NaN
    dtype: object

    Specifying `na` to be `False` instead of `NaN`.

    >>> s.str.startswith('b', na=False)
    0     True
    1    False
    2    False
    3    False
    dtype: bool
    """
    f = lambda x: x.startswith(pat)
    return _na_map(f, arr, na, dtype=bool)


def str_endswith(arr, pat, na=np.nan):
    """
    Test if the end of each string element matches a pattern.

    Equivalent to :meth:`str.endswith`.

    Parameters
    ----------
    pat : str
        Character sequence. Regular expressions are not accepted.
    na : object, default NaN
        Object shown if element tested is not a string.

    Returns
    -------
    Series or Index of bool
        A Series of booleans indicating whether the given pattern matches
        the end of each string element.

    See Also
    --------
    str.endswith : Python standard library string method.
    Series.str.startswith : Same as endswith, but tests the start of string.
    Series.str.contains : Tests if string element contains a pattern.

    Examples
    --------
    >>> s = pd.Series(['bat', 'bear', 'caT', np.nan])
    >>> s
    0     bat
    1    bear
    2     caT
    3     NaN
    dtype: object

    >>> s.str.endswith('t')
    0     True
    1    False
    2    False
    3      NaN
    dtype: object

    Specifying `na` to be `False` instead of `NaN`.

    >>> s.str.endswith('t', na=False)
    0     True
    1    False
    2    False
    3    False
    dtype: bool
    """
    f = lambda x: x.endswith(pat)
    return _na_map(f, arr, na, dtype=bool)


def str_replace(arr, pat, repl, n=-1, case=None, flags=0, regex=True):
    r"""
    Replace occurrences of pattern/regex in the Series/Index with
    some other string. Equivalent to :meth:`str.replace` or
    :func:`re.sub`.

    Parameters
    ----------
    pat : string or compiled regex
        String can be a character sequence or regular expression.

        .. versionadded:: 0.20.0
            `pat` also accepts a compiled regex.

    repl : string or callable
        Replacement string or a callable. The callable is passed the regex
        match object and must return a replacement string to be used.
        See :func:`re.sub`.

        .. versionadded:: 0.20.0
            `repl` also accepts a callable.

    n : int, default -1 (all)
        Number of replacements to make from start
    case : boolean, default None
        - If True, case sensitive (the default if `pat` is a string)
        - Set to False for case insensitive
        - Cannot be set if `pat` is a compiled regex
    flags : int, default 0 (no flags)
        - re module flags, e.g. re.IGNORECASE
        - Cannot be set if `pat` is a compiled regex
    regex : boolean, default True
        - If True, assumes the passed-in pattern is a regular expression.
        - If False, treats the pattern as a literal string
        - Cannot be set to False if `pat` is a compiled regex or `repl` is
          a callable.

        .. versionadded:: 0.23.0

    Returns
    -------
    replaced : Series/Index of objects

    Raises
    ------
    ValueError
        * if `regex` is False and `repl` is a callable or `pat` is a compiled
          regex
        * if `pat` is a compiled regex and `case` or `flags` is set

    Notes
    -----
    When `pat` is a compiled regex, all flags should be included in the
    compiled regex. Use of `case`, `flags`, or `regex=False` with a compiled
    regex will raise an error.

    Examples
    --------
    When `pat` is a string and `regex` is True (the default), the given `pat`
    is compiled as a regex. When `repl` is a string, it replaces matching
    regex patterns as with :meth:`re.sub`. NaN value(s) in the Series are
    left as is:

    >>> pd.Series(['foo', 'fuz', np.nan]).str.replace('f.', 'ba', regex=True)
    0    bao
    1    baz
    2    NaN
    dtype: object

    When `pat` is a string and `regex` is False, every `pat` is replaced with
    `repl` as with :meth:`str.replace`:

    >>> pd.Series(['f.o', 'fuz', np.nan]).str.replace('f.', 'ba', regex=False)
    0    bao
    1    fuz
    2    NaN
    dtype: object

    When `repl` is a callable, it is called on every `pat` using
    :func:`re.sub`. The callable should expect one positional argument
    (a regex object) and return a string.

    To get the idea:

    >>> pd.Series(['foo', 'fuz', np.nan]).str.replace('f', repr)
    0    <_sre.SRE_Match object; span=(0, 1), match='f'>oo
    1    <_sre.SRE_Match object; span=(0, 1), match='f'>uz
    2                                                  NaN
    dtype: object

    Reverse every lowercase alphabetic word:

    >>> repl = lambda m: m.group(0)[::-1]
    >>> pd.Series(['foo 123', 'bar baz', np.nan]).str.replace(r'[a-z]+', repl)
    0    oof 123
    1    rab zab
    2        NaN
    dtype: object

    Using regex groups (extract second group and swap case):

    >>> pat = r"(?P<one>\w+) (?P<two>\w+) (?P<three>\w+)"
    >>> repl = lambda m: m.group('two').swapcase()
    >>> pd.Series(['One Two Three', 'Foo Bar Baz']).str.replace(pat, repl)
    0    tWO
    1    bAR
    dtype: object

    Using a compiled regex with flags

    >>> regex_pat = re.compile(r'FUZ', flags=re.IGNORECASE)
    >>> pd.Series(['foo', 'fuz', np.nan]).str.replace(regex_pat, 'bar')
    0    foo
    1    bar
    2    NaN
    dtype: object

    """

    # Check whether repl is valid (GH 13438, GH 15055)
    if not (is_string_like(repl) or callable(repl)):
        raise TypeError("repl must be a string or callable")

    is_compiled_re = is_re(pat)
    if regex:
        if is_compiled_re:
            if (case is not None) or (flags != 0):
                raise ValueError("case and flags cannot be set"
                                 " when pat is a compiled regex")
        else:
            # not a compiled regex
            # set default case
            if case is None:
                case = True

            # add case flag, if provided
            if case is False:
                flags |= re.IGNORECASE
        if is_compiled_re or len(pat) > 1 or flags or callable(repl):
            n = n if n >= 0 else 0
            compiled = re.compile(pat, flags=flags)
            f = lambda x: compiled.sub(repl=repl, string=x, count=n)
        else:
            f = lambda x: x.replace(pat, repl, n)
    else:
        if is_compiled_re:
            raise ValueError("Cannot use a compiled regex as replacement "
                             "pattern with regex=False")
        if callable(repl):
            raise ValueError("Cannot use a callable replacement when "
                             "regex=False")
        f = lambda x: x.replace(pat, repl, n)

    return _na_map(f, arr)


def str_repeat(arr, repeats):
    """
    Duplicate each string in the Series/Index by indicated number
    of times.

    Parameters
    ----------
    repeats : int or array
        Same value for all (int) or different value per (array)

    Returns
    -------
    repeated : Series/Index of objects
    """
    if is_scalar(repeats):

        def rep(x):
            try:
                return compat.binary_type.__mul__(x, repeats)
            except TypeError:
                return compat.text_type.__mul__(x, repeats)

        return _na_map(rep, arr)
    else:

        def rep(x, r):
            try:
                return compat.binary_type.__mul__(x, r)
            except TypeError:
                return compat.text_type.__mul__(x, r)

        repeats = np.asarray(repeats, dtype=object)
        result = libops.vec_binop(com._values_from_object(arr), repeats, rep)
        return result


def str_match(arr, pat, case=True, flags=0, na=np.nan, as_indexer=None):
    """
    Determine if each string matches a regular expression.

    Parameters
    ----------
    pat : string
        Character sequence or regular expression
    case : boolean, default True
        If True, case sensitive
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE
    na : default NaN, fill value for missing values.
    as_indexer
        .. deprecated:: 0.21.0

    Returns
    -------
    Series/array of boolean values

    See Also
    --------
    contains : analogous, but less strict, relying on re.search instead of
        re.match
    extract : extract matched groups

    """
    if not case:
        flags |= re.IGNORECASE

    regex = re.compile(pat, flags=flags)

    if (as_indexer is False) and (regex.groups > 0):
        raise ValueError("as_indexer=False with a pattern with groups is no "
                         "longer supported. Use '.str.extract(pat)' instead")
    elif as_indexer is not None:
        # Previously, this keyword was used for changing the default but
        # deprecated behaviour. This keyword is now no longer needed.
        warnings.warn("'as_indexer' keyword was specified but is ignored "
                      "(match now returns a boolean indexer by default), "
                      "and will be removed in a future version.",
                      FutureWarning, stacklevel=3)

    dtype = bool
    f = lambda x: bool(regex.match(x))

    return _na_map(f, arr, na, dtype=dtype)


def _get_single_group_name(rx):
    try:
        return list(rx.groupindex.keys()).pop()
    except IndexError:
        return None


def _groups_or_na_fun(regex):
    """Used in both extract_noexpand and extract_frame"""
    if regex.groups == 0:
        raise ValueError("pattern contains no capture groups")
    empty_row = [np.nan] * regex.groups

    def f(x):
        if not isinstance(x, compat.string_types):
            return empty_row
        m = regex.search(x)
        if m:
            return [np.nan if item is None else item for item in m.groups()]
        else:
            return empty_row
    return f


def _str_extract_noexpand(arr, pat, flags=0):
    """
    Find groups in each string in the Series using passed regular
    expression. This function is called from
    str_extract(expand=False), and can return Series, DataFrame, or
    Index.

    """
    from pandas import DataFrame, Index

    regex = re.compile(pat, flags=flags)
    groups_or_na = _groups_or_na_fun(regex)

    if regex.groups == 1:
        result = np.array([groups_or_na(val)[0] for val in arr], dtype=object)
        name = _get_single_group_name(regex)
    else:
        if isinstance(arr, Index):
            raise ValueError("only one regex group is supported with Index")
        name = None
        names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
        columns = [names.get(1 + i, i) for i in range(regex.groups)]
        if arr.empty:
            result = DataFrame(columns=columns, dtype=object)
        else:
            result = DataFrame(
                [groups_or_na(val) for val in arr],
                columns=columns,
                index=arr.index,
                dtype=object)
    return result, name


def _str_extract_frame(arr, pat, flags=0):
    """
    For each subject string in the Series, extract groups from the
    first match of regular expression pat. This function is called from
    str_extract(expand=True), and always returns a DataFrame.

    """
    from pandas import DataFrame

    regex = re.compile(pat, flags=flags)
    groups_or_na = _groups_or_na_fun(regex)
    names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
    columns = [names.get(1 + i, i) for i in range(regex.groups)]

    if len(arr) == 0:
        return DataFrame(columns=columns, dtype=object)
    try:
        result_index = arr.index
    except AttributeError:
        result_index = None
    return DataFrame(
        [groups_or_na(val) for val in arr],
        columns=columns,
        index=result_index,
        dtype=object)


def str_extract(arr, pat, flags=0, expand=True):
    r"""
    For each subject string in the Series, extract groups from the
    first match of regular expression pat.

    Parameters
    ----------
    pat : string
        Regular expression pattern with capturing groups
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    expand : bool, default True
        * If True, return DataFrame.
        * If False, return Series/Index/DataFrame.

        .. versionadded:: 0.18.0

    Returns
    -------
    DataFrame with one row for each subject string, and one column for
    each group. Any capture group names in regular expression pat will
    be used for column names; otherwise capture group numbers will be
    used. The dtype of each result column is always object, even when
    no match is found. If expand=False and pat has only one capture group,
    then return a Series (if subject is a Series) or Index (if subject
    is an Index).

    See Also
    --------
    extractall : returns all matches (not just the first match)

    Examples
    --------
    A pattern with two groups will return a DataFrame with two columns.
    Non-matches will be NaN.

    >>> s = Series(['a1', 'b2', 'c3'])
    >>> s.str.extract(r'([ab])(\d)')
         0    1
    0    a    1
    1    b    2
    2  NaN  NaN

    A pattern may contain optional groups.

    >>> s.str.extract(r'([ab])?(\d)')
         0  1
    0    a  1
    1    b  2
    2  NaN  3

    Named groups will become column names in the result.

    >>> s.str.extract(r'(?P<letter>[ab])(?P<digit>\d)')
      letter digit
    0      a     1
    1      b     2
    2    NaN   NaN

    A pattern with one group will return a DataFrame with one column
    if expand=True.

    >>> s.str.extract(r'[ab](\d)', expand=True)
         0
    0    1
    1    2
    2  NaN

    A pattern with one group will return a Series if expand=False.

    >>> s.str.extract(r'[ab](\d)', expand=False)
    0      1
    1      2
    2    NaN
    dtype: object

    """
    if not isinstance(expand, bool):
        raise ValueError("expand must be True or False")
    if expand:
        return _str_extract_frame(arr._orig, pat, flags=flags)
    else:
        result, name = _str_extract_noexpand(arr._data, pat, flags=flags)
        return arr._wrap_result(result, name=name, expand=expand)


def str_extractall(arr, pat, flags=0):
    r"""
    For each subject string in the Series, extract groups from all
    matches of regular expression pat. When each subject string in the
    Series has exactly one match, extractall(pat).xs(0, level='match')
    is the same as extract(pat).

    .. versionadded:: 0.18.0

    Parameters
    ----------
    pat : string
        Regular expression pattern with capturing groups
    flags : int, default 0 (no flags)
        re module flags, e.g. re.IGNORECASE

    Returns
    -------
    A DataFrame with one row for each match, and one column for each
    group. Its rows have a MultiIndex with first levels that come from
    the subject Series. The last level is named 'match' and indicates
    the order in the subject. Any capture group names in regular
    expression pat will be used for column names; otherwise capture
    group numbers will be used.

    See Also
    --------
    extract : returns first match only (not all matches)

    Examples
    --------
    A pattern with one group will return a DataFrame with one column.
    Indices with no matches will not appear in the result.

    >>> s = Series(["a1a2", "b1", "c1"], index=["A", "B", "C"])
    >>> s.str.extractall(r"[ab](\d)")
             0
      match
    A 0      1
      1      2
    B 0      1

    Capture group names are used for column names of the result.

    >>> s.str.extractall(r"[ab](?P<digit>\d)")
            digit
      match
    A 0         1
      1         2
    B 0         1

    A pattern with two groups will return a DataFrame with two columns.

    >>> s.str.extractall(r"(?P<letter>[ab])(?P<digit>\d)")
            letter digit
      match
    A 0          a     1
      1          a     2
    B 0          b     1

    Optional groups that do not match are NaN in the result.

    >>> s.str.extractall(r"(?P<letter>[ab])?(?P<digit>\d)")
            letter digit
      match
    A 0          a     1
      1          a     2
    B 0          b     1
    C 0        NaN     1

    """

    regex = re.compile(pat, flags=flags)
    # the regex must contain capture groups.
    if regex.groups == 0:
        raise ValueError("pattern contains no capture groups")

    if isinstance(arr, ABCIndex):
        arr = arr.to_series().reset_index(drop=True)

    names = dict(zip(regex.groupindex.values(), regex.groupindex.keys()))
    columns = [names.get(1 + i, i) for i in range(regex.groups)]
    match_list = []
    index_list = []
    is_mi = arr.index.nlevels > 1

    for subject_key, subject in arr.iteritems():
        if isinstance(subject, compat.string_types):

            if not is_mi:
                subject_key = (subject_key, )

            for match_i, match_tuple in enumerate(regex.findall(subject)):
                if isinstance(match_tuple, compat.string_types):
                    match_tuple = (match_tuple,)
                na_tuple = [np.NaN if group == "" else group
                            for group in match_tuple]
                match_list.append(na_tuple)
                result_key = tuple(subject_key + (match_i, ))
                index_list.append(result_key)

    from pandas import MultiIndex
    index = MultiIndex.from_tuples(
        index_list, names=arr.index.names + ["match"])

    result = arr._constructor_expanddim(match_list, index=index,
                                        columns=columns)
    return result


def str_get_dummies(arr, sep='|'):
    """
    Split each string in the Series by sep and return a frame of
    dummy/indicator variables.

    Parameters
    ----------
    sep : string, default "|"
        String to split on.

    Returns
    -------
    dummies : DataFrame

    Examples
    --------
    >>> Series(['a|b', 'a', 'a|c']).str.get_dummies()
       a  b  c
    0  1  1  0
    1  1  0  0
    2  1  0  1

    >>> Series(['a|b', np.nan, 'a|c']).str.get_dummies()
       a  b  c
    0  1  1  0
    1  0  0  0
    2  1  0  1

    See Also
    --------
    pandas.get_dummies
    """
    arr = arr.fillna('')
    try:
        arr = sep + arr + sep
    except TypeError:
        arr = sep + arr.astype(str) + sep

    tags = set()
    for ts in arr.str.split(sep):
        tags.update(ts)
    tags = sorted(tags - set([""]))

    dummies = np.empty((len(arr), len(tags)), dtype=np.int64)

    for i, t in enumerate(tags):
        pat = sep + t + sep
        dummies[:, i] = lib.map_infer(arr.values, lambda x: pat in x)
    return dummies, tags


def str_join(arr, sep):
    """
    Join lists contained as elements in the Series/Index with passed delimiter.

    If the elements of a Series are lists themselves, join the content of these
    lists using the delimiter passed to the function.
    This function is an equivalent to :meth:`str.join`.

    Parameters
    ----------
    sep : str
        Delimiter to use between list entries.

    Returns
    -------
    Series/Index: object

    Notes
    -----
    If any of the lists does not contain string objects the result of the join
    will be `NaN`.

    See Also
    --------
    str.join : Standard library version of this method.
    Series.str.split : Split strings around given separator/delimiter.

    Examples
    --------

    Example with a list that contains non-string elements.

    >>> s = pd.Series([['lion', 'elephant', 'zebra'],
    ...                [1.1, 2.2, 3.3],
    ...                ['cat', np.nan, 'dog'],
    ...                ['cow', 4.5, 'goat']
    ...                ['duck', ['swan', 'fish'], 'guppy']])
    >>> s
    0        [lion, elephant, zebra]
    1                [1.1, 2.2, 3.3]
    2                [cat, nan, dog]
    3               [cow, 4.5, goat]
    4    [duck, [swan, fish], guppy]
    dtype: object

    Join all lists using an '-', the lists containing object(s) of types other
    than str will become a NaN.

    >>> s.str.join('-')
    0    lion-elephant-zebra
    1                    NaN
    2                    NaN
    3                    NaN
    4                    NaN
    dtype: object
    """
    return _na_map(sep.join, arr)


def str_findall(arr, pat, flags=0):
    """
    Find all occurrences of pattern or regular expression in the Series/Index.

    Equivalent to applying :func:`re.findall` to all the elements in the
    Series/Index.

    Parameters
    ----------
    pat : string
        Pattern or regular expression.
    flags : int, default 0
        ``re`` module flags, e.g. `re.IGNORECASE` (default is 0, which means
        no flags).

    Returns
    -------
    Series/Index of lists of strings
        All non-overlapping matches of pattern or regular expression in each
        string of this Series/Index.

    See Also
    --------
    count : Count occurrences of pattern or regular expression in each string
        of the Series/Index.
    extractall : For each string in the Series, extract groups from all matches
        of regular expression and return a DataFrame with one row for each
        match and one column for each group.
    re.findall : The equivalent ``re`` function to all non-overlapping matches
        of pattern or regular expression in string, as a list of strings.

    Examples
    --------

    >>> s = pd.Series(['Lion', 'Monkey', 'Rabbit'])

    The search for the pattern 'Monkey' returns one match:

    >>> s.str.findall('Monkey')
    0          []
    1    [Monkey]
    2          []
    dtype: object

    On the other hand, the search for the pattern 'MONKEY' doesn't return any
    match:

    >>> s.str.findall('MONKEY')
    0    []
    1    []
    2    []
    dtype: object

    Flags can be added to the pattern or regular expression. For instance,
    to find the pattern 'MONKEY' ignoring the case:

    >>> import re
    >>> s.str.findall('MONKEY', flags=re.IGNORECASE)
    0          []
    1    [Monkey]
    2          []
    dtype: object

    When the pattern matches more than one string in the Series, all matches
    are returned:

    >>> s.str.findall('on')
    0    [on]
    1    [on]
    2      []
    dtype: object

    Regular expressions are supported too. For instance, the search for all the
    strings ending with the word 'on' is shown next:

    >>> s.str.findall('on$')
    0    [on]
    1      []
    2      []
    dtype: object

    If the pattern is found more than once in the same string, then a list of
    multiple strings is returned:

    >>> s.str.findall('b')
    0        []
    1        []
    2    [b, b]
    dtype: object

    """
    regex = re.compile(pat, flags=flags)
    return _na_map(regex.findall, arr)


def str_find(arr, sub, start=0, end=None, side='left'):
    """
    Return indexes in each strings in the Series/Index where the
    substring is fully contained between [start:end]. Return -1 on failure.

    Parameters
    ----------
    sub : str
        Substring being searched
    start : int
        Left edge index
    end : int
        Right edge index
    side : {'left', 'right'}, default 'left'
        Specifies a starting side, equivalent to ``find`` or ``rfind``

    Returns
    -------
    found : Series/Index of integer values
    """

    if not isinstance(sub, compat.string_types):
        msg = 'expected a string object, not {0}'
        raise TypeError(msg.format(type(sub).__name__))

    if side == 'left':
        method = 'find'
    elif side == 'right':
        method = 'rfind'
    else:  # pragma: no cover
        raise ValueError('Invalid side')

    if end is None:
        f = lambda x: getattr(x, method)(sub, start)
    else:
        f = lambda x: getattr(x, method)(sub, start, end)

    return _na_map(f, arr, dtype=int)


def str_index(arr, sub, start=0, end=None, side='left'):
    if not isinstance(sub, compat.string_types):
        msg = 'expected a string object, not {0}'
        raise TypeError(msg.format(type(sub).__name__))

    if side == 'left':
        method = 'index'
    elif side == 'right':
        method = 'rindex'
    else:  # pragma: no cover
        raise ValueError('Invalid side')

    if end is None:
        f = lambda x: getattr(x, method)(sub, start)
    else:
        f = lambda x: getattr(x, method)(sub, start, end)

    return _na_map(f, arr, dtype=int)


def str_pad(arr, width, side='left', fillchar=' '):
    """
    Pad strings in the Series/Index with an additional character to
    specified side.

    Parameters
    ----------
    width : int
        Minimum width of resulting string; additional characters will be filled
        with spaces
    side : {'left', 'right', 'both'}, default 'left'
    fillchar : str
        Additional character for filling, default is whitespace

    Returns
    -------
    padded : Series/Index of objects
    """

    if not isinstance(fillchar, compat.string_types):
        msg = 'fillchar must be a character, not {0}'
        raise TypeError(msg.format(type(fillchar).__name__))

    if len(fillchar) != 1:
        raise TypeError('fillchar must be a character, not str')

    if not is_integer(width):
        msg = 'width must be of integer type, not {0}'
        raise TypeError(msg.format(type(width).__name__))

    if side == 'left':
        f = lambda x: x.rjust(width, fillchar)
    elif side == 'right':
        f = lambda x: x.ljust(width, fillchar)
    elif side == 'both':
        f = lambda x: x.center(width, fillchar)
    else:  # pragma: no cover
        raise ValueError('Invalid side')

    return _na_map(f, arr)


def str_split(arr, pat=None, n=None):
    """
    Split strings around given separator/delimiter.

    Split each string in the caller's values by given
    pattern, propagating NaN values. Equivalent to :meth:`str.split`.

    Parameters
    ----------
    pat : str, optional
        String or regular expression to split on.
        If not specified, split on whitespace.
    n : int, default -1 (all)
        Limit number of splits in output.
        ``None``, 0 and -1 will be interpreted as return all splits.
    expand : bool, default False
        Expand the splitted strings into separate columns.

        * If ``True``, return DataFrame/MultiIndex expanding dimensionality.
        * If ``False``, return Series/Index, containing lists of strings.

    Returns
    -------
    Series, Index, DataFrame or MultiIndex
        Type matches caller unless ``expand=True`` (see Notes).

    Notes
    -----
    The handling of the `n` keyword depends on the number of found splits:

    - If found splits > `n`,  make first `n` splits only
    - If found splits <= `n`, make all splits
    - If for a certain row the number of found splits < `n`,
      append `None` for padding up to `n` if ``expand=True``

    If using ``expand=True``, Series and Index callers return DataFrame and
    MultiIndex objects, respectively.

    See Also
    --------
    str.split : Standard library version of this method.
    Series.str.get_dummies : Split each string into dummy variables.
    Series.str.partition : Split string on a separator, returning
        the before, separator, and after components.

    Examples
    --------
    >>> s = pd.Series(["this is good text", "but this is even better"])

    By default, split will return an object of the same size
    having lists containing the split elements

    >>> s.str.split()
    0           [this, is, good, text]
    1    [but, this, is, even, better]
    dtype: object
    >>> s.str.split("random")
    0          [this is good text]
    1    [but this is even better]
    dtype: object

    When using ``expand=True``, the split elements will expand out into
    separate columns.

    For Series object, output return type is DataFrame.

    >>> s.str.split(expand=True)
          0     1     2     3       4
    0  this    is  good  text    None
    1   but  this    is  even  better
    >>> s.str.split(" is ", expand=True)
              0            1
    0      this    good text
    1  but this  even better

    For Index object, output return type is MultiIndex.

    >>> i = pd.Index(["ba 100 001", "ba 101 002", "ba 102 003"])
    >>> i.str.split(expand=True)
    MultiIndex(levels=[['ba'], ['100', '101', '102'], ['001', '002', '003']],
           labels=[[0, 0, 0], [0, 1, 2], [0, 1, 2]])

    Parameter `n` can be used to limit the number of splits in the output.

    >>> s.str.split("is", n=1)
    0          [th,  is good text]
    1    [but th,  is even better]
    dtype: object
    >>> s.str.split("is", n=1, expand=True)
            0                1
    0      th     is good text
    1  but th   is even better

    If NaN is present, it is propagated throughout the columns
    during the split.

    >>> s = pd.Series(["this is good text", "but this is even better", np.nan])
    >>> s.str.split(n=3, expand=True)
          0     1     2            3
    0  this    is  good         text
    1   but  this    is  even better
    2   NaN   NaN   NaN          NaN
    """
    if pat is None:
        if n is None or n == 0:
            n = -1
        f = lambda x: x.split(pat, n)
    else:
        if len(pat) == 1:
            if n is None or n == 0:
                n = -1
            f = lambda x: x.split(pat, n)
        else:
            if n is None or n == -1:
                n = 0
            regex = re.compile(pat)
            f = lambda x: regex.split(x, maxsplit=n)
    res = _na_map(f, arr)
    return res


def str_rsplit(arr, pat=None, n=None):
    """
    Split each string in the Series/Index by the given delimiter
    string, starting at the end of the string and working to the front.
    Equivalent to :meth:`str.rsplit`.

    Parameters
    ----------
    pat : string, default None
        Separator to split on. If None, splits on whitespace
    n : int, default -1 (all)
        None, 0 and -1 will be interpreted as return all splits
    expand : bool, default False
        * If True, return DataFrame/MultiIndex expanding dimensionality.
        * If False, return Series/Index.

    Returns
    -------
    split : Series/Index or DataFrame/MultiIndex of objects
    """
    if n is None or n == 0:
        n = -1
    f = lambda x: x.rsplit(pat, n)
    res = _na_map(f, arr)
    return res


def str_slice(arr, start=None, stop=None, step=None):
    """
    Slice substrings from each element in the Series/Index

    Parameters
    ----------
    start : int or None
    stop : int or None
    step : int or None

    Returns
    -------
    sliced : Series/Index of objects
    """
    obj = slice(start, stop, step)
    f = lambda x: x[obj]
    return _na_map(f, arr)


def str_slice_replace(arr, start=None, stop=None, repl=None):
    """
    Replace a positional slice of a string with another value.

    Parameters
    ----------
    start : int, optional
        Left index position to use for the slice. If not specified (None),
        the slice is unbounded on the left, i.e. slice from the start
        of the string.
    stop : int, optional
        Right index position to use for the slice. If not specified (None),
        the slice is unbounded on the right, i.e. slice until the
        end of the string.
    repl : str, optional
        String for replacement. If not specified (None), the sliced region
        is replaced with an empty string.

    Returns
    -------
    replaced : Series or Index
        Same type as the original object.

    See Also
    --------
    Series.str.slice : Just slicing without replacement.

    Examples
    --------
    >>> s = pd.Series(['a', 'ab', 'abc', 'abdc', 'abcde'])
    >>> s
    0        a
    1       ab
    2      abc
    3     abdc
    4    abcde
    dtype: object

    Specify just `start`, meaning replace `start` until the end of the
    string with `repl`.

    >>> s.str.slice_replace(1, repl='X')
    0    aX
    1    aX
    2    aX
    3    aX
    4    aX
    dtype: object

    Specify just `stop`, meaning the start of the string to `stop` is replaced
    with `repl`, and the rest of the string is included.

    >>> s.str.slice_replace(stop=2, repl='X')
    0       X
    1       X
    2      Xc
    3     Xdc
    4    Xcde
    dtype: object

    Specify `start` and `stop`, meaning the slice from `start` to `stop` is
    replaced with `repl`. Everything before or after `start` and `stop` is
    included as is.

    >>> s.str.slice_replace(start=1, stop=3, repl='X')
    0      aX
    1      aX
    2      aX
    3     aXc
    4    aXde
    dtype: object
    """
    if repl is None:
        repl = ''

    def f(x):
        if x[start:stop] == '':
            local_stop = start
        else:
            local_stop = stop
        y = ''
        if start is not None:
            y += x[:start]
        y += repl
        if stop is not None:
            y += x[local_stop:]
        return y

    return _na_map(f, arr)


def str_strip(arr, to_strip=None, side='both'):
    """
    Strip whitespace (including newlines) from each string in the
    Series/Index.

    Parameters
    ----------
    to_strip : str or unicode
    side : {'left', 'right', 'both'}, default 'both'

    Returns
    -------
    stripped : Series/Index of objects
    """
    if side == 'both':
        f = lambda x: x.strip(to_strip)
    elif side == 'left':
        f = lambda x: x.lstrip(to_strip)
    elif side == 'right':
        f = lambda x: x.rstrip(to_strip)
    else:  # pragma: no cover
        raise ValueError('Invalid side')
    return _na_map(f, arr)


def str_wrap(arr, width, **kwargs):
    r"""
    Wrap long strings in the Series/Index to be formatted in
    paragraphs with length less than a given width.

    This method has the same keyword parameters and defaults as
    :class:`textwrap.TextWrapper`.

    Parameters
    ----------
    width : int
        Maximum line-width
    expand_tabs : bool, optional
        If true, tab characters will be expanded to spaces (default: True)
    replace_whitespace : bool, optional
        If true, each whitespace character (as defined by string.whitespace)
        remaining after tab expansion will be replaced by a single space
        (default: True)
    drop_whitespace : bool, optional
        If true, whitespace that, after wrapping, happens to end up at the
        beginning or end of a line is dropped (default: True)
    break_long_words : bool, optional
        If true, then words longer than width will be broken in order to ensure
        that no lines are longer than width. If it is false, long words will
        not be broken, and some lines may be longer than width. (default: True)
    break_on_hyphens : bool, optional
        If true, wrapping will occur preferably on whitespace and right after
        hyphens in compound words, as it is customary in English. If false,
        only whitespaces will be considered as potentially good places for line
        breaks, but you need to set break_long_words to false if you want truly
        insecable words. (default: True)

    Returns
    -------
    wrapped : Series/Index of objects

    Notes
    -----
    Internally, this method uses a :class:`textwrap.TextWrapper` instance with
    default settings. To achieve behavior matching R's stringr library str_wrap
    function, use the arguments:

    - expand_tabs = False
    - replace_whitespace = True
    - drop_whitespace = True
    - break_long_words = False
    - break_on_hyphens = False

    Examples
    --------

    >>> s = pd.Series(['line to be wrapped', 'another line to be wrapped'])
    >>> s.str.wrap(12)
    0             line to be\nwrapped
    1    another line\nto be\nwrapped
    """
    kwargs['width'] = width

    tw = textwrap.TextWrapper(**kwargs)

    return _na_map(lambda s: '\n'.join(tw.wrap(s)), arr)


def str_translate(arr, table, deletechars=None):
    """
    Map all characters in the string through the given mapping table.
    Equivalent to standard :meth:`str.translate`. Note that the optional
    argument deletechars is only valid if you are using python 2. For python 3,
    character deletion should be specified via the table argument.

    Parameters
    ----------
    table : dict (python 3), str or None (python 2)
        In python 3, table is a mapping of Unicode ordinals to Unicode
        ordinals, strings, or None. Unmapped characters are left untouched.
        Characters mapped to None are deleted. :meth:`str.maketrans` is a
        helper function for making translation tables.
        In python 2, table is either a string of length 256 or None. If the
        table argument is None, no translation is applied and the operation
        simply removes the characters in deletechars. :func:`string.maketrans`
        is a helper function for making translation tables.
    deletechars : str, optional (python 2)
        A string of characters to delete. This argument is only valid
        in python 2.

    Returns
    -------
    translated : Series/Index of objects
    """
    if deletechars is None:
        f = lambda x: x.translate(table)
    else:
        if compat.PY3:
            raise ValueError("deletechars is not a valid argument for "
                             "str.translate in python 3. You should simply "
                             "specify character deletions in the table "
                             "argument")
        f = lambda x: x.translate(table, deletechars)
    return _na_map(f, arr)


def str_get(arr, i):
    """
    Extract element from each component at specified position.

    Extract element from lists, tuples, or strings in each element in the
    Series/Index.

    Parameters
    ----------
    i : int
        Position of element to extract.

    Returns
    -------
    items : Series/Index of objects

    Examples
    --------
    >>> s = pd.Series(["String",
               (1, 2, 3),
               ["a", "b", "c"],
               123, -456,
               {1:"Hello", "2":"World"}])
    >>> s
    0                        String
    1                     (1, 2, 3)
    2                     [a, b, c]
    3                           123
    4                          -456
    5    {1: 'Hello', '2': 'World'}
    dtype: object

    >>> s.str.get(1)
    0        t
    1        2
    2        b
    3      NaN
    4      NaN
    5    Hello
    dtype: object

    >>> s.str.get(-1)
    0      g
    1      3
    2      c
    3    NaN
    4    NaN
    5    NaN
    dtype: object
    """
    def f(x):
        if isinstance(x, dict):
            return x.get(i)
        elif len(x) > i >= -len(x):
            return x[i]
        return np.nan
    return _na_map(f, arr)


def str_decode(arr, encoding, errors="strict"):
    """
    Decode character string in the Series/Index using indicated encoding.
    Equivalent to :meth:`str.decode` in python2 and :meth:`bytes.decode` in
    python3.

    Parameters
    ----------
    encoding : str
    errors : str, optional

    Returns
    -------
    decoded : Series/Index of objects
    """
    if encoding in _cpython_optimized_decoders:
        # CPython optimized implementation
        f = lambda x: x.decode(encoding, errors)
    else:
        decoder = codecs.getdecoder(encoding)
        f = lambda x: decoder(x, errors)[0]
    return _na_map(f, arr)


def str_encode(arr, encoding, errors="strict"):
    """
    Encode character string in the Series/Index using indicated encoding.
    Equivalent to :meth:`str.encode`.

    Parameters
    ----------
    encoding : str
    errors : str, optional

    Returns
    -------
    encoded : Series/Index of objects
    """
    if encoding in _cpython_optimized_encoders:
        # CPython optimized implementation
        f = lambda x: x.encode(encoding, errors)
    else:
        encoder = codecs.getencoder(encoding)
        f = lambda x: encoder(x, errors)[0]
    return _na_map(f, arr)


def _noarg_wrapper(f, docstring=None, **kargs):
    def wrapper(self):
        result = _na_map(f, self._data, **kargs)
        return self._wrap_result(result)

    wrapper.__name__ = f.__name__
    if docstring is not None:
        wrapper.__doc__ = docstring
    else:
        raise ValueError('Provide docstring')

    return wrapper


def _pat_wrapper(f, flags=False, na=False, **kwargs):
    def wrapper1(self, pat):
        result = f(self._data, pat)
        return self._wrap_result(result)

    def wrapper2(self, pat, flags=0, **kwargs):
        result = f(self._data, pat, flags=flags, **kwargs)
        return self._wrap_result(result)

    def wrapper3(self, pat, na=np.nan):
        result = f(self._data, pat, na=na)
        return self._wrap_result(result)

    wrapper = wrapper3 if na else wrapper2 if flags else wrapper1

    wrapper.__name__ = f.__name__
    if f.__doc__:
        wrapper.__doc__ = f.__doc__

    return wrapper


def copy(source):
    "Copy a docstring from another source function (if present)"

    def do_copy(target):
        if source.__doc__:
            target.__doc__ = source.__doc__
        return target

    return do_copy


class StringMethods(NoNewAttributesMixin):
    """
    Vectorized string functions for Series and Index. NAs stay NA unless
    handled otherwise by a particular method. Patterned after Python's string
    methods, with some inspiration from R's stringr package.

    Examples
    --------
    >>> s.str.split('_')
    >>> s.str.replace('_', '')
    """

    def __init__(self, data):
        self._validate(data)
        self._is_categorical = is_categorical_dtype(data)

        # .values.categories works for both Series/Index
        self._data = data.values.categories if self._is_categorical else data
        # save orig to blow up categoricals to the right type
        self._orig = data
        self._freeze()

    @staticmethod
    def _validate(data):
        from pandas.core.index import Index

        if (isinstance(data, ABCSeries) and
                not ((is_categorical_dtype(data.dtype) and
                      is_object_dtype(data.values.categories)) or
                     (is_object_dtype(data.dtype)))):
            # it's neither a string series not a categorical series with
            # strings inside the categories.
            # this really should exclude all series with any non-string values
            # (instead of test for object dtype), but that isn't practical for
            # performance reasons until we have a str dtype (GH 9343)
            raise AttributeError("Can only use .str accessor with string "
                                 "values, which use np.object_ dtype in "
                                 "pandas")
        elif isinstance(data, Index):
            # can't use ABCIndex to exclude non-str

            # see src/inference.pyx which can contain string values
            allowed_types = ('string', 'unicode', 'mixed', 'mixed-integer')
            if is_categorical_dtype(data.dtype):
                inf_type = data.categories.inferred_type
            else:
                inf_type = data.inferred_type
            if inf_type not in allowed_types:
                message = ("Can only use .str accessor with string values "
                           "(i.e. inferred_type is 'string', 'unicode' or "
                           "'mixed')")
                raise AttributeError(message)
            if data.nlevels > 1:
                message = ("Can only use .str accessor with Index, not "
                           "MultiIndex")
                raise AttributeError(message)

    def __getitem__(self, key):
        if isinstance(key, slice):
            return self.slice(start=key.start, stop=key.stop, step=key.step)
        else:
            return self.get(key)

    def __iter__(self):
        i = 0
        g = self.get(i)
        while g.notna().any():
            yield g
            i += 1
            g = self.get(i)

    def _wrap_result(self, result, use_codes=True,
                     name=None, expand=None):

        from pandas.core.index import Index, MultiIndex

        # for category, we do the stuff on the categories, so blow it up
        # to the full series again
        # But for some operations, we have to do the stuff on the full values,
        # so make it possible to skip this step as the method already did this
        # before the transformation...
        if use_codes and self._is_categorical:
            result = take_1d(result, self._orig.cat.codes)

        if not hasattr(result, 'ndim') or not hasattr(result, 'dtype'):
            return result
        assert result.ndim < 3

        if expand is None:
            # infer from ndim if expand is not specified
            expand = False if result.ndim == 1 else True

        elif expand is True and not isinstance(self._orig, Index):
            # required when expand=True is explicitly specified
            # not needed when inferred

            def cons_row(x):
                if is_list_like(x):
                    return x
                else:
                    return [x]

            result = [cons_row(x) for x in result]
            if result:
                # propagate nan values to match longest sequence (GH 18450)
                max_len = max(len(x) for x in result)
                result = [x * max_len if len(x) == 0 or x[0] is np.nan
                          else x for x in result]

        if not isinstance(expand, bool):
            raise ValueError("expand must be True or False")

        if expand is False:
            # if expand is False, result should have the same name
            # as the original otherwise specified
            if name is None:
                name = getattr(result, 'name', None)
            if name is None:
                # do not use logical or, _orig may be a DataFrame
                # which has "name" column
                name = self._orig.name

        # Wait until we are sure result is a Series or Index before
        # checking attributes (GH 12180)
        if isinstance(self._orig, Index):
            # if result is a boolean np.array, return the np.array
            # instead of wrapping it into a boolean Index (GH 8875)
            if is_bool_dtype(result):
                return result

            if expand:
                result = list(result)
                out = MultiIndex.from_tuples(result, names=name)
                if out.nlevels == 1:
                    # We had all tuples of length-one, which are
                    # better represented as a regular Index.
                    out = out.get_level_values(0)
                return out
            else:
                return Index(result, name=name)
        else:
            index = self._orig.index
            if expand:
                cons = self._orig._constructor_expanddim
                return cons(result, columns=name, index=index)
            else:
                # Must be a Series
                cons = self._orig._constructor
                return cons(result, name=name, index=index)

    def _get_series_list(self, others, ignore_index=False):
        """
        Auxiliary function for :meth:`str.cat`. Turn potentially mixed input
        into a list of Series (elements without an index must match the length
        of the calling Series/Index).

        Parameters
        ----------
        others : Series, DataFrame, np.ndarray, list-like or list-like of
            objects that are either Series, np.ndarray (1-dim) or list-like
        ignore_index : boolean, default False
            Determines whether to forcefully align others with index of caller

        Returns
        -------
        tuple : (others transformed into list of Series,
                 boolean whether FutureWarning should be raised)
        """

        # once str.cat defaults to alignment, this function can be simplified;
        # will not need `ignore_index` and the second boolean output anymore

        from pandas import Index, Series, DataFrame

        # self._orig is either Series or Index
        idx = self._orig if isinstance(self._orig, Index) else self._orig.index

        err_msg = ('others must be Series, Index, DataFrame, np.ndarrary or '
                   'list-like (either containing only strings or containing '
                   'only objects of type Series/Index/list-like/np.ndarray)')

        # Generally speaking, all objects without an index inherit the index
        # `idx` of the calling Series/Index - i.e. must have matching length.
        # Objects with an index (i.e. Series/Index/DataFrame) keep their own
        # index, *unless* ignore_index is set to True.
        if isinstance(others, Series):
            warn = not others.index.equals(idx)
            # only reconstruct Series when absolutely necessary
            los = [Series(others.values, index=idx)
                   if ignore_index and warn else others]
            return (los, warn)
        elif isinstance(others, Index):
            warn = not others.equals(idx)
            los = [Series(others.values,
                          index=(idx if ignore_index else others))]
            return (los, warn)
        elif isinstance(others, DataFrame):
            warn = not others.index.equals(idx)
            if ignore_index and warn:
                # without copy, this could change "others"
                # that was passed to str.cat
                others = others.copy()
                others.index = idx
            return ([others[x] for x in others], warn)
        elif isinstance(others, np.ndarray) and others.ndim == 2:
            others = DataFrame(others, index=idx)
            return ([others[x] for x in others], False)
        elif is_list_like(others):
            others = list(others)  # ensure iterators do not get read twice etc

            # in case of list-like `others`, all elements must be
            # either one-dimensional list-likes or scalars
            if all(is_list_like(x) for x in others):
                los = []
                warn = False
                # iterate through list and append list of series for each
                # element (which we check to be one-dimensional and non-nested)
                while others:
                    nxt = others.pop(0)  # nxt is guaranteed list-like by above
                    if not isinstance(nxt, (DataFrame, Series,
                                            Index, np.ndarray)):
                        # safety for non-persistent list-likes (e.g. iterators)
                        # do not map indexed/typed objects; info needed below
                        nxt = list(nxt)

                    # known types for which we can avoid deep inspection
                    no_deep = ((isinstance(nxt, np.ndarray) and nxt.ndim == 1)
                               or isinstance(nxt, (Series, Index)))
                    # nested list-likes are forbidden:
                    # -> elements of nxt must not be list-like
                    is_legal = ((no_deep and nxt.dtype == object)
                                or all(not is_list_like(x) for x in nxt))

                    # DataFrame is false positive of is_legal
                    # because "x in df" returns column names
                    if not is_legal or isinstance(nxt, DataFrame):
                        raise TypeError(err_msg)

                    nxt, wnx = self._get_series_list(nxt,
                                                     ignore_index=ignore_index)
                    los = los + nxt
                    warn = warn or wnx
                return (los, warn)
            elif all(not is_list_like(x) for x in others):
                return ([Series(others, index=idx)], False)
        raise TypeError(err_msg)

    def cat(self, others=None, sep=None, na_rep=None, join=None):
        """
        Concatenate strings in the Series/Index with given separator.

        If `others` is specified, this function concatenates the Series/Index
        and elements of `others` element-wise.
        If `others` is not passed, then all values in the Series/Index are
        concatenated into a single string with a given `sep`.

        Parameters
        ----------
        others : Series, Index, DataFrame, np.ndarrary or list-like
            Series, Index, DataFrame, np.ndarray (one- or two-dimensional) and
            other list-likes of strings must have the same length as the
            calling Series/Index, with the exception of indexed objects (i.e.
            Series/Index/DataFrame) if `join` is not None.

            If others is a list-like that contains a combination of Series,
            np.ndarray (1-dim) or list-like, then all elements will be unpacked
            and must satisfy the above criteria individually.

            If others is None, the method returns the concatenation of all
            strings in the calling Series/Index.
        sep : string or None, default None
            If None, concatenates without any separator.
        na_rep : string or None, default None
            Representation that is inserted for all missing values:

            - If `na_rep` is None, and `others` is None, missing values in the
              Series/Index are omitted from the result.
            - If `na_rep` is None, and `others` is not None, a row containing a
              missing value in any of the columns (before concatenation) will
              have a missing value in the result.
        join : {'left', 'right', 'outer', 'inner'}, default None
            Determines the join-style between the calling Series/Index and any
            Series/Index/DataFrame in `others` (objects without an index need
            to match the length of the calling Series/Index). If None,
            alignment is disabled, but this option will be removed in a future
            version of pandas and replaced with a default of `'left'`. To
            disable alignment, use `.values` on any Series/Index/DataFrame in
            `others`.

            .. versionadded:: 0.23.0

        Returns
        -------
        concat : str if `other is None`, Series/Index of objects if `others is
            not None`. In the latter case, the result will remain categorical
            if the calling Series/Index is categorical.

        See Also
        --------
        split : Split each string in the Series/Index

        Examples
        --------
        When not passing `others`, all values are concatenated into a single
        string:

        >>> s = pd.Series(['a', 'b', np.nan, 'd'])
        >>> s.str.cat(sep=' ')
        'a b d'

        By default, NA values in the Series are ignored. Using `na_rep`, they
        can be given a representation:

        >>> s.str.cat(sep=' ', na_rep='?')
        'a b ? d'

        If `others` is specified, corresponding values are concatenated with
        the separator. Result will be a Series of strings.

        >>> s.str.cat(['A', 'B', 'C', 'D'], sep=',')
        0    a,A
        1    b,B
        2    NaN
        3    d,D
        dtype: object

        Missing values will remain missing in the result, but can again be
        represented using `na_rep`

        >>> s.str.cat(['A', 'B', 'C', 'D'], sep=',', na_rep='-')
        0    a,A
        1    b,B
        2    -,C
        3    d,D
        dtype: object

        If `sep` is not specified, the values are concatenated without
        separation.

        >>> s.str.cat(['A', 'B', 'C', 'D'], na_rep='-')
        0    aA
        1    bB
        2    -C
        3    dD
        dtype: object

        Series with different indexes can be aligned before concatenation. The
        `join`-keyword works as in other methods.

        >>> t = pd.Series(['d', 'a', 'e', 'c'], index=[3, 0, 4, 2])
        >>> s.str.cat(t, join=None, na_rep='-')
        0    ad
        1    ba
        2    -e
        3    dc
        dtype: object
        >>>
        >>> s.str.cat(t, join='left', na_rep='-')
        0    aa
        1    b-
        2    -c
        3    dd
        dtype: object
        >>>
        >>> s.str.cat(t, join='outer', na_rep='-')
        0    aa
        1    b-
        2    -c
        3    dd
        4    -e
        dtype: object
        >>>
        >>> s.str.cat(t, join='inner', na_rep='-')
        0    aa
        2    -c
        3    dd
        dtype: object
        >>>
        >>> s.str.cat(t, join='right', na_rep='-')
        3    dd
        0    aa
        4    -e
        2    -c
        dtype: object

        For more examples, see :ref:`here <text.concatenate>`.
        """
        from pandas import Index, Series, concat

        if isinstance(others, compat.string_types):
            raise ValueError("Did you mean to supply a `sep` keyword?")

        if isinstance(self._orig, Index):
            data = Series(self._orig, index=self._orig)
        else:  # Series
            data = self._orig

        # concatenate Series/Index with itself if no "others"
        if others is None:
            result = str_cat(data, others=others, sep=sep, na_rep=na_rep)
            return self._wrap_result(result,
                                     use_codes=(not self._is_categorical))

        try:
            # turn anything in "others" into lists of Series
            others, warn = self._get_series_list(others,
                                                 ignore_index=(join is None))
        except ValueError:  # do not catch TypeError raised by _get_series_list
            if join is None:
                raise ValueError('All arrays must be same length, except '
                                 'those having an index if `join` is not None')
            else:
                raise ValueError('If `others` contains arrays or lists (or '
                                 'other list-likes without an index), these '
                                 'must all be of the same length as the '
                                 'calling Series/Index.')

        if join is None and warn:
            warnings.warn("A future version of pandas will perform index "
                          "alignment when `others` is a Series/Index/"
                          "DataFrame (or a list-like containing one). To "
                          "disable alignment (the behavior before v.0.23) and "
                          "silence this warning, use `.values` on any Series/"
                          "Index/DataFrame in `others`. To enable alignment "
                          "and silence this warning, pass `join='left'|"
                          "'outer'|'inner'|'right'`. The future default will "
                          "be `join='left'`.", FutureWarning, stacklevel=2)

        # align if required
        if join is not None:
            # Need to add keys for uniqueness in case of duplicate columns
            others = concat(others, axis=1,
                            join=(join if join == 'inner' else 'outer'),
                            keys=range(len(others)))
            data, others = data.align(others, join=join)
            others = [others[x] for x in others]  # again list of Series

        # str_cat discards index
        res = str_cat(data, others=others, sep=sep, na_rep=na_rep)

        if isinstance(self._orig, Index):
            res = Index(res, name=self._orig.name)
        else:  # Series
            res = Series(res, index=data.index, name=self._orig.name)
        return res

    @copy(str_split)
    def split(self, pat=None, n=-1, expand=False):
        result = str_split(self._data, pat, n=n)
        return self._wrap_result(result, expand=expand)

    @copy(str_rsplit)
    def rsplit(self, pat=None, n=-1, expand=False):
        result = str_rsplit(self._data, pat, n=n)
        return self._wrap_result(result, expand=expand)

    _shared_docs['str_partition'] = ("""
    Split the string at the %(side)s occurrence of `sep`, and return 3 elements
    containing the part before the separator, the separator itself,
    and the part after the separator.
    If the separator is not found, return %(return)s.

    Parameters
    ----------
    pat : string, default whitespace
        String to split on.
    expand : bool, default True
        * If True, return DataFrame/MultiIndex expanding dimensionality.
        * If False, return Series/Index.

    Returns
    -------
    split : DataFrame/MultiIndex or Series/Index of objects

    See Also
    --------
    %(also)s

    Examples
    --------

    >>> s = Series(['A_B_C', 'D_E_F', 'X'])
    0    A_B_C
    1    D_E_F
    2        X
    dtype: object

    >>> s.str.partition('_')
       0  1    2
    0  A  _  B_C
    1  D  _  E_F
    2  X

    >>> s.str.rpartition('_')
         0  1  2
    0  A_B  _  C
    1  D_E  _  F
    2          X
    """)

    @Appender(_shared_docs['str_partition'] % {
        'side': 'first',
        'return': '3 elements containing the string itself, followed by two '
                  'empty strings',
        'also': 'rpartition : Split the string at the last occurrence of `sep`'
    })
    def partition(self, pat=' ', expand=True):
        f = lambda x: x.partition(pat)
        result = _na_map(f, self._data)
        return self._wrap_result(result, expand=expand)

    @Appender(_shared_docs['str_partition'] % {
        'side': 'last',
        'return': '3 elements containing two empty strings, followed by the '
                  'string itself',
        'also': 'partition : Split the string at the first occurrence of `sep`'
    })
    def rpartition(self, pat=' ', expand=True):
        f = lambda x: x.rpartition(pat)
        result = _na_map(f, self._data)
        return self._wrap_result(result, expand=expand)

    @copy(str_get)
    def get(self, i):
        result = str_get(self._data, i)
        return self._wrap_result(result)

    @copy(str_join)
    def join(self, sep):
        result = str_join(self._data, sep)
        return self._wrap_result(result)

    @copy(str_contains)
    def contains(self, pat, case=True, flags=0, na=np.nan, regex=True):
        result = str_contains(self._data, pat, case=case, flags=flags, na=na,
                              regex=regex)
        return self._wrap_result(result)

    @copy(str_match)
    def match(self, pat, case=True, flags=0, na=np.nan, as_indexer=None):
        result = str_match(self._data, pat, case=case, flags=flags, na=na,
                           as_indexer=as_indexer)
        return self._wrap_result(result)

    @copy(str_replace)
    def replace(self, pat, repl, n=-1, case=None, flags=0, regex=True):
        result = str_replace(self._data, pat, repl, n=n, case=case,
                             flags=flags, regex=regex)
        return self._wrap_result(result)

    @copy(str_repeat)
    def repeat(self, repeats):
        result = str_repeat(self._data, repeats)
        return self._wrap_result(result)

    @copy(str_pad)
    def pad(self, width, side='left', fillchar=' '):
        result = str_pad(self._data, width, side=side, fillchar=fillchar)
        return self._wrap_result(result)

    _shared_docs['str_pad'] = ("""
    Filling %(side)s side of strings in the Series/Index with an
    additional character. Equivalent to :meth:`str.%(method)s`.

    Parameters
    ----------
    width : int
        Minimum width of resulting string; additional characters will be filled
        with ``fillchar``
    fillchar : str
        Additional character for filling, default is whitespace

    Returns
    -------
    filled : Series/Index of objects
    """)

    @Appender(_shared_docs['str_pad'] % dict(side='left and right',
                                             method='center'))
    def center(self, width, fillchar=' '):
        return self.pad(width, side='both', fillchar=fillchar)

    @Appender(_shared_docs['str_pad'] % dict(side='right', method='ljust'))
    def ljust(self, width, fillchar=' '):
        return self.pad(width, side='right', fillchar=fillchar)

    @Appender(_shared_docs['str_pad'] % dict(side='left', method='rjust'))
    def rjust(self, width, fillchar=' '):
        return self.pad(width, side='left', fillchar=fillchar)

    def zfill(self, width):
        """
        Filling left side of strings in the Series/Index with 0.
        Equivalent to :meth:`str.zfill`.

        Parameters
        ----------
        width : int
            Minimum width of resulting string; additional characters will be
            filled with 0

        Returns
        -------
        filled : Series/Index of objects
        """
        result = str_pad(self._data, width, side='left', fillchar='0')
        return self._wrap_result(result)

    @copy(str_slice)
    def slice(self, start=None, stop=None, step=None):
        result = str_slice(self._data, start, stop, step)
        return self._wrap_result(result)

    @copy(str_slice_replace)
    def slice_replace(self, start=None, stop=None, repl=None):
        result = str_slice_replace(self._data, start, stop, repl)
        return self._wrap_result(result)

    @copy(str_decode)
    def decode(self, encoding, errors="strict"):
        result = str_decode(self._data, encoding, errors)
        return self._wrap_result(result)

    @copy(str_encode)
    def encode(self, encoding, errors="strict"):
        result = str_encode(self._data, encoding, errors)
        return self._wrap_result(result)

    _shared_docs['str_strip'] = ("""
    Strip whitespace (including newlines) from each string in the
    Series/Index from %(side)s. Equivalent to :meth:`str.%(method)s`.

    Returns
    -------
    stripped : Series/Index of objects
    """)

    @Appender(_shared_docs['str_strip'] % dict(side='left and right sides',
                                               method='strip'))
    def strip(self, to_strip=None):
        result = str_strip(self._data, to_strip, side='both')
        return self._wrap_result(result)

    @Appender(_shared_docs['str_strip'] % dict(side='left side',
                                               method='lstrip'))
    def lstrip(self, to_strip=None):
        result = str_strip(self._data, to_strip, side='left')
        return self._wrap_result(result)

    @Appender(_shared_docs['str_strip'] % dict(side='right side',
                                               method='rstrip'))
    def rstrip(self, to_strip=None):
        result = str_strip(self._data, to_strip, side='right')
        return self._wrap_result(result)

    @copy(str_wrap)
    def wrap(self, width, **kwargs):
        result = str_wrap(self._data, width, **kwargs)
        return self._wrap_result(result)

    @copy(str_get_dummies)
    def get_dummies(self, sep='|'):
        # we need to cast to Series of strings as only that has all
        # methods available for making the dummies...
        data = self._orig.astype(str) if self._is_categorical else self._data
        result, name = str_get_dummies(data, sep)
        return self._wrap_result(result, use_codes=(not self._is_categorical),
                                 name=name, expand=True)

    @copy(str_translate)
    def translate(self, table, deletechars=None):
        result = str_translate(self._data, table, deletechars)
        return self._wrap_result(result)

    count = _pat_wrapper(str_count, flags=True)
    startswith = _pat_wrapper(str_startswith, na=True)
    endswith = _pat_wrapper(str_endswith, na=True)
    findall = _pat_wrapper(str_findall, flags=True)

    @copy(str_extract)
    def extract(self, pat, flags=0, expand=True):
        return str_extract(self, pat, flags=flags, expand=expand)

    @copy(str_extractall)
    def extractall(self, pat, flags=0):
        return str_extractall(self._orig, pat, flags=flags)

    _shared_docs['find'] = ("""
    Return %(side)s indexes in each strings in the Series/Index
    where the substring is fully contained between [start:end].
    Return -1 on failure. Equivalent to standard :meth:`str.%(method)s`.

    Parameters
    ----------
    sub : str
        Substring being searched
    start : int
        Left edge index
    end : int
        Right edge index

    Returns
    -------
    found : Series/Index of integer values

    See Also
    --------
    %(also)s
    """)

    @Appender(_shared_docs['find'] %
              dict(side='lowest', method='find',
                   also='rfind : Return highest indexes in each strings'))
    def find(self, sub, start=0, end=None):
        result = str_find(self._data, sub, start=start, end=end, side='left')
        return self._wrap_result(result)

    @Appender(_shared_docs['find'] %
              dict(side='highest', method='rfind',
                   also='find : Return lowest indexes in each strings'))
    def rfind(self, sub, start=0, end=None):
        result = str_find(self._data, sub, start=start, end=end, side='right')
        return self._wrap_result(result)

    def normalize(self, form):
        """Return the Unicode normal form for the strings in the Series/Index.
        For more information on the forms, see the
        :func:`unicodedata.normalize`.

        Parameters
        ----------
        form : {'NFC', 'NFKC', 'NFD', 'NFKD'}
            Unicode form

        Returns
        -------
        normalized : Series/Index of objects
        """
        import unicodedata
        f = lambda x: unicodedata.normalize(form, compat.u_safe(x))
        result = _na_map(f, self._data)
        return self._wrap_result(result)

    _shared_docs['index'] = ("""
    Return %(side)s indexes in each strings where the substring is
    fully contained between [start:end]. This is the same as
    ``str.%(similar)s`` except instead of returning -1, it raises a ValueError
    when the substring is not found. Equivalent to standard ``str.%(method)s``.

    Parameters
    ----------
    sub : str
        Substring being searched
    start : int
        Left edge index
    end : int
        Right edge index

    Returns
    -------
    found : Series/Index of objects

    See Also
    --------
    %(also)s
    """)

    @Appender(_shared_docs['index'] %
              dict(side='lowest', similar='find', method='index',
                   also='rindex : Return highest indexes in each strings'))
    def index(self, sub, start=0, end=None):
        result = str_index(self._data, sub, start=start, end=end, side='left')
        return self._wrap_result(result)

    @Appender(_shared_docs['index'] %
              dict(side='highest', similar='rfind', method='rindex',
                   also='index : Return lowest indexes in each strings'))
    def rindex(self, sub, start=0, end=None):
        result = str_index(self._data, sub, start=start, end=end, side='right')
        return self._wrap_result(result)

    _shared_docs['len'] = ("""
    Compute length of each string in the Series/Index.

    Returns
    -------
    lengths : Series/Index of integer values
    """)
    len = _noarg_wrapper(len, docstring=_shared_docs['len'], dtype=int)

    _shared_docs['casemethods'] = ("""
    Convert strings in the Series/Index to %(type)s.

    Equivalent to :meth:`str.%(method)s`.

    Returns
    -------
    Series/Index of objects

    See Also
    --------
    Series.str.lower : Converts all characters to lowercase.
    Series.str.upper : Converts all characters to uppercase.
    Series.str.title : Converts first character of each word to uppercase and
        remaining to lowercase.
    Series.str.capitalize : Converts first character to uppercase and
        remaining to lowercase.
    Series.str.swapcase : Converts uppercase to lowercase and lowercase to
        uppercase.

    Examples
    --------
    >>> s = pd.Series(['lower', 'CAPITALS', 'this is a sentence', 'SwApCaSe'])
    >>> s
    0                 lower
    1              CAPITALS
    2    this is a sentence
    3              SwApCaSe
    dtype: object

    >>> s.str.lower()
    0                 lower
    1              capitals
    2    this is a sentence
    3              swapcase
    dtype: object

    >>> s.str.upper()
    0                 LOWER
    1              CAPITALS
    2    THIS IS A SENTENCE
    3              SWAPCASE
    dtype: object

    >>> s.str.title()
    0                 Lower
    1              Capitals
    2    This Is A Sentence
    3              Swapcase
    dtype: object

    >>> s.str.capitalize()
    0                 Lower
    1              Capitals
    2    This is a sentence
    3              Swapcase
    dtype: object

    >>> s.str.swapcase()
    0                 LOWER
    1              capitals
    2    THIS IS A SENTENCE
    3              sWaPcAsE
    dtype: object
    """)
    _shared_docs['lower'] = dict(type='lowercase', method='lower')
    _shared_docs['upper'] = dict(type='uppercase', method='upper')
    _shared_docs['title'] = dict(type='titlecase', method='title')
    _shared_docs['capitalize'] = dict(type='be capitalized',
                                      method='capitalize')
    _shared_docs['swapcase'] = dict(type='be swapcased', method='swapcase')
    lower = _noarg_wrapper(lambda x: x.lower(),
                           docstring=_shared_docs['casemethods'] %
                           _shared_docs['lower'])
    upper = _noarg_wrapper(lambda x: x.upper(),
                           docstring=_shared_docs['casemethods'] %
                           _shared_docs['upper'])
    title = _noarg_wrapper(lambda x: x.title(),
                           docstring=_shared_docs['casemethods'] %
                           _shared_docs['title'])
    capitalize = _noarg_wrapper(lambda x: x.capitalize(),
                                docstring=_shared_docs['casemethods'] %
                                _shared_docs['capitalize'])
    swapcase = _noarg_wrapper(lambda x: x.swapcase(),
                              docstring=_shared_docs['casemethods'] %
                              _shared_docs['swapcase'])

    _shared_docs['ismethods'] = ("""
    Check whether all characters in each string in the Series/Index
    are %(type)s. Equivalent to :meth:`str.%(method)s`.

    Returns
    -------
    is : Series/array of boolean values
    """)
    _shared_docs['isalnum'] = dict(type='alphanumeric', method='isalnum')
    _shared_docs['isalpha'] = dict(type='alphabetic', method='isalpha')
    _shared_docs['isdigit'] = dict(type='digits', method='isdigit')
    _shared_docs['isspace'] = dict(type='whitespace', method='isspace')
    _shared_docs['islower'] = dict(type='lowercase', method='islower')
    _shared_docs['isupper'] = dict(type='uppercase', method='isupper')
    _shared_docs['istitle'] = dict(type='titlecase', method='istitle')
    _shared_docs['isnumeric'] = dict(type='numeric', method='isnumeric')
    _shared_docs['isdecimal'] = dict(type='decimal', method='isdecimal')
    isalnum = _noarg_wrapper(lambda x: x.isalnum(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['isalnum'])
    isalpha = _noarg_wrapper(lambda x: x.isalpha(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['isalpha'])
    isdigit = _noarg_wrapper(lambda x: x.isdigit(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['isdigit'])
    isspace = _noarg_wrapper(lambda x: x.isspace(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['isspace'])
    islower = _noarg_wrapper(lambda x: x.islower(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['islower'])
    isupper = _noarg_wrapper(lambda x: x.isupper(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['isupper'])
    istitle = _noarg_wrapper(lambda x: x.istitle(),
                             docstring=_shared_docs['ismethods'] %
                             _shared_docs['istitle'])
    isnumeric = _noarg_wrapper(lambda x: compat.u_safe(x).isnumeric(),
                               docstring=_shared_docs['ismethods'] %
                               _shared_docs['isnumeric'])
    isdecimal = _noarg_wrapper(lambda x: compat.u_safe(x).isdecimal(),
                               docstring=_shared_docs['ismethods'] %
                               _shared_docs['isdecimal'])

    @classmethod
    def _make_accessor(cls, data):
        cls._validate(data)
        return cls(data)
