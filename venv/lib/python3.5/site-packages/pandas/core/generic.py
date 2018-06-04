# pylint: disable=W0231,E1101
import collections
import functools
import warnings
import operator
import weakref
import gc
import json

import numpy as np
import pandas as pd

from pandas._libs import tslib, properties
from pandas.core.dtypes.common import (
    _ensure_int64,
    _ensure_object,
    is_scalar,
    is_number,
    is_integer, is_bool,
    is_bool_dtype,
    is_categorical_dtype,
    is_numeric_dtype,
    is_datetime64_dtype,
    is_timedelta64_dtype,
    is_datetime64tz_dtype,
    is_list_like,
    is_dict_like,
    is_re_compilable,
    is_period_arraylike,
    pandas_dtype)
from pandas.core.dtypes.cast import maybe_promote, maybe_upcast_putmask
from pandas.core.dtypes.inference import is_hashable
from pandas.core.dtypes.missing import isna, notna
from pandas.core.dtypes.generic import ABCSeries, ABCPanel, ABCDataFrame

from pandas.core.base import PandasObject, SelectionMixin
from pandas.core.index import (Index, MultiIndex, _ensure_index,
                               InvalidIndexError, RangeIndex)
import pandas.core.indexing as indexing
from pandas.core.indexes.datetimes import DatetimeIndex
from pandas.core.indexes.period import PeriodIndex, Period
from pandas.core.internals import BlockManager
import pandas.core.algorithms as algos
import pandas.core.common as com
import pandas.core.missing as missing
from pandas.io.formats.printing import pprint_thing
from pandas.io.formats.format import format_percentiles, DataFrameFormatter
from pandas.tseries.frequencies import to_offset
from pandas import compat
from pandas.compat.numpy import function as nv
from pandas.compat import (map, zip, lzip, lrange, string_types, to_str,
                           isidentifier, set_function_name, cPickle as pkl)
from pandas.core.ops import _align_method_FRAME
import pandas.core.nanops as nanops
from pandas.util._decorators import (Appender, Substitution,
                                     deprecate_kwarg)
from pandas.util._validators import validate_bool_kwarg, validate_fillna_kwargs
from pandas.core import config

# goal is to be able to define the docs close to function, while still being
# able to share
_shared_docs = dict()
_shared_doc_kwargs = dict(
    axes='keywords for axes', klass='NDFrame',
    axes_single_arg='int or labels for object',
    args_transpose='axes to permute (int or label for object)',
    optional_by="""
        by : str or list of str
            Name or list of names to sort by""")


def _single_replace(self, to_replace, method, inplace, limit):
    """
    Replaces values in a Series using the fill method specified when no
    replacement value is given in the replace method
    """
    if self.ndim != 1:
        raise TypeError('cannot replace {0} with method {1} on a {2}'
                        .format(to_replace, method, type(self).__name__))

    orig_dtype = self.dtype
    result = self if inplace else self.copy()
    fill_f = missing.get_fill_func(method)

    mask = missing.mask_missing(result.values, to_replace)
    values = fill_f(result.values, limit=limit, mask=mask)

    if values.dtype == orig_dtype and inplace:
        return

    result = pd.Series(values, index=self.index,
                       dtype=self.dtype).__finalize__(self)

    if inplace:
        self._update_inplace(result._data)
        return

    return result


class NDFrame(PandasObject, SelectionMixin):
    """
    N-dimensional analogue of DataFrame. Store multi-dimensional in a
    size-mutable, labeled data structure

    Parameters
    ----------
    data : BlockManager
    axes : list
    copy : boolean, default False
    """
    _internal_names = ['_data', '_cacher', '_item_cache', '_cache', '_is_copy',
                       '_subtyp', '_name', '_index', '_default_kind',
                       '_default_fill_value', '_metadata', '__array_struct__',
                       '__array_interface__']
    _internal_names_set = set(_internal_names)
    _accessors = frozenset([])
    _deprecations = frozenset(['as_blocks', 'blocks',
                               'consolidate', 'convert_objects', 'is_copy'])
    _metadata = []
    _is_copy = None

    def __init__(self, data, axes=None, copy=False, dtype=None,
                 fastpath=False):

        if not fastpath:
            if dtype is not None:
                data = data.astype(dtype)
            elif copy:
                data = data.copy()

            if axes is not None:
                for i, ax in enumerate(axes):
                    data = data.reindex_axis(ax, axis=i)

        object.__setattr__(self, '_is_copy', None)
        object.__setattr__(self, '_data', data)
        object.__setattr__(self, '_item_cache', {})

    @property
    def is_copy(self):
        warnings.warn("Attribute 'is_copy' is deprecated and will be removed "
                      "in a future version.", FutureWarning, stacklevel=2)
        return self._is_copy

    @is_copy.setter
    def is_copy(self, msg):
        warnings.warn("Attribute 'is_copy' is deprecated and will be removed "
                      "in a future version.", FutureWarning, stacklevel=2)
        self._is_copy = msg

    def _repr_data_resource_(self):
        """
        Not a real Jupyter special repr method, but we use the same
        naming convention.
        """
        if config.get_option("display.html.table_schema"):
            data = self.head(config.get_option('display.max_rows'))
            payload = json.loads(data.to_json(orient='table'),
                                 object_pairs_hook=collections.OrderedDict)
            return payload

    def _validate_dtype(self, dtype):
        """ validate the passed dtype """

        if dtype is not None:
            dtype = pandas_dtype(dtype)

            # a compound dtype
            if dtype.kind == 'V':
                raise NotImplementedError("compound dtypes are not implemented"
                                          " in the {0} constructor"
                                          .format(self.__class__.__name__))

        return dtype

    def _init_mgr(self, mgr, axes=None, dtype=None, copy=False):
        """ passed a manager and a axes dict """
        for a, axe in axes.items():
            if axe is not None:
                mgr = mgr.reindex_axis(axe,
                                       axis=self._get_block_manager_axis(a),
                                       copy=False)

        # make a copy if explicitly requested
        if copy:
            mgr = mgr.copy()
        if dtype is not None:
            # avoid further copies if we can
            if len(mgr.blocks) > 1 or mgr.blocks[0].values.dtype != dtype:
                mgr = mgr.astype(dtype=dtype)
        return mgr

    # ----------------------------------------------------------------------
    # Construction

    @property
    def _constructor(self):
        """Used when a manipulation result has the same dimensions as the
        original.
        """
        raise com.AbstractMethodError(self)

    def __unicode__(self):
        # unicode representation based upon iterating over self
        # (since, by definition, `PandasContainers` are iterable)
        prepr = '[%s]' % ','.join(map(pprint_thing, self))
        return '%s(%s)' % (self.__class__.__name__, prepr)

    def _dir_additions(self):
        """ add the string-like attributes from the info_axis.
        If info_axis is a MultiIndex, it's first level values are used.
        """
        additions = {c for c in self._info_axis.unique(level=0)[:100]
                     if isinstance(c, string_types) and isidentifier(c)}
        return super(NDFrame, self)._dir_additions().union(additions)

    @property
    def _constructor_sliced(self):
        """Used when a manipulation result has one lower dimension(s) as the
        original, such as DataFrame single columns slicing.
        """
        raise com.AbstractMethodError(self)

    @property
    def _constructor_expanddim(self):
        """Used when a manipulation result has one higher dimension as the
        original, such as Series.to_frame() and DataFrame.to_panel()
        """
        raise NotImplementedError

    # ----------------------------------------------------------------------
    # Axis

    @classmethod
    def _setup_axes(cls, axes, info_axis=None, stat_axis=None, aliases=None,
                    slicers=None, axes_are_reversed=False, build_axes=True,
                    ns=None, docs=None):
        """Provide axes setup for the major PandasObjects.

        Parameters
        ----------
        axes : the names of the axes in order (lowest to highest)
        info_axis_num : the axis of the selector dimension (int)
        stat_axis_num : the number of axis for the default stats (int)
        aliases : other names for a single axis (dict)
        slicers : how axes slice to others (dict)
        axes_are_reversed : boolean whether to treat passed axes as
            reversed (DataFrame)
        build_axes : setup the axis properties (default True)
        """

        cls._AXIS_ORDERS = axes
        cls._AXIS_NUMBERS = {a: i for i, a in enumerate(axes)}
        cls._AXIS_LEN = len(axes)
        cls._AXIS_ALIASES = aliases or dict()
        cls._AXIS_IALIASES = {v: k for k, v in cls._AXIS_ALIASES.items()}
        cls._AXIS_NAMES = dict(enumerate(axes))
        cls._AXIS_SLICEMAP = slicers or None
        cls._AXIS_REVERSED = axes_are_reversed

        # typ
        setattr(cls, '_typ', cls.__name__.lower())

        # indexing support
        cls._ix = None

        if info_axis is not None:
            cls._info_axis_number = info_axis
            cls._info_axis_name = axes[info_axis]

        if stat_axis is not None:
            cls._stat_axis_number = stat_axis
            cls._stat_axis_name = axes[stat_axis]

        # setup the actual axis
        if build_axes:

            def set_axis(a, i):
                setattr(cls, a, properties.AxisProperty(i, docs.get(a, a)))
                cls._internal_names_set.add(a)

            if axes_are_reversed:
                m = cls._AXIS_LEN - 1
                for i, a in cls._AXIS_NAMES.items():
                    set_axis(a, m - i)
            else:
                for i, a in cls._AXIS_NAMES.items():
                    set_axis(a, i)

        # addtl parms
        if isinstance(ns, dict):
            for k, v in ns.items():
                setattr(cls, k, v)

    def _construct_axes_dict(self, axes=None, **kwargs):
        """Return an axes dictionary for myself."""
        d = {a: self._get_axis(a) for a in (axes or self._AXIS_ORDERS)}
        d.update(kwargs)
        return d

    @staticmethod
    def _construct_axes_dict_from(self, axes, **kwargs):
        """Return an axes dictionary for the passed axes."""
        d = {a: ax for a, ax in zip(self._AXIS_ORDERS, axes)}
        d.update(kwargs)
        return d

    def _construct_axes_dict_for_slice(self, axes=None, **kwargs):
        """Return an axes dictionary for myself."""
        d = {self._AXIS_SLICEMAP[a]: self._get_axis(a)
             for a in (axes or self._AXIS_ORDERS)}
        d.update(kwargs)
        return d

    def _construct_axes_from_arguments(self, args, kwargs, require_all=False):
        """Construct and returns axes if supplied in args/kwargs.

        If require_all, raise if all axis arguments are not supplied
        return a tuple of (axes, kwargs).
        """

        # construct the args
        args = list(args)
        for a in self._AXIS_ORDERS:

            # if we have an alias for this axis
            alias = self._AXIS_IALIASES.get(a)
            if alias is not None:
                if a in kwargs:
                    if alias in kwargs:
                        raise TypeError("arguments are mutually exclusive "
                                        "for [%s,%s]" % (a, alias))
                    continue
                if alias in kwargs:
                    kwargs[a] = kwargs.pop(alias)
                    continue

            # look for a argument by position
            if a not in kwargs:
                try:
                    kwargs[a] = args.pop(0)
                except IndexError:
                    if require_all:
                        raise TypeError("not enough/duplicate arguments "
                                        "specified!")

        axes = {a: kwargs.pop(a, None) for a in self._AXIS_ORDERS}
        return axes, kwargs

    @classmethod
    def _from_axes(cls, data, axes, **kwargs):
        # for construction from BlockManager
        if isinstance(data, BlockManager):
            return cls(data, **kwargs)
        else:
            if cls._AXIS_REVERSED:
                axes = axes[::-1]
            d = cls._construct_axes_dict_from(cls, axes, copy=False)
            d.update(kwargs)
            return cls(data, **d)

    def _get_axis_number(self, axis):
        axis = self._AXIS_ALIASES.get(axis, axis)
        if is_integer(axis):
            if axis in self._AXIS_NAMES:
                return axis
        else:
            try:
                return self._AXIS_NUMBERS[axis]
            except KeyError:
                pass
        raise ValueError('No axis named {0} for object type {1}'
                         .format(axis, type(self)))

    def _get_axis_name(self, axis):
        axis = self._AXIS_ALIASES.get(axis, axis)
        if isinstance(axis, string_types):
            if axis in self._AXIS_NUMBERS:
                return axis
        else:
            try:
                return self._AXIS_NAMES[axis]
            except KeyError:
                pass
        raise ValueError('No axis named {0} for object type {1}'
                         .format(axis, type(self)))

    def _get_axis(self, axis):
        name = self._get_axis_name(axis)
        return getattr(self, name)

    def _get_block_manager_axis(self, axis):
        """Map the axis to the block_manager axis."""
        axis = self._get_axis_number(axis)
        if self._AXIS_REVERSED:
            m = self._AXIS_LEN - 1
            return m - axis
        return axis

    def _get_axis_resolvers(self, axis):
        # index or columns
        axis_index = getattr(self, axis)
        d = dict()
        prefix = axis[0]

        for i, name in enumerate(axis_index.names):
            if name is not None:
                key = level = name
            else:
                # prefix with 'i' or 'c' depending on the input axis
                # e.g., you must do ilevel_0 for the 0th level of an unnamed
                # multiiindex
                key = '{prefix}level_{i}'.format(prefix=prefix, i=i)
                level = i

            level_values = axis_index.get_level_values(level)
            s = level_values.to_series()
            s.index = axis_index
            d[key] = s

        # put the index/columns itself in the dict
        if isinstance(axis_index, MultiIndex):
            dindex = axis_index
        else:
            dindex = axis_index.to_series()

        d[axis] = dindex
        return d

    def _get_index_resolvers(self):
        d = {}
        for axis_name in self._AXIS_ORDERS:
            d.update(self._get_axis_resolvers(axis_name))
        return d

    @property
    def _info_axis(self):
        return getattr(self, self._info_axis_name)

    @property
    def _stat_axis(self):
        return getattr(self, self._stat_axis_name)

    @property
    def shape(self):
        """Return a tuple of axis dimensions"""
        return tuple(len(self._get_axis(a)) for a in self._AXIS_ORDERS)

    @property
    def axes(self):
        """Return index label(s) of the internal NDFrame"""
        # we do it this way because if we have reversed axes, then
        # the block manager shows then reversed
        return [self._get_axis(a) for a in self._AXIS_ORDERS]

    @property
    def ndim(self):
        """
        Return an int representing the number of axes / array dimensions.

        Return 1 if Series. Otherwise return 2 if DataFrame.

        See Also
        --------
        ndarray.ndim

        Examples
        --------
        >>> s = pd.Series({'a': 1, 'b': 2, 'c': 3})
        >>> s.ndim
        1

        >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
        >>> df.ndim
        2
        """
        return self._data.ndim

    @property
    def size(self):
        """
        Return an int representing the number of elements in this object.

        Return the number of rows if Series. Otherwise return the number of
        rows times number of columns if DataFrame.

        See Also
        --------
        ndarray.size

        Examples
        --------
        >>> s = pd.Series({'a': 1, 'b': 2, 'c': 3})
        >>> s.size
        3

        >>> df = pd.DataFrame({'col1': [1, 2], 'col2': [3, 4]})
        >>> df.size
        4
        """
        return np.prod(self.shape)

    @property
    def _selected_obj(self):
        """ internal compat with SelectionMixin """
        return self

    @property
    def _obj_with_exclusions(self):
        """ internal compat with SelectionMixin """
        return self

    def _expand_axes(self, key):
        new_axes = []
        for k, ax in zip(key, self.axes):
            if k not in ax:
                if type(k) != ax.dtype.type:
                    ax = ax.astype('O')
                new_axes.append(ax.insert(len(ax), k))
            else:
                new_axes.append(ax)

        return new_axes

    def set_axis(self, labels, axis=0, inplace=None):
        """
        Assign desired index to given axis.

        Indexes for column or row labels can be changed by assigning
        a list-like or Index.

        .. versionchanged:: 0.21.0

           The signature is now `labels` and `axis`, consistent with
           the rest of pandas API. Previously, the `axis` and `labels`
           arguments were respectively the first and second positional
           arguments.

        Parameters
        ----------
        labels : list-like, Index
            The values for the new index.

        axis : {0 or 'index', 1 or 'columns'}, default 0
            The axis to update. The value 0 identifies the rows, and 1
            identifies the columns.

        inplace : boolean, default None
            Whether to return a new %(klass)s instance.

            .. warning::

               ``inplace=None`` currently falls back to to True, but in a
               future version, will default to False. Use inplace=True
               explicitly rather than relying on the default.

        Returns
        -------
        renamed : %(klass)s or None
            An object of same type as caller if inplace=False, None otherwise.

        See Also
        --------
        pandas.DataFrame.rename_axis : Alter the name of the index or columns.

        Examples
        --------
        **Series**

        >>> s = pd.Series([1, 2, 3])
        >>> s
        0    1
        1    2
        2    3
        dtype: int64

        >>> s.set_axis(['a', 'b', 'c'], axis=0, inplace=False)
        a    1
        b    2
        c    3
        dtype: int64

        The original object is not modified.

        >>> s
        0    1
        1    2
        2    3
        dtype: int64

        **DataFrame**

        >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})

        Change the row labels.

        >>> df.set_axis(['a', 'b', 'c'], axis='index', inplace=False)
           A  B
        a  1  4
        b  2  5
        c  3  6

        Change the column labels.

        >>> df.set_axis(['I', 'II'], axis='columns', inplace=False)
           I  II
        0  1   4
        1  2   5
        2  3   6

        Now, update the labels inplace.

        >>> df.set_axis(['i', 'ii'], axis='columns', inplace=True)
        >>> df
           i  ii
        0  1   4
        1  2   5
        2  3   6
        """
        if is_scalar(labels):
            warnings.warn(
                'set_axis now takes "labels" as first argument, and '
                '"axis" as named parameter. The old form, with "axis" as '
                'first parameter and \"labels\" as second, is still supported '
                'but will be deprecated in a future version of pandas.',
                FutureWarning, stacklevel=2)
            labels, axis = axis, labels

        if inplace is None:
            warnings.warn(
                'set_axis currently defaults to operating inplace.\nThis '
                'will change in a future version of pandas, use '
                'inplace=True to avoid this warning.',
                FutureWarning, stacklevel=2)
            inplace = True
        if inplace:
            setattr(self, self._get_axis_name(axis), labels)
        else:
            obj = self.copy()
            obj.set_axis(labels, axis=axis, inplace=True)
            return obj

    def _set_axis(self, axis, labels):
        self._data.set_axis(axis, labels)
        self._clear_item_cache()

    _shared_docs['transpose'] = """
        Permute the dimensions of the %(klass)s

        Parameters
        ----------
        args : %(args_transpose)s
        copy : boolean, default False
            Make a copy of the underlying data. Mixed-dtype data will
            always result in a copy

        Examples
        --------
        >>> p.transpose(2, 0, 1)
        >>> p.transpose(2, 0, 1, copy=True)

        Returns
        -------
        y : same as input
        """

    @Appender(_shared_docs['transpose'] % _shared_doc_kwargs)
    def transpose(self, *args, **kwargs):

        # construct the args
        axes, kwargs = self._construct_axes_from_arguments(args, kwargs,
                                                           require_all=True)
        axes_names = tuple(self._get_axis_name(axes[a])
                           for a in self._AXIS_ORDERS)
        axes_numbers = tuple(self._get_axis_number(axes[a])
                             for a in self._AXIS_ORDERS)

        # we must have unique axes
        if len(axes) != len(set(axes)):
            raise ValueError('Must specify %s unique axes' % self._AXIS_LEN)

        new_axes = self._construct_axes_dict_from(self, [self._get_axis(x)
                                                         for x in axes_names])
        new_values = self.values.transpose(axes_numbers)
        if kwargs.pop('copy', None) or (len(args) and args[-1]):
            new_values = new_values.copy()

        nv.validate_transpose_for_generic(self, kwargs)
        return self._constructor(new_values, **new_axes).__finalize__(self)

    def swapaxes(self, axis1, axis2, copy=True):
        """
        Interchange axes and swap values axes appropriately

        Returns
        -------
        y : same as input
        """
        i = self._get_axis_number(axis1)
        j = self._get_axis_number(axis2)

        if i == j:
            if copy:
                return self.copy()
            return self

        mapping = {i: j, j: i}

        new_axes = (self._get_axis(mapping.get(k, k))
                    for k in range(self._AXIS_LEN))
        new_values = self.values.swapaxes(i, j)
        if copy:
            new_values = new_values.copy()

        return self._constructor(new_values, *new_axes).__finalize__(self)

    def pop(self, item):
        """
        Return item and drop from frame. Raise KeyError if not found.

        Parameters
        ----------
        item : str
            Column label to be popped

        Returns
        -------
        popped : Series

        Examples
        --------
        >>> df = pd.DataFrame([('falcon', 'bird',    389.0),
        ...                    ('parrot', 'bird',     24.0),
        ...                    ('lion',   'mammal',   80.5),
        ...                    ('monkey', 'mammal', np.nan)],
        ...                   columns=('name', 'class', 'max_speed'))
        >>> df
             name   class  max_speed
        0  falcon    bird      389.0
        1  parrot    bird       24.0
        2    lion  mammal       80.5
        3  monkey  mammal        NaN

        >>> df.pop('class')
        0      bird
        1      bird
        2    mammal
        3    mammal
        Name: class, dtype: object

        >>> df
             name  max_speed
        0  falcon      389.0
        1  parrot       24.0
        2    lion       80.5
        3  monkey        NaN
        """
        result = self[item]
        del self[item]
        try:
            result._reset_cacher()
        except AttributeError:
            pass

        return result

    def squeeze(self, axis=None):
        """
        Squeeze length 1 dimensions.

        Parameters
        ----------
        axis : None, integer or string axis name, optional
            The axis to squeeze if 1-sized.

            .. versionadded:: 0.20.0

        Returns
        -------
        scalar if 1-sized, else original object
        """
        axis = (self._AXIS_NAMES if axis is None else
                (self._get_axis_number(axis),))
        try:
            return self.iloc[
                tuple(0 if i in axis and len(a) == 1 else slice(None)
                      for i, a in enumerate(self.axes))]
        except Exception:
            return self

    def swaplevel(self, i=-2, j=-1, axis=0):
        """
        Swap levels i and j in a MultiIndex on a particular axis

        Parameters
        ----------
        i, j : int, string (can be mixed)
            Level of index to be swapped. Can pass level name as string.

        Returns
        -------
        swapped : type of caller (new object)

        .. versionchanged:: 0.18.1

           The indexes ``i`` and ``j`` are now optional, and default to
           the two innermost levels of the index.

        """
        axis = self._get_axis_number(axis)
        result = self.copy()
        labels = result._data.axes[axis]
        result._data.set_axis(axis, labels.swaplevel(i, j))
        return result

    # ----------------------------------------------------------------------
    # Rename

    # TODO: define separate funcs for DataFrame, Series and Panel so you can
    # get completion on keyword arguments.
    _shared_docs['rename'] = """
        Alter axes input function or functions. Function / dict values must be
        unique (1-to-1). Labels not contained in a dict / Series will be left
        as-is. Extra labels listed don't throw an error. Alternatively, change
        ``Series.name`` with a scalar value (Series only).

        Parameters
        ----------
        %(optional_mapper)s
        %(axes)s : scalar, list-like, dict-like or function, optional
            Scalar or list-like will alter the ``Series.name`` attribute,
            and raise on DataFrame or Panel.
            dict-like or functions are transformations to apply to
            that axis' values
        %(optional_axis)s
        copy : boolean, default True
            Also copy underlying data
        inplace : boolean, default False
            Whether to return a new %(klass)s. If True then value of copy is
            ignored.
        level : int or level name, default None
            In case of a MultiIndex, only rename labels in the specified
            level.

        Returns
        -------
        renamed : %(klass)s (new object)

        See Also
        --------
        pandas.NDFrame.rename_axis

        Examples
        --------

        >>> s = pd.Series([1, 2, 3])
        >>> s
        0    1
        1    2
        2    3
        dtype: int64
        >>> s.rename("my_name") # scalar, changes Series.name
        0    1
        1    2
        2    3
        Name: my_name, dtype: int64
        >>> s.rename(lambda x: x ** 2)  # function, changes labels
        0    1
        1    2
        4    3
        dtype: int64
        >>> s.rename({1: 3, 2: 5})  # mapping, changes labels
        0    1
        3    2
        5    3
        dtype: int64

        Since ``DataFrame`` doesn't have a ``.name`` attribute,
        only mapping-type arguments are allowed.

        >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        >>> df.rename(2)
        Traceback (most recent call last):
        ...
        TypeError: 'int' object is not callable

        ``DataFrame.rename`` supports two calling conventions

        * ``(index=index_mapper, columns=columns_mapper, ...)``
        * ``(mapper, axis={'index', 'columns'}, ...)``

        We *highly* recommend using keyword arguments to clarify your
        intent.

        >>> df.rename(index=str, columns={"A": "a", "B": "c"})
           a  c
        0  1  4
        1  2  5
        2  3  6

        >>> df.rename(index=str, columns={"A": "a", "C": "c"})
           a  B
        0  1  4
        1  2  5
        2  3  6

        Using axis-style parameters

        >>> df.rename(str.lower, axis='columns')
           a  b
        0  1  4
        1  2  5
        2  3  6

        >>> df.rename({1: 2, 2: 4}, axis='index')
           A  B
        0  1  4
        2  2  5
        4  3  6

        See the :ref:`user guide <basics.rename>` for more.
        """

    @Appender(_shared_docs['rename'] % dict(axes='axes keywords for this'
                                            ' object', klass='NDFrame',
                                            optional_mapper='',
                                            optional_axis=''))
    def rename(self, *args, **kwargs):
        axes, kwargs = self._construct_axes_from_arguments(args, kwargs)
        copy = kwargs.pop('copy', True)
        inplace = kwargs.pop('inplace', False)
        level = kwargs.pop('level', None)
        axis = kwargs.pop('axis', None)
        if axis is not None:
            axis = self._get_axis_number(axis)

        if kwargs:
            raise TypeError('rename() got an unexpected keyword '
                            'argument "{0}"'.format(list(kwargs.keys())[0]))

        if com._count_not_none(*axes.values()) == 0:
            raise TypeError('must pass an index to rename')

        # renamer function if passed a dict
        def _get_rename_function(mapper):
            if isinstance(mapper, (dict, ABCSeries)):

                def f(x):
                    if x in mapper:
                        return mapper[x]
                    else:
                        return x
            else:
                f = mapper

            return f

        self._consolidate_inplace()
        result = self if inplace else self.copy(deep=copy)

        # start in the axis order to eliminate too many copies
        for axis in lrange(self._AXIS_LEN):
            v = axes.get(self._AXIS_NAMES[axis])
            if v is None:
                continue
            f = _get_rename_function(v)

            baxis = self._get_block_manager_axis(axis)
            if level is not None:
                level = self.axes[axis]._get_level_number(level)
            result._data = result._data.rename_axis(f, axis=baxis, copy=copy,
                                                    level=level)
            result._clear_item_cache()

        if inplace:
            self._update_inplace(result._data)
        else:
            return result.__finalize__(self)

    rename.__doc__ = _shared_docs['rename']

    def rename_axis(self, mapper, axis=0, copy=True, inplace=False):
        """
        Alter the name of the index or columns.

        Parameters
        ----------
        mapper : scalar, list-like, optional
            Value to set as the axis name attribute.
        axis : {0 or 'index', 1 or 'columns'}, default 0
            The index or the name of the axis.
        copy : boolean, default True
            Also copy underlying data.
        inplace : boolean, default False
            Modifies the object directly, instead of creating a new Series
            or DataFrame.

        Returns
        -------
        renamed : Series, DataFrame, or None
            The same type as the caller or None if `inplace` is True.

        Notes
        -----
        Prior to version 0.21.0, ``rename_axis`` could also be used to change
        the axis *labels* by passing a mapping or scalar. This behavior is
        deprecated and will be removed in a future version. Use ``rename``
        instead.

        See Also
        --------
        pandas.Series.rename : Alter Series index labels or name
        pandas.DataFrame.rename : Alter DataFrame index labels or name
        pandas.Index.rename : Set new names on index

        Examples
        --------
        **Series**

        >>> s = pd.Series([1, 2, 3])
        >>> s.rename_axis("foo")
        foo
        0    1
        1    2
        2    3
        dtype: int64

        **DataFrame**

        >>> df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]})
        >>> df.rename_axis("foo")
             A  B
        foo
        0    1  4
        1    2  5
        2    3  6

        >>> df.rename_axis("bar", axis="columns")
        bar  A  B
        0    1  4
        1    2  5
        2    3  6
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        non_mapper = is_scalar(mapper) or (is_list_like(mapper) and not
                                           is_dict_like(mapper))
        if non_mapper:
            return self._set_axis_name(mapper, axis=axis, inplace=inplace)
        else:
            msg = ("Using 'rename_axis' to alter labels is deprecated. "
                   "Use '.rename' instead")
            warnings.warn(msg, FutureWarning, stacklevel=2)
            axis = self._get_axis_name(axis)
            d = {'copy': copy, 'inplace': inplace}
            d[axis] = mapper
            return self.rename(**d)

    def _set_axis_name(self, name, axis=0, inplace=False):
        """
        Alter the name or names of the axis.

        Parameters
        ----------
        name : str or list of str
            Name for the Index, or list of names for the MultiIndex
        axis : int or str
           0 or 'index' for the index; 1 or 'columns' for the columns
        inplace : bool
            whether to modify `self` directly or return a copy

            .. versionadded:: 0.21.0

        Returns
        -------
        renamed : type of caller or None if inplace=True

        See Also
        --------
        pandas.DataFrame.rename
        pandas.Series.rename
        pandas.Index.rename

        Examples
        --------
        >>> df._set_axis_name("foo")
             A
        foo
        0    1
        1    2
        2    3
        >>> df.index = pd.MultiIndex.from_product([['A'], ['a', 'b', 'c']])
        >>> df._set_axis_name(["bar", "baz"])
                 A
        bar baz
        A   a    1
            b    2
            c    3
        """
        axis = self._get_axis_number(axis)
        idx = self._get_axis(axis).set_names(name)

        inplace = validate_bool_kwarg(inplace, 'inplace')
        renamed = self if inplace else self.copy()
        renamed.set_axis(idx, axis=axis, inplace=True)
        if not inplace:
            return renamed

    # ----------------------------------------------------------------------
    # Comparisons

    def _indexed_same(self, other):
        return all(self._get_axis(a).equals(other._get_axis(a))
                   for a in self._AXIS_ORDERS)

    def __neg__(self):
        values = com._values_from_object(self)
        if is_bool_dtype(values):
            arr = operator.inv(values)
        elif (is_numeric_dtype(values) or is_timedelta64_dtype(values)):
            arr = operator.neg(values)
        else:
            raise TypeError("Unary negative expects numeric dtype, not {}"
                            .format(values.dtype))
        return self.__array_wrap__(arr)

    def __pos__(self):
        values = com._values_from_object(self)
        if (is_bool_dtype(values) or is_period_arraylike(values)):
            arr = values
        elif (is_numeric_dtype(values) or is_timedelta64_dtype(values)):
            arr = operator.pos(values)
        else:
            raise TypeError("Unary plus expects numeric dtype, not {}"
                            .format(values.dtype))
        return self.__array_wrap__(arr)

    def __invert__(self):
        try:
            arr = operator.inv(com._values_from_object(self))
            return self.__array_wrap__(arr)
        except Exception:

            # inv fails with 0 len
            if not np.prod(self.shape):
                return self

            raise

    def equals(self, other):
        """
        Determines if two NDFrame objects contain the same elements. NaNs in
        the same location are considered equal.
        """
        if not isinstance(other, self._constructor):
            return False
        return self._data.equals(other._data)

    # -------------------------------------------------------------------------
    # Label or Level Combination Helpers
    #
    # A collection of helper methods for DataFrame/Series operations that
    # accept a combination of column/index labels and levels.  All such
    # operations should utilize/extend these methods when possible so that we
    # have consistent precedence and validation logic throughout the library.

    def _is_level_reference(self, key, axis=0):
        """
        Test whether a key is a level reference for a given axis.

        To be considered a level reference, `key` must be a string that:
          - (axis=0): Matches the name of an index level and does NOT match
            a column label.
          - (axis=1): Matches the name of a column level and does NOT match
            an index label.

        Parameters
        ----------
        key: str
            Potential level name for the given axis
        axis: int, default 0
            Axis that levels are associated with (0 for index, 1 for columns)

        Returns
        -------
        is_level: bool
        """
        axis = self._get_axis_number(axis)

        if self.ndim > 2:
            raise NotImplementedError(
                "_is_level_reference is not implemented for {type}"
                .format(type=type(self)))

        return (key is not None and
                is_hashable(key) and
                key in self.axes[axis].names and
                not self._is_label_reference(key, axis=axis))

    def _is_label_reference(self, key, axis=0):
        """
        Test whether a key is a label reference for a given axis.

        To be considered a label reference, `key` must be a string that:
          - (axis=0): Matches a column label
          - (axis=1): Matches an index label

        Parameters
        ----------
        key: str
            Potential label name
        axis: int, default 0
            Axis perpendicular to the axis that labels are associated with
            (0 means search for column labels, 1 means search for index labels)

        Returns
        -------
        is_label: bool
        """
        axis = self._get_axis_number(axis)
        other_axes = [ax for ax in range(self._AXIS_LEN) if ax != axis]

        if self.ndim > 2:
            raise NotImplementedError(
                "_is_label_reference is not implemented for {type}"
                .format(type=type(self)))

        return (key is not None and
                is_hashable(key) and
                any(key in self.axes[ax] for ax in other_axes))

    def _is_label_or_level_reference(self, key, axis=0):
        """
        Test whether a key is a label or level reference for a given axis.

        To be considered either a label or a level reference, `key` must be a
        string that:
          - (axis=0): Matches a column label or an index level
          - (axis=1): Matches an index label or a column level

        Parameters
        ----------
        key: str
            Potential label or level name
        axis: int, default 0
            Axis that levels are associated with (0 for index, 1 for columns)

        Returns
        -------
        is_label_or_level: bool
        """

        if self.ndim > 2:
            raise NotImplementedError(
                "_is_label_or_level_reference is not implemented for {type}"
                .format(type=type(self)))

        return (self._is_level_reference(key, axis=axis) or
                self._is_label_reference(key, axis=axis))

    def _check_label_or_level_ambiguity(self, key, axis=0, stacklevel=1):
        """
        Check whether `key` matches both a level of the input `axis` and a
        label of the other axis and raise a ``FutureWarning`` if this is the
        case.

        Note: This method will be altered to raise an ambiguity exception in
        a future version.

        Parameters
        ----------
        key: str or object
            label or level name
        axis: int, default 0
            Axis that levels are associated with (0 for index, 1 for columns)
        stacklevel: int, default 1
            Stack level used when a FutureWarning is raised (see below).

        Returns
        -------
        ambiguous: bool

        Raises
        ------
        FutureWarning
            if `key` is ambiguous. This will become an ambiguity error in a
            future version
        """

        axis = self._get_axis_number(axis)
        other_axes = [ax for ax in range(self._AXIS_LEN) if ax != axis]

        if self.ndim > 2:
            raise NotImplementedError(
                "_check_label_or_level_ambiguity is not implemented for {type}"
                .format(type=type(self)))

        if (key is not None and
                is_hashable(key) and
                key in self.axes[axis].names and
                any(key in self.axes[ax] for ax in other_axes)):

            # Build an informative and grammatical warning
            level_article, level_type = (('an', 'index')
                                         if axis == 0 else
                                         ('a', 'column'))

            label_article, label_type = (('a', 'column')
                                         if axis == 0 else
                                         ('an', 'index'))

            msg = ("'{key}' is both {level_article} {level_type} level and "
                   "{label_article} {label_type} label.\n"
                   "Defaulting to {label_type}, but this will raise an "
                   "ambiguity error in a future version"
                   ).format(key=key,
                            level_article=level_article,
                            level_type=level_type,
                            label_article=label_article,
                            label_type=label_type)

            warnings.warn(msg, FutureWarning, stacklevel=stacklevel + 1)
            return True
        else:
            return False

    def _get_label_or_level_values(self, key, axis=0, stacklevel=1):
        """
        Return a 1-D array of values associated with `key`, a label or level
        from the given `axis`.

        Retrieval logic:
          - (axis=0): Return column values if `key` matches a column label.
            Otherwise return index level values if `key` matches an index
            level.
          - (axis=1): Return row values if `key` matches an index label.
            Otherwise return column level values if 'key' matches a column
            level

        Parameters
        ----------
        key: str
            Label or level name.
        axis: int, default 0
            Axis that levels are associated with (0 for index, 1 for columns)
        stacklevel: int, default 1
            Stack level used when a FutureWarning is raised (see below).

        Returns
        -------
        values: np.ndarray

        Raises
        ------
        KeyError
            if `key` matches neither a label nor a level
        ValueError
            if `key` matches multiple labels
        FutureWarning
            if `key` is ambiguous. This will become an ambiguity error in a
            future version
        """

        axis = self._get_axis_number(axis)
        other_axes = [ax for ax in range(self._AXIS_LEN) if ax != axis]

        if self.ndim > 2:
            raise NotImplementedError(
                "_get_label_or_level_values is not implemented for {type}"
                .format(type=type(self)))

        if self._is_label_reference(key, axis=axis):
            self._check_label_or_level_ambiguity(key, axis=axis,
                                                 stacklevel=stacklevel + 1)
            values = self.xs(key, axis=other_axes[0])._values
        elif self._is_level_reference(key, axis=axis):
            values = self.axes[axis].get_level_values(key)._values
        else:
            raise KeyError(key)

        # Check for duplicates
        if values.ndim > 1:

            if other_axes and isinstance(
                    self._get_axis(other_axes[0]), MultiIndex):
                multi_message = ('\n'
                                 'For a multi-index, the label must be a '
                                 'tuple with elements corresponding to '
                                 'each level.')
            else:
                multi_message = ''

            label_axis_name = 'column' if axis == 0 else 'index'
            raise ValueError(("The {label_axis_name} label '{key}' "
                              "is not unique.{multi_message}")
                             .format(key=key,
                                     label_axis_name=label_axis_name,
                                     multi_message=multi_message))

        return values

    def _drop_labels_or_levels(self, keys, axis=0):
        """
        Drop labels and/or levels for the given `axis`.

        For each key in `keys`:
          - (axis=0): If key matches a column label then drop the column.
            Otherwise if key matches an index level then drop the level.
          - (axis=1): If key matches an index label then drop the row.
            Otherwise if key matches a column level then drop the level.

        Parameters
        ----------
        keys: str or list of str
            labels or levels to drop
        axis: int, default 0
            Axis that levels are associated with (0 for index, 1 for columns)

        Returns
        -------
        dropped: DataFrame

        Raises
        ------
        ValueError
            if any `keys` match neither a label nor a level
        """

        axis = self._get_axis_number(axis)

        if self.ndim > 2:
            raise NotImplementedError(
                "_drop_labels_or_levels is not implemented for {type}"
                .format(type=type(self)))

        # Validate keys
        keys = com._maybe_make_list(keys)
        invalid_keys = [k for k in keys if not
                        self._is_label_or_level_reference(k, axis=axis)]

        if invalid_keys:
            raise ValueError(("The following keys are not valid labels or "
                              "levels for axis {axis}: {invalid_keys}")
                             .format(axis=axis,
                                     invalid_keys=invalid_keys))

        # Compute levels and labels to drop
        levels_to_drop = [k for k in keys
                          if self._is_level_reference(k, axis=axis)]

        labels_to_drop = [k for k in keys
                          if not self._is_level_reference(k, axis=axis)]

        # Perform copy upfront and then use inplace operations below.
        # This ensures that we always perform exactly one copy.
        # ``copy`` and/or ``inplace`` options could be added in the future.
        dropped = self.copy()

        if axis == 0:
            # Handle dropping index levels
            if levels_to_drop:
                dropped.reset_index(levels_to_drop, drop=True, inplace=True)

            # Handle dropping columns labels
            if labels_to_drop:
                dropped.drop(labels_to_drop, axis=1, inplace=True)
        else:
            # Handle dropping column levels
            if levels_to_drop:
                if isinstance(dropped.columns, MultiIndex):
                    # Drop the specified levels from the MultiIndex
                    dropped.columns = dropped.columns.droplevel(levels_to_drop)
                else:
                    # Drop the last level of Index by replacing with
                    # a RangeIndex
                    dropped.columns = RangeIndex(dropped.columns.size)

            # Handle dropping index labels
            if labels_to_drop:
                dropped.drop(labels_to_drop, axis=0, inplace=True)

        return dropped

    # ----------------------------------------------------------------------
    # Iteration

    def __hash__(self):
        raise TypeError('{0!r} objects are mutable, thus they cannot be'
                        ' hashed'.format(self.__class__.__name__))

    def __iter__(self):
        """Iterate over infor axis"""
        return iter(self._info_axis)

    # can we get a better explanation of this?
    def keys(self):
        """Get the 'info axis' (see Indexing for more)

        This is index for Series, columns for DataFrame and major_axis for
        Panel.
        """
        return self._info_axis

    def iteritems(self):
        """Iterate over (label, values) on info axis

        This is index for Series, columns for DataFrame, major_axis for Panel,
        and so on.
        """
        for h in self._info_axis:
            yield h, self[h]

    def __len__(self):
        """Returns length of info axis"""
        return len(self._info_axis)

    def __contains__(self, key):
        """True if the key is in the info axis"""
        return key in self._info_axis

    @property
    def empty(self):
        """
        Indicator whether DataFrame is empty.

        True if DataFrame is entirely empty (no items), meaning any of the
        axes are of length 0.

        Returns
        -------
        bool
            If DataFrame is empty, return True, if not return False.

        Notes
        -----
        If DataFrame contains only NaNs, it is still not considered empty. See
        the example below.

        Examples
        --------
        An example of an actual empty DataFrame. Notice the index is empty:

        >>> df_empty = pd.DataFrame({'A' : []})
        >>> df_empty
        Empty DataFrame
        Columns: [A]
        Index: []
        >>> df_empty.empty
        True

        If we only have NaNs in our DataFrame, it is not considered empty! We
        will need to drop the NaNs to make the DataFrame empty:

        >>> df = pd.DataFrame({'A' : [np.nan]})
        >>> df
            A
        0 NaN
        >>> df.empty
        False
        >>> df.dropna().empty
        True

        See also
        --------
        pandas.Series.dropna
        pandas.DataFrame.dropna
        """
        return any(len(self._get_axis(a)) == 0 for a in self._AXIS_ORDERS)

    def __nonzero__(self):
        raise ValueError("The truth value of a {0} is ambiguous. "
                         "Use a.empty, a.bool(), a.item(), a.any() or a.all()."
                         .format(self.__class__.__name__))

    __bool__ = __nonzero__

    def bool(self):
        """Return the bool of a single element PandasObject.

        This must be a boolean scalar value, either True or False.  Raise a
        ValueError if the PandasObject does not have exactly 1 element, or that
        element is not boolean
        """
        v = self.squeeze()
        if isinstance(v, (bool, np.bool_)):
            return bool(v)
        elif is_scalar(v):
            raise ValueError("bool cannot act on a non-boolean single element "
                             "{0}".format(self.__class__.__name__))

        self.__nonzero__()

    def __abs__(self):
        return self.abs()

    def __round__(self, decimals=0):
        return self.round(decimals)

    # ----------------------------------------------------------------------
    # Array Interface

    def __array__(self, dtype=None):
        return com._values_from_object(self)

    def __array_wrap__(self, result, context=None):
        d = self._construct_axes_dict(self._AXIS_ORDERS, copy=False)
        return self._constructor(result, **d).__finalize__(self)

    # ideally we would define this to avoid the getattr checks, but
    # is slower
    # @property
    # def __array_interface__(self):
    #    """ provide numpy array interface method """
    #    values = self.values
    #    return dict(typestr=values.dtype.str,shape=values.shape,data=values)

    def to_dense(self):
        """Return dense representation of NDFrame (as opposed to sparse)"""
        # compat
        return self

    # ----------------------------------------------------------------------
    # Picklability

    def __getstate__(self):
        meta = {k: getattr(self, k, None) for k in self._metadata}
        return dict(_data=self._data, _typ=self._typ, _metadata=self._metadata,
                    **meta)

    def __setstate__(self, state):

        if isinstance(state, BlockManager):
            self._data = state
        elif isinstance(state, dict):
            typ = state.get('_typ')
            if typ is not None:

                # set in the order of internal names
                # to avoid definitional recursion
                # e.g. say fill_value needing _data to be
                # defined
                meta = set(self._internal_names + self._metadata)
                for k in list(meta):
                    if k in state:
                        v = state[k]
                        object.__setattr__(self, k, v)

                for k, v in state.items():
                    if k not in meta:
                        object.__setattr__(self, k, v)

            else:
                self._unpickle_series_compat(state)
        elif isinstance(state[0], dict):
            if len(state) == 5:
                self._unpickle_sparse_frame_compat(state)
            else:
                self._unpickle_frame_compat(state)
        elif len(state) == 4:
            self._unpickle_panel_compat(state)
        elif len(state) == 2:
            self._unpickle_series_compat(state)
        else:  # pragma: no cover
            # old pickling format, for compatibility
            self._unpickle_matrix_compat(state)

        self._item_cache = {}

    # ----------------------------------------------------------------------
    # IO

    def _repr_latex_(self):
        """
        Returns a LaTeX representation for a particular object.
        Mainly for use with nbconvert (jupyter notebook conversion to pdf).
        """
        if config.get_option('display.latex.repr'):
            return self.to_latex()
        else:
            return None

    # ----------------------------------------------------------------------
    # I/O Methods

    _shared_docs['to_excel'] = """
    Write %(klass)s to an excel sheet
    %(versionadded_to_excel)s

    Parameters
    ----------
    excel_writer : string or ExcelWriter object
        File path or existing ExcelWriter
    sheet_name : string, default 'Sheet1'
        Name of sheet which will contain DataFrame
    na_rep : string, default ''
        Missing data representation
    float_format : string, default None
        Format string for floating point numbers
    columns : sequence, optional
        Columns to write
    header : boolean or list of string, default True
        Write out the column names. If a list of strings is given it is
        assumed to be aliases for the column names
    index : boolean, default True
        Write row names (index)
    index_label : string or sequence, default None
        Column label for index column(s) if desired. If None is given, and
        `header` and `index` are True, then the index names are used. A
        sequence should be given if the DataFrame uses MultiIndex.
    startrow :
        upper left cell row to dump data frame
    startcol :
        upper left cell column to dump data frame
    engine : string, default None
        write engine to use - you can also set this via the options
        ``io.excel.xlsx.writer``, ``io.excel.xls.writer``, and
        ``io.excel.xlsm.writer``.
    merge_cells : boolean, default True
        Write MultiIndex and Hierarchical Rows as merged cells.
    encoding: string, default None
        encoding of the resulting excel file. Only necessary for xlwt,
        other writers support unicode natively.
    inf_rep : string, default 'inf'
        Representation for infinity (there is no native representation for
        infinity in Excel)
    freeze_panes : tuple of integer (length 2), default None
        Specifies the one-based bottommost row and rightmost column that
        is to be frozen

        .. versionadded:: 0.20.0

    Notes
    -----
    If passing an existing ExcelWriter object, then the sheet will be added
    to the existing workbook.  This can be used to save different
    DataFrames to one workbook:

    >>> writer = pd.ExcelWriter('output.xlsx')
    >>> df1.to_excel(writer,'Sheet1')
    >>> df2.to_excel(writer,'Sheet2')
    >>> writer.save()

    For compatibility with to_csv, to_excel serializes lists and dicts to
    strings before writing.
    """

    def to_json(self, path_or_buf=None, orient=None, date_format=None,
                double_precision=10, force_ascii=True, date_unit='ms',
                default_handler=None, lines=False, compression=None,
                index=True):
        """
        Convert the object to a JSON string.

        Note NaN's and None will be converted to null and datetime objects
        will be converted to UNIX timestamps.

        Parameters
        ----------
        path_or_buf : string or file handle, optional
            File path or object. If not specified, the result is returned as
            a string.
        orient : string
            Indication of expected JSON string format.

            * Series

              - default is 'index'
              - allowed values are: {'split','records','index'}

            * DataFrame

              - default is 'columns'
              - allowed values are:
                {'split','records','index','columns','values'}

            * The format of the JSON string

              - 'split' : dict like {'index' -> [index],
                'columns' -> [columns], 'data' -> [values]}
              - 'records' : list like
                [{column -> value}, ... , {column -> value}]
              - 'index' : dict like {index -> {column -> value}}
              - 'columns' : dict like {column -> {index -> value}}
              - 'values' : just the values array
              - 'table' : dict like {'schema': {schema}, 'data': {data}}
                describing the data, and the data component is
                like ``orient='records'``.

                .. versionchanged:: 0.20.0

        date_format : {None, 'epoch', 'iso'}
            Type of date conversion. 'epoch' = epoch milliseconds,
            'iso' = ISO8601. The default depends on the `orient`. For
            ``orient='table'``, the default is 'iso'. For all other orients,
            the default is 'epoch'.
        double_precision : int, default 10
            The number of decimal places to use when encoding
            floating point values.
        force_ascii : boolean, default True
            Force encoded string to be ASCII.
        date_unit : string, default 'ms' (milliseconds)
            The time unit to encode to, governs timestamp and ISO8601
            precision.  One of 's', 'ms', 'us', 'ns' for second, millisecond,
            microsecond, and nanosecond respectively.
        default_handler : callable, default None
            Handler to call if object cannot otherwise be converted to a
            suitable format for JSON. Should receive a single argument which is
            the object to convert and return a serialisable object.
        lines : boolean, default False
            If 'orient' is 'records' write out line delimited json format. Will
            throw ValueError if incorrect 'orient' since others are not list
            like.

            .. versionadded:: 0.19.0

        compression : {None, 'gzip', 'bz2', 'zip', 'xz'}
            A string representing the compression to use in the output file,
            only used when the first argument is a filename.

            .. versionadded:: 0.21.0

        index : boolean, default True
            Whether to include the index values in the JSON string. Not
            including the index (``index=False``) is only supported when
            orient is 'split' or 'table'.

            .. versionadded:: 0.23.0

        See Also
        --------
        pandas.read_json

        Examples
        --------

        >>> df = pd.DataFrame([['a', 'b'], ['c', 'd']],
        ...                   index=['row 1', 'row 2'],
        ...                   columns=['col 1', 'col 2'])
        >>> df.to_json(orient='split')
        '{"columns":["col 1","col 2"],
          "index":["row 1","row 2"],
          "data":[["a","b"],["c","d"]]}'

        Encoding/decoding a Dataframe using ``'records'`` formatted JSON.
        Note that index labels are not preserved with this encoding.

        >>> df.to_json(orient='records')
        '[{"col 1":"a","col 2":"b"},{"col 1":"c","col 2":"d"}]'

        Encoding/decoding a Dataframe using ``'index'`` formatted JSON:

        >>> df.to_json(orient='index')
        '{"row 1":{"col 1":"a","col 2":"b"},"row 2":{"col 1":"c","col 2":"d"}}'

        Encoding/decoding a Dataframe using ``'columns'`` formatted JSON:

        >>> df.to_json(orient='columns')
        '{"col 1":{"row 1":"a","row 2":"c"},"col 2":{"row 1":"b","row 2":"d"}}'

        Encoding/decoding a Dataframe using ``'values'`` formatted JSON:

        >>> df.to_json(orient='values')
        '[["a","b"],["c","d"]]'

        Encoding with Table Schema

        >>> df.to_json(orient='table')
        '{"schema": {"fields": [{"name": "index", "type": "string"},
                                {"name": "col 1", "type": "string"},
                                {"name": "col 2", "type": "string"}],
                     "primaryKey": "index",
                     "pandas_version": "0.20.0"},
          "data": [{"index": "row 1", "col 1": "a", "col 2": "b"},
                   {"index": "row 2", "col 1": "c", "col 2": "d"}]}'
        """

        from pandas.io import json
        if date_format is None and orient == 'table':
            date_format = 'iso'
        elif date_format is None:
            date_format = 'epoch'
        return json.to_json(path_or_buf=path_or_buf, obj=self, orient=orient,
                            date_format=date_format,
                            double_precision=double_precision,
                            force_ascii=force_ascii, date_unit=date_unit,
                            default_handler=default_handler,
                            lines=lines, compression=compression,
                            index=index)

    def to_hdf(self, path_or_buf, key, **kwargs):
        """
        Write the contained data to an HDF5 file using HDFStore.

        Hierarchical Data Format (HDF) is self-describing, allowing an
        application to interpret the structure and contents of a file with
        no outside information. One HDF file can hold a mix of related objects
        which can be accessed as a group or as individual objects.

        In order to add another DataFrame or Series to an existing HDF file
        please use append mode and a different a key.

        For more information see the :ref:`user guide <io.hdf5>`.

        Parameters
        ----------
        path_or_buf : str or pandas.HDFStore
            File path or HDFStore object.
        key : str
            Identifier for the group in the store.
        mode : {'a', 'w', 'r+'}, default 'a'
            Mode to open file:

            - 'w': write, a new file is created (an existing file with
              the same name would be deleted).
            - 'a': append, an existing file is opened for reading and
              writing, and if the file does not exist it is created.
            - 'r+': similar to 'a', but the file must already exist.
        format : {'fixed', 'table'}, default 'fixed'
            Possible values:

            - 'fixed': Fixed format. Fast writing/reading. Not-appendable,
              nor searchable.
            - 'table': Table format. Write as a PyTables Table structure
              which may perform worse but allow more flexible operations
              like searching / selecting subsets of the data.
        append : bool, default False
            For Table formats, append the input data to the existing.
        data_columns :  list of columns or True, optional
            List of columns to create as indexed data columns for on-disk
            queries, or True to use all columns. By default only the axes
            of the object are indexed. See :ref:`io.hdf5-query-data-columns`.
            Applicable only to format='table'.
        complevel : {0-9}, optional
            Specifies a compression level for data.
            A value of 0 disables compression.
        complib : {'zlib', 'lzo', 'bzip2', 'blosc'}, default 'zlib'
            Specifies the compression library to be used.
            As of v0.20.2 these additional compressors for Blosc are supported
            (default if no compressor specified: 'blosc:blosclz'):
            {'blosc:blosclz', 'blosc:lz4', 'blosc:lz4hc', 'blosc:snappy',
            'blosc:zlib', 'blosc:zstd'}.
            Specifying a compression library which is not available issues
            a ValueError.
        fletcher32 : bool, default False
            If applying compression use the fletcher32 checksum.
        dropna : bool, default False
            If true, ALL nan rows will not be written to store.
        errors : str, default 'strict'
            Specifies how encoding and decoding errors are to be handled.
            See the errors argument for :func:`open` for a full list
            of options.

        See Also
        --------
        DataFrame.read_hdf : Read from HDF file.
        DataFrame.to_parquet : Write a DataFrame to the binary parquet format.
        DataFrame.to_sql : Write to a sql table.
        DataFrame.to_feather : Write out feather-format for DataFrames.
        DataFrame.to_csv : Write out to a csv file.

        Examples
        --------
        >>> df = pd.DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]},
        ...                   index=['a', 'b', 'c'])
        >>> df.to_hdf('data.h5', key='df', mode='w')

        We can add another object to the same file:

        >>> s = pd.Series([1, 2, 3, 4])
        >>> s.to_hdf('data.h5', key='s')

        Reading from HDF file:

        >>> pd.read_hdf('data.h5', 'df')
        A  B
        a  1  4
        b  2  5
        c  3  6
        >>> pd.read_hdf('data.h5', 's')
        0    1
        1    2
        2    3
        3    4
        dtype: int64

        Deleting file with data:

        >>> import os
        >>> os.remove('data.h5')

        """
        from pandas.io import pytables
        return pytables.to_hdf(path_or_buf, key, self, **kwargs)

    def to_msgpack(self, path_or_buf=None, encoding='utf-8', **kwargs):
        """
        msgpack (serialize) object to input file path

        THIS IS AN EXPERIMENTAL LIBRARY and the storage format
        may not be stable until a future release.

        Parameters
        ----------
        path : string File path, buffer-like, or None
            if None, return generated string
        append : boolean whether to append to an existing msgpack
            (default is False)
        compress : type of compressor (zlib or blosc), default to None (no
            compression)
        """

        from pandas.io import packers
        return packers.to_msgpack(path_or_buf, self, encoding=encoding,
                                  **kwargs)

    def to_sql(self, name, con, schema=None, if_exists='fail', index=True,
               index_label=None, chunksize=None, dtype=None):
        """
        Write records stored in a DataFrame to a SQL database.

        Databases supported by SQLAlchemy [1]_ are supported. Tables can be
        newly created, appended to, or overwritten.

        Parameters
        ----------
        name : string
            Name of SQL table.
        con : sqlalchemy.engine.Engine or sqlite3.Connection
            Using SQLAlchemy makes it possible to use any DB supported by that
            library. Legacy support is provided for sqlite3.Connection objects.
        schema : string, optional
            Specify the schema (if database flavor supports this). If None, use
            default schema.
        if_exists : {'fail', 'replace', 'append'}, default 'fail'
            How to behave if the table already exists.

            * fail: Raise a ValueError.
            * replace: Drop the table before inserting new values.
            * append: Insert new values to the existing table.

        index : boolean, default True
            Write DataFrame index as a column. Uses `index_label` as the column
            name in the table.
        index_label : string or sequence, default None
            Column label for index column(s). If None is given (default) and
            `index` is True, then the index names are used.
            A sequence should be given if the DataFrame uses MultiIndex.
        chunksize : int, optional
            Rows will be written in batches of this size at a time. By default,
            all rows will be written at once.
        dtype : dict, optional
            Specifying the datatype for columns. The keys should be the column
            names and the values should be the SQLAlchemy types or strings for
            the sqlite3 legacy mode.

        Raises
        ------
        ValueError
            When the table already exists and `if_exists` is 'fail' (the
            default).

        See Also
        --------
        pandas.read_sql : read a DataFrame from a table

        References
        ----------
        .. [1] http://docs.sqlalchemy.org
        .. [2] https://www.python.org/dev/peps/pep-0249/

        Examples
        --------

        Create an in-memory SQLite database.

        >>> from sqlalchemy import create_engine
        >>> engine = create_engine('sqlite://', echo=False)

        Create a table from scratch with 3 rows.

        >>> df = pd.DataFrame({'name' : ['User 1', 'User 2', 'User 3']})
        >>> df
             name
        0  User 1
        1  User 2
        2  User 3

        >>> df.to_sql('users', con=engine)
        >>> engine.execute("SELECT * FROM users").fetchall()
        [(0, 'User 1'), (1, 'User 2'), (2, 'User 3')]

        >>> df1 = pd.DataFrame({'name' : ['User 4', 'User 5']})
        >>> df1.to_sql('users', con=engine, if_exists='append')
        >>> engine.execute("SELECT * FROM users").fetchall()
        [(0, 'User 1'), (1, 'User 2'), (2, 'User 3'),
         (0, 'User 4'), (1, 'User 5')]

        Overwrite the table with just ``df1``.

        >>> df1.to_sql('users', con=engine, if_exists='replace',
        ...            index_label='id')
        >>> engine.execute("SELECT * FROM users").fetchall()
        [(0, 'User 4'), (1, 'User 5')]

        Specify the dtype (especially useful for integers with missing values).
        Notice that while pandas is forced to store the data as floating point,
        the database supports nullable integers. When fetching the data with
        Python, we get back integer scalars.

        >>> df = pd.DataFrame({"A": [1, None, 2]})
        >>> df
             A
        0  1.0
        1  NaN
        2  2.0

        >>> from sqlalchemy.types import Integer
        >>> df.to_sql('integers', con=engine, index=False,
        ...           dtype={"A": Integer()})

        >>> engine.execute("SELECT * FROM integers").fetchall()
        [(1,), (None,), (2,)]
        """
        from pandas.io import sql
        sql.to_sql(self, name, con, schema=schema, if_exists=if_exists,
                   index=index, index_label=index_label, chunksize=chunksize,
                   dtype=dtype)

    def to_pickle(self, path, compression='infer',
                  protocol=pkl.HIGHEST_PROTOCOL):
        """
        Pickle (serialize) object to file.

        Parameters
        ----------
        path : str
            File path where the pickled object will be stored.
        compression : {'infer', 'gzip', 'bz2', 'zip', 'xz', None}, \
        default 'infer'
            A string representing the compression to use in the output file. By
            default, infers from the file extension in specified path.

            .. versionadded:: 0.20.0
        protocol : int
            Int which indicates which protocol should be used by the pickler,
            default HIGHEST_PROTOCOL (see [1]_ paragraph 12.1.2). The possible
            values for this parameter depend on the version of Python. For
            Python 2.x, possible values are 0, 1, 2. For Python>=3.0, 3 is a
            valid value. For Python >= 3.4, 4 is a valid value. A negative
            value for the protocol parameter is equivalent to setting its value
            to HIGHEST_PROTOCOL.

            .. [1] https://docs.python.org/3/library/pickle.html
            .. versionadded:: 0.21.0

        See Also
        --------
        read_pickle : Load pickled pandas object (or any object) from file.
        DataFrame.to_hdf : Write DataFrame to an HDF5 file.
        DataFrame.to_sql : Write DataFrame to a SQL database.
        DataFrame.to_parquet : Write a DataFrame to the binary parquet format.

        Examples
        --------
        >>> original_df = pd.DataFrame({"foo": range(5), "bar": range(5, 10)})
        >>> original_df
           foo  bar
        0    0    5
        1    1    6
        2    2    7
        3    3    8
        4    4    9
        >>> original_df.to_pickle("./dummy.pkl")

        >>> unpickled_df = pd.read_pickle("./dummy.pkl")
        >>> unpickled_df
           foo  bar
        0    0    5
        1    1    6
        2    2    7
        3    3    8
        4    4    9

        >>> import os
        >>> os.remove("./dummy.pkl")
        """
        from pandas.io.pickle import to_pickle
        return to_pickle(self, path, compression=compression,
                         protocol=protocol)

    def to_clipboard(self, excel=True, sep=None, **kwargs):
        r"""
        Copy object to the system clipboard.

        Write a text representation of object to the system clipboard.
        This can be pasted into Excel, for example.

        Parameters
        ----------
        excel : bool, default True
            - True, use the provided separator, writing in a csv format for
              allowing easy pasting into excel.
            - False, write a string representation of the object to the
              clipboard.

        sep : str, default ``'\t'``
            Field delimiter.
        **kwargs
            These parameters will be passed to DataFrame.to_csv.

        See Also
        --------
        DataFrame.to_csv : Write a DataFrame to a comma-separated values
            (csv) file.
        read_clipboard : Read text from clipboard and pass to read_table.

        Notes
        -----
        Requirements for your platform.

          - Linux : `xclip`, or `xsel` (with `gtk` or `PyQt4` modules)
          - Windows : none
          - OS X : none

        Examples
        --------
        Copy the contents of a DataFrame to the clipboard.

        >>> df = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=['A', 'B', 'C'])
        >>> df.to_clipboard(sep=',')
        ... # Wrote the following to the system clipboard:
        ... # ,A,B,C
        ... # 0,1,2,3
        ... # 1,4,5,6

        We can omit the the index by passing the keyword `index` and setting
        it to false.

        >>> df.to_clipboard(sep=',', index=False)
        ... # Wrote the following to the system clipboard:
        ... # A,B,C
        ... # 1,2,3
        ... # 4,5,6
        """
        from pandas.io import clipboards
        clipboards.to_clipboard(self, excel=excel, sep=sep, **kwargs)

    def to_xarray(self):
        """
        Return an xarray object from the pandas object.

        Returns
        -------
        a DataArray for a Series
        a Dataset for a DataFrame
        a DataArray for higher dims

        Examples
        --------
        >>> df = pd.DataFrame({'A' : [1, 1, 2],
                               'B' : ['foo', 'bar', 'foo'],
                               'C' : np.arange(4.,7)})
        >>> df
           A    B    C
        0  1  foo  4.0
        1  1  bar  5.0
        2  2  foo  6.0

        >>> df.to_xarray()
        <xarray.Dataset>
        Dimensions:  (index: 3)
        Coordinates:
          * index    (index) int64 0 1 2
        Data variables:
            A        (index) int64 1 1 2
            B        (index) object 'foo' 'bar' 'foo'
            C        (index) float64 4.0 5.0 6.0

        >>> df = pd.DataFrame({'A' : [1, 1, 2],
                               'B' : ['foo', 'bar', 'foo'],
                               'C' : np.arange(4.,7)}
                             ).set_index(['B','A'])
        >>> df
                 C
        B   A
        foo 1  4.0
        bar 1  5.0
        foo 2  6.0

        >>> df.to_xarray()
        <xarray.Dataset>
        Dimensions:  (A: 2, B: 2)
        Coordinates:
          * B        (B) object 'bar' 'foo'
          * A        (A) int64 1 2
        Data variables:
            C        (B, A) float64 5.0 nan 4.0 6.0

        >>> p = pd.Panel(np.arange(24).reshape(4,3,2),
                         items=list('ABCD'),
                         major_axis=pd.date_range('20130101', periods=3),
                         minor_axis=['first', 'second'])
        >>> p
        <class 'pandas.core.panel.Panel'>
        Dimensions: 4 (items) x 3 (major_axis) x 2 (minor_axis)
        Items axis: A to D
        Major_axis axis: 2013-01-01 00:00:00 to 2013-01-03 00:00:00
        Minor_axis axis: first to second

        >>> p.to_xarray()
        <xarray.DataArray (items: 4, major_axis: 3, minor_axis: 2)>
        array([[[ 0,  1],
                [ 2,  3],
                [ 4,  5]],
               [[ 6,  7],
                [ 8,  9],
                [10, 11]],
               [[12, 13],
                [14, 15],
                [16, 17]],
               [[18, 19],
                [20, 21],
                [22, 23]]])
        Coordinates:
          * items       (items) object 'A' 'B' 'C' 'D'
          * major_axis  (major_axis) datetime64[ns] 2013-01-01 2013-01-02 2013-01-03  # noqa
          * minor_axis  (minor_axis) object 'first' 'second'

        Notes
        -----
        See the `xarray docs <http://xarray.pydata.org/en/stable/>`__
        """

        try:
            import xarray
        except ImportError:
            # Give a nice error message
            raise ImportError("the xarray library is not installed\n"
                              "you can install via conda\n"
                              "conda install xarray\n"
                              "or via pip\n"
                              "pip install xarray\n")

        if self.ndim == 1:
            return xarray.DataArray.from_series(self)
        elif self.ndim == 2:
            return xarray.Dataset.from_dataframe(self)

        # > 2 dims
        coords = [(a, self._get_axis(a)) for a in self._AXIS_ORDERS]
        return xarray.DataArray(self,
                                coords=coords,
                                )

    _shared_docs['to_latex'] = r"""
        Render an object to a tabular environment table. You can splice
        this into a LaTeX document. Requires \\usepackage{booktabs}.

        .. versionchanged:: 0.20.2
           Added to Series

        `to_latex`-specific options:

        bold_rows : boolean, default False
            Make the row labels bold in the output
        column_format : str, default None
            The columns format as specified in `LaTeX table format
            <https://en.wikibooks.org/wiki/LaTeX/Tables>`__ e.g 'rcl' for 3
            columns
        longtable : boolean, default will be read from the pandas config module
            Default: False.
            Use a longtable environment instead of tabular. Requires adding
            a \\usepackage{longtable} to your LaTeX preamble.
        escape : boolean, default will be read from the pandas config module
            Default: True.
            When set to False prevents from escaping latex special
            characters in column names.
        encoding : str, default None
            A string representing the encoding to use in the output file,
            defaults to 'ascii' on Python 2 and 'utf-8' on Python 3.
        decimal : string, default '.'
            Character recognized as decimal separator, e.g. ',' in Europe.

            .. versionadded:: 0.18.0

        multicolumn : boolean, default True
            Use \multicolumn to enhance MultiIndex columns.
            The default will be read from the config module.

            .. versionadded:: 0.20.0

        multicolumn_format : str, default 'l'
            The alignment for multicolumns, similar to `column_format`
            The default will be read from the config module.

            .. versionadded:: 0.20.0

        multirow : boolean, default False
            Use \multirow to enhance MultiIndex rows.
            Requires adding a \\usepackage{multirow} to your LaTeX preamble.
            Will print centered labels (instead of top-aligned)
            across the contained rows, separating groups via clines.
            The default will be read from the pandas config module.

            .. versionadded:: 0.20.0
            """

    @Substitution(header='Write out the column names. If a list of strings '
                         'is given, it is assumed to be aliases for the '
                         'column names.')
    @Appender(_shared_docs['to_latex'] % _shared_doc_kwargs)
    def to_latex(self, buf=None, columns=None, col_space=None, header=True,
                 index=True, na_rep='NaN', formatters=None, float_format=None,
                 sparsify=None, index_names=True, bold_rows=False,
                 column_format=None, longtable=None, escape=None,
                 encoding=None, decimal='.', multicolumn=None,
                 multicolumn_format=None, multirow=None):
        # Get defaults from the pandas config
        if self.ndim == 1:
            self = self.to_frame()
        if longtable is None:
            longtable = config.get_option("display.latex.longtable")
        if escape is None:
            escape = config.get_option("display.latex.escape")
        if multicolumn is None:
            multicolumn = config.get_option("display.latex.multicolumn")
        if multicolumn_format is None:
            multicolumn_format = config.get_option(
                "display.latex.multicolumn_format")
        if multirow is None:
            multirow = config.get_option("display.latex.multirow")

        formatter = DataFrameFormatter(self, buf=buf, columns=columns,
                                       col_space=col_space, na_rep=na_rep,
                                       header=header, index=index,
                                       formatters=formatters,
                                       float_format=float_format,
                                       bold_rows=bold_rows,
                                       sparsify=sparsify,
                                       index_names=index_names,
                                       escape=escape, decimal=decimal)
        formatter.to_latex(column_format=column_format, longtable=longtable,
                           encoding=encoding, multicolumn=multicolumn,
                           multicolumn_format=multicolumn_format,
                           multirow=multirow)

        if buf is None:
            return formatter.buf.getvalue()

    # ----------------------------------------------------------------------
    # Fancy Indexing

    @classmethod
    def _create_indexer(cls, name, indexer):
        """Create an indexer like _name in the class."""
        if getattr(cls, name, None) is None:
            _indexer = functools.partial(indexer, name)
            setattr(cls, name, property(_indexer, doc=indexer.__doc__))

    def get(self, key, default=None):
        """
        Get item from object for given key (DataFrame column, Panel slice,
        etc.). Returns default value if not found.

        Parameters
        ----------
        key : object

        Returns
        -------
        value : type of items contained in object
        """
        try:
            return self[key]
        except (KeyError, ValueError, IndexError):
            return default

    def __getitem__(self, item):
        return self._get_item_cache(item)

    def _get_item_cache(self, item):
        """Return the cached item, item represents a label indexer."""
        cache = self._item_cache
        res = cache.get(item)
        if res is None:
            values = self._data.get(item)
            res = self._box_item_values(item, values)
            cache[item] = res
            res._set_as_cached(item, self)

            # for a chain
            res._is_copy = self._is_copy
        return res

    def _set_as_cached(self, item, cacher):
        """Set the _cacher attribute on the calling object with a weakref to
        cacher.
        """
        self._cacher = (item, weakref.ref(cacher))

    def _reset_cacher(self):
        """Reset the cacher."""
        if hasattr(self, '_cacher'):
            del self._cacher

    def _iget_item_cache(self, item):
        """Return the cached item, item represents a positional indexer."""
        ax = self._info_axis
        if ax.is_unique:
            lower = self._get_item_cache(ax[item])
        else:
            lower = self._take(item, axis=self._info_axis_number)
        return lower

    def _box_item_values(self, key, values):
        raise com.AbstractMethodError(self)

    def _maybe_cache_changed(self, item, value):
        """The object has called back to us saying maybe it has changed.
        """
        self._data.set(item, value, check=False)

    @property
    def _is_cached(self):
        """Return boolean indicating if self is cached or not."""
        return getattr(self, '_cacher', None) is not None

    def _get_cacher(self):
        """return my cacher or None"""
        cacher = getattr(self, '_cacher', None)
        if cacher is not None:
            cacher = cacher[1]()
        return cacher

    @property
    def _is_view(self):
        """Return boolean indicating if self is view of another array """
        return self._data.is_view

    def _maybe_update_cacher(self, clear=False, verify_is_copy=True):
        """
        See if we need to update our parent cacher if clear, then clear our
        cache.

        Parameters
        ----------
        clear : boolean, default False
            clear the item cache
        verify_is_copy : boolean, default True
            provide is_copy checks

        """

        cacher = getattr(self, '_cacher', None)
        if cacher is not None:
            ref = cacher[1]()

            # we are trying to reference a dead referant, hence
            # a copy
            if ref is None:
                del self._cacher
            else:
                try:
                    ref._maybe_cache_changed(cacher[0], self)
                except Exception:
                    pass

        if verify_is_copy:
            self._check_setitem_copy(stacklevel=5, t='referant')

        if clear:
            self._clear_item_cache()

    def _clear_item_cache(self, i=None):
        if i is not None:
            self._item_cache.pop(i, None)
        else:
            self._item_cache.clear()

    def _slice(self, slobj, axis=0, kind=None):
        """
        Construct a slice of this container.

        kind parameter is maintained for compatibility with Series slicing.
        """
        axis = self._get_block_manager_axis(axis)
        result = self._constructor(self._data.get_slice(slobj, axis=axis))
        result = result.__finalize__(self)

        # this could be a view
        # but only in a single-dtyped view slicable case
        is_copy = axis != 0 or result._is_view
        result._set_is_copy(self, copy=is_copy)
        return result

    def _set_item(self, key, value):
        self._data.set(key, value)
        self._clear_item_cache()

    def _set_is_copy(self, ref=None, copy=True):
        if not copy:
            self._is_copy = None
        else:
            if ref is not None:
                self._is_copy = weakref.ref(ref)
            else:
                self._is_copy = None

    def _check_is_chained_assignment_possible(self):
        """
        Check if we are a view, have a cacher, and are of mixed type.
        If so, then force a setitem_copy check.

        Should be called just near setting a value

        Will return a boolean if it we are a view and are cached, but a
        single-dtype meaning that the cacher should be updated following
        setting.
        """
        if self._is_view and self._is_cached:
            ref = self._get_cacher()
            if ref is not None and ref._is_mixed_type:
                self._check_setitem_copy(stacklevel=4, t='referant',
                                         force=True)
            return True
        elif self._is_copy:
            self._check_setitem_copy(stacklevel=4, t='referant')
        return False

    def _check_setitem_copy(self, stacklevel=4, t='setting', force=False):
        """

        Parameters
        ----------
        stacklevel : integer, default 4
           the level to show of the stack when the error is output
        t : string, the type of setting error
        force : boolean, default False
           if True, then force showing an error

        validate if we are doing a settitem on a chained copy.

        If you call this function, be sure to set the stacklevel such that the
        user will see the error *at the level of setting*

        It is technically possible to figure out that we are setting on
        a copy even WITH a multi-dtyped pandas object. In other words, some
        blocks may be views while other are not. Currently _is_view will ALWAYS
        return False for multi-blocks to avoid having to handle this case.

        df = DataFrame(np.arange(0,9), columns=['count'])
        df['group'] = 'b'

        # This technically need not raise SettingWithCopy if both are view
        # (which is not # generally guaranteed but is usually True.  However,
        # this is in general not a good practice and we recommend using .loc.
        df.iloc[0:5]['group'] = 'a'

        """

        if force or self._is_copy:

            value = config.get_option('mode.chained_assignment')
            if value is None:
                return

            # see if the copy is not actually referred; if so, then dissolve
            # the copy weakref
            try:
                gc.collect(2)
                if not gc.get_referents(self._is_copy()):
                    self._is_copy = None
                    return
            except Exception:
                pass

            # we might be a false positive
            try:
                if self._is_copy().shape == self.shape:
                    self._is_copy = None
                    return
            except Exception:
                pass

            # a custom message
            if isinstance(self._is_copy, string_types):
                t = self._is_copy

            elif t == 'referant':
                t = ("\n"
                     "A value is trying to be set on a copy of a slice from a "
                     "DataFrame\n\n"
                     "See the caveats in the documentation: "
                     "http://pandas.pydata.org/pandas-docs/stable/"
                     "indexing.html#indexing-view-versus-copy"
                     )

            else:
                t = ("\n"
                     "A value is trying to be set on a copy of a slice from a "
                     "DataFrame.\n"
                     "Try using .loc[row_indexer,col_indexer] = value "
                     "instead\n\nSee the caveats in the documentation: "
                     "http://pandas.pydata.org/pandas-docs/stable/"
                     "indexing.html#indexing-view-versus-copy"
                     )

            if value == 'raise':
                raise com.SettingWithCopyError(t)
            elif value == 'warn':
                warnings.warn(t, com.SettingWithCopyWarning,
                              stacklevel=stacklevel)

    def __delitem__(self, key):
        """
        Delete item
        """
        deleted = False

        maybe_shortcut = False
        if hasattr(self, 'columns') and isinstance(self.columns, MultiIndex):
            try:
                maybe_shortcut = key not in self.columns._engine
            except TypeError:
                pass

        if maybe_shortcut:
            # Allow shorthand to delete all columns whose first len(key)
            # elements match key:
            if not isinstance(key, tuple):
                key = (key, )
            for col in self.columns:
                if isinstance(col, tuple) and col[:len(key)] == key:
                    del self[col]
                    deleted = True
        if not deleted:
            # If the above loop ran and didn't delete anything because
            # there was no match, this call should raise the appropriate
            # exception:
            self._data.delete(key)

        # delete from the caches
        try:
            del self._item_cache[key]
        except KeyError:
            pass

    _shared_docs['_take'] = """
        Return the elements in the given *positional* indices along an axis.

        This means that we are not indexing according to actual values in
        the index attribute of the object. We are indexing according to the
        actual position of the element in the object.

        This is the internal version of ``.take()`` and will contain a wider
        selection of parameters useful for internal use but not as suitable
        for public usage.

        Parameters
        ----------
        indices : array-like
            An array of ints indicating which positions to take.
        axis : int, default 0
            The axis on which to select elements. "0" means that we are
            selecting rows, "1" means that we are selecting columns, etc.
        is_copy : bool, default True
            Whether to return a copy of the original object or not.

        Returns
        -------
        taken : type of caller
            An array-like containing the elements taken from the object.

        See Also
        --------
        numpy.ndarray.take
        numpy.take
        """

    @Appender(_shared_docs['_take'])
    def _take(self, indices, axis=0, is_copy=True):
        self._consolidate_inplace()

        new_data = self._data.take(indices,
                                   axis=self._get_block_manager_axis(axis),
                                   verify=True)
        result = self._constructor(new_data).__finalize__(self)

        # Maybe set copy if we didn't actually change the index.
        if is_copy:
            if not result._get_axis(axis).equals(self._get_axis(axis)):
                result._set_is_copy(self)

        return result

    _shared_docs['take'] = """
        Return the elements in the given *positional* indices along an axis.

        This means that we are not indexing according to actual values in
        the index attribute of the object. We are indexing according to the
        actual position of the element in the object.

        Parameters
        ----------
        indices : array-like
            An array of ints indicating which positions to take.
        axis : {0 or 'index', 1 or 'columns', None}, default 0
            The axis on which to select elements. ``0`` means that we are
            selecting rows, ``1`` means that we are selecting columns.
        convert : bool, default True
            Whether to convert negative indices into positive ones.
            For example, ``-1`` would map to the ``len(axis) - 1``.
            The conversions are similar to the behavior of indexing a
            regular Python list.

            .. deprecated:: 0.21.0
               In the future, negative indices will always be converted.

        is_copy : bool, default True
            Whether to return a copy of the original object or not.
        **kwargs
            For compatibility with :meth:`numpy.take`. Has no effect on the
            output.

        Returns
        -------
        taken : type of caller
            An array-like containing the elements taken from the object.

        See Also
        --------
        DataFrame.loc : Select a subset of a DataFrame by labels.
        DataFrame.iloc : Select a subset of a DataFrame by positions.
        numpy.take : Take elements from an array along an axis.

        Examples
        --------
        >>> df = pd.DataFrame([('falcon', 'bird',    389.0),
        ...                    ('parrot', 'bird',     24.0),
        ...                    ('lion',   'mammal',   80.5),
        ...                    ('monkey', 'mammal', np.nan)],
        ...                    columns=['name', 'class', 'max_speed'],
        ...                    index=[0, 2, 3, 1])
        >>> df
             name   class  max_speed
        0  falcon    bird      389.0
        2  parrot    bird       24.0
        3    lion  mammal       80.5
        1  monkey  mammal        NaN

        Take elements at positions 0 and 3 along the axis 0 (default).

        Note how the actual indices selected (0 and 1) do not correspond to
        our selected indices 0 and 3. That's because we are selecting the 0th
        and 3rd rows, not rows whose indices equal 0 and 3.

        >>> df.take([0, 3])
             name   class  max_speed
        0  falcon    bird      389.0
        1  monkey  mammal        NaN

        Take elements at indices 1 and 2 along the axis 1 (column selection).

        >>> df.take([1, 2], axis=1)
            class  max_speed
        0    bird      389.0
        2    bird       24.0
        3  mammal       80.5
        1  mammal        NaN

        We may take elements using negative integers for positive indices,
        starting from the end of the object, just like with Python lists.

        >>> df.take([-1, -2])
             name   class  max_speed
        1  monkey  mammal        NaN
        3    lion  mammal       80.5
        """

    @Appender(_shared_docs['take'])
    def take(self, indices, axis=0, convert=None, is_copy=True, **kwargs):
        if convert is not None:
            msg = ("The 'convert' parameter is deprecated "
                   "and will be removed in a future version.")
            warnings.warn(msg, FutureWarning, stacklevel=2)

        nv.validate_take(tuple(), kwargs)
        return self._take(indices, axis=axis, is_copy=is_copy)

    def xs(self, key, axis=0, level=None, drop_level=True):
        """
        Returns a cross-section (row(s) or column(s)) from the
        Series/DataFrame. Defaults to cross-section on the rows (axis=0).

        Parameters
        ----------
        key : object
            Some label contained in the index, or partially in a MultiIndex
        axis : int, default 0
            Axis to retrieve cross-section on
        level : object, defaults to first n levels (n=1 or len(key))
            In case of a key partially contained in a MultiIndex, indicate
            which levels are used. Levels can be referred by label or position.
        drop_level : boolean, default True
            If False, returns object with same levels as self.

        Examples
        --------
        >>> df
           A  B  C
        a  4  5  2
        b  4  0  9
        c  9  7  3
        >>> df.xs('a')
        A    4
        B    5
        C    2
        Name: a
        >>> df.xs('C', axis=1)
        a    2
        b    9
        c    3
        Name: C

        >>> df
                            A  B  C  D
        first second third
        bar   one    1      4  1  8  9
              two    1      7  5  5  0
        baz   one    1      6  6  8  0
              three  2      5  3  5  3
        >>> df.xs(('baz', 'three'))
               A  B  C  D
        third
        2      5  3  5  3
        >>> df.xs('one', level=1)
                     A  B  C  D
        first third
        bar   1      4  1  8  9
        baz   1      6  6  8  0
        >>> df.xs(('baz', 2), level=[0, 'third'])
                A  B  C  D
        second
        three   5  3  5  3

        Returns
        -------
        xs : Series or DataFrame

        Notes
        -----
        xs is only for getting, not setting values.

        MultiIndex Slicers is a generic way to get/set values on any level or
        levels.  It is a superset of xs functionality, see
        :ref:`MultiIndex Slicers <advanced.mi_slicers>`

        """
        axis = self._get_axis_number(axis)
        labels = self._get_axis(axis)
        if level is not None:
            loc, new_ax = labels.get_loc_level(key, level=level,
                                               drop_level=drop_level)

            # create the tuple of the indexer
            indexer = [slice(None)] * self.ndim
            indexer[axis] = loc
            indexer = tuple(indexer)

            result = self.iloc[indexer]
            setattr(result, result._get_axis_name(axis), new_ax)
            return result

        if axis == 1:
            return self[key]

        self._consolidate_inplace()

        index = self.index
        if isinstance(index, MultiIndex):
            loc, new_index = self.index.get_loc_level(key,
                                                      drop_level=drop_level)
        else:
            loc = self.index.get_loc(key)

            if isinstance(loc, np.ndarray):
                if loc.dtype == np.bool_:
                    inds, = loc.nonzero()
                    return self._take(inds, axis=axis)
                else:
                    return self._take(loc, axis=axis)

            if not is_scalar(loc):
                new_index = self.index[loc]

        if is_scalar(loc):
            new_values = self._data.fast_xs(loc)

            # may need to box a datelike-scalar
            #
            # if we encounter an array-like and we only have 1 dim
            # that means that their are list/ndarrays inside the Series!
            # so just return them (GH 6394)
            if not is_list_like(new_values) or self.ndim == 1:
                return com._maybe_box_datetimelike(new_values)

            result = self._constructor_sliced(
                new_values, index=self.columns,
                name=self.index[loc], dtype=new_values.dtype)

        else:
            result = self.iloc[loc]
            result.index = new_index

        # this could be a view
        # but only in a single-dtyped view slicable case
        result._set_is_copy(self, copy=not result._is_view)
        return result

    _xs = xs

    def select(self, crit, axis=0):
        """Return data corresponding to axis labels matching criteria

        .. deprecated:: 0.21.0
            Use df.loc[df.index.map(crit)] to select via labels

        Parameters
        ----------
        crit : function
            To be called on each index (label). Should return True or False
        axis : int

        Returns
        -------
        selection : type of caller
        """
        warnings.warn("'select' is deprecated and will be removed in a "
                      "future release. You can use "
                      ".loc[labels.map(crit)] as a replacement",
                      FutureWarning, stacklevel=2)

        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        axis_values = self._get_axis(axis)

        if len(axis_values) > 0:
            new_axis = axis_values[
                np.asarray([bool(crit(label)) for label in axis_values])]
        else:
            new_axis = axis_values

        return self.reindex(**{axis_name: new_axis})

    def reindex_like(self, other, method=None, copy=True, limit=None,
                     tolerance=None):
        """Return an object with matching indices to myself.

        Parameters
        ----------
        other : Object
        method : string or None
        copy : boolean, default True
        limit : int, default None
            Maximum number of consecutive labels to fill for inexact matches.
        tolerance : optional
            Maximum distance between labels of the other object and this
            object for inexact matches. Can be list-like.

            .. versionadded:: 0.21.0 (list-like tolerance)

        Notes
        -----
        Like calling s.reindex(index=other.index, columns=other.columns,
                               method=...)

        Returns
        -------
        reindexed : same as input
        """
        d = other._construct_axes_dict(axes=self._AXIS_ORDERS, method=method,
                                       copy=copy, limit=limit,
                                       tolerance=tolerance)

        return self.reindex(**d)

    def drop(self, labels=None, axis=0, index=None, columns=None, level=None,
             inplace=False, errors='raise'):

        inplace = validate_bool_kwarg(inplace, 'inplace')

        if labels is not None:
            if index is not None or columns is not None:
                raise ValueError("Cannot specify both 'labels' and "
                                 "'index'/'columns'")
            axis_name = self._get_axis_name(axis)
            axes = {axis_name: labels}
        elif index is not None or columns is not None:
            axes, _ = self._construct_axes_from_arguments((index, columns), {})
        else:
            raise ValueError("Need to specify at least one of 'labels', "
                             "'index' or 'columns'")

        obj = self

        for axis, labels in axes.items():
            if labels is not None:
                obj = obj._drop_axis(labels, axis, level=level, errors=errors)

        if inplace:
            self._update_inplace(obj)
        else:
            return obj

    def _drop_axis(self, labels, axis, level=None, errors='raise'):
        """
        Drop labels from specified axis. Used in the ``drop`` method
        internally.

        Parameters
        ----------
        labels : single label or list-like
        axis : int or axis name
        level : int or level name, default None
            For MultiIndex
        errors : {'ignore', 'raise'}, default 'raise'
            If 'ignore', suppress error and existing labels are dropped.

        """
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        axis, axis_ = self._get_axis(axis), axis

        if axis.is_unique:
            if level is not None:
                if not isinstance(axis, MultiIndex):
                    raise AssertionError('axis must be a MultiIndex')
                new_axis = axis.drop(labels, level=level, errors=errors)
            else:
                new_axis = axis.drop(labels, errors=errors)
            dropped = self.reindex(**{axis_name: new_axis})
            try:
                dropped.axes[axis_].set_names(axis.names, inplace=True)
            except AttributeError:
                pass
            result = dropped

        else:
            labels = _ensure_object(com._index_labels_to_array(labels))
            if level is not None:
                if not isinstance(axis, MultiIndex):
                    raise AssertionError('axis must be a MultiIndex')
                indexer = ~axis.get_level_values(level).isin(labels)
            else:
                indexer = ~axis.isin(labels)

            if errors == 'raise' and indexer.all():
                raise KeyError('{} not found in axis'.format(labels))

            slicer = [slice(None)] * self.ndim
            slicer[self._get_axis_number(axis_name)] = indexer

            result = self.loc[tuple(slicer)]

        return result

    def _update_inplace(self, result, verify_is_copy=True):
        """
        Replace self internals with result.

        Parameters
        ----------
        verify_is_copy : boolean, default True
            provide is_copy checks

        """
        # NOTE: This does *not* call __finalize__ and that's an explicit
        # decision that we may revisit in the future.

        self._reset_cache()
        self._clear_item_cache()
        self._data = getattr(result, '_data', result)
        self._maybe_update_cacher(verify_is_copy=verify_is_copy)

    def add_prefix(self, prefix):
        """
        Prefix labels with string `prefix`.

        For Series, the row labels are prefixed.
        For DataFrame, the column labels are prefixed.

        Parameters
        ----------
        prefix : str
            The string to add before each label.

        Returns
        -------
        Series or DataFrame
            New Series or DataFrame with updated labels.

        See Also
        --------
        Series.add_suffix: Suffix row labels with string `suffix`.
        DataFrame.add_suffix: Suffix column labels with string `suffix`.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s
        0    1
        1    2
        2    3
        3    4
        dtype: int64

        >>> s.add_prefix('item_')
        item_0    1
        item_1    2
        item_2    3
        item_3    4
        dtype: int64

        >>> df = pd.DataFrame({'A': [1, 2, 3, 4],  'B': [3, 4, 5, 6]})
        >>> df
           A  B
        0  1  3
        1  2  4
        2  3  5
        3  4  6

        >>> df.add_prefix('col_')
             col_A  col_B
        0       1       3
        1       2       4
        2       3       5
        3       4       6
        """
        new_data = self._data.add_prefix(prefix)
        return self._constructor(new_data).__finalize__(self)

    def add_suffix(self, suffix):
        """
        Suffix labels with string `suffix`.

        For Series, the row labels are suffixed.
        For DataFrame, the column labels are suffixed.

        Parameters
        ----------
        suffix : str
            The string to add after each label.

        Returns
        -------
        Series or DataFrame
            New Series or DataFrame with updated labels.

        See Also
        --------
        Series.add_prefix: Prefix row labels with string `prefix`.
        DataFrame.add_prefix: Prefix column labels with string `prefix`.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])
        >>> s
        0    1
        1    2
        2    3
        3    4
        dtype: int64

        >>> s.add_suffix('_item')
        0_item    1
        1_item    2
        2_item    3
        3_item    4
        dtype: int64

        >>> df = pd.DataFrame({'A': [1, 2, 3, 4],  'B': [3, 4, 5, 6]})
        >>> df
           A  B
        0  1  3
        1  2  4
        2  3  5
        3  4  6

        >>> df.add_suffix('_col')
             A_col  B_col
        0       1       3
        1       2       4
        2       3       5
        3       4       6
        """
        new_data = self._data.add_suffix(suffix)
        return self._constructor(new_data).__finalize__(self)

    _shared_docs['sort_values'] = """
        Sort by the values along either axis

        Parameters
        ----------%(optional_by)s
        axis : %(axes_single_arg)s, default 0
             Axis to be sorted
        ascending : bool or list of bool, default True
             Sort ascending vs. descending. Specify list for multiple sort
             orders.  If this is a list of bools, must match the length of
             the by.
        inplace : bool, default False
             if True, perform operation in-place
        kind : {'quicksort', 'mergesort', 'heapsort'}, default 'quicksort'
             Choice of sorting algorithm. See also ndarray.np.sort for more
             information.  `mergesort` is the only stable algorithm. For
             DataFrames, this option is only applied when sorting on a single
             column or label.
        na_position : {'first', 'last'}, default 'last'
             `first` puts NaNs at the beginning, `last` puts NaNs at the end

        Returns
        -------
        sorted_obj : %(klass)s

        Examples
        --------
        >>> df = pd.DataFrame({
        ...     'col1' : ['A', 'A', 'B', np.nan, 'D', 'C'],
        ...     'col2' : [2, 1, 9, 8, 7, 4],
        ...     'col3': [0, 1, 9, 4, 2, 3],
        ... })
        >>> df
            col1 col2 col3
        0   A    2    0
        1   A    1    1
        2   B    9    9
        3   NaN  8    4
        4   D    7    2
        5   C    4    3

        Sort by col1

        >>> df.sort_values(by=['col1'])
            col1 col2 col3
        0   A    2    0
        1   A    1    1
        2   B    9    9
        5   C    4    3
        4   D    7    2
        3   NaN  8    4

        Sort by multiple columns

        >>> df.sort_values(by=['col1', 'col2'])
            col1 col2 col3
        1   A    1    1
        0   A    2    0
        2   B    9    9
        5   C    4    3
        4   D    7    2
        3   NaN  8    4

        Sort Descending

        >>> df.sort_values(by='col1', ascending=False)
            col1 col2 col3
        4   D    7    2
        5   C    4    3
        2   B    9    9
        0   A    2    0
        1   A    1    1
        3   NaN  8    4

        Putting NAs first

        >>> df.sort_values(by='col1', ascending=False, na_position='first')
            col1 col2 col3
        3   NaN  8    4
        4   D    7    2
        5   C    4    3
        2   B    9    9
        0   A    2    0
        1   A    1    1
        """

    def sort_values(self, by=None, axis=0, ascending=True, inplace=False,
                    kind='quicksort', na_position='last'):
        """
        NOT IMPLEMENTED: do not call this method, as sorting values is not
        supported for Panel objects and will raise an error.
        """
        raise NotImplementedError("sort_values has not been implemented "
                                  "on Panel or Panel4D objects.")

    _shared_docs['sort_index'] = """
        Sort object by labels (along an axis)

        Parameters
        ----------
        axis : %(axes)s to direct sorting
        level : int or level name or list of ints or list of level names
            if not None, sort on values in specified index level(s)
        ascending : boolean, default True
            Sort ascending vs. descending
        inplace : bool, default False
            if True, perform operation in-place
        kind : {'quicksort', 'mergesort', 'heapsort'}, default 'quicksort'
             Choice of sorting algorithm. See also ndarray.np.sort for more
             information.  `mergesort` is the only stable algorithm. For
             DataFrames, this option is only applied when sorting on a single
             column or label.
        na_position : {'first', 'last'}, default 'last'
             `first` puts NaNs at the beginning, `last` puts NaNs at the end.
             Not implemented for MultiIndex.
        sort_remaining : bool, default True
            if true and sorting by level and index is multilevel, sort by other
            levels too (in order) after sorting by specified level

        Returns
        -------
        sorted_obj : %(klass)s
        """

    @Appender(_shared_docs['sort_index'] % dict(axes="axes", klass="NDFrame"))
    def sort_index(self, axis=0, level=None, ascending=True, inplace=False,
                   kind='quicksort', na_position='last', sort_remaining=True):
        inplace = validate_bool_kwarg(inplace, 'inplace')
        axis = self._get_axis_number(axis)
        axis_name = self._get_axis_name(axis)
        labels = self._get_axis(axis)

        if level is not None:
            raise NotImplementedError("level is not implemented")
        if inplace:
            raise NotImplementedError("inplace is not implemented")

        sort_index = labels.argsort()
        if not ascending:
            sort_index = sort_index[::-1]

        new_axis = labels.take(sort_index)
        return self.reindex(**{axis_name: new_axis})

    _shared_docs['reindex'] = """
        Conform %(klass)s to new index with optional filling logic, placing
        NA/NaN in locations having no value in the previous index. A new object
        is produced unless the new index is equivalent to the current one and
        copy=False

        Parameters
        ----------
        %(optional_labels)s
        %(axes)s : array-like, optional (should be specified using keywords)
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        %(optional_axis)s
        method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}, optional
            method to use for filling holes in reindexed DataFrame.
            Please note: this is only applicable to DataFrames/Series with a
            monotonically increasing/decreasing index.

            * default: don't fill gaps
            * pad / ffill: propagate last valid observation forward to next
              valid
            * backfill / bfill: use next valid observation to fill gap
            * nearest: use nearest valid observations to fill gap

        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        limit : int, default None
            Maximum number of consecutive elements to forward or backward fill
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations most
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.

            Tolerance may be a scalar value, which applies the same tolerance
            to all values, or list-like, which applies variable tolerance per
            element. List-like includes list, tuple, array, Series, and must be
            the same size as the index and its dtype must exactly match the
            index's type.

            .. versionadded:: 0.21.0 (list-like tolerance)

        Examples
        --------

        ``DataFrame.reindex`` supports two calling conventions

        * ``(index=index_labels, columns=column_labels, ...)``
        * ``(labels, axis={'index', 'columns'}, ...)``

        We *highly* recommend using keyword arguments to clarify your
        intent.

        Create a dataframe with some fictional data.

        >>> index = ['Firefox', 'Chrome', 'Safari', 'IE10', 'Konqueror']
        >>> df = pd.DataFrame({
        ...      'http_status': [200,200,404,404,301],
        ...      'response_time': [0.04, 0.02, 0.07, 0.08, 1.0]},
        ...       index=index)
        >>> df
                   http_status  response_time
        Firefox            200           0.04
        Chrome             200           0.02
        Safari             404           0.07
        IE10               404           0.08
        Konqueror          301           1.00

        Create a new index and reindex the dataframe. By default
        values in the new index that do not have corresponding
        records in the dataframe are assigned ``NaN``.

        >>> new_index= ['Safari', 'Iceweasel', 'Comodo Dragon', 'IE10',
        ...             'Chrome']
        >>> df.reindex(new_index)
                       http_status  response_time
        Safari               404.0           0.07
        Iceweasel              NaN            NaN
        Comodo Dragon          NaN            NaN
        IE10                 404.0           0.08
        Chrome               200.0           0.02

        We can fill in the missing values by passing a value to
        the keyword ``fill_value``. Because the index is not monotonically
        increasing or decreasing, we cannot use arguments to the keyword
        ``method`` to fill the ``NaN`` values.

        >>> df.reindex(new_index, fill_value=0)
                       http_status  response_time
        Safari                 404           0.07
        Iceweasel                0           0.00
        Comodo Dragon            0           0.00
        IE10                   404           0.08
        Chrome                 200           0.02

        >>> df.reindex(new_index, fill_value='missing')
                      http_status response_time
        Safari                404          0.07
        Iceweasel         missing       missing
        Comodo Dragon     missing       missing
        IE10                  404          0.08
        Chrome                200          0.02

        We can also reindex the columns.

        >>> df.reindex(columns=['http_status', 'user_agent'])
                   http_status  user_agent
        Firefox            200         NaN
        Chrome             200         NaN
        Safari             404         NaN
        IE10               404         NaN
        Konqueror          301         NaN

        Or we can use "axis-style" keyword arguments

        >>> df.reindex(['http_status', 'user_agent'], axis="columns")
                   http_status  user_agent
        Firefox            200         NaN
        Chrome             200         NaN
        Safari             404         NaN
        IE10               404         NaN
        Konqueror          301         NaN

        To further illustrate the filling functionality in
        ``reindex``, we will create a dataframe with a
        monotonically increasing index (for example, a sequence
        of dates).

        >>> date_index = pd.date_range('1/1/2010', periods=6, freq='D')
        >>> df2 = pd.DataFrame({"prices": [100, 101, np.nan, 100, 89, 88]},
        ...                    index=date_index)
        >>> df2
                    prices
        2010-01-01     100
        2010-01-02     101
        2010-01-03     NaN
        2010-01-04     100
        2010-01-05      89
        2010-01-06      88

        Suppose we decide to expand the dataframe to cover a wider
        date range.

        >>> date_index2 = pd.date_range('12/29/2009', periods=10, freq='D')
        >>> df2.reindex(date_index2)
                    prices
        2009-12-29     NaN
        2009-12-30     NaN
        2009-12-31     NaN
        2010-01-01     100
        2010-01-02     101
        2010-01-03     NaN
        2010-01-04     100
        2010-01-05      89
        2010-01-06      88
        2010-01-07     NaN

        The index entries that did not have a value in the original data frame
        (for example, '2009-12-29') are by default filled with ``NaN``.
        If desired, we can fill in the missing values using one of several
        options.

        For example, to backpropagate the last valid value to fill the ``NaN``
        values, pass ``bfill`` as an argument to the ``method`` keyword.

        >>> df2.reindex(date_index2, method='bfill')
                    prices
        2009-12-29     100
        2009-12-30     100
        2009-12-31     100
        2010-01-01     100
        2010-01-02     101
        2010-01-03     NaN
        2010-01-04     100
        2010-01-05      89
        2010-01-06      88
        2010-01-07     NaN

        Please note that the ``NaN`` value present in the original dataframe
        (at index value 2010-01-03) will not be filled by any of the
        value propagation schemes. This is because filling while reindexing
        does not look at dataframe values, but only compares the original and
        desired indexes. If you do want to fill in the ``NaN`` values present
        in the original dataframe, use the ``fillna()`` method.

        See the :ref:`user guide <basics.reindexing>` for more.

        Returns
        -------
        reindexed : %(klass)s
        """

    # TODO: Decide if we care about having different examples for different
    #       kinds

    @Appender(_shared_docs['reindex'] % dict(axes="axes", klass="NDFrame",
                                             optional_labels="",
                                             optional_axis=""))
    def reindex(self, *args, **kwargs):

        # construct the args
        axes, kwargs = self._construct_axes_from_arguments(args, kwargs)
        method = missing.clean_reindex_fill_method(kwargs.pop('method', None))
        level = kwargs.pop('level', None)
        copy = kwargs.pop('copy', True)
        limit = kwargs.pop('limit', None)
        tolerance = kwargs.pop('tolerance', None)
        fill_value = kwargs.pop('fill_value', None)

        # Series.reindex doesn't use / need the axis kwarg
        # We pop and ignore it here, to make writing Series/Frame generic code
        # easier
        kwargs.pop("axis", None)

        if kwargs:
            raise TypeError('reindex() got an unexpected keyword '
                            'argument "{0}"'.format(list(kwargs.keys())[0]))

        self._consolidate_inplace()

        # if all axes that are requested to reindex are equal, then only copy
        # if indicated must have index names equal here as well as values
        if all(self._get_axis(axis).identical(ax)
               for axis, ax in axes.items() if ax is not None):
            if copy:
                return self.copy()
            return self

        # check if we are a multi reindex
        if self._needs_reindex_multi(axes, method, level):
            try:
                return self._reindex_multi(axes, copy, fill_value)
            except Exception:
                pass

        # perform the reindex on the axes
        return self._reindex_axes(axes, level, limit, tolerance, method,
                                  fill_value, copy).__finalize__(self)

    def _reindex_axes(self, axes, level, limit, tolerance, method, fill_value,
                      copy):
        """Perform the reindex for all the axes."""
        obj = self
        for a in self._AXIS_ORDERS:
            labels = axes[a]
            if labels is None:
                continue

            ax = self._get_axis(a)
            new_index, indexer = ax.reindex(labels, level=level, limit=limit,
                                            tolerance=tolerance, method=method)

            axis = self._get_axis_number(a)
            obj = obj._reindex_with_indexers({axis: [new_index, indexer]},
                                             fill_value=fill_value,
                                             copy=copy, allow_dups=False)

        return obj

    def _needs_reindex_multi(self, axes, method, level):
        """Check if we do need a multi reindex."""
        return ((com._count_not_none(*axes.values()) == self._AXIS_LEN) and
                method is None and level is None and not self._is_mixed_type)

    def _reindex_multi(self, axes, copy, fill_value):
        return NotImplemented

    _shared_docs[
        'reindex_axis'] = ("""Conform input object to new index with optional
        filling logic, placing NA/NaN in locations having no value in the
        previous index. A new object is produced unless the new index is
        equivalent to the current one and copy=False

        Parameters
        ----------
        labels : array-like
            New labels / index to conform to. Preferably an Index object to
            avoid duplicating data
        axis : %(axes_single_arg)s
        method : {None, 'backfill'/'bfill', 'pad'/'ffill', 'nearest'}, optional
            Method to use for filling holes in reindexed DataFrame:

            * default: don't fill gaps
            * pad / ffill: propagate last valid observation forward to next
              valid
            * backfill / bfill: use next valid observation to fill gap
            * nearest: use nearest valid observations to fill gap

        copy : boolean, default True
            Return a new object, even if the passed indexes are the same
        level : int or name
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        limit : int, default None
            Maximum number of consecutive elements to forward or backward fill
        tolerance : optional
            Maximum distance between original and new labels for inexact
            matches. The values of the index at the matching locations most
            satisfy the equation ``abs(index[indexer] - target) <= tolerance``.

            Tolerance may be a scalar value, which applies the same tolerance
            to all values, or list-like, which applies variable tolerance per
            element. List-like includes list, tuple, array, Series, and must be
            the same size as the index and its dtype must exactly match the
            index's type.

            .. versionadded:: 0.21.0 (list-like tolerance)

        Examples
        --------
        >>> df.reindex_axis(['A', 'B', 'C'], axis=1)

        See Also
        --------
        reindex, reindex_like

        Returns
        -------
        reindexed : %(klass)s
        """)

    @Appender(_shared_docs['reindex_axis'] % _shared_doc_kwargs)
    def reindex_axis(self, labels, axis=0, method=None, level=None, copy=True,
                     limit=None, fill_value=None):
        msg = ("'.reindex_axis' is deprecated and will be removed in a future "
               "version. Use '.reindex' instead.")
        self._consolidate_inplace()

        axis_name = self._get_axis_name(axis)
        axis_values = self._get_axis(axis_name)
        method = missing.clean_reindex_fill_method(method)
        warnings.warn(msg, FutureWarning, stacklevel=3)
        new_index, indexer = axis_values.reindex(labels, method, level,
                                                 limit=limit)
        return self._reindex_with_indexers({axis: [new_index, indexer]},
                                           fill_value=fill_value, copy=copy)

    def _reindex_with_indexers(self, reindexers, fill_value=None, copy=False,
                               allow_dups=False):
        """allow_dups indicates an internal call here """

        # reindex doing multiple operations on different axes if indicated
        new_data = self._data
        for axis in sorted(reindexers.keys()):
            index, indexer = reindexers[axis]
            baxis = self._get_block_manager_axis(axis)

            if index is None:
                continue

            index = _ensure_index(index)
            if indexer is not None:
                indexer = _ensure_int64(indexer)

            # TODO: speed up on homogeneous DataFrame objects
            new_data = new_data.reindex_indexer(index, indexer, axis=baxis,
                                                fill_value=fill_value,
                                                allow_dups=allow_dups,
                                                copy=copy)

        if copy and new_data is self._data:
            new_data = new_data.copy()

        return self._constructor(new_data).__finalize__(self)

    def _reindex_axis(self, new_index, fill_method, axis, copy):
        new_data = self._data.reindex_axis(new_index, axis=axis,
                                           method=fill_method, copy=copy)

        if new_data is self._data and not copy:
            return self
        else:
            return self._constructor(new_data).__finalize__(self)

    def filter(self, items=None, like=None, regex=None, axis=None):
        """
        Subset rows or columns of dataframe according to labels in
        the specified index.

        Note that this routine does not filter a dataframe on its
        contents. The filter is applied to the labels of the index.

        Parameters
        ----------
        items : list-like
            List of info axis to restrict to (must not all be present)
        like : string
            Keep info axis where "arg in col == True"
        regex : string (regular expression)
            Keep info axis with re.search(regex, col) == True
        axis : int or string axis name
            The axis to filter on.  By default this is the info axis,
            'index' for Series, 'columns' for DataFrame

        Returns
        -------
        same type as input object

        Examples
        --------
        >>> df
        one  two  three
        mouse     1    2      3
        rabbit    4    5      6

        >>> # select columns by name
        >>> df.filter(items=['one', 'three'])
        one  three
        mouse     1      3
        rabbit    4      6

        >>> # select columns by regular expression
        >>> df.filter(regex='e$', axis=1)
        one  three
        mouse     1      3
        rabbit    4      6

        >>> # select rows containing 'bbi'
        >>> df.filter(like='bbi', axis=0)
        one  two  three
        rabbit    4    5      6

        See Also
        --------
        pandas.DataFrame.loc

        Notes
        -----
        The ``items``, ``like``, and ``regex`` parameters are
        enforced to be mutually exclusive.

        ``axis`` defaults to the info axis that is used when indexing
        with ``[]``.
        """
        import re

        nkw = com._count_not_none(items, like, regex)
        if nkw > 1:
            raise TypeError('Keyword arguments `items`, `like`, or `regex` '
                            'are mutually exclusive')

        if axis is None:
            axis = self._info_axis_name
        labels = self._get_axis(axis)

        if items is not None:
            name = self._get_axis_name(axis)
            return self.reindex(
                **{name: [r for r in items if r in labels]})
        elif like:
            def f(x):
                return like in to_str(x)
            values = labels.map(f)
            return self.loc(axis=axis)[values]
        elif regex:
            def f(x):
                return matcher.search(to_str(x)) is not None
            matcher = re.compile(regex)
            values = labels.map(f)
            return self.loc(axis=axis)[values]
        else:
            raise TypeError('Must pass either `items`, `like`, or `regex`')

    def head(self, n=5):
        """
        Return the first `n` rows.

        This function returns the first `n` rows for the object based
        on position. It is useful for quickly testing if your object
        has the right type of data in it.

        Parameters
        ----------
        n : int, default 5
            Number of rows to select.

        Returns
        -------
        obj_head : type of caller
            The first `n` rows of the caller object.

        See Also
        --------
        pandas.DataFrame.tail: Returns the last `n` rows.

        Examples
        --------
        >>> df = pd.DataFrame({'animal':['alligator', 'bee', 'falcon', 'lion',
        ...                    'monkey', 'parrot', 'shark', 'whale', 'zebra']})
        >>> df
              animal
        0  alligator
        1        bee
        2     falcon
        3       lion
        4     monkey
        5     parrot
        6      shark
        7      whale
        8      zebra

        Viewing the first 5 lines

        >>> df.head()
              animal
        0  alligator
        1        bee
        2     falcon
        3       lion
        4     monkey

        Viewing the first `n` lines (three in this case)

        >>> df.head(3)
              animal
        0  alligator
        1        bee
        2     falcon
        """

        return self.iloc[:n]

    def tail(self, n=5):
        """
        Return the last `n` rows.

        This function returns last `n` rows from the object based on
        position. It is useful for quickly verifying data, for example,
        after sorting or appending rows.

        Parameters
        ----------
        n : int, default 5
            Number of rows to select.

        Returns
        -------
        type of caller
            The last `n` rows of the caller object.

        See Also
        --------
        pandas.DataFrame.head : The first `n` rows of the caller object.

        Examples
        --------
        >>> df = pd.DataFrame({'animal':['alligator', 'bee', 'falcon', 'lion',
        ...                    'monkey', 'parrot', 'shark', 'whale', 'zebra']})
        >>> df
              animal
        0  alligator
        1        bee
        2     falcon
        3       lion
        4     monkey
        5     parrot
        6      shark
        7      whale
        8      zebra

        Viewing the last 5 lines

        >>> df.tail()
           animal
        4  monkey
        5  parrot
        6   shark
        7   whale
        8   zebra

        Viewing the last `n` lines (three in this case)

        >>> df.tail(3)
          animal
        6  shark
        7  whale
        8  zebra
        """

        if n == 0:
            return self.iloc[0:0]
        return self.iloc[-n:]

    def sample(self, n=None, frac=None, replace=False, weights=None,
               random_state=None, axis=None):
        """
        Return a random sample of items from an axis of object.

        You can use `random_state` for reproducibility.

        Parameters
        ----------
        n : int, optional
            Number of items from axis to return. Cannot be used with `frac`.
            Default = 1 if `frac` = None.
        frac : float, optional
            Fraction of axis items to return. Cannot be used with `n`.
        replace : boolean, optional
            Sample with or without replacement. Default = False.
        weights : str or ndarray-like, optional
            Default 'None' results in equal probability weighting.
            If passed a Series, will align with target object on index. Index
            values in weights not found in sampled object will be ignored and
            index values in sampled object not in weights will be assigned
            weights of zero.
            If called on a DataFrame, will accept the name of a column
            when axis = 0.
            Unless weights are a Series, weights must be same length as axis
            being sampled.
            If weights do not sum to 1, they will be normalized to sum to 1.
            Missing values in the weights column will be treated as zero.
            inf and -inf values not allowed.
        random_state : int or numpy.random.RandomState, optional
            Seed for the random number generator (if int), or numpy RandomState
            object.
        axis : int or string, optional
            Axis to sample. Accepts axis number or name. Default is stat axis
            for given data type (0 for Series and DataFrames, 1 for Panels).

        Returns
        -------
        A new object of same type as caller.

        Examples
        --------
        Generate an example ``Series`` and ``DataFrame``:

        >>> s = pd.Series(np.random.randn(50))
        >>> s.head()
        0   -0.038497
        1    1.820773
        2   -0.972766
        3   -1.598270
        4   -1.095526
        dtype: float64
        >>> df = pd.DataFrame(np.random.randn(50, 4), columns=list('ABCD'))
        >>> df.head()
                  A         B         C         D
        0  0.016443 -2.318952 -0.566372 -1.028078
        1 -1.051921  0.438836  0.658280 -0.175797
        2 -1.243569 -0.364626 -0.215065  0.057736
        3  1.768216  0.404512 -0.385604 -1.457834
        4  1.072446 -1.137172  0.314194 -0.046661

        Next extract a random sample from both of these objects...

        3 random elements from the ``Series``:

        >>> s.sample(n=3)
        27   -0.994689
        55   -1.049016
        67   -0.224565
        dtype: float64

        And a random 10% of the ``DataFrame`` with replacement:

        >>> df.sample(frac=0.1, replace=True)
                   A         B         C         D
        35  1.981780  0.142106  1.817165 -0.290805
        49 -1.336199 -0.448634 -0.789640  0.217116
        40  0.823173 -0.078816  1.009536  1.015108
        15  1.421154 -0.055301 -1.922594 -0.019696
        6  -0.148339  0.832938  1.787600 -1.383767

        You can use `random state` for reproducibility:

        >>> df.sample(random_state=1)
        A         B         C         D
        37 -2.027662  0.103611  0.237496 -0.165867
        43 -0.259323 -0.583426  1.516140 -0.479118
        12 -1.686325 -0.579510  0.985195 -0.460286
        8   1.167946  0.429082  1.215742 -1.636041
        9   1.197475 -0.864188  1.554031 -1.505264
        """

        if axis is None:
            axis = self._stat_axis_number

        axis = self._get_axis_number(axis)
        axis_length = self.shape[axis]

        # Process random_state argument
        rs = com._random_state(random_state)

        # Check weights for compliance
        if weights is not None:

            # If a series, align with frame
            if isinstance(weights, pd.Series):
                weights = weights.reindex(self.axes[axis])

            # Strings acceptable if a dataframe and axis = 0
            if isinstance(weights, string_types):
                if isinstance(self, pd.DataFrame):
                    if axis == 0:
                        try:
                            weights = self[weights]
                        except KeyError:
                            raise KeyError("String passed to weights not a "
                                           "valid column")
                    else:
                        raise ValueError("Strings can only be passed to "
                                         "weights when sampling from rows on "
                                         "a DataFrame")
                else:
                    raise ValueError("Strings cannot be passed as weights "
                                     "when sampling from a Series or Panel.")

            weights = pd.Series(weights, dtype='float64')

            if len(weights) != axis_length:
                raise ValueError("Weights and axis to be sampled must be of "
                                 "same length")

            if (weights == np.inf).any() or (weights == -np.inf).any():
                raise ValueError("weight vector may not include `inf` values")

            if (weights < 0).any():
                raise ValueError("weight vector many not include negative "
                                 "values")

            # If has nan, set to zero.
            weights = weights.fillna(0)

            # Renormalize if don't sum to 1
            if weights.sum() != 1:
                if weights.sum() != 0:
                    weights = weights / weights.sum()
                else:
                    raise ValueError("Invalid weights: weights sum to zero")

            weights = weights.values

        # If no frac or n, default to n=1.
        if n is None and frac is None:
            n = 1
        elif n is not None and frac is None and n % 1 != 0:
            raise ValueError("Only integers accepted as `n` values")
        elif n is None and frac is not None:
            n = int(round(frac * axis_length))
        elif n is not None and frac is not None:
            raise ValueError('Please enter a value for `frac` OR `n`, not '
                             'both')

        # Check for negative sizes
        if n < 0:
            raise ValueError("A negative number of rows requested. Please "
                             "provide positive value.")

        locs = rs.choice(axis_length, size=n, replace=replace, p=weights)
        return self.take(locs, axis=axis, is_copy=False)

    _shared_docs['pipe'] = (r"""
        Apply func(self, \*args, \*\*kwargs)

        Parameters
        ----------
        func : function
            function to apply to the %(klass)s.
            ``args``, and ``kwargs`` are passed into ``func``.
            Alternatively a ``(callable, data_keyword)`` tuple where
            ``data_keyword`` is a string indicating the keyword of
            ``callable`` that expects the %(klass)s.
        args : iterable, optional
            positional arguments passed into ``func``.
        kwargs : mapping, optional
            a dictionary of keyword arguments passed into ``func``.

        Returns
        -------
        object : the return type of ``func``.

        Notes
        -----

        Use ``.pipe`` when chaining together functions that expect
        Series, DataFrames or GroupBy objects. Instead of writing

        >>> f(g(h(df), arg1=a), arg2=b, arg3=c)

        You can write

        >>> (df.pipe(h)
        ...    .pipe(g, arg1=a)
        ...    .pipe(f, arg2=b, arg3=c)
        ... )

        If you have a function that takes the data as (say) the second
        argument, pass a tuple indicating which keyword expects the
        data. For example, suppose ``f`` takes its data as ``arg2``:

        >>> (df.pipe(h)
        ...    .pipe(g, arg1=a)
        ...    .pipe((f, 'arg2'), arg1=a, arg3=c)
        ...  )

        See Also
        --------
        pandas.DataFrame.apply
        pandas.DataFrame.applymap
        pandas.Series.map
    """)

    @Appender(_shared_docs['pipe'] % _shared_doc_kwargs)
    def pipe(self, func, *args, **kwargs):
        return com._pipe(self, func, *args, **kwargs)

    _shared_docs['aggregate'] = ("""
    Aggregate using one or more operations over the specified axis.

    %(versionadded)s

    Parameters
    ----------
    func : function, string, dictionary, or list of string/functions
        Function to use for aggregating the data. If a function, must either
        work when passed a %(klass)s or when passed to %(klass)s.apply. For
        a DataFrame, can pass a dict, if the keys are DataFrame column names.

        Accepted combinations are:

        - string function name.
        - function.
        - list of functions.
        - dict of column names -> functions (or list of functions).

    %(axis)s
    *args
        Positional arguments to pass to `func`.
    **kwargs
        Keyword arguments to pass to `func`.

    Returns
    -------
    aggregated : %(klass)s

    Notes
    -----
    `agg` is an alias for `aggregate`. Use the alias.

    A passed user-defined-function will be passed a Series for evaluation.
    """)

    _shared_docs['transform'] = ("""
    Call function producing a like-indexed %(klass)s
    and return a %(klass)s with the transformed values

    .. versionadded:: 0.20.0

    Parameters
    ----------
    func : callable, string, dictionary, or list of string/callables
        To apply to column

        Accepted Combinations are:

        - string function name
        - function
        - list of functions
        - dict of column names -> functions (or list of functions)

    Returns
    -------
    transformed : %(klass)s

    Examples
    --------
    >>> df = pd.DataFrame(np.random.randn(10, 3), columns=['A', 'B', 'C'],
    ...                   index=pd.date_range('1/1/2000', periods=10))
    df.iloc[3:7] = np.nan

    >>> df.transform(lambda x: (x - x.mean()) / x.std())
                       A         B         C
    2000-01-01  0.579457  1.236184  0.123424
    2000-01-02  0.370357 -0.605875 -1.231325
    2000-01-03  1.455756 -0.277446  0.288967
    2000-01-04       NaN       NaN       NaN
    2000-01-05       NaN       NaN       NaN
    2000-01-06       NaN       NaN       NaN
    2000-01-07       NaN       NaN       NaN
    2000-01-08 -0.498658  1.274522  1.642524
    2000-01-09 -0.540524 -1.012676 -0.828968
    2000-01-10 -1.366388 -0.614710  0.005378

    See also
    --------
    pandas.%(klass)s.aggregate
    pandas.%(klass)s.apply
    """)

    # ----------------------------------------------------------------------
    # Attribute access

    def __finalize__(self, other, method=None, **kwargs):
        """
        Propagate metadata from other to self.

        Parameters
        ----------
        other : the object from which to get the attributes that we are going
            to propagate
        method : optional, a passed method name ; possibly to take different
            types of propagation actions based on this

        """
        if isinstance(other, NDFrame):
            for name in self._metadata:
                object.__setattr__(self, name, getattr(other, name, None))
        return self

    def __getattr__(self, name):
        """After regular attribute access, try looking up the name
        This allows simpler access to columns for interactive use.
        """

        # Note: obj.x will always call obj.__getattribute__('x') prior to
        # calling obj.__getattr__('x').

        if (name in self._internal_names_set or name in self._metadata or
                name in self._accessors):
            return object.__getattribute__(self, name)
        else:
            if self._info_axis._can_hold_identifiers_and_holds_name(name):
                return self[name]
            return object.__getattribute__(self, name)

    def __setattr__(self, name, value):
        """After regular attribute access, try setting the name
        This allows simpler access to columns for interactive use.
        """

        # first try regular attribute access via __getattribute__, so that
        # e.g. ``obj.x`` and ``obj.x = 4`` will always reference/modify
        # the same attribute.

        try:
            object.__getattribute__(self, name)
            return object.__setattr__(self, name, value)
        except AttributeError:
            pass

        # if this fails, go on to more involved attribute setting
        # (note that this matches __getattr__, above).
        if name in self._internal_names_set:
            object.__setattr__(self, name, value)
        elif name in self._metadata:
            object.__setattr__(self, name, value)
        else:
            try:
                existing = getattr(self, name)
                if isinstance(existing, Index):
                    object.__setattr__(self, name, value)
                elif name in self._info_axis:
                    self[name] = value
                else:
                    object.__setattr__(self, name, value)
            except (AttributeError, TypeError):
                if isinstance(self, ABCDataFrame) and (is_list_like(value)):
                    warnings.warn("Pandas doesn't allow columns to be "
                                  "created via a new attribute name - see "
                                  "https://pandas.pydata.org/pandas-docs/"
                                  "stable/indexing.html#attribute-access",
                                  stacklevel=2)
                object.__setattr__(self, name, value)

    # ----------------------------------------------------------------------
    # Getting and setting elements

    # ----------------------------------------------------------------------
    # Consolidation of internals

    def _protect_consolidate(self, f):
        """Consolidate _data -- if the blocks have changed, then clear the
        cache
        """
        blocks_before = len(self._data.blocks)
        result = f()
        if len(self._data.blocks) != blocks_before:
            self._clear_item_cache()
        return result

    def _consolidate_inplace(self):
        """Consolidate data in place and return None"""

        def f():
            self._data = self._data.consolidate()

        self._protect_consolidate(f)

    def _consolidate(self, inplace=False):
        """
        Compute NDFrame with "consolidated" internals (data of each dtype
        grouped together in a single ndarray).

        Parameters
        ----------
        inplace : boolean, default False
            If False return new object, otherwise modify existing object

        Returns
        -------
        consolidated : type of caller
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        if inplace:
            self._consolidate_inplace()
        else:
            f = lambda: self._data.consolidate()
            cons_data = self._protect_consolidate(f)
            return self._constructor(cons_data).__finalize__(self)

    def consolidate(self, inplace=False):
        """Compute NDFrame with "consolidated" internals (data of each dtype
        grouped together in a single ndarray).

        .. deprecated:: 0.20.0
            Consolidate will be an internal implementation only.
        """
        # 15483
        warnings.warn("consolidate is deprecated and will be removed in a "
                      "future release.", FutureWarning, stacklevel=2)
        return self._consolidate(inplace)

    @property
    def _is_mixed_type(self):
        f = lambda: self._data.is_mixed_type
        return self._protect_consolidate(f)

    @property
    def _is_numeric_mixed_type(self):
        f = lambda: self._data.is_numeric_mixed_type
        return self._protect_consolidate(f)

    @property
    def _is_datelike_mixed_type(self):
        f = lambda: self._data.is_datelike_mixed_type
        return self._protect_consolidate(f)

    def _check_inplace_setting(self, value):
        """ check whether we allow in-place setting with this type of value """

        if self._is_mixed_type:
            if not self._is_numeric_mixed_type:

                # allow an actual np.nan thru
                try:
                    if np.isnan(value):
                        return True
                except Exception:
                    pass

                raise TypeError('Cannot do inplace boolean setting on '
                                'mixed-types with a non np.nan value')

        return True

    def _get_numeric_data(self):
        return self._constructor(
            self._data.get_numeric_data()).__finalize__(self)

    def _get_bool_data(self):
        return self._constructor(self._data.get_bool_data()).__finalize__(self)

    # ----------------------------------------------------------------------
    # Internal Interface Methods

    def as_matrix(self, columns=None):
        """Convert the frame to its Numpy-array representation.

        .. deprecated:: 0.23.0
            Use :meth:`DataFrame.values` instead.

        Parameters
        ----------
        columns: list, optional, default:None
            If None, return all columns, otherwise, returns specified columns.

        Returns
        -------
        values : ndarray
            If the caller is heterogeneous and contains booleans or objects,
            the result will be of dtype=object. See Notes.


        Notes
        -----
        Return is NOT a Numpy-matrix, rather, a Numpy-array.

        The dtype will be a lower-common-denominator dtype (implicit
        upcasting); that is to say if the dtypes (even of numeric types)
        are mixed, the one that accommodates all will be chosen. Use this
        with care if you are not dealing with the blocks.

        e.g. If the dtypes are float16 and float32, dtype will be upcast to
        float32.  If dtypes are int32 and uint8, dtype will be upcase to
        int32. By numpy.find_common_type convention, mixing int64 and uint64
        will result in a flot64 dtype.

        This method is provided for backwards compatibility. Generally,
        it is recommended to use '.values'.

        See Also
        --------
        pandas.DataFrame.values
        """
        warnings.warn("Method .as_matrix will be removed in a future version. "
                      "Use .values instead.", FutureWarning, stacklevel=2)
        self._consolidate_inplace()
        return self._data.as_array(transpose=self._AXIS_REVERSED,
                                   items=columns)

    @property
    def values(self):
        """
        Return a Numpy representation of the DataFrame.

        Only the values in the DataFrame will be returned, the axes labels
        will be removed.

        Returns
        -------
        numpy.ndarray
            The values of the DataFrame.

        Examples
        --------
        A DataFrame where all columns are the same type (e.g., int64) results
        in an array of the same type.

        >>> df = pd.DataFrame({'age':    [ 3,  29],
        ...                    'height': [94, 170],
        ...                    'weight': [31, 115]})
        >>> df
           age  height  weight
        0    3      94      31
        1   29     170     115
        >>> df.dtypes
        age       int64
        height    int64
        weight    int64
        dtype: object
        >>> df.values
        array([[  3,  94,  31],
               [ 29, 170, 115]], dtype=int64)

        A DataFrame with mixed type columns(e.g., str/object, int64, float32)
        results in an ndarray of the broadest type that accommodates these
        mixed types (e.g., object).

        >>> df2 = pd.DataFrame([('parrot',   24.0, 'second'),
        ...                     ('lion',     80.5, 1),
        ...                     ('monkey', np.nan, None)],
        ...                   columns=('name', 'max_speed', 'rank'))
        >>> df2.dtypes
        name          object
        max_speed    float64
        rank          object
        dtype: object
        >>> df2.values
        array([['parrot', 24.0, 'second'],
               ['lion', 80.5, 1],
               ['monkey', nan, None]], dtype=object)

        Notes
        -----
        The dtype will be a lower-common-denominator dtype (implicit
        upcasting); that is to say if the dtypes (even of numeric types)
        are mixed, the one that accommodates all will be chosen. Use this
        with care if you are not dealing with the blocks.

        e.g. If the dtypes are float16 and float32, dtype will be upcast to
        float32.  If dtypes are int32 and uint8, dtype will be upcast to
        int32. By :func:`numpy.find_common_type` convention, mixing int64
        and uint64 will result in a float64 dtype.

        See Also
        --------
        pandas.DataFrame.index : Retrievie the index labels
        pandas.DataFrame.columns : Retrieving the column names
        """
        self._consolidate_inplace()
        return self._data.as_array(transpose=self._AXIS_REVERSED)

    @property
    def _values(self):
        """internal implementation"""
        return self.values

    @property
    def _get_values(self):
        # compat
        return self.values

    def get_values(self):
        """
        Return an ndarray after converting sparse values to dense.

        This is the same as ``.values`` for non-sparse data. For sparse
        data contained in a `pandas.SparseArray`, the data are first
        converted to a dense representation.

        Returns
        -------
        numpy.ndarray
            Numpy representation of DataFrame

        See Also
        --------
        values : Numpy representation of DataFrame.
        pandas.SparseArray : Container for sparse data.

        Examples
        --------
        >>> df = pd.DataFrame({'a': [1, 2], 'b': [True, False],
        ...                    'c': [1.0, 2.0]})
        >>> df
           a      b    c
        0  1   True  1.0
        1  2  False  2.0

        >>> df.get_values()
        array([[1, True, 1.0], [2, False, 2.0]], dtype=object)

        >>> df = pd.DataFrame({"a": pd.SparseArray([1, None, None]),
        ...                    "c": [1.0, 2.0, 3.0]})
        >>> df
             a    c
        0  1.0  1.0
        1  NaN  2.0
        2  NaN  3.0

        >>> df.get_values()
        array([[ 1.,  1.],
               [nan,  2.],
               [nan,  3.]])
        """
        return self.values

    def get_dtype_counts(self):
        """
        Return counts of unique dtypes in this object.

        Returns
        -------
        dtype : Series
            Series with the count of columns with each dtype.

        See Also
        --------
        dtypes : Return the dtypes in this object.

        Examples
        --------
        >>> a = [['a', 1, 1.0], ['b', 2, 2.0], ['c', 3, 3.0]]
        >>> df = pd.DataFrame(a, columns=['str', 'int', 'float'])
        >>> df
          str  int  float
        0   a    1    1.0
        1   b    2    2.0
        2   c    3    3.0

        >>> df.get_dtype_counts()
        float64    1
        int64      1
        object     1
        dtype: int64
        """
        from pandas import Series
        return Series(self._data.get_dtype_counts())

    def get_ftype_counts(self):
        """
        Return counts of unique ftypes in this object.

        .. deprecated:: 0.23.0

        This is useful for SparseDataFrame or for DataFrames containing
        sparse arrays.

        Returns
        -------
        dtype : Series
            Series with the count of columns with each type and
            sparsity (dense/sparse)

        See Also
        --------
        ftypes : Return ftypes (indication of sparse/dense and dtype) in
            this object.

        Examples
        --------
        >>> a = [['a', 1, 1.0], ['b', 2, 2.0], ['c', 3, 3.0]]
        >>> df = pd.DataFrame(a, columns=['str', 'int', 'float'])
        >>> df
          str  int  float
        0   a    1    1.0
        1   b    2    2.0
        2   c    3    3.0

        >>> df.get_ftype_counts()
        float64:dense    1
        int64:dense      1
        object:dense     1
        dtype: int64
        """
        warnings.warn("get_ftype_counts is deprecated and will "
                      "be removed in a future version",
                      FutureWarning, stacklevel=2)

        from pandas import Series
        return Series(self._data.get_ftype_counts())

    @property
    def dtypes(self):
        """
        Return the dtypes in the DataFrame.

        This returns a Series with the data type of each column.
        The result's index is the original DataFrame's columns. Columns
        with mixed types are stored with the ``object`` dtype. See
        :ref:`the User Guide <basics.dtypes>` for more.

        Returns
        -------
        pandas.Series
            The data type of each column.

        See Also
        --------
        pandas.DataFrame.ftypes : dtype and sparsity information.

        Examples
        --------
        >>> df = pd.DataFrame({'float': [1.0],
        ...                    'int': [1],
        ...                    'datetime': [pd.Timestamp('20180310')],
        ...                    'string': ['foo']})
        >>> df.dtypes
        float              float64
        int                  int64
        datetime    datetime64[ns]
        string              object
        dtype: object
        """
        from pandas import Series
        return Series(self._data.get_dtypes(), index=self._info_axis,
                      dtype=np.object_)

    @property
    def ftypes(self):
        """
        Return the ftypes (indication of sparse/dense and dtype) in DataFrame.

        This returns a Series with the data type of each column.
        The result's index is the original DataFrame's columns. Columns
        with mixed types are stored with the ``object`` dtype.  See
        :ref:`the User Guide <basics.dtypes>` for more.

        Returns
        -------
        pandas.Series
            The data type and indication of sparse/dense of each column.

        See Also
        --------
        pandas.DataFrame.dtypes: Series with just dtype information.
        pandas.SparseDataFrame : Container for sparse tabular data.

        Notes
        -----
        Sparse data should have the same dtypes as its dense representation.

        Examples
        --------
        >>> import numpy as np
        >>> arr = np.random.RandomState(0).randn(100, 4)
        >>> arr[arr < .8] = np.nan
        >>> pd.DataFrame(arr).ftypes
        0    float64:dense
        1    float64:dense
        2    float64:dense
        3    float64:dense
        dtype: object

        >>> pd.SparseDataFrame(arr).ftypes
        0    float64:sparse
        1    float64:sparse
        2    float64:sparse
        3    float64:sparse
        dtype: object
        """
        from pandas import Series
        return Series(self._data.get_ftypes(), index=self._info_axis,
                      dtype=np.object_)

    def as_blocks(self, copy=True):
        """
        Convert the frame to a dict of dtype -> Constructor Types that each has
        a homogeneous dtype.

        .. deprecated:: 0.21.0

        NOTE: the dtypes of the blocks WILL BE PRESERVED HERE (unlike in
              as_matrix)

        Parameters
        ----------
        copy : boolean, default True

        Returns
        -------
        values : a dict of dtype -> Constructor Types
        """
        warnings.warn("as_blocks is deprecated and will "
                      "be removed in a future version",
                      FutureWarning, stacklevel=2)
        return self._to_dict_of_blocks(copy=copy)

    @property
    def blocks(self):
        """
        Internal property, property synonym for as_blocks()

        .. deprecated:: 0.21.0
        """
        return self.as_blocks()

    def _to_dict_of_blocks(self, copy=True):
        """
        Return a dict of dtype -> Constructor Types that
        each is a homogeneous dtype.

        Internal ONLY
        """
        return {k: self._constructor(v).__finalize__(self)
                for k, v, in self._data.to_dict(copy=copy).items()}

    @deprecate_kwarg(old_arg_name='raise_on_error', new_arg_name='errors',
                     mapping={True: 'raise', False: 'ignore'})
    def astype(self, dtype, copy=True, errors='raise', **kwargs):
        """
        Cast a pandas object to a specified dtype ``dtype``.

        Parameters
        ----------
        dtype : data type, or dict of column name -> data type
            Use a numpy.dtype or Python type to cast entire pandas object to
            the same type. Alternatively, use {col: dtype, ...}, where col is a
            column label and dtype is a numpy.dtype or Python type to cast one
            or more of the DataFrame's columns to column-specific types.
        copy : bool, default True.
            Return a copy when ``copy=True`` (be very careful setting
            ``copy=False`` as changes to values then may propagate to other
            pandas objects).
        errors : {'raise', 'ignore'}, default 'raise'.
            Control raising of exceptions on invalid data for provided dtype.

            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object

            .. versionadded:: 0.20.0

        raise_on_error : raise on invalid input
            .. deprecated:: 0.20.0
               Use ``errors`` instead
        kwargs : keyword arguments to pass on to the constructor

        Returns
        -------
        casted : type of caller

        Examples
        --------
        >>> ser = pd.Series([1, 2], dtype='int32')
        >>> ser
        0    1
        1    2
        dtype: int32
        >>> ser.astype('int64')
        0    1
        1    2
        dtype: int64

        Convert to categorical type:

        >>> ser.astype('category')
        0    1
        1    2
        dtype: category
        Categories (2, int64): [1, 2]

        Convert to ordered categorical type with custom ordering:

        >>> ser.astype('category', ordered=True, categories=[2, 1])
        0    1
        1    2
        dtype: category
        Categories (2, int64): [2 < 1]

        Note that using ``copy=False`` and changing data on a new
        pandas object may propagate changes:

        >>> s1 = pd.Series([1,2])
        >>> s2 = s1.astype('int64', copy=False)
        >>> s2[0] = 10
        >>> s1  # note that s1[0] has changed too
        0    10
        1     2
        dtype: int64

        See also
        --------
        pandas.to_datetime : Convert argument to datetime.
        pandas.to_timedelta : Convert argument to timedelta.
        pandas.to_numeric : Convert argument to a numeric type.
        numpy.ndarray.astype : Cast a numpy array to a specified type.
        """
        if is_dict_like(dtype):
            if self.ndim == 1:  # i.e. Series
                if len(dtype) > 1 or self.name not in dtype:
                    raise KeyError('Only the Series name can be used for '
                                   'the key in Series dtype mappings.')
                new_type = dtype[self.name]
                return self.astype(new_type, copy, errors, **kwargs)
            elif self.ndim > 2:
                raise NotImplementedError(
                    'astype() only accepts a dtype arg of type dict when '
                    'invoked on Series and DataFrames. A single dtype must be '
                    'specified when invoked on a Panel.'
                )
            for col_name in dtype.keys():
                if col_name not in self:
                    raise KeyError('Only a column name can be used for the '
                                   'key in a dtype mappings argument.')
            results = []
            for col_name, col in self.iteritems():
                if col_name in dtype:
                    results.append(col.astype(dtype[col_name], copy=copy))
                else:
                    results.append(results.append(col.copy() if copy else col))

        elif is_categorical_dtype(dtype) and self.ndim > 1:
            # GH 18099: columnwise conversion to categorical
            results = (self[col].astype(dtype, copy=copy) for col in self)

        else:
            # else, only a single dtype is given
            new_data = self._data.astype(dtype=dtype, copy=copy, errors=errors,
                                         **kwargs)
            return self._constructor(new_data).__finalize__(self)

        # GH 19920: retain column metadata after concat
        result = pd.concat(results, axis=1, copy=False)
        result.columns = self.columns
        return result

    def copy(self, deep=True):
        """
        Make a copy of this object's indices and data.

        When ``deep=True`` (default), a new object will be created with a
        copy of the calling object's data and indices. Modifications to
        the data or indices of the copy will not be reflected in the
        original object (see notes below).

        When ``deep=False``, a new object will be created without copying
        the calling object's data or index (only references to the data
        and index are copied). Any changes to the data of the original
        will be reflected in the shallow copy (and vice versa).

        Parameters
        ----------
        deep : bool, default True
            Make a deep copy, including a copy of the data and the indices.
            With ``deep=False`` neither the indices nor the data are copied.

        Returns
        -------
        copy : Series, DataFrame or Panel
            Object type matches caller.

        Notes
        -----
        When ``deep=True``, data is copied but actual Python objects
        will not be copied recursively, only the reference to the object.
        This is in contrast to `copy.deepcopy` in the Standard Library,
        which recursively copies object data (see examples below).

        While ``Index`` objects are copied when ``deep=True``, the underlying
        numpy array is not copied for performance reasons. Since ``Index`` is
        immutable, the underlying data can be safely shared and a copy
        is not needed.

        Examples
        --------
        >>> s = pd.Series([1, 2], index=["a", "b"])
        >>> s
        a    1
        b    2
        dtype: int64

        >>> s_copy = s.copy()
        >>> s_copy
        a    1
        b    2
        dtype: int64

        **Shallow copy versus default (deep) copy:**

        >>> s = pd.Series([1, 2], index=["a", "b"])
        >>> deep = s.copy()
        >>> shallow = s.copy(deep=False)

        Shallow copy shares data and index with original.

        >>> s is shallow
        False
        >>> s.values is shallow.values and s.index is shallow.index
        True

        Deep copy has own copy of data and index.

        >>> s is deep
        False
        >>> s.values is deep.values or s.index is deep.index
        False

        Updates to the data shared by shallow copy and original is reflected
        in both; deep copy remains unchanged.

        >>> s[0] = 3
        >>> shallow[1] = 4
        >>> s
        a    3
        b    4
        dtype: int64
        >>> shallow
        a    3
        b    4
        dtype: int64
        >>> deep
        a    1
        b    2
        dtype: int64

        Note that when copying an object containing Python objects, a deep copy
        will copy the data, but will not do so recursively. Updating a nested
        data object will be reflected in the deep copy.

        >>> s = pd.Series([[1, 2], [3, 4]])
        >>> deep = s.copy()
        >>> s[0][0] = 10
        >>> s
        0    [10, 2]
        1     [3, 4]
        dtype: object
        >>> deep
        0    [10, 2]
        1     [3, 4]
        dtype: object
        """
        data = self._data.copy(deep=deep)
        return self._constructor(data).__finalize__(self)

    def __copy__(self, deep=True):
        return self.copy(deep=deep)

    def __deepcopy__(self, memo=None):
        if memo is None:
            memo = {}
        return self.copy(deep=True)

    def _convert(self, datetime=False, numeric=False, timedelta=False,
                 coerce=False, copy=True):
        """
        Attempt to infer better dtype for object columns

        Parameters
        ----------
        datetime : boolean, default False
            If True, convert to date where possible.
        numeric : boolean, default False
            If True, attempt to convert to numbers (including strings), with
            unconvertible values becoming NaN.
        timedelta : boolean, default False
            If True, convert to timedelta where possible.
        coerce : boolean, default False
            If True, force conversion with unconvertible values converted to
            nulls (NaN or NaT)
        copy : boolean, default True
            If True, return a copy even if no copy is necessary (e.g. no
            conversion was done). Note: This is meant for internal use, and
            should not be confused with inplace.

        Returns
        -------
        converted : same as input object
        """
        return self._constructor(
            self._data.convert(datetime=datetime, numeric=numeric,
                               timedelta=timedelta, coerce=coerce,
                               copy=copy)).__finalize__(self)

    def convert_objects(self, convert_dates=True, convert_numeric=False,
                        convert_timedeltas=True, copy=True):
        """Attempt to infer better dtype for object columns.

        .. deprecated:: 0.21.0

        Parameters
        ----------
        convert_dates : boolean, default True
            If True, convert to date where possible. If 'coerce', force
            conversion, with unconvertible values becoming NaT.
        convert_numeric : boolean, default False
            If True, attempt to coerce to numbers (including strings), with
            unconvertible values becoming NaN.
        convert_timedeltas : boolean, default True
            If True, convert to timedelta where possible. If 'coerce', force
            conversion, with unconvertible values becoming NaT.
        copy : boolean, default True
            If True, return a copy even if no copy is necessary (e.g. no
            conversion was done). Note: This is meant for internal use, and
            should not be confused with inplace.

        See Also
        --------
        pandas.to_datetime : Convert argument to datetime.
        pandas.to_timedelta : Convert argument to timedelta.
        pandas.to_numeric : Return a fixed frequency timedelta index,
            with day as the default.

        Returns
        -------
        converted : same as input object
        """
        msg = ("convert_objects is deprecated.  To re-infer data dtypes for "
               "object columns, use {klass}.infer_objects()\nFor all "
               "other conversions use the data-type specific converters "
               "pd.to_datetime, pd.to_timedelta and pd.to_numeric."
               ).format(klass=self.__class__.__name__)
        warnings.warn(msg, FutureWarning, stacklevel=2)

        return self._constructor(
            self._data.convert(convert_dates=convert_dates,
                               convert_numeric=convert_numeric,
                               convert_timedeltas=convert_timedeltas,
                               copy=copy)).__finalize__(self)

    def infer_objects(self):
        """
        Attempt to infer better dtypes for object columns.

        Attempts soft conversion of object-dtyped
        columns, leaving non-object and unconvertible
        columns unchanged. The inference rules are the
        same as during normal Series/DataFrame construction.

        .. versionadded:: 0.21.0

        See Also
        --------
        pandas.to_datetime : Convert argument to datetime.
        pandas.to_timedelta : Convert argument to timedelta.
        pandas.to_numeric : Convert argument to numeric typeR

        Returns
        -------
        converted : same type as input object

        Examples
        --------
        >>> df = pd.DataFrame({"A": ["a", 1, 2, 3]})
        >>> df = df.iloc[1:]
        >>> df
           A
        1  1
        2  2
        3  3

        >>> df.dtypes
        A    object
        dtype: object

        >>> df.infer_objects().dtypes
        A    int64
        dtype: object
        """
        # numeric=False necessary to only soft convert;
        # python objects will still be converted to
        # native numpy numeric types
        return self._constructor(
            self._data.convert(datetime=True, numeric=False,
                               timedelta=True, coerce=False,
                               copy=True)).__finalize__(self)

    # ----------------------------------------------------------------------
    # Filling NA's

    def fillna(self, value=None, method=None, axis=None, inplace=False,
               limit=None, downcast=None):
        """
        Fill NA/NaN values using the specified method

        Parameters
        ----------
        value : scalar, dict, Series, or DataFrame
            Value to use to fill holes (e.g. 0), alternately a
            dict/Series/DataFrame of values specifying which value to use for
            each index (for a Series) or column (for a DataFrame). (values not
            in the dict/Series/DataFrame will not be filled). This value cannot
            be a list.
        method : {'backfill', 'bfill', 'pad', 'ffill', None}, default None
            Method to use for filling holes in reindexed Series
            pad / ffill: propagate last valid observation forward to next valid
            backfill / bfill: use NEXT valid observation to fill gap
        axis : %(axes_single_arg)s
        inplace : boolean, default False
            If True, fill in place. Note: this will modify any
            other views on this object, (e.g. a no-copy slice for a column in a
            DataFrame).
        limit : int, default None
            If method is specified, this is the maximum number of consecutive
            NaN values to forward/backward fill. In other words, if there is
            a gap with more than this number of consecutive NaNs, it will only
            be partially filled. If method is not specified, this is the
            maximum number of entries along the entire axis where NaNs will be
            filled. Must be greater than 0 if not None.
        downcast : dict, default is None
            a dict of item->dtype of what to downcast if possible,
            or the string 'infer' which will try to downcast to an appropriate
            equal type (e.g. float64 to int64 if possible)

        See Also
        --------
        interpolate : Fill NaN values using interpolation.
        reindex, asfreq

        Returns
        -------
        filled : %(klass)s

        Examples
        --------
        >>> df = pd.DataFrame([[np.nan, 2, np.nan, 0],
        ...                    [3, 4, np.nan, 1],
        ...                    [np.nan, np.nan, np.nan, 5],
        ...                    [np.nan, 3, np.nan, 4]],
        ...                    columns=list('ABCD'))
        >>> df
             A    B   C  D
        0  NaN  2.0 NaN  0
        1  3.0  4.0 NaN  1
        2  NaN  NaN NaN  5
        3  NaN  3.0 NaN  4

        Replace all NaN elements with 0s.

        >>> df.fillna(0)
            A   B   C   D
        0   0.0 2.0 0.0 0
        1   3.0 4.0 0.0 1
        2   0.0 0.0 0.0 5
        3   0.0 3.0 0.0 4

        We can also propagate non-null values forward or backward.

        >>> df.fillna(method='ffill')
            A   B   C   D
        0   NaN 2.0 NaN 0
        1   3.0 4.0 NaN 1
        2   3.0 4.0 NaN 5
        3   3.0 3.0 NaN 4

        Replace all NaN elements in column 'A', 'B', 'C', and 'D', with 0, 1,
        2, and 3 respectively.

        >>> values = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
        >>> df.fillna(value=values)
            A   B   C   D
        0   0.0 2.0 2.0 0
        1   3.0 4.0 2.0 1
        2   0.0 1.0 2.0 5
        3   0.0 3.0 2.0 4

        Only replace the first NaN element.

        >>> df.fillna(value=values, limit=1)
            A   B   C   D
        0   0.0 2.0 2.0 0
        1   3.0 4.0 NaN 1
        2   NaN 1.0 NaN 5
        3   NaN 3.0 NaN 4
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')
        value, method = validate_fillna_kwargs(value, method)

        self._consolidate_inplace()

        # set the default here, so functions examining the signaure
        # can detect if something was set (e.g. in groupby) (GH9221)
        if axis is None:
            axis = 0
        axis = self._get_axis_number(axis)

        from pandas import DataFrame
        if value is None:

            if self._is_mixed_type and axis == 1:
                if inplace:
                    raise NotImplementedError()
                result = self.T.fillna(method=method, limit=limit).T

                # need to downcast here because of all of the transposes
                result._data = result._data.downcast()

                return result

            # > 3d
            if self.ndim > 3:
                raise NotImplementedError('Cannot fillna with a method for > '
                                          '3dims')

            # 3d
            elif self.ndim == 3:
                # fill in 2d chunks
                result = {col: s.fillna(method=method, value=value)
                          for col, s in self.iteritems()}
                new_obj = self._constructor.\
                    from_dict(result).__finalize__(self)
                new_data = new_obj._data

            else:
                # 2d or less
                new_data = self._data.interpolate(method=method, axis=axis,
                                                  limit=limit, inplace=inplace,
                                                  coerce=True,
                                                  downcast=downcast)
        else:
            if len(self._get_axis(axis)) == 0:
                return self

            if self.ndim == 1:
                if isinstance(value, (dict, ABCSeries)):
                    from pandas import Series
                    value = Series(value)
                elif not is_list_like(value):
                    pass
                else:
                    raise TypeError('"value" parameter must be a scalar, dict '
                                    'or Series, but you passed a '
                                    '"{0}"'.format(type(value).__name__))

                new_data = self._data.fillna(value=value, limit=limit,
                                             inplace=inplace,
                                             downcast=downcast)

            elif isinstance(value, (dict, ABCSeries)):
                if axis == 1:
                    raise NotImplementedError('Currently only can fill '
                                              'with dict/Series column '
                                              'by column')

                result = self if inplace else self.copy()
                for k, v in compat.iteritems(value):
                    if k not in result:
                        continue
                    obj = result[k]
                    obj.fillna(v, limit=limit, inplace=True, downcast=downcast)
                return result if not inplace else None

            elif not is_list_like(value):
                new_data = self._data.fillna(value=value, limit=limit,
                                             inplace=inplace,
                                             downcast=downcast)
            elif isinstance(value, DataFrame) and self.ndim == 2:
                new_data = self.where(self.notna(), value)
            else:
                raise ValueError("invalid fill value with a %s" % type(value))

        if inplace:
            self._update_inplace(new_data)
        else:
            return self._constructor(new_data).__finalize__(self)

    def ffill(self, axis=None, inplace=False, limit=None, downcast=None):
        """
        Synonym for :meth:`DataFrame.fillna(method='ffill') <DataFrame.fillna>`
        """
        return self.fillna(method='ffill', axis=axis, inplace=inplace,
                           limit=limit, downcast=downcast)

    def bfill(self, axis=None, inplace=False, limit=None, downcast=None):
        """
        Synonym for :meth:`DataFrame.fillna(method='bfill') <DataFrame.fillna>`
        """
        return self.fillna(method='bfill', axis=axis, inplace=inplace,
                           limit=limit, downcast=downcast)

    _shared_docs['replace'] = ("""
        Replace values given in `to_replace` with `value`.

        Values of the %(klass)s are replaced with other values dynamically.
        This differs from updating with ``.loc`` or ``.iloc``, which require
        you to specify a location to update with some value.

        Parameters
        ----------
        to_replace : str, regex, list, dict, Series, int, float, or None
            How to find the values that will be replaced.

            * numeric, str or regex:

                - numeric: numeric values equal to `to_replace` will be
                  replaced with `value`
                - str: string exactly matching `to_replace` will be replaced
                  with `value`
                - regex: regexs matching `to_replace` will be replaced with
                  `value`

            * list of str, regex, or numeric:

                - First, if `to_replace` and `value` are both lists, they
                  **must** be the same length.
                - Second, if ``regex=True`` then all of the strings in **both**
                  lists will be interpreted as regexs otherwise they will match
                  directly. This doesn't matter much for `value` since there
                  are only a few possible substitution regexes you can use.
                - str, regex and numeric rules apply as above.

            * dict:

                - Dicts can be used to specify different replacement values
                  for different existing values. For example,
                  ``{'a': 'b', 'y': 'z'}`` replaces the value 'a' with 'b' and
                  'y' with 'z'. To use a dict in this way the `value`
                  parameter should be `None`.
                - For a DataFrame a dict can specify that different values
                  should be replaced in different columns. For example,
                  ``{'a': 1, 'b': 'z'}`` looks for the value 1 in column 'a'
                  and the value 'z' in column 'b' and replaces these values
                  with whatever is specified in `value`. The `value` parameter
                  should not be ``None`` in this case. You can treat this as a
                  special case of passing two lists except that you are
                  specifying the column to search in.
                - For a DataFrame nested dictionaries, e.g.,
                  ``{'a': {'b': np.nan}}``, are read as follows: look in column
                  'a' for the value 'b' and replace it with NaN. The `value`
                  parameter should be ``None`` to use a nested dict in this
                  way. You can nest regular expressions as well. Note that
                  column names (the top-level dictionary keys in a nested
                  dictionary) **cannot** be regular expressions.

            * None:

                - This means that the `regex` argument must be a string,
                  compiled regular expression, or list, dict, ndarray or
                  Series of such elements. If `value` is also ``None`` then
                  this **must** be a nested dictionary or Series.

            See the examples section for examples of each of these.
        value : scalar, dict, list, str, regex, default None
            Value to replace any values matching `to_replace` with.
            For a DataFrame a dict of values can be used to specify which
            value to use for each column (columns not in the dict will not be
            filled). Regular expressions, strings and lists or dicts of such
            objects are also allowed.
        inplace : boolean, default False
            If True, in place. Note: this will modify any
            other views on this object (e.g. a column from a DataFrame).
            Returns the caller if this is True.
        limit : int, default None
            Maximum size gap to forward or backward fill.
        regex : bool or same types as `to_replace`, default False
            Whether to interpret `to_replace` and/or `value` as regular
            expressions. If this is ``True`` then `to_replace` *must* be a
            string. Alternatively, this could be a regular expression or a
            list, dict, or array of regular expressions in which case
            `to_replace` must be ``None``.
        method : {'pad', 'ffill', 'bfill', `None`}
            The method to use when for replacement, when `to_replace` is a
            scalar, list or tuple and `value` is ``None``.

            .. versionchanged:: 0.23.0
                Added to DataFrame.

        See Also
        --------
        %(klass)s.fillna : Fill NA values
        %(klass)s.where : Replace values based on boolean condition
        Series.str.replace : Simple string replacement.

        Returns
        -------
        %(klass)s
            Object after replacement.

        Raises
        ------
        AssertionError
            * If `regex` is not a ``bool`` and `to_replace` is not
              ``None``.
        TypeError
            * If `to_replace` is a ``dict`` and `value` is not a ``list``,
              ``dict``, ``ndarray``, or ``Series``
            * If `to_replace` is ``None`` and `regex` is not compilable
              into a regular expression or is a list, dict, ndarray, or
              Series.
            * When replacing multiple ``bool`` or ``datetime64`` objects and
              the arguments to `to_replace` does not match the type of the
              value being replaced
        ValueError
            * If a ``list`` or an ``ndarray`` is passed to `to_replace` and
              `value` but they are not the same length.

        Notes
        -----
        * Regex substitution is performed under the hood with ``re.sub``. The
          rules for substitution for ``re.sub`` are the same.
        * Regular expressions will only substitute on strings, meaning you
          cannot provide, for example, a regular expression matching floating
          point numbers and expect the columns in your frame that have a
          numeric dtype to be matched. However, if those floating point
          numbers *are* strings, then you can do this.
        * This method has *a lot* of options. You are encouraged to experiment
          and play with this method to gain intuition about how it works.
        * When dict is used as the `to_replace` value, it is like
          key(s) in the dict are the to_replace part and
          value(s) in the dict are the value parameter.

        Examples
        --------

        **Scalar `to_replace` and `value`**

        >>> s = pd.Series([0, 1, 2, 3, 4])
        >>> s.replace(0, 5)
        0    5
        1    1
        2    2
        3    3
        4    4
        dtype: int64

        >>> df = pd.DataFrame({'A': [0, 1, 2, 3, 4],
        ...                    'B': [5, 6, 7, 8, 9],
        ...                    'C': ['a', 'b', 'c', 'd', 'e']})
        >>> df.replace(0, 5)
           A  B  C
        0  5  5  a
        1  1  6  b
        2  2  7  c
        3  3  8  d
        4  4  9  e

        **List-like `to_replace`**

        >>> df.replace([0, 1, 2, 3], 4)
           A  B  C
        0  4  5  a
        1  4  6  b
        2  4  7  c
        3  4  8  d
        4  4  9  e

        >>> df.replace([0, 1, 2, 3], [4, 3, 2, 1])
           A  B  C
        0  4  5  a
        1  3  6  b
        2  2  7  c
        3  1  8  d
        4  4  9  e

        >>> s.replace([1, 2], method='bfill')
        0    0
        1    3
        2    3
        3    3
        4    4
        dtype: int64

        **dict-like `to_replace`**

        >>> df.replace({0: 10, 1: 100})
             A  B  C
        0   10  5  a
        1  100  6  b
        2    2  7  c
        3    3  8  d
        4    4  9  e

        >>> df.replace({'A': 0, 'B': 5}, 100)
             A    B  C
        0  100  100  a
        1    1    6  b
        2    2    7  c
        3    3    8  d
        4    4    9  e

        >>> df.replace({'A': {0: 100, 4: 400}})
             A  B  C
        0  100  5  a
        1    1  6  b
        2    2  7  c
        3    3  8  d
        4  400  9  e

        **Regular expression `to_replace`**

        >>> df = pd.DataFrame({'A': ['bat', 'foo', 'bait'],
        ...                    'B': ['abc', 'bar', 'xyz']})
        >>> df.replace(to_replace=r'^ba.$', value='new', regex=True)
              A    B
        0   new  abc
        1   foo  new
        2  bait  xyz

        >>> df.replace({'A': r'^ba.$'}, {'A': 'new'}, regex=True)
              A    B
        0   new  abc
        1   foo  bar
        2  bait  xyz

        >>> df.replace(regex=r'^ba.$', value='new')
              A    B
        0   new  abc
        1   foo  new
        2  bait  xyz

        >>> df.replace(regex={r'^ba.$':'new', 'foo':'xyz'})
              A    B
        0   new  abc
        1   xyz  new
        2  bait  xyz

        >>> df.replace(regex=[r'^ba.$', 'foo'], value='new')
              A    B
        0   new  abc
        1   new  new
        2  bait  xyz

        Note that when replacing multiple ``bool`` or ``datetime64`` objects,
        the data types in the `to_replace` parameter must match the data
        type of the value being replaced:

        >>> df = pd.DataFrame({'A': [True, False, True],
        ...                    'B': [False, True, False]})
        >>> df.replace({'a string': 'new value', True: False})  # raises
        Traceback (most recent call last):
            ...
        TypeError: Cannot compare types 'ndarray(dtype=bool)' and 'str'

        This raises a ``TypeError`` because one of the ``dict`` keys is not of
        the correct type for replacement.

        Compare the behavior of ``s.replace({'a': None})`` and
        ``s.replace('a', None)`` to understand the pecularities
        of the `to_replace` parameter:

        >>> s = pd.Series([10, 'a', 'a', 'b', 'a'])

        When one uses a dict as the `to_replace` value, it is like the
        value(s) in the dict are equal to the `value` parameter.
        ``s.replace({'a': None})`` is equivalent to
        ``s.replace(to_replace={'a': None}, value=None, method=None)``:

        >>> s.replace({'a': None})
        0      10
        1    None
        2    None
        3       b
        4    None
        dtype: object

        When ``value=None`` and `to_replace` is a scalar, list or
        tuple, `replace` uses the method parameter (default 'pad') to do the
        replacement. So this is why the 'a' values are being replaced by 10
        in rows 1 and 2 and 'b' in row 4 in this case.
        The command ``s.replace('a', None)`` is actually equivalent to
        ``s.replace(to_replace='a', value=None, method='pad')``:

        >>> s.replace('a', None)
        0    10
        1    10
        2    10
        3     b
        4     b
        dtype: object
    """)

    @Appender(_shared_docs['replace'] % _shared_doc_kwargs)
    def replace(self, to_replace=None, value=None, inplace=False, limit=None,
                regex=False, method='pad'):
        inplace = validate_bool_kwarg(inplace, 'inplace')
        if not is_bool(regex) and to_replace is not None:
            raise AssertionError("'to_replace' must be 'None' if 'regex' is "
                                 "not a bool")

        self._consolidate_inplace()

        if value is None:
            # passing a single value that is scalar like
            # when value is None (GH5319), for compat
            if not is_dict_like(to_replace) and not is_dict_like(regex):
                to_replace = [to_replace]

            if isinstance(to_replace, (tuple, list)):
                if isinstance(self, pd.DataFrame):
                    return self.apply(_single_replace,
                                      args=(to_replace, method, inplace,
                                            limit))
                return _single_replace(self, to_replace, method, inplace,
                                       limit)

            if not is_dict_like(to_replace):
                if not is_dict_like(regex):
                    raise TypeError('If "to_replace" and "value" are both None'
                                    ' and "to_replace" is not a list, then '
                                    'regex must be a mapping')
                to_replace = regex
                regex = True

            items = list(compat.iteritems(to_replace))
            keys, values = lzip(*items) or ([], [])

            are_mappings = [is_dict_like(v) for v in values]

            if any(are_mappings):
                if not all(are_mappings):
                    raise TypeError("If a nested mapping is passed, all values"
                                    " of the top level mapping must be "
                                    "mappings")
                # passed a nested dict/Series
                to_rep_dict = {}
                value_dict = {}

                for k, v in items:
                    keys, values = lzip(*v.items()) or ([], [])
                    if set(keys) & set(values):
                        raise ValueError("Replacement not allowed with "
                                         "overlapping keys and values")
                    to_rep_dict[k] = list(keys)
                    value_dict[k] = list(values)

                to_replace, value = to_rep_dict, value_dict
            else:
                to_replace, value = keys, values

            return self.replace(to_replace, value, inplace=inplace,
                                limit=limit, regex=regex)
        else:

            # need a non-zero len on all axes
            for a in self._AXIS_ORDERS:
                if not len(self._get_axis(a)):
                    return self

            new_data = self._data
            if is_dict_like(to_replace):
                if is_dict_like(value):  # {'A' : NA} -> {'A' : 0}
                    res = self if inplace else self.copy()
                    for c, src in compat.iteritems(to_replace):
                        if c in value and c in self:
                            # object conversion is handled in
                            # series.replace which is called recursivelly
                            res[c] = res[c].replace(to_replace=src,
                                                    value=value[c],
                                                    inplace=False,
                                                    regex=regex)
                    return None if inplace else res

                # {'A': NA} -> 0
                elif not is_list_like(value):
                    keys = [(k, src) for k, src in compat.iteritems(to_replace)
                            if k in self]
                    keys_len = len(keys) - 1
                    for i, (k, src) in enumerate(keys):
                        convert = i == keys_len
                        new_data = new_data.replace(to_replace=src,
                                                    value=value,
                                                    filter=[k],
                                                    inplace=inplace,
                                                    regex=regex,
                                                    convert=convert)
                else:
                    raise TypeError('value argument must be scalar, dict, or '
                                    'Series')

            elif is_list_like(to_replace):  # [NA, ''] -> [0, 'missing']
                if is_list_like(value):
                    if len(to_replace) != len(value):
                        raise ValueError('Replacement lists must match '
                                         'in length. Expecting %d got %d ' %
                                         (len(to_replace), len(value)))

                    new_data = self._data.replace_list(src_list=to_replace,
                                                       dest_list=value,
                                                       inplace=inplace,
                                                       regex=regex)

                else:  # [NA, ''] -> 0
                    new_data = self._data.replace(to_replace=to_replace,
                                                  value=value, inplace=inplace,
                                                  regex=regex)
            elif to_replace is None:
                if not (is_re_compilable(regex) or
                        is_list_like(regex) or is_dict_like(regex)):
                    raise TypeError("'regex' must be a string or a compiled "
                                    "regular expression or a list or dict of "
                                    "strings or regular expressions, you "
                                    "passed a"
                                    " {0!r}".format(type(regex).__name__))
                return self.replace(regex, value, inplace=inplace, limit=limit,
                                    regex=True)
            else:

                # dest iterable dict-like
                if is_dict_like(value):  # NA -> {'A' : 0, 'B' : -1}
                    new_data = self._data

                    for k, v in compat.iteritems(value):
                        if k in self:
                            new_data = new_data.replace(to_replace=to_replace,
                                                        value=v, filter=[k],
                                                        inplace=inplace,
                                                        regex=regex)

                elif not is_list_like(value):  # NA -> 0
                    new_data = self._data.replace(to_replace=to_replace,
                                                  value=value, inplace=inplace,
                                                  regex=regex)
                else:
                    msg = ('Invalid "to_replace" type: '
                           '{0!r}').format(type(to_replace).__name__)
                    raise TypeError(msg)  # pragma: no cover

        if inplace:
            self._update_inplace(new_data)
        else:
            return self._constructor(new_data).__finalize__(self)

    _shared_docs['interpolate'] = """
        Please note that only ``method='linear'`` is supported for
        DataFrames/Series with a MultiIndex.

        Parameters
        ----------
        method : {'linear', 'time', 'index', 'values', 'nearest', 'zero',
                  'slinear', 'quadratic', 'cubic', 'barycentric', 'krogh',
                  'polynomial', 'spline', 'piecewise_polynomial',
                  'from_derivatives', 'pchip', 'akima'}

            * 'linear': ignore the index and treat the values as equally
              spaced. This is the only method supported on MultiIndexes.
              default
            * 'time': interpolation works on daily and higher resolution
              data to interpolate given length of interval
            * 'index', 'values': use the actual numerical values of the index
            * 'nearest', 'zero', 'slinear', 'quadratic', 'cubic',
              'barycentric', 'polynomial' is passed to
              ``scipy.interpolate.interp1d``. Both 'polynomial' and 'spline'
              require that you also specify an `order` (int),
              e.g. df.interpolate(method='polynomial', order=4).
              These use the actual numerical values of the index.
            * 'krogh', 'piecewise_polynomial', 'spline', 'pchip' and 'akima'
              are all wrappers around the scipy interpolation methods of
              similar names. These use the actual numerical values of the
              index. For more information on their behavior, see the
              `scipy documentation
              <http://docs.scipy.org/doc/scipy/reference/interpolate.html#univariate-interpolation>`__
              and `tutorial documentation
              <http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html>`__
            * 'from_derivatives' refers to BPoly.from_derivatives which
              replaces 'piecewise_polynomial' interpolation method in
              scipy 0.18

            .. versionadded:: 0.18.1

               Added support for the 'akima' method
               Added interpolate method 'from_derivatives' which replaces
               'piecewise_polynomial' in scipy 0.18; backwards-compatible with
               scipy < 0.18

        axis : {0, 1}, default 0
            * 0: fill column-by-column
            * 1: fill row-by-row
        limit : int, default None.
            Maximum number of consecutive NaNs to fill. Must be greater than 0.
        limit_direction : {'forward', 'backward', 'both'}, default 'forward'
        limit_area : {'inside', 'outside'}, default None
            * None: (default) no fill restriction
            * 'inside' Only fill NaNs surrounded by valid values (interpolate).
            * 'outside' Only fill NaNs outside valid values (extrapolate).

            If limit is specified, consecutive NaNs will be filled in this
            direction.

            .. versionadded:: 0.21.0
        inplace : bool, default False
            Update the NDFrame in place if possible.
        downcast : optional, 'infer' or None, defaults to None
            Downcast dtypes if possible.
        kwargs : keyword arguments to pass on to the interpolating function.

        Returns
        -------
        Series or DataFrame of same shape interpolated at the NaNs

        See Also
        --------
        reindex, replace, fillna

        Examples
        --------

        Filling in NaNs

        >>> s = pd.Series([0, 1, np.nan, 3])
        >>> s.interpolate()
        0    0
        1    1
        2    2
        3    3
        dtype: float64

        """

    @Appender(_shared_docs['interpolate'] % _shared_doc_kwargs)
    def interpolate(self, method='linear', axis=0, limit=None, inplace=False,
                    limit_direction='forward', limit_area=None,
                    downcast=None, **kwargs):
        """
        Interpolate values according to different methods.
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')

        if self.ndim > 2:
            raise NotImplementedError("Interpolate has not been implemented "
                                      "on Panel and Panel 4D objects.")

        if axis == 0:
            ax = self._info_axis_name
            _maybe_transposed_self = self
        elif axis == 1:
            _maybe_transposed_self = self.T
            ax = 1
        else:
            _maybe_transposed_self = self
        ax = _maybe_transposed_self._get_axis_number(ax)

        if _maybe_transposed_self.ndim == 2:
            alt_ax = 1 - ax
        else:
            alt_ax = ax

        if (isinstance(_maybe_transposed_self.index, MultiIndex) and
                method != 'linear'):
            raise ValueError("Only `method=linear` interpolation is supported "
                             "on MultiIndexes.")

        if _maybe_transposed_self._data.get_dtype_counts().get(
                'object') == len(_maybe_transposed_self.T):
            raise TypeError("Cannot interpolate with all NaNs.")

        # create/use the index
        if method == 'linear':
            # prior default
            index = np.arange(len(_maybe_transposed_self._get_axis(alt_ax)))
        else:
            index = _maybe_transposed_self._get_axis(alt_ax)

        if isna(index).any():
            raise NotImplementedError("Interpolation with NaNs in the index "
                                      "has not been implemented. Try filling "
                                      "those NaNs before interpolating.")
        data = _maybe_transposed_self._data
        new_data = data.interpolate(method=method, axis=ax, index=index,
                                    values=_maybe_transposed_self, limit=limit,
                                    limit_direction=limit_direction,
                                    limit_area=limit_area,
                                    inplace=inplace, downcast=downcast,
                                    **kwargs)

        if inplace:
            if axis == 1:
                new_data = self._constructor(new_data).T._data
            self._update_inplace(new_data)
        else:
            res = self._constructor(new_data).__finalize__(self)
            if axis == 1:
                res = res.T
            return res

    # ----------------------------------------------------------------------
    # Timeseries methods Methods

    def asof(self, where, subset=None):
        """
        The last row without any NaN is taken (or the last row without
        NaN considering only the subset of columns in the case of a DataFrame)

        .. versionadded:: 0.19.0 For DataFrame

        If there is no good value, NaN is returned for a Series
        a Series of NaN values for a DataFrame

        Parameters
        ----------
        where : date or array of dates
        subset : string or list of strings, default None
           if not None use these columns for NaN propagation

        Notes
        -----
        Dates are assumed to be sorted
        Raises if this is not the case

        Returns
        -------
        where is scalar

          - value or NaN if input is Series
          - Series if input is DataFrame

        where is Index: same shape object as input

        See Also
        --------
        merge_asof

        """

        if isinstance(where, compat.string_types):
            from pandas import to_datetime
            where = to_datetime(where)

        if not self.index.is_monotonic:
            raise ValueError("asof requires a sorted index")

        is_series = isinstance(self, ABCSeries)
        if is_series:
            if subset is not None:
                raise ValueError("subset is not valid for Series")
        elif self.ndim > 2:
            raise NotImplementedError("asof is not implemented "
                                      "for {type}".format(type=type(self)))
        else:
            if subset is None:
                subset = self.columns
            if not is_list_like(subset):
                subset = [subset]

        is_list = is_list_like(where)
        if not is_list:
            start = self.index[0]
            if isinstance(self.index, PeriodIndex):
                where = Period(where, freq=self.index.freq).ordinal
                start = start.ordinal

            if where < start:
                if not is_series:
                    from pandas import Series
                    return Series(index=self.columns, name=where)
                return np.nan

            # It's always much faster to use a *while* loop here for
            # Series than pre-computing all the NAs. However a
            # *while* loop is extremely expensive for DataFrame
            # so we later pre-compute all the NAs and use the same
            # code path whether *where* is a scalar or list.
            # See PR: https://github.com/pandas-dev/pandas/pull/14476
            if is_series:
                loc = self.index.searchsorted(where, side='right')
                if loc > 0:
                    loc -= 1

                values = self._values
                while loc > 0 and isna(values[loc]):
                    loc -= 1
                return values[loc]

        if not isinstance(where, Index):
            where = Index(where) if is_list else Index([where])

        nulls = self.isna() if is_series else self[subset].isna().any(1)
        if nulls.all():
            if is_series:
                return self._constructor(np.nan, index=where, name=self.name)
            elif is_list:
                from pandas import DataFrame
                return DataFrame(np.nan, index=where, columns=self.columns)
            else:
                from pandas import Series
                return Series(np.nan, index=self.columns, name=where[0])

        locs = self.index.asof_locs(where, ~(nulls.values))

        # mask the missing
        missing = locs == -1
        data = self.take(locs, is_copy=False)
        data.index = where
        data.loc[missing] = np.nan
        return data if is_list else data.iloc[-1]

    # ----------------------------------------------------------------------
    # Action Methods

    _shared_docs['isna'] = """
        Detect missing values.

        Return a boolean same-sized object indicating if the values are NA.
        NA values, such as None or :attr:`numpy.NaN`, gets mapped to True
        values.
        Everything else gets mapped to False values. Characters such as empty
        strings ``''`` or :attr:`numpy.inf` are not considered NA values
        (unless you set ``pandas.options.mode.use_inf_as_na = True``).

        Returns
        -------
        %(klass)s
            Mask of bool values for each element in %(klass)s that
            indicates whether an element is not an NA value.

        See Also
        --------
        %(klass)s.isnull : alias of isna
        %(klass)s.notna : boolean inverse of isna
        %(klass)s.dropna : omit axes labels with missing values
        isna : top-level isna

        Examples
        --------
        Show which entries in a DataFrame are NA.

        >>> df = pd.DataFrame({'age': [5, 6, np.NaN],
        ...                    'born': [pd.NaT, pd.Timestamp('1939-05-27'),
        ...                             pd.Timestamp('1940-04-25')],
        ...                    'name': ['Alfred', 'Batman', ''],
        ...                    'toy': [None, 'Batmobile', 'Joker']})
        >>> df
           age       born    name        toy
        0  5.0        NaT  Alfred       None
        1  6.0 1939-05-27  Batman  Batmobile
        2  NaN 1940-04-25              Joker

        >>> df.isna()
             age   born   name    toy
        0  False   True  False   True
        1  False  False  False  False
        2   True  False  False  False

        Show which entries in a Series are NA.

        >>> ser = pd.Series([5, 6, np.NaN])
        >>> ser
        0    5.0
        1    6.0
        2    NaN
        dtype: float64

        >>> ser.isna()
        0    False
        1    False
        2     True
        dtype: bool
        """

    @Appender(_shared_docs['isna'] % _shared_doc_kwargs)
    def isna(self):
        return isna(self).__finalize__(self)

    @Appender(_shared_docs['isna'] % _shared_doc_kwargs)
    def isnull(self):
        return isna(self).__finalize__(self)

    _shared_docs['notna'] = """
        Detect existing (non-missing) values.

        Return a boolean same-sized object indicating if the values are not NA.
        Non-missing values get mapped to True. Characters such as empty
        strings ``''`` or :attr:`numpy.inf` are not considered NA values
        (unless you set ``pandas.options.mode.use_inf_as_na = True``).
        NA values, such as None or :attr:`numpy.NaN`, get mapped to False
        values.

        Returns
        -------
        %(klass)s
            Mask of bool values for each element in %(klass)s that
            indicates whether an element is not an NA value.

        See Also
        --------
        %(klass)s.notnull : alias of notna
        %(klass)s.isna : boolean inverse of notna
        %(klass)s.dropna : omit axes labels with missing values
        notna : top-level notna

        Examples
        --------
        Show which entries in a DataFrame are not NA.

        >>> df = pd.DataFrame({'age': [5, 6, np.NaN],
        ...                    'born': [pd.NaT, pd.Timestamp('1939-05-27'),
        ...                             pd.Timestamp('1940-04-25')],
        ...                    'name': ['Alfred', 'Batman', ''],
        ...                    'toy': [None, 'Batmobile', 'Joker']})
        >>> df
           age       born    name        toy
        0  5.0        NaT  Alfred       None
        1  6.0 1939-05-27  Batman  Batmobile
        2  NaN 1940-04-25              Joker

        >>> df.notna()
             age   born  name    toy
        0   True  False  True  False
        1   True   True  True   True
        2  False   True  True   True

        Show which entries in a Series are not NA.

        >>> ser = pd.Series([5, 6, np.NaN])
        >>> ser
        0    5.0
        1    6.0
        2    NaN
        dtype: float64

        >>> ser.notna()
        0     True
        1     True
        2    False
        dtype: bool
        """

    @Appender(_shared_docs['notna'] % _shared_doc_kwargs)
    def notna(self):
        return notna(self).__finalize__(self)

    @Appender(_shared_docs['notna'] % _shared_doc_kwargs)
    def notnull(self):
        return notna(self).__finalize__(self)

    def _clip_with_scalar(self, lower, upper, inplace=False):
        if ((lower is not None and np.any(isna(lower))) or
                (upper is not None and np.any(isna(upper)))):
            raise ValueError("Cannot use an NA value as a clip threshold")

        result = self.values
        mask = isna(result)

        with np.errstate(all='ignore'):
            if upper is not None:
                result = np.where(result >= upper, upper, result)
            if lower is not None:
                result = np.where(result <= lower, lower, result)
        if np.any(mask):
            result[mask] = np.nan

        axes_dict = self._construct_axes_dict()
        result = self._constructor(result, **axes_dict).__finalize__(self)

        if inplace:
            self._update_inplace(result)
        else:
            return result

    def _clip_with_one_bound(self, threshold, method, axis, inplace):

        inplace = validate_bool_kwarg(inplace, 'inplace')
        if axis is not None:
            axis = self._get_axis_number(axis)

        # method is self.le for upper bound and self.ge for lower bound
        if is_scalar(threshold) and is_number(threshold):
            if method.__name__ == 'le':
                return self._clip_with_scalar(None, threshold, inplace=inplace)
            return self._clip_with_scalar(threshold, None, inplace=inplace)

        subset = method(threshold, axis=axis) | isna(self)

        # GH #15390
        # In order for where method to work, the threshold must
        # be transformed to NDFrame from other array like structure.
        if (not isinstance(threshold, ABCSeries)) and is_list_like(threshold):
            if isinstance(self, ABCSeries):
                threshold = pd.Series(threshold, index=self.index)
            else:
                threshold = _align_method_FRAME(self, np.asarray(threshold),
                                                axis)
        return self.where(subset, threshold, axis=axis, inplace=inplace)

    def clip(self, lower=None, upper=None, axis=None, inplace=False,
             *args, **kwargs):
        """
        Trim values at input threshold(s).

        Assigns values outside boundary to boundary values. Thresholds
        can be singular values or array like, and in the latter case
        the clipping is performed element-wise in the specified axis.

        Parameters
        ----------
        lower : float or array_like, default None
            Minimum threshold value. All values below this
            threshold will be set to it.
        upper : float or array_like, default None
            Maximum threshold value. All values above this
            threshold will be set to it.
        axis : int or string axis name, optional
            Align object with lower and upper along the given axis.
        inplace : boolean, default False
            Whether to perform the operation in place on the data.

            .. versionadded:: 0.21.0
        *args, **kwargs
            Additional keywords have no effect but might be accepted
            for compatibility with numpy.

        See Also
        --------
        clip_lower : Clip values below specified threshold(s).
        clip_upper : Clip values above specified threshold(s).

        Returns
        -------
        Series or DataFrame
            Same type as calling object with the values outside the
            clip boundaries replaced

        Examples
        --------
        >>> data = {'col_0': [9, -3, 0, -1, 5], 'col_1': [-2, -7, 6, 8, -5]}
        >>> df = pd.DataFrame(data)
        >>> df
           col_0  col_1
        0      9     -2
        1     -3     -7
        2      0      6
        3     -1      8
        4      5     -5

        Clips per column using lower and upper thresholds:

        >>> df.clip(-4, 6)
           col_0  col_1
        0      6     -2
        1     -3     -4
        2      0      6
        3     -1      6
        4      5     -4

        Clips using specific lower and upper thresholds per column element:

        >>> t = pd.Series([2, -4, -1, 6, 3])
        >>> t
        0    2
        1   -4
        2   -1
        3    6
        4    3
        dtype: int64

        >>> df.clip(t, t + 4, axis=0)
           col_0  col_1
        0      6      2
        1     -3     -4
        2      0      3
        3      6      8
        4      5      3
        """
        if isinstance(self, ABCPanel):
            raise NotImplementedError("clip is not supported yet for panels")

        inplace = validate_bool_kwarg(inplace, 'inplace')

        axis = nv.validate_clip_with_axis(axis, args, kwargs)
        if axis is not None:
            axis = self._get_axis_number(axis)

        # GH 17276
        # numpy doesn't like NaN as a clip value
        # so ignore
        if np.any(pd.isnull(lower)):
            lower = None
        if np.any(pd.isnull(upper)):
            upper = None

        # GH 2747 (arguments were reversed)
        if lower is not None and upper is not None:
            if is_scalar(lower) and is_scalar(upper):
                lower, upper = min(lower, upper), max(lower, upper)

        # fast-path for scalars
        if ((lower is None or (is_scalar(lower) and is_number(lower))) and
                (upper is None or (is_scalar(upper) and is_number(upper)))):
            return self._clip_with_scalar(lower, upper, inplace=inplace)

        result = self
        if lower is not None:
            result = result.clip_lower(lower, axis, inplace=inplace)
        if upper is not None:
            if inplace:
                result = self
            result = result.clip_upper(upper, axis, inplace=inplace)

        return result

    def clip_upper(self, threshold, axis=None, inplace=False):
        """
        Return copy of input with values above given value(s) truncated.

        Parameters
        ----------
        threshold : float or array_like
        axis : int or string axis name, optional
            Align object with threshold along the given axis.
        inplace : boolean, default False
            Whether to perform the operation in place on the data

            .. versionadded:: 0.21.0

        See Also
        --------
        clip

        Returns
        -------
        clipped : same type as input
        """
        return self._clip_with_one_bound(threshold, method=self.le,
                                         axis=axis, inplace=inplace)

    def clip_lower(self, threshold, axis=None, inplace=False):
        """
        Return copy of the input with values below a threshold truncated.

        Parameters
        ----------
        threshold : numeric or array-like
            Minimum value allowed. All values below threshold will be set to
            this value.

            * float : every value is compared to `threshold`.
            * array-like : The shape of `threshold` should match the object
              it's compared to. When `self` is a Series, `threshold` should be
              the length. When `self` is a DataFrame, `threshold` should 2-D
              and the same shape as `self` for ``axis=None``, or 1-D and the
              same length as the axis being compared.

        axis : {0 or 'index', 1 or 'columns'}, default 0
            Align `self` with `threshold` along the given axis.

        inplace : boolean, default False
            Whether to perform the operation in place on the data.

            .. versionadded:: 0.21.0

        See Also
        --------
        Series.clip : Return copy of input with values below and above
            thresholds truncated.
        Series.clip_upper : Return copy of input with values above
            threshold truncated.

        Returns
        -------
        clipped : same type as input

        Examples
        --------
        Series single threshold clipping:

        >>> s = pd.Series([5, 6, 7, 8, 9])
        >>> s.clip_lower(8)
        0    8
        1    8
        2    8
        3    8
        4    9
        dtype: int64

        Series clipping element-wise using an array of thresholds. `threshold`
        should be the same length as the Series.

        >>> elemwise_thresholds = [4, 8, 7, 2, 5]
        >>> s.clip_lower(elemwise_thresholds)
        0    5
        1    8
        2    7
        3    8
        4    9
        dtype: int64

        DataFrames can be compared to a scalar.

        >>> df = pd.DataFrame({"A": [1, 3, 5], "B": [2, 4, 6]})
        >>> df
           A  B
        0  1  2
        1  3  4
        2  5  6

        >>> df.clip_lower(3)
           A  B
        0  3  3
        1  3  4
        2  5  6

        Or to an array of values. By default, `threshold` should be the same
        shape as the DataFrame.

        >>> df.clip_lower(np.array([[3, 4], [2, 2], [6, 2]]))
           A  B
        0  3  4
        1  3  4
        2  6  6

        Control how `threshold` is broadcast with `axis`. In this case
        `threshold` should be the same length as the axis specified by
        `axis`.

        >>> df.clip_lower(np.array([3, 3, 5]), axis='index')
           A  B
        0  3  3
        1  3  4
        2  5  6

        >>> df.clip_lower(np.array([4, 5]), axis='columns')
           A  B
        0  4  5
        1  4  5
        2  5  6
        """
        return self._clip_with_one_bound(threshold, method=self.ge,
                                         axis=axis, inplace=inplace)

    def groupby(self, by=None, axis=0, level=None, as_index=True, sort=True,
                group_keys=True, squeeze=False, observed=False, **kwargs):
        """
        Group series using mapper (dict or key function, apply given function
        to group, return result as series) or by a series of columns.

        Parameters
        ----------
        by : mapping, function, label, or list of labels
            Used to determine the groups for the groupby.
            If ``by`` is a function, it's called on each value of the object's
            index. If a dict or Series is passed, the Series or dict VALUES
            will be used to determine the groups (the Series' values are first
            aligned; see ``.align()`` method). If an ndarray is passed, the
            values are used as-is determine the groups. A label or list of
            labels may be passed to group by the columns in ``self``. Notice
            that a tuple is interpreted a (single) key.
        axis : int, default 0
        level : int, level name, or sequence of such, default None
            If the axis is a MultiIndex (hierarchical), group by a particular
            level or levels
        as_index : boolean, default True
            For aggregated output, return object with group labels as the
            index. Only relevant for DataFrame input. as_index=False is
            effectively "SQL-style" grouped output
        sort : boolean, default True
            Sort group keys. Get better performance by turning this off.
            Note this does not influence the order of observations within each
            group.  groupby preserves the order of rows within each group.
        group_keys : boolean, default True
            When calling apply, add group keys to index to identify pieces
        squeeze : boolean, default False
            reduce the dimensionality of the return type if possible,
            otherwise return a consistent type
        observed : boolean, default False
            This only applies if any of the groupers are Categoricals
            If True: only show observed values for categorical groupers.
            If False: show all values for categorical groupers.

            .. versionadded:: 0.23.0

        Returns
        -------
        GroupBy object

        Examples
        --------
        DataFrame results

        >>> data.groupby(func, axis=0).mean()
        >>> data.groupby(['col1', 'col2'])['col3'].mean()

        DataFrame with hierarchical index

        >>> data.groupby(['col1', 'col2']).mean()

        Notes
        -----
        See the `user guide
        <http://pandas.pydata.org/pandas-docs/stable/groupby.html>`_ for more.

        See also
        --------
        resample : Convenience method for frequency conversion and resampling
            of time series.
        """
        from pandas.core.groupby.groupby import groupby

        if level is None and by is None:
            raise TypeError("You have to supply one of 'by' and 'level'")
        axis = self._get_axis_number(axis)
        return groupby(self, by=by, axis=axis, level=level, as_index=as_index,
                       sort=sort, group_keys=group_keys, squeeze=squeeze,
                       observed=observed, **kwargs)

    def asfreq(self, freq, method=None, how=None, normalize=False,
               fill_value=None):
        """
        Convert TimeSeries to specified frequency.

        Optionally provide filling method to pad/backfill missing values.

        Returns the original data conformed to a new index with the specified
        frequency. ``resample`` is more appropriate if an operation, such as
        summarization, is necessary to represent the data at the new frequency.

        Parameters
        ----------
        freq : DateOffset object, or string
        method : {'backfill'/'bfill', 'pad'/'ffill'}, default None
            Method to use for filling holes in reindexed Series (note this
            does not fill NaNs that already were present):

            * 'pad' / 'ffill': propagate last valid observation forward to next
              valid
            * 'backfill' / 'bfill': use NEXT valid observation to fill
        how : {'start', 'end'}, default end
            For PeriodIndex only, see PeriodIndex.asfreq
        normalize : bool, default False
            Whether to reset output index to midnight
        fill_value: scalar, optional
            Value to use for missing values, applied during upsampling (note
            this does not fill NaNs that already were present).

            .. versionadded:: 0.20.0

        Returns
        -------
        converted : type of caller

        Examples
        --------

        Start by creating a series with 4 one minute timestamps.

        >>> index = pd.date_range('1/1/2000', periods=4, freq='T')
        >>> series = pd.Series([0.0, None, 2.0, 3.0], index=index)
        >>> df = pd.DataFrame({'s':series})
        >>> df
                               s
        2000-01-01 00:00:00    0.0
        2000-01-01 00:01:00    NaN
        2000-01-01 00:02:00    2.0
        2000-01-01 00:03:00    3.0

        Upsample the series into 30 second bins.

        >>> df.asfreq(freq='30S')
                               s
        2000-01-01 00:00:00    0.0
        2000-01-01 00:00:30    NaN
        2000-01-01 00:01:00    NaN
        2000-01-01 00:01:30    NaN
        2000-01-01 00:02:00    2.0
        2000-01-01 00:02:30    NaN
        2000-01-01 00:03:00    3.0

        Upsample again, providing a ``fill value``.

        >>> df.asfreq(freq='30S', fill_value=9.0)
                               s
        2000-01-01 00:00:00    0.0
        2000-01-01 00:00:30    9.0
        2000-01-01 00:01:00    NaN
        2000-01-01 00:01:30    9.0
        2000-01-01 00:02:00    2.0
        2000-01-01 00:02:30    9.0
        2000-01-01 00:03:00    3.0

        Upsample again, providing a ``method``.

        >>> df.asfreq(freq='30S', method='bfill')
                               s
        2000-01-01 00:00:00    0.0
        2000-01-01 00:00:30    NaN
        2000-01-01 00:01:00    NaN
        2000-01-01 00:01:30    2.0
        2000-01-01 00:02:00    2.0
        2000-01-01 00:02:30    3.0
        2000-01-01 00:03:00    3.0

        See Also
        --------
        reindex

        Notes
        -----
        To learn more about the frequency strings, please see `this link
        <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.
        """
        from pandas.core.resample import asfreq
        return asfreq(self, freq, method=method, how=how, normalize=normalize,
                      fill_value=fill_value)

    def at_time(self, time, asof=False):
        """
        Select values at particular time of day (e.g. 9:30AM).

        Raises
        ------
        TypeError
            If the index is not  a :class:`DatetimeIndex`

        Parameters
        ----------
        time : datetime.time or string

        Returns
        -------
        values_at_time : type of caller

        Examples
        --------
        >>> i = pd.date_range('2018-04-09', periods=4, freq='12H')
        >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
        >>> ts
                             A
        2018-04-09 00:00:00  1
        2018-04-09 12:00:00  2
        2018-04-10 00:00:00  3
        2018-04-10 12:00:00  4

        >>> ts.at_time('12:00')
                             A
        2018-04-09 12:00:00  2
        2018-04-10 12:00:00  4

        See Also
        --------
        between_time : Select values between particular times of the day
        first : Select initial periods of time series based on a date offset
        last : Select final periods of time series based on a date offset
        DatetimeIndex.indexer_at_time : Get just the index locations for
            values at particular time of the day
        """
        try:
            indexer = self.index.indexer_at_time(time, asof=asof)
            return self._take(indexer)
        except AttributeError:
            raise TypeError('Index must be DatetimeIndex')

    def between_time(self, start_time, end_time, include_start=True,
                     include_end=True):
        """
        Select values between particular times of the day (e.g., 9:00-9:30 AM).

        By setting ``start_time`` to be later than ``end_time``,
        you can get the times that are *not* between the two times.

        Raises
        ------
        TypeError
            If the index is not  a :class:`DatetimeIndex`

        Parameters
        ----------
        start_time : datetime.time or string
        end_time : datetime.time or string
        include_start : boolean, default True
        include_end : boolean, default True

        Returns
        -------
        values_between_time : type of caller

        Examples
        --------
        >>> i = pd.date_range('2018-04-09', periods=4, freq='1D20min')
        >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
        >>> ts
                             A
        2018-04-09 00:00:00  1
        2018-04-10 00:20:00  2
        2018-04-11 00:40:00  3
        2018-04-12 01:00:00  4

        >>> ts.between_time('0:15', '0:45')
                             A
        2018-04-10 00:20:00  2
        2018-04-11 00:40:00  3

        You get the times that are *not* between two times by setting
        ``start_time`` later than ``end_time``:

        >>> ts.between_time('0:45', '0:15')
                             A
        2018-04-09 00:00:00  1
        2018-04-12 01:00:00  4

        See Also
        --------
        at_time : Select values at a particular time of the day
        first : Select initial periods of time series based on a date offset
        last : Select final periods of time series based on a date offset
        DatetimeIndex.indexer_between_time : Get just the index locations for
            values between particular times of the day
        """
        try:
            indexer = self.index.indexer_between_time(
                start_time, end_time, include_start=include_start,
                include_end=include_end)
            return self._take(indexer)
        except AttributeError:
            raise TypeError('Index must be DatetimeIndex')

    def resample(self, rule, how=None, axis=0, fill_method=None, closed=None,
                 label=None, convention='start', kind=None, loffset=None,
                 limit=None, base=0, on=None, level=None):
        """
        Convenience method for frequency conversion and resampling of time
        series.  Object must have a datetime-like index (DatetimeIndex,
        PeriodIndex, or TimedeltaIndex), or pass datetime-like values
        to the on or level keyword.

        Parameters
        ----------
        rule : string
            the offset string or object representing target conversion
        axis : int, optional, default 0
        closed : {'right', 'left'}
            Which side of bin interval is closed. The default is 'left'
            for all frequency offsets except for 'M', 'A', 'Q', 'BM',
            'BA', 'BQ', and 'W' which all have a default of 'right'.
        label : {'right', 'left'}
            Which bin edge label to label bucket with. The default is 'left'
            for all frequency offsets except for 'M', 'A', 'Q', 'BM',
            'BA', 'BQ', and 'W' which all have a default of 'right'.
        convention : {'start', 'end', 's', 'e'}
            For PeriodIndex only, controls whether to use the start or end of
            `rule`
        kind: {'timestamp', 'period'}, optional
            Pass 'timestamp' to convert the resulting index to a
            ``DateTimeIndex`` or 'period' to convert it to a ``PeriodIndex``.
            By default the input representation is retained.
        loffset : timedelta
            Adjust the resampled time labels
        base : int, default 0
            For frequencies that evenly subdivide 1 day, the "origin" of the
            aggregated intervals. For example, for '5min' frequency, base could
            range from 0 through 4. Defaults to 0
        on : string, optional
            For a DataFrame, column to use instead of index for resampling.
            Column must be datetime-like.

            .. versionadded:: 0.19.0

        level : string or int, optional
            For a MultiIndex, level (name or number) to use for
            resampling.  Level must be datetime-like.

            .. versionadded:: 0.19.0

        Returns
        -------
        Resampler object

        Notes
        -----
        See the `user guide
        <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#resampling>`_
        for more.

        To learn more about the offset strings, please see `this link
        <http://pandas.pydata.org/pandas-docs/stable/timeseries.html#offset-aliases>`__.

        Examples
        --------

        Start by creating a series with 9 one minute timestamps.

        >>> index = pd.date_range('1/1/2000', periods=9, freq='T')
        >>> series = pd.Series(range(9), index=index)
        >>> series
        2000-01-01 00:00:00    0
        2000-01-01 00:01:00    1
        2000-01-01 00:02:00    2
        2000-01-01 00:03:00    3
        2000-01-01 00:04:00    4
        2000-01-01 00:05:00    5
        2000-01-01 00:06:00    6
        2000-01-01 00:07:00    7
        2000-01-01 00:08:00    8
        Freq: T, dtype: int64

        Downsample the series into 3 minute bins and sum the values
        of the timestamps falling into a bin.

        >>> series.resample('3T').sum()
        2000-01-01 00:00:00     3
        2000-01-01 00:03:00    12
        2000-01-01 00:06:00    21
        Freq: 3T, dtype: int64

        Downsample the series into 3 minute bins as above, but label each
        bin using the right edge instead of the left. Please note that the
        value in the bucket used as the label is not included in the bucket,
        which it labels. For example, in the original series the
        bucket ``2000-01-01 00:03:00`` contains the value 3, but the summed
        value in the resampled bucket with the label ``2000-01-01 00:03:00``
        does not include 3 (if it did, the summed value would be 6, not 3).
        To include this value close the right side of the bin interval as
        illustrated in the example below this one.

        >>> series.resample('3T', label='right').sum()
        2000-01-01 00:03:00     3
        2000-01-01 00:06:00    12
        2000-01-01 00:09:00    21
        Freq: 3T, dtype: int64

        Downsample the series into 3 minute bins as above, but close the right
        side of the bin interval.

        >>> series.resample('3T', label='right', closed='right').sum()
        2000-01-01 00:00:00     0
        2000-01-01 00:03:00     6
        2000-01-01 00:06:00    15
        2000-01-01 00:09:00    15
        Freq: 3T, dtype: int64

        Upsample the series into 30 second bins.

        >>> series.resample('30S').asfreq()[0:5] #select first 5 rows
        2000-01-01 00:00:00   0.0
        2000-01-01 00:00:30   NaN
        2000-01-01 00:01:00   1.0
        2000-01-01 00:01:30   NaN
        2000-01-01 00:02:00   2.0
        Freq: 30S, dtype: float64

        Upsample the series into 30 second bins and fill the ``NaN``
        values using the ``pad`` method.

        >>> series.resample('30S').pad()[0:5]
        2000-01-01 00:00:00    0
        2000-01-01 00:00:30    0
        2000-01-01 00:01:00    1
        2000-01-01 00:01:30    1
        2000-01-01 00:02:00    2
        Freq: 30S, dtype: int64

        Upsample the series into 30 second bins and fill the
        ``NaN`` values using the ``bfill`` method.

        >>> series.resample('30S').bfill()[0:5]
        2000-01-01 00:00:00    0
        2000-01-01 00:00:30    1
        2000-01-01 00:01:00    1
        2000-01-01 00:01:30    2
        2000-01-01 00:02:00    2
        Freq: 30S, dtype: int64

        Pass a custom function via ``apply``

        >>> def custom_resampler(array_like):
        ...     return np.sum(array_like)+5

        >>> series.resample('3T').apply(custom_resampler)
        2000-01-01 00:00:00     8
        2000-01-01 00:03:00    17
        2000-01-01 00:06:00    26
        Freq: 3T, dtype: int64

        For a Series with a PeriodIndex, the keyword `convention` can be
        used to control whether to use the start or end of `rule`.

        >>> s = pd.Series([1, 2], index=pd.period_range('2012-01-01',
                                                        freq='A',
                                                        periods=2))
        >>> s
        2012    1
        2013    2
        Freq: A-DEC, dtype: int64

        Resample by month using 'start' `convention`. Values are assigned to
        the first month of the period.

        >>> s.resample('M', convention='start').asfreq().head()
        2012-01    1.0
        2012-02    NaN
        2012-03    NaN
        2012-04    NaN
        2012-05    NaN
        Freq: M, dtype: float64

        Resample by month using 'end' `convention`. Values are assigned to
        the last month of the period.

        >>> s.resample('M', convention='end').asfreq()
        2012-12    1.0
        2013-01    NaN
        2013-02    NaN
        2013-03    NaN
        2013-04    NaN
        2013-05    NaN
        2013-06    NaN
        2013-07    NaN
        2013-08    NaN
        2013-09    NaN
        2013-10    NaN
        2013-11    NaN
        2013-12    2.0
        Freq: M, dtype: float64

        For DataFrame objects, the keyword ``on`` can be used to specify the
        column instead of the index for resampling.

        >>> df = pd.DataFrame(data=9*[range(4)], columns=['a', 'b', 'c', 'd'])
        >>> df['time'] = pd.date_range('1/1/2000', periods=9, freq='T')
        >>> df.resample('3T', on='time').sum()
                             a  b  c  d
        time
        2000-01-01 00:00:00  0  3  6  9
        2000-01-01 00:03:00  0  3  6  9
        2000-01-01 00:06:00  0  3  6  9

        For a DataFrame with MultiIndex, the keyword ``level`` can be used to
        specify on level the resampling needs to take place.

        >>> time = pd.date_range('1/1/2000', periods=5, freq='T')
        >>> df2 = pd.DataFrame(data=10*[range(4)],
                               columns=['a', 'b', 'c', 'd'],
                               index=pd.MultiIndex.from_product([time, [1, 2]])
                               )
        >>> df2.resample('3T', level=0).sum()
                             a  b   c   d
        2000-01-01 00:00:00  0  6  12  18
        2000-01-01 00:03:00  0  4   8  12

        See also
        --------
        groupby : Group by mapping, function, label, or list of labels.
        """
        from pandas.core.resample import (resample,
                                          _maybe_process_deprecations)
        axis = self._get_axis_number(axis)
        r = resample(self, freq=rule, label=label, closed=closed,
                     axis=axis, kind=kind, loffset=loffset,
                     convention=convention,
                     base=base, key=on, level=level)
        return _maybe_process_deprecations(r,
                                           how=how,
                                           fill_method=fill_method,
                                           limit=limit)

    def first(self, offset):
        """
        Convenience method for subsetting initial periods of time series data
        based on a date offset.

        Raises
        ------
        TypeError
            If the index is not  a :class:`DatetimeIndex`

        Parameters
        ----------
        offset : string, DateOffset, dateutil.relativedelta

        Examples
        --------
        >>> i = pd.date_range('2018-04-09', periods=4, freq='2D')
        >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
        >>> ts
                    A
        2018-04-09  1
        2018-04-11  2
        2018-04-13  3
        2018-04-15  4

        Get the rows for the first 3 days:

        >>> ts.first('3D')
                    A
        2018-04-09  1
        2018-04-11  2

        Notice the data for 3 first calender days were returned, not the first
        3 days observed in the dataset, and therefore data for 2018-04-13 was
        not returned.

        Returns
        -------
        subset : type of caller

        See Also
        --------
        last : Select final periods of time series based on a date offset
        at_time : Select values at a particular time of the day
        between_time : Select values between particular times of the day
        """
        from pandas.tseries.frequencies import to_offset
        if not isinstance(self.index, DatetimeIndex):
            raise TypeError("'first' only supports a DatetimeIndex index")

        if len(self.index) == 0:
            return self

        offset = to_offset(offset)
        end_date = end = self.index[0] + offset

        # Tick-like, e.g. 3 weeks
        if not offset.isAnchored() and hasattr(offset, '_inc'):
            if end_date in self.index:
                end = self.index.searchsorted(end_date, side='left')
                return self.iloc[:end]

        return self.loc[:end]

    def last(self, offset):
        """
        Convenience method for subsetting final periods of time series data
        based on a date offset.

        Raises
        ------
        TypeError
            If the index is not  a :class:`DatetimeIndex`

        Parameters
        ----------
        offset : string, DateOffset, dateutil.relativedelta

        Examples
        --------
        >>> i = pd.date_range('2018-04-09', periods=4, freq='2D')
        >>> ts = pd.DataFrame({'A': [1,2,3,4]}, index=i)
        >>> ts
                    A
        2018-04-09  1
        2018-04-11  2
        2018-04-13  3
        2018-04-15  4

        Get the rows for the last 3 days:

        >>> ts.last('3D')
                    A
        2018-04-13  3
        2018-04-15  4

        Notice the data for 3 last calender days were returned, not the last
        3 observed days in the dataset, and therefore data for 2018-04-11 was
        not returned.

        Returns
        -------
        subset : type of caller

        See Also
        --------
        first : Select initial periods of time series based on a date offset
        at_time : Select values at a particular time of the day
        between_time : Select values between particular times of the day
        """
        from pandas.tseries.frequencies import to_offset
        if not isinstance(self.index, DatetimeIndex):
            raise TypeError("'last' only supports a DatetimeIndex index")

        if len(self.index) == 0:
            return self

        offset = to_offset(offset)

        start_date = self.index[-1] - offset
        start = self.index.searchsorted(start_date, side='right')
        return self.iloc[start:]

    def rank(self, axis=0, method='average', numeric_only=None,
             na_option='keep', ascending=True, pct=False):
        """
        Compute numerical data ranks (1 through n) along axis. Equal values are
        assigned a rank that is the average of the ranks of those values

        Parameters
        ----------
        axis : {0 or 'index', 1 or 'columns'}, default 0
            index to direct ranking
        method : {'average', 'min', 'max', 'first', 'dense'}
            * average: average rank of group
            * min: lowest rank in group
            * max: highest rank in group
            * first: ranks assigned in order they appear in the array
            * dense: like 'min', but rank always increases by 1 between groups
        numeric_only : boolean, default None
            Include only float, int, boolean data. Valid only for DataFrame or
            Panel objects
        na_option : {'keep', 'top', 'bottom'}
            * keep: leave NA values where they are
            * top: smallest rank if ascending
            * bottom: smallest rank if descending
        ascending : boolean, default True
            False for ranks by high (1) to low (N)
        pct : boolean, default False
            Computes percentage rank of data

        Returns
        -------
        ranks : same type as caller
        """
        axis = self._get_axis_number(axis)

        if self.ndim > 2:
            msg = "rank does not make sense when ndim > 2"
            raise NotImplementedError(msg)

        def ranker(data):
            ranks = algos.rank(data.values, axis=axis, method=method,
                               ascending=ascending, na_option=na_option,
                               pct=pct)
            ranks = self._constructor(ranks, **data._construct_axes_dict())
            return ranks.__finalize__(self)

        # if numeric_only is None, and we can't get anything, we try with
        # numeric_only=True
        if numeric_only is None:
            try:
                return ranker(self)
            except TypeError:
                numeric_only = True

        if numeric_only:
            data = self._get_numeric_data()
        else:
            data = self

        return ranker(data)

    _shared_docs['align'] = ("""
        Align two objects on their axes with the
        specified join method for each axis Index

        Parameters
        ----------
        other : DataFrame or Series
        join : {'outer', 'inner', 'left', 'right'}, default 'outer'
        axis : allowed axis of the other object, default None
            Align on index (0), columns (1), or both (None)
        level : int or level name, default None
            Broadcast across a level, matching Index values on the
            passed MultiIndex level
        copy : boolean, default True
            Always returns new objects. If copy=False and no reindexing is
            required then original objects are returned.
        fill_value : scalar, default np.NaN
            Value to use for missing values. Defaults to NaN, but can be any
            "compatible" value
        method : str, default None
        limit : int, default None
        fill_axis : %(axes_single_arg)s, default 0
            Filling axis, method and limit
        broadcast_axis : %(axes_single_arg)s, default None
            Broadcast values along this axis, if aligning two objects of
            different dimensions

        Returns
        -------
        (left, right) : (%(klass)s, type of other)
            Aligned objects
        """)

    @Appender(_shared_docs['align'] % _shared_doc_kwargs)
    def align(self, other, join='outer', axis=None, level=None, copy=True,
              fill_value=None, method=None, limit=None, fill_axis=0,
              broadcast_axis=None):
        from pandas import DataFrame, Series
        method = missing.clean_fill_method(method)

        if broadcast_axis == 1 and self.ndim != other.ndim:
            if isinstance(self, Series):
                # this means other is a DataFrame, and we need to broadcast
                # self
                cons = self._constructor_expanddim
                df = cons({c: self for c in other.columns},
                          **other._construct_axes_dict())
                return df._align_frame(other, join=join, axis=axis,
                                       level=level, copy=copy,
                                       fill_value=fill_value, method=method,
                                       limit=limit, fill_axis=fill_axis)
            elif isinstance(other, Series):
                # this means self is a DataFrame, and we need to broadcast
                # other
                cons = other._constructor_expanddim
                df = cons({c: other for c in self.columns},
                          **self._construct_axes_dict())
                return self._align_frame(df, join=join, axis=axis, level=level,
                                         copy=copy, fill_value=fill_value,
                                         method=method, limit=limit,
                                         fill_axis=fill_axis)

        if axis is not None:
            axis = self._get_axis_number(axis)
        if isinstance(other, DataFrame):
            return self._align_frame(other, join=join, axis=axis, level=level,
                                     copy=copy, fill_value=fill_value,
                                     method=method, limit=limit,
                                     fill_axis=fill_axis)
        elif isinstance(other, Series):
            return self._align_series(other, join=join, axis=axis, level=level,
                                      copy=copy, fill_value=fill_value,
                                      method=method, limit=limit,
                                      fill_axis=fill_axis)
        else:  # pragma: no cover
            raise TypeError('unsupported type: %s' % type(other))

    def _align_frame(self, other, join='outer', axis=None, level=None,
                     copy=True, fill_value=None, method=None, limit=None,
                     fill_axis=0):
        # defaults
        join_index, join_columns = None, None
        ilidx, iridx = None, None
        clidx, cridx = None, None

        is_series = isinstance(self, ABCSeries)

        if axis is None or axis == 0:
            if not self.index.equals(other.index):
                join_index, ilidx, iridx = self.index.join(
                    other.index, how=join, level=level, return_indexers=True)

        if axis is None or axis == 1:
            if not is_series and not self.columns.equals(other.columns):
                join_columns, clidx, cridx = self.columns.join(
                    other.columns, how=join, level=level, return_indexers=True)

        if is_series:
            reindexers = {0: [join_index, ilidx]}
        else:
            reindexers = {0: [join_index, ilidx], 1: [join_columns, clidx]}

        left = self._reindex_with_indexers(reindexers, copy=copy,
                                           fill_value=fill_value,
                                           allow_dups=True)
        # other must be always DataFrame
        right = other._reindex_with_indexers({0: [join_index, iridx],
                                              1: [join_columns, cridx]},
                                             copy=copy, fill_value=fill_value,
                                             allow_dups=True)

        if method is not None:
            left = left.fillna(axis=fill_axis, method=method, limit=limit)
            right = right.fillna(axis=fill_axis, method=method, limit=limit)

        # if DatetimeIndex have different tz, convert to UTC
        if is_datetime64tz_dtype(left.index):
            if left.index.tz != right.index.tz:
                if join_index is not None:
                    left.index = join_index
                    right.index = join_index

        return left.__finalize__(self), right.__finalize__(other)

    def _align_series(self, other, join='outer', axis=None, level=None,
                      copy=True, fill_value=None, method=None, limit=None,
                      fill_axis=0):

        is_series = isinstance(self, ABCSeries)

        # series/series compat, other must always be a Series
        if is_series:
            if axis:
                raise ValueError('cannot align series to a series other than '
                                 'axis 0')

            # equal
            if self.index.equals(other.index):
                join_index, lidx, ridx = None, None, None
            else:
                join_index, lidx, ridx = self.index.join(other.index, how=join,
                                                         level=level,
                                                         return_indexers=True)

            left = self._reindex_indexer(join_index, lidx, copy)
            right = other._reindex_indexer(join_index, ridx, copy)

        else:
            # one has > 1 ndim
            fdata = self._data
            if axis == 0:
                join_index = self.index
                lidx, ridx = None, None
                if not self.index.equals(other.index):
                    join_index, lidx, ridx = self.index.join(
                        other.index, how=join, level=level,
                        return_indexers=True)

                if lidx is not None:
                    fdata = fdata.reindex_indexer(join_index, lidx, axis=1)

            elif axis == 1:
                join_index = self.columns
                lidx, ridx = None, None
                if not self.columns.equals(other.index):
                    join_index, lidx, ridx = self.columns.join(
                        other.index, how=join, level=level,
                        return_indexers=True)

                if lidx is not None:
                    fdata = fdata.reindex_indexer(join_index, lidx, axis=0)
            else:
                raise ValueError('Must specify axis=0 or 1')

            if copy and fdata is self._data:
                fdata = fdata.copy()

            left = self._constructor(fdata)

            if ridx is None:
                right = other
            else:
                right = other.reindex(join_index, level=level)

        # fill
        fill_na = notna(fill_value) or (method is not None)
        if fill_na:
            left = left.fillna(fill_value, method=method, limit=limit,
                               axis=fill_axis)
            right = right.fillna(fill_value, method=method, limit=limit)

        # if DatetimeIndex have different tz, convert to UTC
        if is_series or (not is_series and axis == 0):
            if is_datetime64tz_dtype(left.index):
                if left.index.tz != right.index.tz:
                    if join_index is not None:
                        left.index = join_index
                        right.index = join_index

        return left.__finalize__(self), right.__finalize__(other)

    def _where(self, cond, other=np.nan, inplace=False, axis=None, level=None,
               errors='raise', try_cast=False):
        """
        Equivalent to public method `where`, except that `other` is not
        applied as a function even if callable. Used in __setitem__.
        """
        inplace = validate_bool_kwarg(inplace, 'inplace')

        # align the cond to same shape as myself
        cond = com._apply_if_callable(cond, self)
        if isinstance(cond, NDFrame):
            cond, _ = cond.align(self, join='right', broadcast_axis=1)
        else:
            if not hasattr(cond, 'shape'):
                cond = np.asanyarray(cond)
            if cond.shape != self.shape:
                raise ValueError('Array conditional must be same shape as '
                                 'self')
            cond = self._constructor(cond, **self._construct_axes_dict())

        # make sure we are boolean
        fill_value = True if inplace else False
        cond = cond.fillna(fill_value)

        msg = "Boolean array expected for the condition, not {dtype}"

        if not isinstance(cond, pd.DataFrame):
            # This is a single-dimensional object.
            if not is_bool_dtype(cond):
                raise ValueError(msg.format(dtype=cond.dtype))
        else:
            for dt in cond.dtypes:
                if not is_bool_dtype(dt):
                    raise ValueError(msg.format(dtype=dt))

        cond = -cond if inplace else cond

        # try to align with other
        try_quick = True
        if hasattr(other, 'align'):

            # align with me
            if other.ndim <= self.ndim:

                _, other = self.align(other, join='left', axis=axis,
                                      level=level, fill_value=np.nan)

                # if we are NOT aligned, raise as we cannot where index
                if (axis is None and
                        not all(other._get_axis(i).equals(ax)
                                for i, ax in enumerate(self.axes))):
                    raise InvalidIndexError

            # slice me out of the other
            else:
                raise NotImplementedError("cannot align with a higher "
                                          "dimensional NDFrame")

        if isinstance(other, np.ndarray):

            if other.shape != self.shape:

                if self.ndim == 1:

                    icond = cond.values

                    # GH 2745 / GH 4192
                    # treat like a scalar
                    if len(other) == 1:
                        other = np.array(other[0])

                    # GH 3235
                    # match True cond to other
                    elif len(cond[icond]) == len(other):

                        # try to not change dtype at first (if try_quick)
                        if try_quick:

                            try:
                                new_other = com._values_from_object(self)
                                new_other = new_other.copy()
                                new_other[icond] = other
                                other = new_other
                            except Exception:
                                try_quick = False

                        # let's create a new (if we failed at the above
                        # or not try_quick
                        if not try_quick:

                            dtype, fill_value = maybe_promote(other.dtype)
                            new_other = np.empty(len(icond), dtype=dtype)
                            new_other.fill(fill_value)
                            maybe_upcast_putmask(new_other, icond, other)
                            other = new_other

                    else:
                        raise ValueError('Length of replacements must equal '
                                         'series length')

                else:
                    raise ValueError('other must be the same shape as self '
                                     'when an ndarray')

            # we are the same shape, so create an actual object for alignment
            else:
                other = self._constructor(other, **self._construct_axes_dict())

        if axis is None:
            axis = 0

        if self.ndim == getattr(other, 'ndim', 0):
            align = True
        else:
            align = (self._get_axis_number(axis) == 1)

        block_axis = self._get_block_manager_axis(axis)

        if inplace:
            # we may have different type blocks come out of putmask, so
            # reconstruct the block manager

            self._check_inplace_setting(other)
            new_data = self._data.putmask(mask=cond, new=other, align=align,
                                          inplace=True, axis=block_axis,
                                          transpose=self._AXIS_REVERSED)
            self._update_inplace(new_data)

        else:
            new_data = self._data.where(other=other, cond=cond, align=align,
                                        errors=errors,
                                        try_cast=try_cast, axis=block_axis,
                                        transpose=self._AXIS_REVERSED)

            return self._constructor(new_data).__finalize__(self)

    _shared_docs['where'] = ("""
        Return an object of same shape as self and whose corresponding
        entries are from self where `cond` is %(cond)s and otherwise are from
        `other`.

        Parameters
        ----------
        cond : boolean %(klass)s, array-like, or callable
            Where `cond` is %(cond)s, keep the original value. Where
            %(cond_rev)s, replace with corresponding value from `other`.
            If `cond` is callable, it is computed on the %(klass)s and
            should return boolean %(klass)s or array. The callable must
            not change input %(klass)s (though pandas doesn't check it).

            .. versionadded:: 0.18.1
                A callable can be used as cond.

        other : scalar, %(klass)s, or callable
            Entries where `cond` is %(cond_rev)s are replaced with
            corresponding value from `other`.
            If other is callable, it is computed on the %(klass)s and
            should return scalar or %(klass)s. The callable must not
            change input %(klass)s (though pandas doesn't check it).

            .. versionadded:: 0.18.1
                A callable can be used as other.

        inplace : boolean, default False
            Whether to perform the operation in place on the data
        axis : alignment axis if needed, default None
        level : alignment level if needed, default None
        errors : str, {'raise', 'ignore'}, default 'raise'
            - ``raise`` : allow exceptions to be raised
            - ``ignore`` : suppress exceptions. On error return original object

            Note that currently this parameter won't affect
            the results and will always coerce to a suitable dtype.

        try_cast : boolean, default False
            try to cast the result back to the input type (if possible),
        raise_on_error : boolean, default True
            Whether to raise on invalid data types (e.g. trying to where on
            strings)

            .. deprecated:: 0.21.0

        Returns
        -------
        wh : same type as caller

        Notes
        -----
        The %(name)s method is an application of the if-then idiom. For each
        element in the calling DataFrame, if ``cond`` is ``%(cond)s`` the
        element is used; otherwise the corresponding element from the DataFrame
        ``other`` is used.

        The signature for :func:`DataFrame.where` differs from
        :func:`numpy.where`. Roughly ``df1.where(m, df2)`` is equivalent to
        ``np.where(m, df1, df2)``.

        For further details and examples see the ``%(name)s`` documentation in
        :ref:`indexing <indexing.where_mask>`.

        Examples
        --------
        >>> s = pd.Series(range(5))
        >>> s.where(s > 0)
        0    NaN
        1    1.0
        2    2.0
        3    3.0
        4    4.0

        >>> s.mask(s > 0)
        0    0.0
        1    NaN
        2    NaN
        3    NaN
        4    NaN

        >>> s.where(s > 1, 10)
        0    10.0
        1    10.0
        2    2.0
        3    3.0
        4    4.0

        >>> df = pd.DataFrame(np.arange(10).reshape(-1, 2), columns=['A', 'B'])
        >>> m = df %% 3 == 0
        >>> df.where(m, -df)
           A  B
        0  0 -1
        1 -2  3
        2 -4 -5
        3  6 -7
        4 -8  9
        >>> df.where(m, -df) == np.where(m, df, -df)
              A     B
        0  True  True
        1  True  True
        2  True  True
        3  True  True
        4  True  True
        >>> df.where(m, -df) == df.mask(~m, -df)
              A     B
        0  True  True
        1  True  True
        2  True  True
        3  True  True
        4  True  True

        See Also
        --------
        :func:`DataFrame.%(name_other)s`
        """)

    @Appender(_shared_docs['where'] % dict(_shared_doc_kwargs, cond="True",
                                           cond_rev="False", name='where',
                                           name_other='mask'))
    def where(self, cond, other=np.nan, inplace=False, axis=None, level=None,
              errors='raise', try_cast=False, raise_on_error=None):

        if raise_on_error is not None:
            warnings.warn(
                "raise_on_error is deprecated in "
                "favor of errors='raise|ignore'",
                FutureWarning, stacklevel=2)

            if raise_on_error:
                errors = 'raise'
            else:
                errors = 'ignore'

        other = com._apply_if_callable(other, self)
        return self._where(cond, other, inplace, axis, level,
                           errors=errors, try_cast=try_cast)

    @Appender(_shared_docs['where'] % dict(_shared_doc_kwargs, cond="False",
                                           cond_rev="True", name='mask',
                                           name_other='where'))
    def mask(self, cond, other=np.nan, inplace=False, axis=None, level=None,
             errors='raise', try_cast=False, raise_on_error=None):

        if raise_on_error is not None:
            warnings.warn(
                "raise_on_error is deprecated in "
                "favor of errors='raise|ignore'",
                FutureWarning, stacklevel=2)

            if raise_on_error:
                errors = 'raise'
            else:
                errors = 'ignore'

        inplace = validate_bool_kwarg(inplace, 'inplace')
        cond = com._apply_if_callable(cond, self)

        return self.where(~cond, other=other, inplace=inplace, axis=axis,
                          level=level, try_cast=try_cast,
                          errors=errors)

    _shared_docs['shift'] = ("""
        Shift index by desired number of periods with an optional time freq

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        freq : DateOffset, timedelta, or time rule string, optional
            Increment to use from the tseries module or time rule (e.g. 'EOM').
            See Notes.
        axis : %(axes_single_arg)s

        Notes
        -----
        If freq is specified then the index values are shifted but the data
        is not realigned. That is, use freq if you would like to extend the
        index when shifting and preserve the original data.

        Returns
        -------
        shifted : %(klass)s
    """)

    @Appender(_shared_docs['shift'] % _shared_doc_kwargs)
    def shift(self, periods=1, freq=None, axis=0):
        if periods == 0:
            return self

        block_axis = self._get_block_manager_axis(axis)
        if freq is None:
            new_data = self._data.shift(periods=periods, axis=block_axis)
        else:
            return self.tshift(periods, freq)

        return self._constructor(new_data).__finalize__(self)

    def slice_shift(self, periods=1, axis=0):
        """
        Equivalent to `shift` without copying data. The shifted data will
        not include the dropped periods and the shifted axis will be smaller
        than the original.

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative

        Notes
        -----
        While the `slice_shift` is faster than `shift`, you may pay for it
        later during alignment.

        Returns
        -------
        shifted : same type as caller
        """
        if periods == 0:
            return self

        if periods > 0:
            vslicer = slice(None, -periods)
            islicer = slice(periods, None)
        else:
            vslicer = slice(-periods, None)
            islicer = slice(None, periods)

        new_obj = self._slice(vslicer, axis=axis)
        shifted_axis = self._get_axis(axis)[islicer]
        new_obj.set_axis(shifted_axis, axis=axis, inplace=True)

        return new_obj.__finalize__(self)

    def tshift(self, periods=1, freq=None, axis=0):
        """
        Shift the time index, using the index's frequency if available.

        Parameters
        ----------
        periods : int
            Number of periods to move, can be positive or negative
        freq : DateOffset, timedelta, or time rule string, default None
            Increment to use from the tseries module or time rule (e.g. 'EOM')
        axis : int or basestring
            Corresponds to the axis that contains the Index

        Notes
        -----
        If freq is not specified then tries to use the freq or inferred_freq
        attributes of the index. If neither of those attributes exist, a
        ValueError is thrown

        Returns
        -------
        shifted : NDFrame
        """

        index = self._get_axis(axis)
        if freq is None:
            freq = getattr(index, 'freq', None)

        if freq is None:
            freq = getattr(index, 'inferred_freq', None)

        if freq is None:
            msg = 'Freq was not given and was not set in the index'
            raise ValueError(msg)

        if periods == 0:
            return self

        if isinstance(freq, string_types):
            freq = to_offset(freq)

        block_axis = self._get_block_manager_axis(axis)
        if isinstance(index, PeriodIndex):
            orig_freq = to_offset(index.freq)
            if freq == orig_freq:
                new_data = self._data.copy()
                new_data.axes[block_axis] = index.shift(periods)
            else:
                msg = ('Given freq %s does not match PeriodIndex freq %s' %
                       (freq.rule_code, orig_freq.rule_code))
                raise ValueError(msg)
        else:
            new_data = self._data.copy()
            new_data.axes[block_axis] = index.shift(periods, freq)

        return self._constructor(new_data).__finalize__(self)

    def truncate(self, before=None, after=None, axis=None, copy=True):
        """
        Truncate a Series or DataFrame before and after some index value.

        This is a useful shorthand for boolean indexing based on index
        values above or below certain thresholds.

        Parameters
        ----------
        before : date, string, int
            Truncate all rows before this index value.
        after : date, string, int
            Truncate all rows after this index value.
        axis : {0 or 'index', 1 or 'columns'}, optional
            Axis to truncate. Truncates the index (rows) by default.
        copy : boolean, default is True,
            Return a copy of the truncated section.

        Returns
        -------
        type of caller
            The truncated Series or DataFrame.

        See Also
        --------
        DataFrame.loc : Select a subset of a DataFrame by label.
        DataFrame.iloc : Select a subset of a DataFrame by position.

        Notes
        -----
        If the index being truncated contains only datetime values,
        `before` and `after` may be specified as strings instead of
        Timestamps.

        Examples
        --------
        >>> df = pd.DataFrame({'A': ['a', 'b', 'c', 'd', 'e'],
        ...                    'B': ['f', 'g', 'h', 'i', 'j'],
        ...                    'C': ['k', 'l', 'm', 'n', 'o']},
        ...                    index=[1, 2, 3, 4, 5])
        >>> df
           A  B  C
        1  a  f  k
        2  b  g  l
        3  c  h  m
        4  d  i  n
        5  e  j  o

        >>> df.truncate(before=2, after=4)
           A  B  C
        2  b  g  l
        3  c  h  m
        4  d  i  n

        The columns of a DataFrame can be truncated.

        >>> df.truncate(before="A", after="B", axis="columns")
           A  B
        1  a  f
        2  b  g
        3  c  h
        4  d  i
        5  e  j

        For Series, only rows can be truncated.

        >>> df['A'].truncate(before=2, after=4)
        2    b
        3    c
        4    d
        Name: A, dtype: object

        The index values in ``truncate`` can be datetimes or string
        dates.

        >>> dates = pd.date_range('2016-01-01', '2016-02-01', freq='s')
        >>> df = pd.DataFrame(index=dates, data={'A': 1})
        >>> df.tail()
                             A
        2016-01-31 23:59:56  1
        2016-01-31 23:59:57  1
        2016-01-31 23:59:58  1
        2016-01-31 23:59:59  1
        2016-02-01 00:00:00  1

        >>> df.truncate(before=pd.Timestamp('2016-01-05'),
        ...             after=pd.Timestamp('2016-01-10')).tail()
                             A
        2016-01-09 23:59:56  1
        2016-01-09 23:59:57  1
        2016-01-09 23:59:58  1
        2016-01-09 23:59:59  1
        2016-01-10 00:00:00  1

        Because the index is a DatetimeIndex containing only dates, we can
        specify `before` and `after` as strings. They will be coerced to
        Timestamps before truncation.

        >>> df.truncate('2016-01-05', '2016-01-10').tail()
                             A
        2016-01-09 23:59:56  1
        2016-01-09 23:59:57  1
        2016-01-09 23:59:58  1
        2016-01-09 23:59:59  1
        2016-01-10 00:00:00  1

        Note that ``truncate`` assumes a 0 value for any unspecified time
        component (midnight). This differs from partial string slicing, which
        returns any partially matching dates.

        >>> df.loc['2016-01-05':'2016-01-10', :].tail()
                             A
        2016-01-10 23:59:55  1
        2016-01-10 23:59:56  1
        2016-01-10 23:59:57  1
        2016-01-10 23:59:58  1
        2016-01-10 23:59:59  1
        """

        if axis is None:
            axis = self._stat_axis_number
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        # GH 17935
        # Check that index is sorted
        if not ax.is_monotonic_increasing and not ax.is_monotonic_decreasing:
            raise ValueError("truncate requires a sorted index")

        # if we have a date index, convert to dates, otherwise
        # treat like a slice
        if ax.is_all_dates:
            from pandas.core.tools.datetimes import to_datetime
            before = to_datetime(before)
            after = to_datetime(after)

        if before is not None and after is not None:
            if before > after:
                raise ValueError('Truncate: %s must be after %s' %
                                 (after, before))

        slicer = [slice(None, None)] * self._AXIS_LEN
        slicer[axis] = slice(before, after)
        result = self.loc[tuple(slicer)]

        if isinstance(ax, MultiIndex):
            setattr(result, self._get_axis_name(axis),
                    ax.truncate(before, after))

        if copy:
            result = result.copy()

        return result

    def tz_convert(self, tz, axis=0, level=None, copy=True):
        """
        Convert tz-aware axis to target time zone.

        Parameters
        ----------
        tz : string or pytz.timezone object
        axis : the axis to convert
        level : int, str, default None
            If axis ia a MultiIndex, convert a specific level. Otherwise
            must be None
        copy : boolean, default True
            Also make a copy of the underlying data

        Returns
        -------

        Raises
        ------
        TypeError
            If the axis is tz-naive.
        """
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        def _tz_convert(ax, tz):
            if not hasattr(ax, 'tz_convert'):
                if len(ax) > 0:
                    ax_name = self._get_axis_name(axis)
                    raise TypeError('%s is not a valid DatetimeIndex or '
                                    'PeriodIndex' % ax_name)
                else:
                    ax = DatetimeIndex([], tz=tz)
            else:
                ax = ax.tz_convert(tz)
            return ax

        # if a level is given it must be a MultiIndex level or
        # equivalent to the axis name
        if isinstance(ax, MultiIndex):
            level = ax._get_level_number(level)
            new_level = _tz_convert(ax.levels[level], tz)
            ax = ax.set_levels(new_level, level=level)
        else:
            if level not in (None, 0, ax.name):
                raise ValueError("The level {0} is not valid".format(level))
            ax = _tz_convert(ax, tz)

        result = self._constructor(self._data, copy=copy)
        result.set_axis(ax, axis=axis, inplace=True)
        return result.__finalize__(self)

    def tz_localize(self, tz, axis=0, level=None, copy=True,
                    ambiguous='raise'):
        """
        Localize tz-naive TimeSeries to target time zone.

        Parameters
        ----------
        tz : string or pytz.timezone object
        axis : the axis to localize
        level : int, str, default None
            If axis ia a MultiIndex, localize a specific level. Otherwise
            must be None
        copy : boolean, default True
            Also make a copy of the underlying data
        ambiguous : 'infer', bool-ndarray, 'NaT', default 'raise'
            - 'infer' will attempt to infer fall dst-transition hours based on
              order
            - bool-ndarray where True signifies a DST time, False designates
              a non-DST time (note that this flag is only applicable for
              ambiguous times)
            - 'NaT' will return NaT where there are ambiguous times
            - 'raise' will raise an AmbiguousTimeError if there are ambiguous
              times

        Returns
        -------

        Raises
        ------
        TypeError
            If the TimeSeries is tz-aware and tz is not None.
        """
        axis = self._get_axis_number(axis)
        ax = self._get_axis(axis)

        def _tz_localize(ax, tz, ambiguous):
            if not hasattr(ax, 'tz_localize'):
                if len(ax) > 0:
                    ax_name = self._get_axis_name(axis)
                    raise TypeError('%s is not a valid DatetimeIndex or '
                                    'PeriodIndex' % ax_name)
                else:
                    ax = DatetimeIndex([], tz=tz)
            else:
                ax = ax.tz_localize(tz, ambiguous=ambiguous)
            return ax

        # if a level is given it must be a MultiIndex level or
        # equivalent to the axis name
        if isinstance(ax, MultiIndex):
            level = ax._get_level_number(level)
            new_level = _tz_localize(ax.levels[level], tz, ambiguous)
            ax = ax.set_levels(new_level, level=level)
        else:
            if level not in (None, 0, ax.name):
                raise ValueError("The level {0} is not valid".format(level))
            ax = _tz_localize(ax, tz, ambiguous)

        result = self._constructor(self._data, copy=copy)
        result.set_axis(ax, axis=axis, inplace=True)
        return result.__finalize__(self)

    # ----------------------------------------------------------------------
    # Numeric Methods
    def abs(self):
        """
        Return a Series/DataFrame with absolute numeric value of each element.

        This function only applies to elements that are all numeric.

        Returns
        -------
        abs
            Series/DataFrame containing the absolute value of each element.

        Notes
        -----
        For ``complex`` inputs, ``1.2 + 1j``, the absolute value is
        :math:`\\sqrt{ a^2 + b^2 }`.

        Examples
        --------
        Absolute numeric values in a Series.

        >>> s = pd.Series([-1.10, 2, -3.33, 4])
        >>> s.abs()
        0    1.10
        1    2.00
        2    3.33
        3    4.00
        dtype: float64

        Absolute numeric values in a Series with complex numbers.

        >>> s = pd.Series([1.2 + 1j])
        >>> s.abs()
        0    1.56205
        dtype: float64

        Absolute numeric values in a Series with a Timedelta element.

        >>> s = pd.Series([pd.Timedelta('1 days')])
        >>> s.abs()
        0   1 days
        dtype: timedelta64[ns]

        Select rows with data closest to certain value using argsort (from
        `StackOverflow <https://stackoverflow.com/a/17758115>`__).

        >>> df = pd.DataFrame({
        ...     'a': [4, 5, 6, 7],
        ...     'b': [10, 20, 30, 40],
        ...     'c': [100, 50, -30, -50]
        ... })
        >>> df
             a    b    c
        0    4   10  100
        1    5   20   50
        2    6   30  -30
        3    7   40  -50
        >>> df.loc[(df.c - 43).abs().argsort()]
             a    b    c
        1    5   20   50
        0    4   10  100
        2    6   30  -30
        3    7   40  -50

        See Also
        --------
        numpy.absolute : calculate the absolute value element-wise.
        """
        return np.abs(self)

    def describe(self, percentiles=None, include=None, exclude=None):
        """
        Generates descriptive statistics that summarize the central tendency,
        dispersion and shape of a dataset's distribution, excluding
        ``NaN`` values.

        Analyzes both numeric and object series, as well
        as ``DataFrame`` column sets of mixed data types. The output
        will vary depending on what is provided. Refer to the notes
        below for more detail.

        Parameters
        ----------
        percentiles : list-like of numbers, optional
            The percentiles to include in the output. All should
            fall between 0 and 1. The default is
            ``[.25, .5, .75]``, which returns the 25th, 50th, and
            75th percentiles.
        include : 'all', list-like of dtypes or None (default), optional
            A white list of data types to include in the result. Ignored
            for ``Series``. Here are the options:

            - 'all' : All columns of the input will be included in the output.
            - A list-like of dtypes : Limits the results to the
              provided data types.
              To limit the result to numeric types submit
              ``numpy.number``. To limit it instead to object columns submit
              the ``numpy.object`` data type. Strings
              can also be used in the style of
              ``select_dtypes`` (e.g. ``df.describe(include=['O'])``). To
              select pandas categorical columns, use ``'category'``
            - None (default) : The result will include all numeric columns.
        exclude : list-like of dtypes or None (default), optional,
            A black list of data types to omit from the result. Ignored
            for ``Series``. Here are the options:

            - A list-like of dtypes : Excludes the provided data types
              from the result. To exclude numeric types submit
              ``numpy.number``. To exclude object columns submit the data
              type ``numpy.object``. Strings can also be used in the style of
              ``select_dtypes`` (e.g. ``df.describe(include=['O'])``). To
              exclude pandas categorical columns, use ``'category'``
            - None (default) : The result will exclude nothing.

        Returns
        -------
        summary:  Series/DataFrame of summary statistics

        Notes
        -----
        For numeric data, the result's index will include ``count``,
        ``mean``, ``std``, ``min``, ``max`` as well as lower, ``50`` and
        upper percentiles. By default the lower percentile is ``25`` and the
        upper percentile is ``75``. The ``50`` percentile is the
        same as the median.

        For object data (e.g. strings or timestamps), the result's index
        will include ``count``, ``unique``, ``top``, and ``freq``. The ``top``
        is the most common value. The ``freq`` is the most common value's
        frequency. Timestamps also include the ``first`` and ``last`` items.

        If multiple object values have the highest count, then the
        ``count`` and ``top`` results will be arbitrarily chosen from
        among those with the highest count.

        For mixed data types provided via a ``DataFrame``, the default is to
        return only an analysis of numeric columns. If the dataframe consists
        only of object and categorical data without any numeric columns, the
        default is to return an analysis of both the object and categorical
        columns. If ``include='all'`` is provided as an option, the result
        will include a union of attributes of each type.

        The `include` and `exclude` parameters can be used to limit
        which columns in a ``DataFrame`` are analyzed for the output.
        The parameters are ignored when analyzing a ``Series``.

        Examples
        --------
        Describing a numeric ``Series``.

        >>> s = pd.Series([1, 2, 3])
        >>> s.describe()
        count    3.0
        mean     2.0
        std      1.0
        min      1.0
        25%      1.5
        50%      2.0
        75%      2.5
        max      3.0

        Describing a categorical ``Series``.

        >>> s = pd.Series(['a', 'a', 'b', 'c'])
        >>> s.describe()
        count     4
        unique    3
        top       a
        freq      2
        dtype: object

        Describing a timestamp ``Series``.

        >>> s = pd.Series([
        ...   np.datetime64("2000-01-01"),
        ...   np.datetime64("2010-01-01"),
        ...   np.datetime64("2010-01-01")
        ... ])
        >>> s.describe()
        count                       3
        unique                      2
        top       2010-01-01 00:00:00
        freq                        2
        first     2000-01-01 00:00:00
        last      2010-01-01 00:00:00
        dtype: object

        Describing a ``DataFrame``. By default only numeric fields
        are returned.

        >>> df = pd.DataFrame({ 'object': ['a', 'b', 'c'],
        ...                     'numeric': [1, 2, 3],
        ...                     'categorical': pd.Categorical(['d','e','f'])
        ...                   })
        >>> df.describe()
               numeric
        count      3.0
        mean       2.0
        std        1.0
        min        1.0
        25%        1.5
        50%        2.0
        75%        2.5
        max        3.0

        Describing all columns of a ``DataFrame`` regardless of data type.

        >>> df.describe(include='all')
                categorical  numeric object
        count            3      3.0      3
        unique           3      NaN      3
        top              f      NaN      c
        freq             1      NaN      1
        mean           NaN      2.0    NaN
        std            NaN      1.0    NaN
        min            NaN      1.0    NaN
        25%            NaN      1.5    NaN
        50%            NaN      2.0    NaN
        75%            NaN      2.5    NaN
        max            NaN      3.0    NaN

        Describing a column from a ``DataFrame`` by accessing it as
        an attribute.

        >>> df.numeric.describe()
        count    3.0
        mean     2.0
        std      1.0
        min      1.0
        25%      1.5
        50%      2.0
        75%      2.5
        max      3.0
        Name: numeric, dtype: float64

        Including only numeric columns in a ``DataFrame`` description.

        >>> df.describe(include=[np.number])
               numeric
        count      3.0
        mean       2.0
        std        1.0
        min        1.0
        25%        1.5
        50%        2.0
        75%        2.5
        max        3.0

        Including only string columns in a ``DataFrame`` description.

        >>> df.describe(include=[np.object])
               object
        count       3
        unique      3
        top         c
        freq        1

        Including only categorical columns from a ``DataFrame`` description.

        >>> df.describe(include=['category'])
               categorical
        count            3
        unique           3
        top              f
        freq             1

        Excluding numeric columns from a ``DataFrame`` description.

        >>> df.describe(exclude=[np.number])
               categorical object
        count            3      3
        unique           3      3
        top              f      c
        freq             1      1

        Excluding object columns from a ``DataFrame`` description.

        >>> df.describe(exclude=[np.object])
                categorical  numeric
        count            3      3.0
        unique           3      NaN
        top              f      NaN
        freq             1      NaN
        mean           NaN      2.0
        std            NaN      1.0
        min            NaN      1.0
        25%            NaN      1.5
        50%            NaN      2.0
        75%            NaN      2.5
        max            NaN      3.0

        See Also
        --------
        DataFrame.count
        DataFrame.max
        DataFrame.min
        DataFrame.mean
        DataFrame.std
        DataFrame.select_dtypes
        """
        if self.ndim >= 3:
            msg = "describe is not implemented on Panel objects."
            raise NotImplementedError(msg)
        elif self.ndim == 2 and self.columns.size == 0:
            raise ValueError("Cannot describe a DataFrame without columns")

        if percentiles is not None:
            # explicit conversion of `percentiles` to list
            percentiles = list(percentiles)

            # get them all to be in [0, 1]
            self._check_percentile(percentiles)

            # median should always be included
            if 0.5 not in percentiles:
                percentiles.append(0.5)
            percentiles = np.asarray(percentiles)
        else:
            percentiles = np.array([0.25, 0.5, 0.75])

        # sort and check for duplicates
        unique_pcts = np.unique(percentiles)
        if len(unique_pcts) < len(percentiles):
            raise ValueError("percentiles cannot contain duplicates")
        percentiles = unique_pcts

        formatted_percentiles = format_percentiles(percentiles)

        def describe_numeric_1d(series):
            stat_index = (['count', 'mean', 'std', 'min'] +
                          formatted_percentiles + ['max'])
            d = ([series.count(), series.mean(), series.std(), series.min()] +
                 [series.quantile(x) for x in percentiles] + [series.max()])
            return pd.Series(d, index=stat_index, name=series.name)

        def describe_categorical_1d(data):
            names = ['count', 'unique']
            objcounts = data.value_counts()
            count_unique = len(objcounts[objcounts != 0])
            result = [data.count(), count_unique]
            if result[1] > 0:
                top, freq = objcounts.index[0], objcounts.iloc[0]

                if is_datetime64_dtype(data):
                    asint = data.dropna().values.view('i8')
                    names += ['top', 'freq', 'first', 'last']
                    result += [tslib.Timestamp(top), freq,
                               tslib.Timestamp(asint.min()),
                               tslib.Timestamp(asint.max())]
                else:
                    names += ['top', 'freq']
                    result += [top, freq]

            return pd.Series(result, index=names, name=data.name)

        def describe_1d(data):
            if is_bool_dtype(data):
                return describe_categorical_1d(data)
            elif is_numeric_dtype(data):
                return describe_numeric_1d(data)
            elif is_timedelta64_dtype(data):
                return describe_numeric_1d(data)
            else:
                return describe_categorical_1d(data)

        if self.ndim == 1:
            return describe_1d(self)
        elif (include is None) and (exclude is None):
            # when some numerics are found, keep only numerics
            data = self.select_dtypes(include=[np.number])
            if len(data.columns) == 0:
                data = self
        elif include == 'all':
            if exclude is not None:
                msg = "exclude must be None when include is 'all'"
                raise ValueError(msg)
            data = self
        else:
            data = self.select_dtypes(include=include, exclude=exclude)

        ldesc = [describe_1d(s) for _, s in data.iteritems()]
        # set a convenient order for rows
        names = []
        ldesc_indexes = sorted([x.index for x in ldesc], key=len)
        for idxnames in ldesc_indexes:
            for name in idxnames:
                if name not in names:
                    names.append(name)

        d = pd.concat(ldesc, join_axes=pd.Index([names]), axis=1)
        d.columns = data.columns.copy()
        return d

    def _check_percentile(self, q):
        """Validate percentiles (used by describe and quantile)."""

        msg = ("percentiles should all be in the interval [0, 1]. "
               "Try {0} instead.")
        q = np.asarray(q)
        if q.ndim == 0:
            if not 0 <= q <= 1:
                raise ValueError(msg.format(q / 100.0))
        else:
            if not all(0 <= qs <= 1 for qs in q):
                raise ValueError(msg.format(q / 100.0))
        return q

    _shared_docs['pct_change'] = """
        Percentage change between the current and a prior element.

        Computes the percentage change from the immediately previous row by
        default. This is useful in comparing the percentage of change in a time
        series of elements.

        Parameters
        ----------
        periods : int, default 1
            Periods to shift for forming percent change.
        fill_method : str, default 'pad'
            How to handle NAs before computing percent changes.
        limit : int, default None
            The number of consecutive NAs to fill before stopping.
        freq : DateOffset, timedelta, or offset alias string, optional
            Increment to use from time series API (e.g. 'M' or BDay()).
        **kwargs
            Additional keyword arguments are passed into
            `DataFrame.shift` or `Series.shift`.

        Returns
        -------
        chg : Series or DataFrame
            The same type as the calling object.

        See Also
        --------
        Series.diff : Compute the difference of two elements in a Series.
        DataFrame.diff : Compute the difference of two elements in a DataFrame.
        Series.shift : Shift the index by some number of periods.
        DataFrame.shift : Shift the index by some number of periods.

        Examples
        --------
        **Series**

        >>> s = pd.Series([90, 91, 85])
        >>> s
        0    90
        1    91
        2    85
        dtype: int64

        >>> s.pct_change()
        0         NaN
        1    0.011111
        2   -0.065934
        dtype: float64

        >>> s.pct_change(periods=2)
        0         NaN
        1         NaN
        2   -0.055556
        dtype: float64

        See the percentage change in a Series where filling NAs with last
        valid observation forward to next valid.

        >>> s = pd.Series([90, 91, None, 85])
        >>> s
        0    90.0
        1    91.0
        2     NaN
        3    85.0
        dtype: float64

        >>> s.pct_change(fill_method='ffill')
        0         NaN
        1    0.011111
        2    0.000000
        3   -0.065934
        dtype: float64

        **DataFrame**

        Percentage change in French franc, Deutsche Mark, and Italian lira from
        1980-01-01 to 1980-03-01.

        >>> df = pd.DataFrame({
        ...     'FR': [4.0405, 4.0963, 4.3149],
        ...     'GR': [1.7246, 1.7482, 1.8519],
        ...     'IT': [804.74, 810.01, 860.13]},
        ...     index=['1980-01-01', '1980-02-01', '1980-03-01'])
        >>> df
                        FR      GR      IT
        1980-01-01  4.0405  1.7246  804.74
        1980-02-01  4.0963  1.7482  810.01
        1980-03-01  4.3149  1.8519  860.13

        >>> df.pct_change()
                          FR        GR        IT
        1980-01-01       NaN       NaN       NaN
        1980-02-01  0.013810  0.013684  0.006549
        1980-03-01  0.053365  0.059318  0.061876

        Percentage of change in GOOG and APPL stock volume. Shows computing
        the percentage change between columns.

        >>> df = pd.DataFrame({
        ...     '2016': [1769950, 30586265],
        ...     '2015': [1500923, 40912316],
        ...     '2014': [1371819, 41403351]},
        ...     index=['GOOG', 'APPL'])
        >>> df
                  2016      2015      2014
        GOOG   1769950   1500923   1371819
        APPL  30586265  40912316  41403351

        >>> df.pct_change(axis='columns')
              2016      2015      2014
        GOOG   NaN -0.151997 -0.086016
        APPL   NaN  0.337604  0.012002
        """

    @Appender(_shared_docs['pct_change'] % _shared_doc_kwargs)
    def pct_change(self, periods=1, fill_method='pad', limit=None, freq=None,
                   **kwargs):
        # TODO: Not sure if above is correct - need someone to confirm.
        axis = self._get_axis_number(kwargs.pop('axis', self._stat_axis_name))
        if fill_method is None:
            data = self
        else:
            data = self.fillna(method=fill_method, limit=limit, axis=axis)

        rs = (data.div(data.shift(periods=periods, freq=freq, axis=axis,
                                  **kwargs)) - 1)
        rs = rs.reindex_like(data)
        if freq is None:
            mask = isna(com._values_from_object(data))
            np.putmask(rs.values, mask, np.nan)
        return rs

    def _agg_by_level(self, name, axis=0, level=0, skipna=True, **kwargs):
        grouped = self.groupby(level=level, axis=axis, sort=False)
        if hasattr(grouped, name) and skipna:
            return getattr(grouped, name)(**kwargs)
        axis = self._get_axis_number(axis)
        method = getattr(type(self), name)
        applyf = lambda x: method(x, axis=axis, skipna=skipna, **kwargs)
        return grouped.aggregate(applyf)

    @classmethod
    def _add_numeric_operations(cls):
        """Add the operations to the cls; evaluate the doc strings again"""

        axis_descr, name, name2 = _doc_parms(cls)

        cls.any = _make_logical_function(
            cls, 'any', name, name2, axis_descr,
            _any_desc, nanops.nanany, _any_examples, _any_see_also)
        cls.all = _make_logical_function(
            cls, 'all', name, name2, axis_descr, _all_doc,
            nanops.nanall, _all_examples, _all_see_also)

        @Substitution(outname='mad',
                      desc="Return the mean absolute deviation of the values "
                           "for the requested axis",
                      name1=name, name2=name2, axis_descr=axis_descr,
                      min_count='', examples='')
        @Appender(_num_doc)
        def mad(self, axis=None, skipna=None, level=None):
            if skipna is None:
                skipna = True
            if axis is None:
                axis = self._stat_axis_number
            if level is not None:
                return self._agg_by_level('mad', axis=axis, level=level,
                                          skipna=skipna)

            data = self._get_numeric_data()
            if axis == 0:
                demeaned = data - data.mean(axis=0)
            else:
                demeaned = data.sub(data.mean(axis=1), axis=0)
            return np.abs(demeaned).mean(axis=axis, skipna=skipna)

        cls.mad = mad

        cls.sem = _make_stat_function_ddof(
            cls, 'sem', name, name2, axis_descr,
            "Return unbiased standard error of the mean over requested "
            "axis.\n\nNormalized by N-1 by default. This can be changed "
            "using the ddof argument",
            nanops.nansem)
        cls.var = _make_stat_function_ddof(
            cls, 'var', name, name2, axis_descr,
            "Return unbiased variance over requested axis.\n\nNormalized by "
            "N-1 by default. This can be changed using the ddof argument",
            nanops.nanvar)
        cls.std = _make_stat_function_ddof(
            cls, 'std', name, name2, axis_descr,
            "Return sample standard deviation over requested axis."
            "\n\nNormalized by N-1 by default. This can be changed using the "
            "ddof argument",
            nanops.nanstd)

        @Substitution(outname='compounded',
                      desc="Return the compound percentage of the values for "
                      "the requested axis", name1=name, name2=name2,
                      axis_descr=axis_descr,
                      min_count='', examples='')
        @Appender(_num_doc)
        def compound(self, axis=None, skipna=None, level=None):
            if skipna is None:
                skipna = True
            return (1 + self).prod(axis=axis, skipna=skipna, level=level) - 1

        cls.compound = compound

        cls.cummin = _make_cum_function(
            cls, 'cummin', name, name2, axis_descr, "minimum",
            lambda y, axis: np.minimum.accumulate(y, axis), "min",
            np.inf, np.nan, _cummin_examples)
        cls.cumsum = _make_cum_function(
            cls, 'cumsum', name, name2, axis_descr, "sum",
            lambda y, axis: y.cumsum(axis), "sum", 0.,
            np.nan, _cumsum_examples)
        cls.cumprod = _make_cum_function(
            cls, 'cumprod', name, name2, axis_descr, "product",
            lambda y, axis: y.cumprod(axis), "prod", 1.,
            np.nan, _cumprod_examples)
        cls.cummax = _make_cum_function(
            cls, 'cummax', name, name2, axis_descr, "maximum",
            lambda y, axis: np.maximum.accumulate(y, axis), "max",
            -np.inf, np.nan, _cummax_examples)

        cls.sum = _make_min_count_stat_function(
            cls, 'sum', name, name2, axis_descr,
            'Return the sum of the values for the requested axis',
            nanops.nansum, _sum_examples)
        cls.mean = _make_stat_function(
            cls, 'mean', name, name2, axis_descr,
            'Return the mean of the values for the requested axis',
            nanops.nanmean)
        cls.skew = _make_stat_function(
            cls, 'skew', name, name2, axis_descr,
            'Return unbiased skew over requested axis\nNormalized by N-1',
            nanops.nanskew)
        cls.kurt = _make_stat_function(
            cls, 'kurt', name, name2, axis_descr,
            "Return unbiased kurtosis over requested axis using Fisher's "
            "definition of\nkurtosis (kurtosis of normal == 0.0). Normalized "
            "by N-1\n",
            nanops.nankurt)
        cls.kurtosis = cls.kurt
        cls.prod = _make_min_count_stat_function(
            cls, 'prod', name, name2, axis_descr,
            'Return the product of the values for the requested axis',
            nanops.nanprod, _prod_examples)
        cls.product = cls.prod
        cls.median = _make_stat_function(
            cls, 'median', name, name2, axis_descr,
            'Return the median of the values for the requested axis',
            nanops.nanmedian)
        cls.max = _make_stat_function(
            cls, 'max', name, name2, axis_descr,
            """This method returns the maximum of the values in the object.
            If you want the *index* of the maximum, use ``idxmax``. This is
            the equivalent of the ``numpy.ndarray`` method ``argmax``.""",
            nanops.nanmax)
        cls.min = _make_stat_function(
            cls, 'min', name, name2, axis_descr,
            """This method returns the minimum of the values in the object.
            If you want the *index* of the minimum, use ``idxmin``. This is
            the equivalent of the ``numpy.ndarray`` method ``argmin``.""",
            nanops.nanmin)

    @classmethod
    def _add_series_only_operations(cls):
        """Add the series only operations to the cls; evaluate the doc
        strings again.
        """

        axis_descr, name, name2 = _doc_parms(cls)

        def nanptp(values, axis=0, skipna=True):
            nmax = nanops.nanmax(values, axis, skipna)
            nmin = nanops.nanmin(values, axis, skipna)
            return nmax - nmin

        cls.ptp = _make_stat_function(
            cls, 'ptp', name, name2, axis_descr,
            """Returns the difference between the maximum value and the
            minimum value in the object. This is the equivalent of the
            ``numpy.ndarray`` method ``ptp``.""",
            nanptp)

    @classmethod
    def _add_series_or_dataframe_operations(cls):
        """Add the series or dataframe only operations to the cls; evaluate
        the doc strings again.
        """

        from pandas.core import window as rwindow

        @Appender(rwindow.rolling.__doc__)
        def rolling(self, window, min_periods=None, center=False,
                    win_type=None, on=None, axis=0, closed=None):
            axis = self._get_axis_number(axis)
            return rwindow.rolling(self, window=window,
                                   min_periods=min_periods,
                                   center=center, win_type=win_type,
                                   on=on, axis=axis, closed=closed)

        cls.rolling = rolling

        @Appender(rwindow.expanding.__doc__)
        def expanding(self, min_periods=1, center=False, axis=0):
            axis = self._get_axis_number(axis)
            return rwindow.expanding(self, min_periods=min_periods,
                                     center=center, axis=axis)

        cls.expanding = expanding

        @Appender(rwindow.ewm.__doc__)
        def ewm(self, com=None, span=None, halflife=None, alpha=None,
                min_periods=0, adjust=True, ignore_na=False,
                axis=0):
            axis = self._get_axis_number(axis)
            return rwindow.ewm(self, com=com, span=span, halflife=halflife,
                               alpha=alpha, min_periods=min_periods,
                               adjust=adjust, ignore_na=ignore_na, axis=axis)

        cls.ewm = ewm

        @Appender(_shared_docs['transform'] % _shared_doc_kwargs)
        def transform(self, func, *args, **kwargs):
            result = self.agg(func, *args, **kwargs)
            if is_scalar(result) or len(result) != len(self):
                raise ValueError("transforms cannot produce "
                                 "aggregated results")

            return result

        cls.transform = transform

    # ----------------------------------------------------------------------
    # Misc methods

    _shared_docs['valid_index'] = """
        Return index for %(position)s non-NA/null value.

        Notes
        --------
        If all elements are non-NA/null, returns None.
        Also returns None for empty %(klass)s.

        Returns
        --------
        scalar : type of index
        """

    def _find_valid_index(self, how):
        """Retrieves the index of the first valid value.

        Parameters
        ----------
        how : {'first', 'last'}
            Use this parameter to change between the first or last valid index.

        Returns
        -------
        idx_first_valid : type of index
        """
        assert how in ['first', 'last']

        if len(self) == 0:  # early stop
            return None
        is_valid = ~self.isna()

        if self.ndim == 2:
            is_valid = is_valid.any(1)  # reduce axis 1

        if how == 'first':
            # First valid value case
            i = is_valid.idxmax()
            if not is_valid[i]:
                return None
            return i

        elif how == 'last':
            # Last valid value case
            i = is_valid.values[::-1].argmax()
            if not is_valid.iat[len(self) - i - 1]:
                return None
            return self.index[len(self) - i - 1]

    @Appender(_shared_docs['valid_index'] % {'position': 'first',
                                             'klass': 'NDFrame'})
    def first_valid_index(self):
        return self._find_valid_index('first')

    @Appender(_shared_docs['valid_index'] % {'position': 'last',
                                             'klass': 'NDFrame'})
    def last_valid_index(self):
        return self._find_valid_index('last')


def _doc_parms(cls):
    """Return a tuple of the doc parms."""
    axis_descr = "{%s}" % ', '.join(["{0} ({1})".format(a, i)
                                     for i, a in enumerate(cls._AXIS_ORDERS)])
    name = (cls._constructor_sliced.__name__
            if cls._AXIS_LEN > 1 else 'scalar')
    name2 = cls.__name__
    return axis_descr, name, name2


_num_doc = """

%(desc)s

Parameters
----------
axis : %(axis_descr)s
skipna : boolean, default True
    Exclude NA/null values when computing the result.
level : int or level name, default None
    If the axis is a MultiIndex (hierarchical), count along a
    particular level, collapsing into a %(name1)s
numeric_only : boolean, default None
    Include only float, int, boolean columns. If None, will attempt to use
    everything, then use only numeric data. Not implemented for Series.
%(min_count)s\

Returns
-------
%(outname)s : %(name1)s or %(name2)s (if level specified)

%(examples)s"""

_num_ddof_doc = """

%(desc)s

Parameters
----------
axis : %(axis_descr)s
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA
level : int or level name, default None
    If the axis is a MultiIndex (hierarchical), count along a
    particular level, collapsing into a %(name1)s
ddof : int, default 1
    Delta Degrees of Freedom. The divisor used in calculations is N - ddof,
    where N represents the number of elements.
numeric_only : boolean, default None
    Include only float, int, boolean columns. If None, will attempt to use
    everything, then use only numeric data. Not implemented for Series.

Returns
-------
%(outname)s : %(name1)s or %(name2)s (if level specified)\n"""

_bool_doc = """
%(desc)s

Parameters
----------
axis : int, default 0
    Select the axis which can be 0 for indices and 1 for columns.
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA.
level : int or level name, default None
    If the axis is a MultiIndex (hierarchical), count along a
    particular level, collapsing into a %(name1)s.
bool_only : boolean, default None
    Include only boolean columns. If None, will attempt to use everything,
    then use only boolean data. Not implemented for Series.
**kwargs : any, default None
    Additional keywords have no effect but might be accepted for
    compatibility with NumPy.

Returns
-------
%(outname)s : %(name1)s or %(name2)s (if level specified)

%(see_also)s
%(examples)s"""

_all_doc = """\
Return whether all elements are True over series or dataframe axis.

Returns True if all elements within a series or along a dataframe
axis are non-zero, not-empty or not-False."""

_all_examples = """\
Examples
--------
Series

>>> pd.Series([True, True]).all()
True
>>> pd.Series([True, False]).all()
False

Dataframes

Create a dataframe from a dictionary.

>>> df = pd.DataFrame({'col1': [True, True], 'col2': [True, False]})
>>> df
   col1   col2
0  True   True
1  True  False

Default behaviour checks if column-wise values all return True.

>>> df.all()
col1     True
col2    False
dtype: bool

Adding axis=1 argument will check if row-wise values all return True.

>>> df.all(axis=1)
0     True
1    False
dtype: bool
"""

_all_see_also = """\
See also
--------
pandas.Series.all : Return True if all elements are True
pandas.DataFrame.any : Return True if one (or more) elements are True
"""

_cnum_doc = """
Return cumulative %(desc)s over a DataFrame or Series axis.

Returns a DataFrame or Series of the same size containing the cumulative
%(desc)s.

Parameters
----------
axis : {0 or 'index', 1 or 'columns'}, default 0
    The index or the name of the axis. 0 is equivalent to None or 'index'.
skipna : boolean, default True
    Exclude NA/null values. If an entire row/column is NA, the result
    will be NA.
*args, **kwargs :
    Additional keywords have no effect but might be accepted for
    compatibility with NumPy.

Returns
-------
%(outname)s : %(name1)s or %(name2)s\n
%(examples)s
See also
--------
pandas.core.window.Expanding.%(accum_func_name)s : Similar functionality
    but ignores ``NaN`` values.
%(name2)s.%(accum_func_name)s : Return the %(desc)s over
    %(name2)s axis.
%(name2)s.cummax : Return cumulative maximum over %(name2)s axis.
%(name2)s.cummin : Return cumulative minimum over %(name2)s axis.
%(name2)s.cumsum : Return cumulative sum over %(name2)s axis.
%(name2)s.cumprod : Return cumulative product over %(name2)s axis.
"""

_cummin_examples = """\
Examples
--------
**Series**

>>> s = pd.Series([2, np.nan, 5, -1, 0])
>>> s
0    2.0
1    NaN
2    5.0
3   -1.0
4    0.0
dtype: float64

By default, NA values are ignored.

>>> s.cummin()
0    2.0
1    NaN
2    2.0
3   -1.0
4   -1.0
dtype: float64

To include NA values in the operation, use ``skipna=False``

>>> s.cummin(skipna=False)
0    2.0
1    NaN
2    NaN
3    NaN
4    NaN
dtype: float64

**DataFrame**

>>> df = pd.DataFrame([[2.0, 1.0],
...                    [3.0, np.nan],
...                    [1.0, 0.0]],
...                    columns=list('AB'))
>>> df
     A    B
0  2.0  1.0
1  3.0  NaN
2  1.0  0.0

By default, iterates over rows and finds the minimum
in each column. This is equivalent to ``axis=None`` or ``axis='index'``.

>>> df.cummin()
     A    B
0  2.0  1.0
1  2.0  NaN
2  1.0  0.0

To iterate over columns and find the minimum in each row,
use ``axis=1``

>>> df.cummin(axis=1)
     A    B
0  2.0  1.0
1  3.0  NaN
2  1.0  0.0
"""

_cumsum_examples = """\
Examples
--------
**Series**

>>> s = pd.Series([2, np.nan, 5, -1, 0])
>>> s
0    2.0
1    NaN
2    5.0
3   -1.0
4    0.0
dtype: float64

By default, NA values are ignored.

>>> s.cumsum()
0    2.0
1    NaN
2    7.0
3    6.0
4    6.0
dtype: float64

To include NA values in the operation, use ``skipna=False``

>>> s.cumsum(skipna=False)
0    2.0
1    NaN
2    NaN
3    NaN
4    NaN
dtype: float64

**DataFrame**

>>> df = pd.DataFrame([[2.0, 1.0],
...                    [3.0, np.nan],
...                    [1.0, 0.0]],
...                    columns=list('AB'))
>>> df
     A    B
0  2.0  1.0
1  3.0  NaN
2  1.0  0.0

By default, iterates over rows and finds the sum
in each column. This is equivalent to ``axis=None`` or ``axis='index'``.

>>> df.cumsum()
     A    B
0  2.0  1.0
1  5.0  NaN
2  6.0  1.0

To iterate over columns and find the sum in each row,
use ``axis=1``

>>> df.cumsum(axis=1)
     A    B
0  2.0  3.0
1  3.0  NaN
2  1.0  1.0
"""

_cumprod_examples = """\
Examples
--------
**Series**

>>> s = pd.Series([2, np.nan, 5, -1, 0])
>>> s
0    2.0
1    NaN
2    5.0
3   -1.0
4    0.0
dtype: float64

By default, NA values are ignored.

>>> s.cumprod()
0     2.0
1     NaN
2    10.0
3   -10.0
4    -0.0
dtype: float64

To include NA values in the operation, use ``skipna=False``

>>> s.cumprod(skipna=False)
0    2.0
1    NaN
2    NaN
3    NaN
4    NaN
dtype: float64

**DataFrame**

>>> df = pd.DataFrame([[2.0, 1.0],
...                    [3.0, np.nan],
...                    [1.0, 0.0]],
...                    columns=list('AB'))
>>> df
     A    B
0  2.0  1.0
1  3.0  NaN
2  1.0  0.0

By default, iterates over rows and finds the product
in each column. This is equivalent to ``axis=None`` or ``axis='index'``.

>>> df.cumprod()
     A    B
0  2.0  1.0
1  6.0  NaN
2  6.0  0.0

To iterate over columns and find the product in each row,
use ``axis=1``

>>> df.cumprod(axis=1)
     A    B
0  2.0  2.0
1  3.0  NaN
2  1.0  0.0
"""

_cummax_examples = """\
Examples
--------
**Series**

>>> s = pd.Series([2, np.nan, 5, -1, 0])
>>> s
0    2.0
1    NaN
2    5.0
3   -1.0
4    0.0
dtype: float64

By default, NA values are ignored.

>>> s.cummax()
0    2.0
1    NaN
2    5.0
3    5.0
4    5.0
dtype: float64

To include NA values in the operation, use ``skipna=False``

>>> s.cummax(skipna=False)
0    2.0
1    NaN
2    NaN
3    NaN
4    NaN
dtype: float64

**DataFrame**

>>> df = pd.DataFrame([[2.0, 1.0],
...                    [3.0, np.nan],
...                    [1.0, 0.0]],
...                    columns=list('AB'))
>>> df
     A    B
0  2.0  1.0
1  3.0  NaN
2  1.0  0.0

By default, iterates over rows and finds the maximum
in each column. This is equivalent to ``axis=None`` or ``axis='index'``.

>>> df.cummax()
     A    B
0  2.0  1.0
1  3.0  NaN
2  3.0  1.0

To iterate over columns and find the maximum in each row,
use ``axis=1``

>>> df.cummax(axis=1)
     A    B
0  2.0  2.0
1  3.0  NaN
2  1.0  1.0
"""

_any_see_also = """\
See Also
--------
pandas.DataFrame.all : Return whether all elements are True.
"""

_any_desc = """\
Return whether any element is True over requested axis.

Unlike :meth:`DataFrame.all`, this performs an *or* operation. If any of the
values along the specified axis is True, this will return True."""

_any_examples = """\
Examples
--------
**Series**

For Series input, the output is a scalar indicating whether any element
is True.

>>> pd.Series([True, False]).any()
True

**DataFrame**

Whether each column contains at least one True element (the default).

>>> df = pd.DataFrame({"A": [1, 2], "B": [0, 2], "C": [0, 0]})
>>> df
   A  B  C
0  1  0  0
1  2  2  0

>>> df.any()
A     True
B     True
C    False
dtype: bool

Aggregating over the columns.

>>> df = pd.DataFrame({"A": [True, False], "B": [1, 2]})
>>> df
       A  B
0   True  1
1  False  2

>>> df.any(axis='columns')
0    True
1    True
dtype: bool

>>> df = pd.DataFrame({"A": [True, False], "B": [1, 0]})
>>> df
       A  B
0   True  1
1  False  0

>>> df.any(axis='columns')
0    True
1    False
dtype: bool

`any` for an empty DataFrame is an empty Series.

>>> pd.DataFrame([]).any()
Series([], dtype: bool)
"""

_sum_examples = """\
Examples
--------
By default, the sum of an empty or all-NA Series is ``0``.

>>> pd.Series([]).sum()  # min_count=0 is the default
0.0

This can be controlled with the ``min_count`` parameter. For example, if
you'd like the sum of an empty series to be NaN, pass ``min_count=1``.

>>> pd.Series([]).sum(min_count=1)
nan

Thanks to the ``skipna`` parameter, ``min_count`` handles all-NA and
empty series identically.

>>> pd.Series([np.nan]).sum()
0.0

>>> pd.Series([np.nan]).sum(min_count=1)
nan
"""

_prod_examples = """\
Examples
--------
By default, the product of an empty or all-NA Series is ``1``

>>> pd.Series([]).prod()
1.0

This can be controlled with the ``min_count`` parameter

>>> pd.Series([]).prod(min_count=1)
nan

Thanks to the ``skipna`` parameter, ``min_count`` handles all-NA and
empty series identically.

>>> pd.Series([np.nan]).prod()
1.0

>>> pd.Series([np.nan]).prod(min_count=1)
nan
"""


_min_count_stub = """\
min_count : int, default 0
    The required number of valid values to perform the operation. If fewer than
    ``min_count`` non-NA values are present the result will be NA.

    .. versionadded :: 0.22.0

       Added with the default being 0. This means the sum of an all-NA
       or empty Series is 0, and the product of an all-NA or empty
       Series is 1.
"""


def _make_min_count_stat_function(cls, name, name1, name2, axis_descr, desc,
                                  f, examples):
    @Substitution(outname=name, desc=desc, name1=name1, name2=name2,
                  axis_descr=axis_descr, min_count=_min_count_stub,
                  examples=examples)
    @Appender(_num_doc)
    def stat_func(self, axis=None, skipna=None, level=None, numeric_only=None,
                  min_count=0,
                  **kwargs):
        nv.validate_stat_func(tuple(), kwargs, fname=name)
        if skipna is None:
            skipna = True
        if axis is None:
            axis = self._stat_axis_number
        if level is not None:
            return self._agg_by_level(name, axis=axis, level=level,
                                      skipna=skipna, min_count=min_count)
        return self._reduce(f, name, axis=axis, skipna=skipna,
                            numeric_only=numeric_only, min_count=min_count)

    return set_function_name(stat_func, name, cls)


def _make_stat_function(cls, name, name1, name2, axis_descr, desc, f):
    @Substitution(outname=name, desc=desc, name1=name1, name2=name2,
                  axis_descr=axis_descr, min_count='', examples='')
    @Appender(_num_doc)
    def stat_func(self, axis=None, skipna=None, level=None, numeric_only=None,
                  **kwargs):
        nv.validate_stat_func(tuple(), kwargs, fname=name)
        if skipna is None:
            skipna = True
        if axis is None:
            axis = self._stat_axis_number
        if level is not None:
            return self._agg_by_level(name, axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(f, name, axis=axis, skipna=skipna,
                            numeric_only=numeric_only)

    return set_function_name(stat_func, name, cls)


def _make_stat_function_ddof(cls, name, name1, name2, axis_descr, desc, f):
    @Substitution(outname=name, desc=desc, name1=name1, name2=name2,
                  axis_descr=axis_descr)
    @Appender(_num_ddof_doc)
    def stat_func(self, axis=None, skipna=None, level=None, ddof=1,
                  numeric_only=None, **kwargs):
        nv.validate_stat_ddof_func(tuple(), kwargs, fname=name)
        if skipna is None:
            skipna = True
        if axis is None:
            axis = self._stat_axis_number
        if level is not None:
            return self._agg_by_level(name, axis=axis, level=level,
                                      skipna=skipna, ddof=ddof)
        return self._reduce(f, name, axis=axis, numeric_only=numeric_only,
                            skipna=skipna, ddof=ddof)

    return set_function_name(stat_func, name, cls)


def _make_cum_function(cls, name, name1, name2, axis_descr, desc,
                       accum_func, accum_func_name, mask_a, mask_b, examples):
    @Substitution(outname=name, desc=desc, name1=name1, name2=name2,
                  axis_descr=axis_descr, accum_func_name=accum_func_name,
                  examples=examples)
    @Appender(_cnum_doc)
    def cum_func(self, axis=None, skipna=True, *args, **kwargs):
        skipna = nv.validate_cum_func_with_skipna(skipna, args, kwargs, name)
        if axis is None:
            axis = self._stat_axis_number
        else:
            axis = self._get_axis_number(axis)

        y = com._values_from_object(self).copy()

        if (skipna and
                issubclass(y.dtype.type, (np.datetime64, np.timedelta64))):
            result = accum_func(y, axis)
            mask = isna(self)
            np.putmask(result, mask, tslib.iNaT)
        elif skipna and not issubclass(y.dtype.type, (np.integer, np.bool_)):
            mask = isna(self)
            np.putmask(y, mask, mask_a)
            result = accum_func(y, axis)
            np.putmask(result, mask, mask_b)
        else:
            result = accum_func(y, axis)

        d = self._construct_axes_dict()
        d['copy'] = False
        return self._constructor(result, **d).__finalize__(self)

    return set_function_name(cum_func, name, cls)


def _make_logical_function(cls, name, name1, name2, axis_descr, desc, f,
                           examples, see_also):
    @Substitution(outname=name, desc=desc, name1=name1, name2=name2,
                  axis_descr=axis_descr, examples=examples, see_also=see_also)
    @Appender(_bool_doc)
    def logical_func(self, axis=None, bool_only=None, skipna=None, level=None,
                     **kwargs):
        nv.validate_logical_func(tuple(), kwargs, fname=name)
        if skipna is None:
            skipna = True
        if axis is None:
            axis = self._stat_axis_number
        if level is not None:
            if bool_only is not None:
                raise NotImplementedError("Option bool_only is not "
                                          "implemented with option level.")
            return self._agg_by_level(name, axis=axis, level=level,
                                      skipna=skipna)
        return self._reduce(f, axis=axis, skipna=skipna,
                            numeric_only=bool_only, filter_type='bool',
                            name=name)

    return set_function_name(logical_func, name, cls)


# install the indexes
for _name, _indexer in indexing.get_indexers_list():
    NDFrame._create_indexer(_name, _indexer)
