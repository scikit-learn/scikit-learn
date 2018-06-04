from __future__ import absolute_import, division, print_function

import warnings

import pandas as pd
import numpy as np

from ..core import DataFrame, Series
from ..utils import PANDAS_VERSION
from ...base import tokenize
from ...utils import derived_from

if PANDAS_VERSION >= '0.20.0':
    from pandas.core.resample import Resampler as pd_Resampler
else:
    from pandas.tseries.resample import Resampler as pd_Resampler


def getnanos(rule):
    try:
        return getattr(rule, 'nanos', None)
    except ValueError:
        return None


def _resample(obj, rule, how, **kwargs):
    resampler = Resampler(obj, rule, **kwargs)
    if how is not None:
        w = FutureWarning(("how in .resample() is deprecated "
                           "the new syntax is .resample(...)"
                           ".{0}()").format(how))
        warnings.warn(w)
        return getattr(resampler, how)()
    return resampler


def _resample_series(series, start, end, reindex_closed, rule,
                     resample_kwargs, how, fill_value, how_args, how_kwargs):
    out = getattr(series.resample(rule, **resample_kwargs), how)(*how_args, **how_kwargs)
    return out.reindex(pd.date_range(start, end, freq=rule,
                                     closed=reindex_closed),
                       fill_value=fill_value)


def _resample_bin_and_out_divs(divisions, rule, closed='left', label='left'):
    rule = pd.tseries.frequencies.to_offset(rule)
    g = pd.Grouper(freq=rule, how='count', closed=closed, label=label)

    # Determine bins to apply `how` to. Disregard labeling scheme.
    divs = pd.Series(range(len(divisions)), index=divisions)
    temp = divs.resample(rule, closed=closed, label='left').count()
    tempdivs = temp.loc[temp > 0].index

    # Cleanup closed == 'right' and label == 'right'
    res = pd.offsets.Nano() if hasattr(rule, 'delta') else pd.offsets.Day()
    if g.closed == 'right':
        newdivs = tempdivs + res
    else:
        newdivs = tempdivs
    if g.label == 'right':
        outdivs = tempdivs + rule
    else:
        outdivs = tempdivs

    newdivs = newdivs.tolist()
    outdivs = outdivs.tolist()

    # Adjust ends
    if newdivs[0] < divisions[0]:
        newdivs[0] = divisions[0]
    if newdivs[-1] < divisions[-1]:
        if len(newdivs) < len(divs):
            setter = lambda a, val: a.append(val)
        else:
            setter = lambda a, val: a.__setitem__(-1, val)
        setter(newdivs, divisions[-1])
        if outdivs[-1] > divisions[-1]:
            setter(outdivs, outdivs[-1])
        elif outdivs[-1] < divisions[-1]:
            setter(outdivs, temp.index[-1])

    return tuple(map(pd.Timestamp, newdivs)), tuple(map(pd.Timestamp, outdivs))


class Resampler(object):
    def __init__(self, obj, rule, **kwargs):
        if not obj.known_divisions:
            msg = ("Can only resample dataframes with known divisions\n"
                   "See dask.pydata.org/en/latest/dataframe-design.html#partitions\n"
                   "for more information.")
            raise ValueError(msg)
        self.obj = obj
        rule = pd.tseries.frequencies.to_offset(rule)
        day_nanos = pd.tseries.frequencies.Day().nanos

        if getnanos(rule) and day_nanos % rule.nanos:
            raise NotImplementedError('Resampling frequency %s that does'
                                      ' not evenly divide a day is not '
                                      'implemented' % rule)
        self._rule = rule
        self._kwargs = kwargs

    def _agg(self, how, meta=None, fill_value=np.nan, how_args=(), how_kwargs={}):
        rule = self._rule
        kwargs = self._kwargs
        name = 'resample-' + tokenize(self.obj, rule, kwargs, how, *how_args,
                                      **how_kwargs)

        # Create a grouper to determine closed and label conventions
        newdivs, outdivs = _resample_bin_and_out_divs(self.obj.divisions, rule,
                                                      **kwargs)

        # Repartition divs into bins. These won't match labels after mapping
        partitioned = self.obj.repartition(newdivs, force=True)

        keys = partitioned.__dask_keys__()
        dsk = partitioned.dask

        args = zip(keys, outdivs, outdivs[1:], ['left'] * (len(keys) - 1) + [None])
        for i, (k, s, e, c) in enumerate(args):
            dsk[(name, i)] = (_resample_series, k, s, e, c,
                              rule, kwargs, how, fill_value, list(how_args),
                              how_kwargs)

        # Infer output metadata
        meta_r = self.obj._meta_nonempty.resample(self._rule, **self._kwargs)
        meta = getattr(meta_r, how)(*how_args, **how_kwargs)

        if isinstance(meta, pd.DataFrame):
            return DataFrame(dsk, name, meta, outdivs)
        return Series(dsk, name, meta, outdivs)

    @derived_from(pd_Resampler)
    def agg(self, agg_funcs, *args, **kwargs):
        return self._agg('agg', how_args=(agg_funcs,) + args, how_kwargs=kwargs)

    @derived_from(pd_Resampler)
    def count(self):
        return self._agg('count', fill_value=0)

    @derived_from(pd_Resampler)
    def first(self):
        return self._agg('first')

    @derived_from(pd_Resampler)
    def last(self):
        return self._agg('last')

    @derived_from(pd_Resampler)
    def mean(self):
        return self._agg('mean')

    @derived_from(pd_Resampler)
    def min(self):
        return self._agg('min')

    @derived_from(pd_Resampler)
    def median(self):
        return self._agg('median')

    @derived_from(pd_Resampler)
    def max(self):
        return self._agg('max')

    @derived_from(pd_Resampler)
    def ohlc(self):
        return self._agg('ohlc')

    @derived_from(pd_Resampler)
    def prod(self):
        return self._agg('prod')

    @derived_from(pd_Resampler)
    def sem(self):
        return self._agg('sem')

    @derived_from(pd_Resampler)
    def std(self):
        return self._agg('std')

    @derived_from(pd_Resampler)
    def sum(self):
        return self._agg('sum')

    @derived_from(pd_Resampler)
    def var(self):
        return self._agg('var')
