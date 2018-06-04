"""
Alternate namespace for toolz such that all functions are curried

Currying provides implicit partial evaluation of all functions

Example:

    Get usually requires two arguments, an index and a collection
    >>> from toolz.curried import get
    >>> get(0, ('a', 'b'))
    'a'

    When we use it in higher order functions we often want to pass a partially
    evaluated form
    >>> data = [(1, 2), (11, 22), (111, 222)]
    >>> list(map(lambda seq: get(0, seq), data))
    [1, 11, 111]

    The curried version allows simple expression of partial evaluation
    >>> list(map(get(0), data))
    [1, 11, 111]

See Also:
    toolz.functoolz.curry
"""
import toolz
from . import operator
from toolz import (
    comp,
    complement,
    compose,
    concat,
    concatv,
    count,
    curry,
    diff,
    dissoc,
    first,
    flip,
    frequencies,
    identity,
    interleave,
    isdistinct,
    isiterable,
    juxt,
    last,
    memoize,
    merge_sorted,
    peek,
    pipe,
    second,
    thread_first,
    thread_last,
)
from .exceptions import merge, merge_with

accumulate = toolz.curry(toolz.accumulate)
assoc = toolz.curry(toolz.assoc)
assoc_in = toolz.curry(toolz.assoc_in)
cons = toolz.curry(toolz.cons)
countby = toolz.curry(toolz.countby)
do = toolz.curry(toolz.do)
drop = toolz.curry(toolz.drop)
excepts = toolz.curry(toolz.excepts)
filter = toolz.curry(toolz.filter)
get = toolz.curry(toolz.get)
get_in = toolz.curry(toolz.get_in)
groupby = toolz.curry(toolz.groupby)
interpose = toolz.curry(toolz.interpose)
itemfilter = toolz.curry(toolz.itemfilter)
itemmap = toolz.curry(toolz.itemmap)
iterate = toolz.curry(toolz.iterate)
join = toolz.curry(toolz.join)
keyfilter = toolz.curry(toolz.keyfilter)
keymap = toolz.curry(toolz.keymap)
map = toolz.curry(toolz.map)
mapcat = toolz.curry(toolz.mapcat)
nth = toolz.curry(toolz.nth)
partial = toolz.curry(toolz.partial)
partition = toolz.curry(toolz.partition)
partition_all = toolz.curry(toolz.partition_all)
partitionby = toolz.curry(toolz.partitionby)
pluck = toolz.curry(toolz.pluck)
random_sample = toolz.curry(toolz.random_sample)
reduce = toolz.curry(toolz.reduce)
reduceby = toolz.curry(toolz.reduceby)
remove = toolz.curry(toolz.remove)
sliding_window = toolz.curry(toolz.sliding_window)
sorted = toolz.curry(toolz.sorted)
tail = toolz.curry(toolz.tail)
take = toolz.curry(toolz.take)
take_nth = toolz.curry(toolz.take_nth)
topk = toolz.curry(toolz.topk)
unique = toolz.curry(toolz.unique)
update_in = toolz.curry(toolz.update_in)
valfilter = toolz.curry(toolz.valfilter)
valmap = toolz.curry(toolz.valmap)

del exceptions
del toolz
