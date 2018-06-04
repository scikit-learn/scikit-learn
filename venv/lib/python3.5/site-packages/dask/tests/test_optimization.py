from operator import getitem
from functools import partial

import pytest

from dask.utils_test import add, inc
from dask.core import get_dependencies
from dask.optimization import (cull, fuse, inline, inline_functions,
                               functions_of, fuse_getitem, fuse_selections,
                               fuse_linear)


def double(x):
    return x * 2


def test_cull():
    # 'out' depends on 'x' and 'y', but not 'z'
    d = {'x': 1, 'y': (inc, 'x'), 'z': (inc, 'x'), 'out': (add, 'y', 10)}
    culled, dependencies = cull(d, 'out')
    assert culled == {'x': 1, 'y': (inc, 'x'), 'out': (add, 'y', 10)}
    assert dependencies == {'x': [], 'y': ['x'], 'out': ['y']}

    assert cull(d, 'out') == cull(d, ['out'])
    assert cull(d, ['out', 'z'])[0] == d
    assert cull(d, [['out'], ['z']]) == cull(d, ['out', 'z'])
    pytest.raises(KeyError, lambda: cull(d, 'badkey'))


def fuse2(*args, **kwargs):
    """Run both ``fuse`` and ``fuse_linear`` and compare results"""
    rv1 = fuse_linear(*args, **kwargs)
    if kwargs.get('rename_keys') is not False:
        return rv1
    rv2 = fuse(*args, **kwargs)
    assert rv1 == rv2
    return rv1


def with_deps(dsk):
    return dsk, {k: get_dependencies(dsk, k) for k in dsk}


def test_fuse():
    fuse = fuse2  # tests both `fuse` and `fuse_linear`
    d = {
        'w': (inc, 'x'),
        'x': (inc, 'y'),
        'y': (inc, 'z'),
        'z': (add, 'a', 'b'),
        'a': 1,
        'b': 2,
    }
    assert fuse(d, rename_keys=False) == with_deps({
        'w': (inc, (inc, (inc, (add, 'a', 'b')))),
        'a': 1,
        'b': 2,
    })
    assert fuse(d, rename_keys=True) == with_deps({
        'z-y-x-w': (inc, (inc, (inc, (add, 'a', 'b')))),
        'a': 1,
        'b': 2,
        'w': 'z-y-x-w',
    })

    d = {
        'NEW': (inc, 'y'),
        'w': (inc, 'x'),
        'x': (inc, 'y'),
        'y': (inc, 'z'),
        'z': (add, 'a', 'b'),
        'a': 1,
        'b': 2,
    }
    assert fuse(d, rename_keys=False) == with_deps({
        'NEW': (inc, 'y'),
        'w': (inc, (inc, 'y')),
        'y': (inc, (add, 'a', 'b')),
        'a': 1,
        'b': 2,
    })
    assert fuse(d, rename_keys=True) == with_deps({
        'NEW': (inc, 'z-y'),
        'x-w': (inc, (inc, 'z-y')),
        'z-y': (inc, (add, 'a', 'b')),
        'a': 1,
        'b': 2,
        'w': 'x-w',
        'y': 'z-y',
    })

    d = {
        'v': (inc, 'y'),
        'u': (inc, 'w'),
        'w': (inc, 'x'),
        'x': (inc, 'y'),
        'y': (inc, 'z'),
        'z': (add, 'a', 'b'),
        'a': (inc, 'c'),
        'b': (inc, 'd'),
        'c': 1,
        'd': 2,
    }
    assert fuse(d, rename_keys=False) == with_deps({
        'u': (inc, (inc, (inc, 'y'))),
        'v': (inc, 'y'),
        'y': (inc, (add, 'a', 'b')),
        'a': (inc, 1),
        'b': (inc, 2),
    })
    assert fuse(d, rename_keys=True) == with_deps({
        'x-w-u': (inc, (inc, (inc, 'z-y'))),
        'v': (inc, 'z-y'),
        'z-y': (inc, (add, 'c-a', 'd-b')),
        'c-a': (inc, 1),
        'd-b': (inc, 2),
        'a': 'c-a',
        'b': 'd-b',
        'u': 'x-w-u',
        'y': 'z-y',
    })

    d = {
        'a': (inc, 'x'),
        'b': (inc, 'x'),
        'c': (inc, 'x'),
        'd': (inc, 'c'),
        'x': (inc, 'y'),
        'y': 0,
    }
    assert fuse(d, rename_keys=False) == with_deps({
        'a': (inc, 'x'),
        'b': (inc, 'x'),
        'd': (inc, (inc, 'x')),
        'x': (inc, 0)
    })
    assert fuse(d, rename_keys=True) == with_deps({
        'a': (inc, 'y-x'),
        'b': (inc, 'y-x'),
        'c-d': (inc, (inc, 'y-x')),
        'y-x': (inc, 0),
        'd': 'c-d',
        'x': 'y-x',
    })

    d = {
        'a': 1,
        'b': (inc, 'a'),
        'c': (add, 'b', 'b'),
    }
    assert fuse(d, rename_keys=False) == with_deps({
        'b': (inc, 1),
        'c': (add, 'b', 'b'),
    })
    assert fuse(d, rename_keys=True) == with_deps({
        'a-b': (inc, 1),
        'c': (add, 'a-b', 'a-b'),
        'b': 'a-b',
    })


def test_fuse_keys():
    fuse = fuse2  # tests both `fuse` and `fuse_linear`
    d = {
        'a': 1,
        'b': (inc, 'a'),
        'c': (inc, 'b'),
    }
    keys = ['b']
    assert fuse(d, keys, rename_keys=False) == with_deps({
        'b': (inc, 1),
        'c': (inc, 'b'),
    })
    assert fuse(d, keys, rename_keys=True) == with_deps({
        'a-b': (inc, 1),
        'c': (inc, 'a-b'),
        'b': 'a-b',
    })

    d = {
        'w': (inc, 'x'),
        'x': (inc, 'y'),
        'y': (inc, 'z'),
        'z': (add, 'a', 'b'),
        'a': 1,
        'b': 2,
    }
    keys = ['x', 'z']
    assert fuse(d, keys, rename_keys=False) == with_deps({
        'w': (inc, 'x'),
        'x': (inc, (inc, 'z')),
        'z': (add, 'a', 'b'),
        'a': 1,
        'b': 2 ,
    })
    assert fuse(d, keys, rename_keys=True) == with_deps({
        'w': (inc, 'y-x'),
        'y-x': (inc, (inc, 'z')),
        'z': (add, 'a', 'b'),
        'a': 1,
        'b': 2 ,
        'x': 'y-x',
    })


def test_inline():
    d = {'a': 1,
         'b': (inc, 'a'),
         'c': (inc, 'b'),
         'd': (add, 'a', 'c')}
    assert inline(d) == {'a': 1,
                         'b': (inc, 1),
                         'c': (inc, 'b'),
                         'd': (add, 1, 'c')}
    assert inline(d, ['a', 'b', 'c']) == {'a': 1,
                                          'b': (inc, 1),
                                          'c': (inc, (inc, 1)),
                                          'd': (add, 1, (inc, (inc, 1)))}
    d = {'x': 1,
         'y': (inc, 'x'),
         'z': (add, 'x', 'y')}
    assert inline(d) == {'x': 1,
                         'y': (inc, 1),
                         'z': (add, 1, 'y')}
    assert inline(d, keys='y') == {'x': 1,
                                   'y': (inc, 1),
                                   'z': (add, 1, (inc, 1))}
    assert inline(d, keys='y',
                  inline_constants=False) == {'x': 1,
                                              'y': (inc, 'x'),
                                              'z': (add, 'x', (inc, 'x'))}

    d = {'a': 1,
         'b': 'a',
         'c': 'b',
         'd': ['a', 'b', 'c'],
         'e': (add, (len, 'd'), 'a')}
    assert inline(d, 'd') == {'a': 1,
                              'b': 1,
                              'c': 1,
                              'd': [1, 1, 1],
                              'e': (add, (len, [1, 1, 1]), 1)}
    assert inline(d, 'a',
                  inline_constants=False) == {'a': 1,
                                              'b': 1,
                                              'c': 'b',
                                              'd': [1, 'b', 'c'],
                                              'e': (add, (len, 'd'), 1)}


def test_inline_functions():
    x, y, i, d = 'xyid'
    dsk = {'out': (add, i, d),
           i: (inc, x),
           d: (double, y),
           x: 1, y: 1}

    result = inline_functions(dsk, [], fast_functions=set([inc]))
    expected = {'out': (add, (inc, x), d),
                d: (double, y),
                x: 1, y: 1}
    assert result == expected


def test_inline_ignores_curries_and_partials():
    dsk = {'x': 1, 'y': 2,
           'a': (partial(add, 1), 'x'),
           'b': (inc, 'a')}

    result = inline_functions(dsk, [], fast_functions=set([add]))
    assert result['b'] == (inc, dsk['a'])
    assert 'a' not in result


def test_inline_doesnt_shrink_fast_functions_at_top():
    dsk = {'x': (inc, 'y'), 'y': 1}
    result = inline_functions(dsk, [], fast_functions=set([inc]))
    assert result == dsk


def test_inline_traverses_lists():
    x, y, i, d = 'xyid'
    dsk = {'out': (sum, [i, d]),
           i: (inc, x),
           d: (double, y),
           x: 1, y: 1}
    expected = {'out': (sum, [(inc, x), d]),
                d: (double, y),
                x: 1, y: 1}
    result = inline_functions(dsk, [], fast_functions=set([inc]))
    assert result == expected


def test_inline_functions_protects_output_keys():
    dsk = {'x': (inc, 1), 'y': (double, 'x')}
    assert inline_functions(dsk, [], [inc]) == {'y': (double, (inc, 1))}
    assert inline_functions(dsk, ['x'], [inc]) == {'y': (double, 'x'),
                                                   'x': (inc, 1)}


def test_functions_of():
    a = lambda x: x
    b = lambda x: x
    assert functions_of((a, 1)) == set([a])
    assert functions_of((a, (b, 1))) == set([a, b])
    assert functions_of((a, [(b, 1)])) == set([a, b])
    assert functions_of((a, [[[(b, 1)]]])) == set([a, b])
    assert functions_of(1) == set()
    assert functions_of(a) == set()
    assert functions_of((a,)) == set([a])


def test_fuse_getitem():
    def load(*args):
        pass
    dsk = {'x': (load, 'store', 'part', ['a', 'b']),
           'y': (getitem, 'x', 'a')}
    dsk2 = fuse_getitem(dsk, load, 3)
    dsk2, dependencies = cull(dsk2, 'y')
    assert dsk2 == {'y': (load, 'store', 'part', 'a')}


def test_fuse_selections():
    def load(*args):
        pass
    dsk = {'x': (load, 'store', 'part', ['a', 'b']),
           'y': (getitem, 'x', 'a')}
    merge = lambda t1, t2: (load, t2[1], t2[2], t1[2])
    dsk2 = fuse_selections(dsk, getitem, load, merge)
    dsk2, dependencies = cull(dsk2, 'y')
    assert dsk2 == {'y': (load, 'store', 'part', 'a')}


def test_inline_cull_dependencies():
    d = {'a': 1,
         'b': 'a',
         'c': 'b',
         'd': ['a', 'b', 'c'],
         'e': (add, (len, 'd'), 'a')}

    d2, dependencies = cull(d, ['d', 'e'])
    inline(d2, {'b'}, dependencies=dependencies)


def test_fuse_reductions_single_input():
    def f(*args):
        return args

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a', 'a'),
        'c': (f, 'b1', 'b2'),
    }
    assert fuse(d, ave_width=1.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1.9, rename_keys=True) == with_deps(d)
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'a': 1,
        'c': (f, (f, 'a'), (f, 'a', 'a')),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-c': (f, (f, 'a'), (f, 'a', 'a')),
        'c': 'b1-b2-c',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a', 'a'),
        'b3': (f, 'a', 'a', 'a'),
        'c': (f, 'b1', 'b2', 'b3'),
    }
    assert fuse(d, ave_width=2.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=2.9, rename_keys=True) == with_deps(d)
    assert fuse(d, ave_width=3, rename_keys=False) == with_deps({
        'a': 1,
        'c': (f, (f, 'a'), (f, 'a', 'a'), (f, 'a', 'a', 'a')),
    })
    assert fuse(d, ave_width=3, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-b3-c': (f, (f, 'a'), (f, 'a', 'a'), (f, 'a', 'a', 'a')),
        'c': 'b1-b2-b3-c',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'c': (f, 'a', 'b1', 'b2'),
    }
    assert fuse(d, ave_width=1.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1.9, rename_keys=True) == with_deps(d)
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'a': 1,
        'c': (f, 'a', (f, 'a'), (f, 'a')),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-c': (f, 'a', (f, 'a'), (f, 'a')),
        'c': 'b1-b2-c',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'c': (f, 'b1', 'b2'),
        'd1': (f, 'c'),
        'd2': (f, 'c'),
        'e': (f, 'd1', 'd2'),
    }
    assert fuse(d, ave_width=1.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1.9, rename_keys=True) == with_deps(d)
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'a': 1,
        'c': (f, (f, 'a'), (f, 'a')),
        'e': (f, (f, 'c'), (f, 'c')),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-c': (f, (f, 'a'), (f, 'a')),
        'd1-d2-e': (f, (f, 'c'), (f, 'c')),
        'c': 'b1-b2-c',
        'e': 'd1-d2-e',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'b3': (f, 'a'),
        'b4': (f, 'a'),
        'c1': (f, 'b1', 'b2'),
        'c2': (f, 'b3', 'b4'),
        'd': (f, 'c1', 'c2'),
    }
    assert fuse(d, ave_width=1.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1.9, rename_keys=True) == with_deps(d)
    expected = with_deps({
        'a': 1,
        'c1': (f, (f, 'a'), (f, 'a')),
        'c2': (f, (f, 'a'), (f, 'a')),
        'd': (f, 'c1', 'c2'),
    })
    assert fuse(d, ave_width=2, rename_keys=False) == expected
    assert fuse(d, ave_width=2.9, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b1-b2-c1': (f, (f, 'a'), (f, 'a')),
        'b3-b4-c2': (f, (f, 'a'), (f, 'a')),
        'd': (f, 'c1', 'c2'),
        'c1': 'b1-b2-c1',
        'c2': 'b3-b4-c2',
    })
    assert fuse(d, ave_width=2, rename_keys=True) == expected
    assert fuse(d, ave_width=2.9, rename_keys=True) == expected
    assert fuse(d, ave_width=3, rename_keys=False) == with_deps({
        'a': 1,
        'd': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
    })
    assert fuse(d, ave_width=3, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-b3-b4-c1-c2-d': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'd': 'b1-b2-b3-b4-c1-c2-d',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'b3': (f, 'a'),
        'b4': (f, 'a'),
        'b5': (f, 'a'),
        'b6': (f, 'a'),
        'b7': (f, 'a'),
        'b8': (f, 'a'),
        'c1': (f, 'b1', 'b2'),
        'c2': (f, 'b3', 'b4'),
        'c3': (f, 'b5', 'b6'),
        'c4': (f, 'b7', 'b8'),
        'd1': (f, 'c1', 'c2'),
        'd2': (f, 'c3', 'c4'),
        'e': (f, 'd1', 'd2'),
    }
    assert fuse(d, ave_width=1.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1.9, rename_keys=True) == with_deps(d)
    expected = with_deps({
        'a': 1,
        'c1': (f, (f, 'a'), (f, 'a')),
        'c2': (f, (f, 'a'), (f, 'a')),
        'c3': (f, (f, 'a'), (f, 'a')),
        'c4': (f, (f, 'a'), (f, 'a')),
        'd1': (f, 'c1', 'c2'),
        'd2': (f, 'c3', 'c4'),
        'e': (f, 'd1', 'd2'),
    })
    assert fuse(d, ave_width=2, rename_keys=False) == expected
    assert fuse(d, ave_width=2.9, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b1-b2-c1': (f, (f, 'a'), (f, 'a')),
        'b3-b4-c2': (f, (f, 'a'), (f, 'a')),
        'b5-b6-c3': (f, (f, 'a'), (f, 'a')),
        'b7-b8-c4': (f, (f, 'a'), (f, 'a')),
        'd1': (f, 'c1', 'c2'),
        'd2': (f, 'c3', 'c4'),
        'e': (f, 'd1', 'd2'),
        'c1': 'b1-b2-c1',
        'c2': 'b3-b4-c2',
        'c3': 'b5-b6-c3',
        'c4': 'b7-b8-c4',
    })
    assert fuse(d, ave_width=2, rename_keys=True) == expected
    assert fuse(d, ave_width=2.9, rename_keys=True) == expected
    expected = with_deps({
        'a': 1,
        'd1': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'd2': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'e': (f, 'd1', 'd2'),
    })
    assert fuse(d, ave_width=3, rename_keys=False) == expected
    assert fuse(d, ave_width=4.6, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b1-b2-b3-b4-c1-c2-d1': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'b5-b6-b7-b8-c3-c4-d2': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'e': (f, 'd1', 'd2'),
        'd1': 'b1-b2-b3-b4-c1-c2-d1',
        'd2': 'b5-b6-b7-b8-c3-c4-d2',

    })
    assert fuse(d, ave_width=3, rename_keys=True) == expected
    assert fuse(d, ave_width=4.6, rename_keys=True) == expected
    assert fuse(d, ave_width=4.7, rename_keys=False) == with_deps({
        'a': 1,
        'e': (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
              (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))))
    })
    assert fuse(d, ave_width=4.7, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-b3-b4-b5-b6-b7-b8-c1-c2-c3-c4-d1-d2-e': (
            f,
            (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
            (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))
        ),
        'e': 'b1-b2-b3-b4-b5-b6-b7-b8-c1-c2-c3-c4-d1-d2-e',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'b3': (f, 'a'),
        'b4': (f, 'a'),
        'b5': (f, 'a'),
        'b6': (f, 'a'),
        'b7': (f, 'a'),
        'b8': (f, 'a'),
        'b9': (f, 'a'),
        'b10': (f, 'a'),
        'b11': (f, 'a'),
        'b12': (f, 'a'),
        'b13': (f, 'a'),
        'b14': (f, 'a'),
        'b15': (f, 'a'),
        'b16': (f, 'a'),
        'c1': (f, 'b1', 'b2'),
        'c2': (f, 'b3', 'b4'),
        'c3': (f, 'b5', 'b6'),
        'c4': (f, 'b7', 'b8'),
        'c5': (f, 'b9', 'b10'),
        'c6': (f, 'b11', 'b12'),
        'c7': (f, 'b13', 'b14'),
        'c8': (f, 'b15', 'b16'),
        'd1': (f, 'c1', 'c2'),
        'd2': (f, 'c3', 'c4'),
        'd3': (f, 'c5', 'c6'),
        'd4': (f, 'c7', 'c8'),
        'e1': (f, 'd1', 'd2'),
        'e2': (f, 'd3', 'd4'),
        'f': (f, 'e1', 'e2'),
    }
    assert fuse(d, ave_width=1.9, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1.9, rename_keys=True) == with_deps(d)
    expected = with_deps({
        'a': 1,
        'c1': (f, (f, 'a'), (f, 'a')),
        'c2': (f, (f, 'a'), (f, 'a')),
        'c3': (f, (f, 'a'), (f, 'a')),
        'c4': (f, (f, 'a'), (f, 'a')),
        'c5': (f, (f, 'a'), (f, 'a')),
        'c6': (f, (f, 'a'), (f, 'a')),
        'c7': (f, (f, 'a'), (f, 'a')),
        'c8': (f, (f, 'a'), (f, 'a')),
        'd1': (f, 'c1', 'c2'),
        'd2': (f, 'c3', 'c4'),
        'd3': (f, 'c5', 'c6'),
        'd4': (f, 'c7', 'c8'),
        'e1': (f, 'd1', 'd2'),
        'e2': (f, 'd3', 'd4'),
        'f': (f, 'e1', 'e2'),
    })
    assert fuse(d, ave_width=2, rename_keys=False) == expected
    assert fuse(d, ave_width=2.9, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b1-b2-c1': (f, (f, 'a'), (f, 'a')),
        'b3-b4-c2': (f, (f, 'a'), (f, 'a')),
        'b5-b6-c3': (f, (f, 'a'), (f, 'a')),
        'b7-b8-c4': (f, (f, 'a'), (f, 'a')),
        'b10-b9-c5': (f, (f, 'a'), (f, 'a')),
        'b11-b12-c6': (f, (f, 'a'), (f, 'a')),
        'b13-b14-c7': (f, (f, 'a'), (f, 'a')),
        'b15-b16-c8': (f, (f, 'a'), (f, 'a')),
        'd1': (f, 'c1', 'c2'),
        'd2': (f, 'c3', 'c4'),
        'd3': (f, 'c5', 'c6'),
        'd4': (f, 'c7', 'c8'),
        'e1': (f, 'd1', 'd2'),
        'e2': (f, 'd3', 'd4'),
        'f': (f, 'e1', 'e2'),
        'c1': 'b1-b2-c1',
        'c2': 'b3-b4-c2',
        'c3': 'b5-b6-c3',
        'c4': 'b7-b8-c4',
        'c5': 'b10-b9-c5',
        'c6': 'b11-b12-c6',
        'c7': 'b13-b14-c7',
        'c8': 'b15-b16-c8',
    })
    assert fuse(d, ave_width=2, rename_keys=True) == expected
    assert fuse(d, ave_width=2.9, rename_keys=True) == expected
    expected = with_deps({
        'a': 1,
        'd1': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'd2': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'd3': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'd4': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'e1': (f, 'd1', 'd2'),
        'e2': (f, 'd3', 'd4'),
        'f': (f, 'e1', 'e2'),
    })
    assert fuse(d, ave_width=3, rename_keys=False) == expected
    assert fuse(d, ave_width=4.6, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b1-b2-b3-b4-c1-c2-d1': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'b5-b6-b7-b8-c3-c4-d2': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'b10-b11-b12-b9-c5-c6-d3': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'b13-b14-b15-b16-c7-c8-d4': (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
        'e1': (f, 'd1', 'd2'),
        'e2': (f, 'd3', 'd4'),
        'f': (f, 'e1', 'e2'),
        'd1': 'b1-b2-b3-b4-c1-c2-d1',
        'd2': 'b5-b6-b7-b8-c3-c4-d2',
        'd3': 'b10-b11-b12-b9-c5-c6-d3',
        'd4': 'b13-b14-b15-b16-c7-c8-d4',
    })
    assert fuse(d, ave_width=3, rename_keys=True) == expected
    assert fuse(d, ave_width=4.6, rename_keys=True) == expected
    expected = with_deps({
        'a': 1,
        'e1': (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
               (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))),
        'e2': (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
               (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))),
        'f': (f, 'e1', 'e2'),
    })
    assert fuse(d, ave_width=4.7, rename_keys=False) == expected
    assert fuse(d, ave_width=7.4, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b1-b2-b3-b4-b5-b6-b7-b8-c1-c2-c3-c4-d1-d2-e1': (
            f,
            (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
            (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))
        ),
        'b10-b11-b12-b13-b14-b15-b16-b9-c5-c6-c7-c8-d3-d4-e2': (
            f,
            (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
            (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))
        ),
        'f': (f, 'e1', 'e2'),
        'e1': 'b1-b2-b3-b4-b5-b6-b7-b8-c1-c2-c3-c4-d1-d2-e1',
        'e2': 'b10-b11-b12-b13-b14-b15-b16-b9-c5-c6-c7-c8-d3-d4-e2',

    })
    assert fuse(d, ave_width=4.7, rename_keys=True) == expected
    assert fuse(d, ave_width=7.4, rename_keys=True) == expected
    assert fuse(d, ave_width=7.5, rename_keys=False) == with_deps({
        'a': 1,
        'f': (f, (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
                  (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))),
              (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
               (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))))),
    })
    assert fuse(d, ave_width=7.5, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b10-b11-b12-b13-b14-b15-b16-b2-b3-b4-b5-b6-b7-b8-b9-c1-c2-c3-c4-c5-c6-c7-c8-d1-d2-d3-d4-e1-e2-f': (
            f,
            (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
             (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a')))),
            (f, (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))),
             (f, (f, (f, 'a'), (f, 'a')), (f, (f, 'a'), (f, 'a'))))
        ),
        'f': 'b1-b10-b11-b12-b13-b14-b15-b16-b2-b3-b4-b5-b6-b7-b8-b9-c1-c2-c3-c4-c5-c6-c7-c8-d1-d2-d3-d4-e1-e2-f',

    })

    d = {
        'a': 1,
        'b': (f, 'a'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps({
        'b': (f, 1)
    })
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps({
        'a-b': (f, 1),
        'b': 'a-b',
    })

    d = {
        'a': 1,
        'b': (f, 'a'),
        'c': (f, 'b'),
        'd': (f, 'c'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps({
        'd': (f, (f, (f, 1)))
    })
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps({
        'a-b-c-d': (f, (f, (f, 1))),
        'd': 'a-b-c-d',
    })

    d = {
        'a': 1,
        'b': (f, 'a'),
        'c': (f, 'a', 'b'),
        'd': (f, 'a', 'c'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps({
        'a': 1,
        'd': (f, 'a', (f, 'a', (f, 'a'))),
    })
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps({
        'a': 1,
        'b-c-d': (f, 'a', (f, 'a', (f, 'a'))),
        'd': 'b-c-d',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'c1': (f, 'b1'),
        'd1': (f, 'c1'),
        'e1': (f, 'd1'),
        'f': (f, 'e1', 'b2'),
    }
    expected = with_deps({
        'a': 1,
        'b2': (f, 'a'),
        'e1': (f, (f, (f, (f, 'a')))),
        'f': (f, 'e1', 'b2'),

    })
    assert fuse(d, ave_width=1, rename_keys=False) == expected
    assert fuse(d, ave_width=1.9, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b2': (f, 'a'),
        'b1-c1-d1-e1': (f, (f, (f, (f, 'a')))),
        'f': (f, 'e1', 'b2'),
        'e1': 'b1-c1-d1-e1',

    })
    assert fuse(d, ave_width=1, rename_keys=True) == expected
    assert fuse(d, ave_width=1.9, rename_keys=True) == expected
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'a': 1,
        'f': (f, (f, (f, (f, (f, 'a')))), (f, 'a')),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-c1-d1-e1-f': (f, (f, (f, (f, (f, 'a')))), (f, 'a')),
        'f': 'b1-b2-c1-d1-e1-f',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'c1': (f, 'a', 'b1'),
        'd1': (f, 'a', 'c1'),
        'e1': (f, 'a', 'd1'),
        'f': (f, 'a', 'e1', 'b2'),
    }
    expected = with_deps({
        'a': 1,
        'b2': (f, 'a'),
        'e1': (f, 'a', (f, 'a', (f, 'a', (f, 'a')))),
        'f': (f, 'a', 'e1', 'b2'),

    })
    assert fuse(d, ave_width=1, rename_keys=False) == expected
    assert fuse(d, ave_width=1.9, rename_keys=False) == expected
    expected = with_deps({
        'a': 1,
        'b2': (f, 'a'),
        'b1-c1-d1-e1': (f, 'a', (f, 'a', (f, 'a', (f, 'a')))),
        'f': (f, 'a', 'e1', 'b2'),
        'e1': 'b1-c1-d1-e1',
    })
    assert fuse(d, ave_width=1, rename_keys=True) == expected
    assert fuse(d, ave_width=1.9, rename_keys=True) == expected
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'a': 1,
        'f': (f, 'a', (f, 'a', (f, 'a', (f, 'a', (f, 'a')))), (f, 'a')),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a': 1,
        'b1-b2-c1-d1-e1-f': (f, 'a', (f, 'a', (f, 'a', (f, 'a', (f, 'a')))), (f, 'a')),
        'f': 'b1-b2-c1-d1-e1-f',
    })

    d = {
        'a': 1,
        'b1': (f, 'a'),
        'b2': (f, 'a'),
        'b3': (f, 'a'),
        'c1': (f, 'b1'),
        'c2': (f, 'b2'),
        'c3': (f, 'b3'),
        'd1': (f, 'c1'),
        'd2': (f, 'c2'),
        'd3': (f, 'c3'),
        'e': (f, 'd1', 'd2', 'd3'),
        'f': (f, 'e'),
        'g': (f, 'f'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps({
        'a': 1,
        'd1': (f, (f, (f, 'a'))),
        'd2': (f, (f, (f, 'a'))),
        'd3': (f, (f, (f, 'a'))),
        'g': (f, (f, (f, 'd1', 'd2', 'd3'))),
    })
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps({
        'a': 1,
        'b1-c1-d1': (f, (f, (f, 'a'))),
        'b2-c2-d2': (f, (f, (f, 'a'))),
        'b3-c3-d3': (f, (f, (f, 'a'))),
        'e-f-g': (f, (f, (f, 'd1', 'd2', 'd3'))),
        'd1': 'b1-c1-d1',
        'd2': 'b2-c2-d2',
        'd3': 'b3-c3-d3',
        'g': 'e-f-g',
    })

    d = {
        'a': 1,
        'b': (f, 'a'),
        'c': (f, 'b'),
        'd': (f, 'b', 'c'),
        'e': (f, 'd'),
        'f': (f, 'e'),
        'g': (f, 'd', 'f'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps({
        'b': (f, 1),
        'd': (f, 'b', (f, 'b')),
        'g': (f, 'd', (f, (f, 'd'))),
    })
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps({
        'a-b': (f, 1),
        'c-d': (f, 'b', (f, 'b')),
        'e-f-g': (f, 'd', (f, (f, 'd'))),
        'b': 'a-b',
        'd': 'c-d',
        'g': 'e-f-g',
    })


def test_fuse_stressed():
    def f(*args):
        return args

    d = {
        'array-original-27b9f9d257a80fa6adae06a98faf71eb': 1,
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 0): (
            f,
            ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 0),
        ),
        ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 0): (
            f,
            ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 1),
        ),
        ('array-27b9f9d257a80fa6adae06a98faf71eb', 0, 0): (
            f,
            'array-original-27b9f9d257a80fa6adae06a98faf71eb',
            (slice(0, 10, None), slice(0, 10, None)),
        ),
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 0): ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 1),
        ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 1): (
            f,
            (f,
             ('array-27b9f9d257a80fa6adae06a98faf71eb', 1, 1),
             (f, [('cholesky-lt-dot-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 0, 1, 0)]))
        ),
        ('cholesky-lt-dot-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 0, 1, 0): (
            f,
            ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 0),
            ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 1),
        ),
        ('array-27b9f9d257a80fa6adae06a98faf71eb', 0, 1): (
            f,
            'array-original-27b9f9d257a80fa6adae06a98faf71eb',
            (slice(0, 10, None), slice(10, 20, None)),
        ),
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 1): (
            f,
            ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 1)
        ),
        ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 1): (
            f,
            (10, 10)
        ),
        ('array-27b9f9d257a80fa6adae06a98faf71eb', 1, 1): (
            f,
            'array-original-27b9f9d257a80fa6adae06a98faf71eb',
            (slice(10, 20, None), slice(10, 20, None)),
        ),
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 1): (
            f,
            ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 0),
            ('array-27b9f9d257a80fa6adae06a98faf71eb', 0, 1),
        ),
        ('cholesky-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 0): (
            f,
            ('array-27b9f9d257a80fa6adae06a98faf71eb', 0, 0),
        ),
    }
    keys = {
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 0),
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 0, 1),
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 0),
        ('cholesky-upper-26a6b670a8aabb7e2f8936db7ccb6a88', 1, 1),
    }
    rv = fuse(d, keys=keys, ave_width=2, rename_keys=True)
    assert rv == with_deps(rv[0])


def test_fuse_reductions_multiple_input():
    def f(*args):
        return args

    d = {
        'a1': 1,
        'a2': 2,
        'b': (f, 'a1', 'a2'),
        'c': (f, 'b'),
    }
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'c': (f, (f, 1, 2)),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a1-a2-b-c': (f, (f, 1, 2)),
        'c': 'a1-a2-b-c',
    })
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps({
        'a1': 1,
        'a2': 2,
        'c': (f, (f, 'a1', 'a2')),
    })
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps({
        'a1': 1,
        'a2': 2,
        'b-c': (f, (f, 'a1', 'a2')),
        'c': 'b-c',
    })

    d = {
        'a1': 1,
        'a2': 2,
        'b1': (f, 'a1'),
        'b2': (f, 'a1', 'a2'),
        'b3': (f, 'a2'),
        'c': (f, 'b1', 'b2', 'b3'),
    }
    expected = with_deps(d)
    assert fuse(d, ave_width=1, rename_keys=False) == expected
    assert fuse(d, ave_width=2.9, rename_keys=False) == expected
    assert fuse(d, ave_width=1, rename_keys=True) == expected
    assert fuse(d, ave_width=2.9, rename_keys=True) == expected
    assert fuse(d, ave_width=3, rename_keys=False) == with_deps({
        'a1': 1,
        'a2': 2,
        'c': (f, (f, 'a1'), (f, 'a1', 'a2'), (f, 'a2')),
    })
    assert fuse(d, ave_width=3, rename_keys=True) == with_deps({
        'a1': 1,
        'a2': 2,
        'b1-b2-b3-c': (f, (f, 'a1'), (f, 'a1', 'a2'), (f, 'a2')),
        'c': 'b1-b2-b3-c',
    })

    d = {
        'a1': 1,
        'a2': 2,
        'b1': (f, 'a1'),
        'b2': (f, 'a1', 'a2'),
        'b3': (f, 'a2'),
        'c1': (f, 'b1', 'b2'),
        'c2': (f, 'b2', 'b3'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps(d)
    assert fuse(d, ave_width=2, rename_keys=False) == with_deps({
        'a1': 1,
        'a2': 2,
        'b2': (f, 'a1', 'a2'),
        'c1': (f, (f, 'a1'), 'b2'),
        'c2': (f, 'b2', (f, 'a2')),
    })
    assert fuse(d, ave_width=2, rename_keys=True) == with_deps({
        'a1': 1,
        'a2': 2,
        'b2': (f, 'a1', 'a2'),
        'b1-c1': (f, (f, 'a1'), 'b2'),
        'b3-c2': (f, 'b2', (f, 'a2')),
        'c1': 'b1-c1',
        'c2': 'b3-c2',
    })

    d = {
        'a1': 1,
        'a2': 2,
        'b1': (f, 'a1'),
        'b2': (f, 'a1', 'a2'),
        'b3': (f, 'a2'),
        'c1': (f, 'b1', 'b2'),
        'c2': (f, 'b2', 'b3'),
        'd': (f, 'c1', 'c2'),
    }
    assert fuse(d, ave_width=1, rename_keys=False) == with_deps(d)
    assert fuse(d, ave_width=1, rename_keys=True) == with_deps(d)

    # A more aggressive heuristic could do this at `ave_width=2`.  Perhaps
    # we can improve this.  Nevertheless, this is behaving as intended.
    assert fuse(d, ave_width=3, rename_keys=False) == with_deps({
        'a1': 1,
        'a2': 2,
        'b2': (f, 'a1', 'a2'),
        'd': (f, (f, (f, 'a1'), 'b2'), (f, 'b2', (f, 'a2'))),
    })
    assert fuse(d, ave_width=3, rename_keys=True) == with_deps({
        'a1': 1,
        'a2': 2,
        'b2': (f, 'a1', 'a2'),
        'b1-b3-c1-c2-d': (f, (f, (f, 'a1'), 'b2'), (f, 'b2', (f, 'a2'))),
        'd': 'b1-b3-c1-c2-d',
    })
