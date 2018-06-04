from collections import namedtuple

import pytest
import pickle

from dask.utils_test import GetFunctionTestMixin, inc, add
from dask import core
from dask.core import (istask, get_dependencies, get_deps, flatten, subs,
                       preorder_traversal, literal, quote, has_tasks)


def contains(a, b):
    """

    >>> contains({'x': 1, 'y': 2}, {'x': 1})
    True
    >>> contains({'x': 1, 'y': 2}, {'z': 3})
    False
    """
    return all(a.get(k) == v for k, v in b.items())


def test_istask():
    assert istask((inc, 1))
    assert not istask(1)
    assert not istask((1, 2))
    f = namedtuple('f', ['x', 'y'])
    assert not istask(f(sum, 2))


def test_has_tasks():
    dsk = {'a': [1, 2, 3],
           'b': 'a',
           'c': [1, (inc, 1)],
           'd': [(sum, 'a')],
           'e': ['a', 'b'],
           'f': [['a', 'b'], 2, 3]}
    assert not has_tasks(dsk, dsk['a'])
    assert has_tasks(dsk, dsk['b'])
    assert has_tasks(dsk, dsk['c'])
    assert has_tasks(dsk, dsk['d'])
    assert has_tasks(dsk, dsk['e'])
    assert has_tasks(dsk, dsk['f'])


def test_preorder_traversal():
    t = (add, 1, 2)
    assert list(preorder_traversal(t)) == [add, 1, 2]
    t = (add, (add, 1, 2), (add, 3, 4))
    assert list(preorder_traversal(t)) == [add, add, 1, 2, add, 3, 4]
    t = (add, (sum, [1, 2]), 3)
    assert list(preorder_traversal(t)) == [add, sum, list, 1, 2, 3]


class TestGet(GetFunctionTestMixin):
    get = staticmethod(core.get)


class TestRecursiveGet(GetFunctionTestMixin):
    get = staticmethod(lambda d, k: core.get(d, k, recursive=True))

    def test_get_stack_limit(self):
        # will blow stack in recursive mode
        pass


def test_GetFunctionTestMixin_class():
    class TestCustomGetFail(GetFunctionTestMixin):
        get = staticmethod(lambda x, y: 1)

    custom_testget = TestCustomGetFail()
    pytest.raises(AssertionError, custom_testget.test_get)

    class TestCustomGetPass(GetFunctionTestMixin):
        get = staticmethod(core.get)

    custom_testget = TestCustomGetPass()
    custom_testget.test_get()


def test_get_dependencies_nested():
    dsk = {'x': 1, 'y': 2,
           'z': (add, (inc, [['x']]), 'y')}

    assert get_dependencies(dsk, 'z') == set(['x', 'y'])
    assert sorted(get_dependencies(dsk, 'z', as_list=True)) == ['x', 'y']


def test_get_dependencies_empty():
    dsk = {'x': (inc,)}
    assert get_dependencies(dsk, 'x') == set()
    assert get_dependencies(dsk, 'x', as_list=True) == []


def test_get_dependencies_list():
    dsk = {'x': 1, 'y': 2, 'z': ['x', [(inc, 'y')]]}
    assert get_dependencies(dsk, 'z') == set(['x', 'y'])
    assert sorted(get_dependencies(dsk, 'z', as_list=True)) == ['x', 'y']


def test_get_dependencies_task():
    dsk = {'x': 1, 'y': 2, 'z': ['x', [(inc, 'y')]]}
    assert get_dependencies(dsk, task=(inc, 'x')) == set(['x'])
    assert get_dependencies(dsk, task=(inc, 'x'), as_list=True) == ['x']


def test_get_dependencies_nothing():
    with pytest.raises(ValueError):
        get_dependencies({})


def test_get_dependencies_many():
    dsk = {'a': [1, 2, 3],
           'b': 'a',
           'c': [1, (inc, 1)],
           'd': [(sum, 'c')],
           'e': ['a', 'b', 'zzz'],
           'f': [['a', 'b'], 2, 3]}

    tasks = [dsk[k] for k in ('d', 'f')]
    s = get_dependencies(dsk, task=tasks)
    assert s == {'a', 'b', 'c'}
    s = get_dependencies(dsk, task=tasks, as_list=True)
    assert sorted(s) == ['a', 'b', 'c']

    s = get_dependencies(dsk, task=[])
    assert s == set()
    s = get_dependencies(dsk, task=[], as_list=True)
    assert s == []


def test_get_deps():
    """
    >>> dsk = {'a': 1, 'b': (inc, 'a'), 'c': (inc, 'b')}
    >>> dependencies, dependents = get_deps(dsk)
    >>> dependencies
    {'a': set([]), 'c': set(['b']), 'b': set(['a'])}
    >>> dependents
    {'a': set(['b']), 'c': set([]), 'b': set(['c'])}
    """
    dsk = {'a': [1, 2, 3],
           'b': 'a',
           'c': [1, (inc, 1)],
           'd': [(sum, 'c')],
           'e': ['b', 'zzz', 'b'],
           'f': [['a', 'b'], 2, 3]}
    dependencies, dependents = get_deps(dsk)
    assert dependencies == {'a': set(),
                            'b': {'a'},
                            'c': set(),
                            'd': {'c'},
                            'e': {'b'},
                            'f': {'a', 'b'},
                            }
    assert dependents == {'a': {'b', 'f'},
                          'b': {'e', 'f'},
                          'c': {'d'},
                          'd': set(),
                          'e': set(),
                          'f': set(),
                          }


def test_flatten():
    assert list(flatten(())) == []
    assert list(flatten('foo')) == ['foo']


def test_subs():
    assert subs((sum, [1, 'x']), 'x', 2) == (sum, [1, 2])
    assert subs((sum, [1, ['x']]), 'x', 2) == (sum, [1, [2]])


class MutateOnEq(object):
    hit_eq = 0

    def __eq__(self, other):
        self.hit_eq += 1
        return False


def test_subs_no_key_data_eq():
    # Numpy throws a deprecation warning on bool(array == scalar), which
    # pollutes the terminal. This test checks that `subs` never tries to
    # compare keys (scalars) with values (which could be arrays)`subs` never
    # tries to compare keys (scalars) with values (which could be arrays).
    a = MutateOnEq()
    subs(a, 'x', 1)
    assert a.hit_eq == 0
    subs((add, a, 'x'), 'x', 1)
    assert a.hit_eq == 0


def test_subs_with_unfriendly_eq():
    try:
        import numpy as np
    except ImportError:
        return
    else:
        task = (np.sum, np.array([1, 2]))
        assert (subs(task, (4, 5), 1) == task) is True

    class MyException(Exception):
        pass

    class F():
        def __eq__(self, other):
            raise MyException()

    task = F()
    assert subs(task, 1, 2) is task


def test_subs_with_surprisingly_friendly_eq():
    try:
        import pandas as pd
    except ImportError:
        return
    else:
        df = pd.DataFrame()
        assert subs(df, 'x', 1) is df


def test_subs_unexpected_hashable_key():
    class UnexpectedButHashable:
        def __init__(self):
            self.name = "a"

        def __hash__(self):
            return hash(self.name)

        def __eq__(self, other):
            return isinstance(other, UnexpectedButHashable)

    assert subs((id, UnexpectedButHashable()), UnexpectedButHashable(), 1) == (id, 1)


def test_quote():
    literals = [[1, 2, 3], (add, 1, 2),
                [1, [2, 3]], (add, 1, (add, 2, 3))]

    for l in literals:
        assert core.get({'x': quote(l)}, 'x') == l


def test_literal_serializable():
    l = literal((add, 1, 2))
    assert pickle.loads(pickle.dumps(l)).data == (add, 1, 2)
