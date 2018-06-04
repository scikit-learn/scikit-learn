from collections import Mapping

import pytest
from toolz import merge

from dask.sharedict import ShareDict


a = {'x': 1, 'y': 2}
b = {'z': 3}
c = {'w': 2}


def test_core():
    s = ShareDict()
    assert isinstance(s, Mapping)

    s.update(a)
    s.update(b)

    assert s['x'] == 1
    with pytest.raises(KeyError):
        s['abc']

    with pytest.raises((NotImplementedError, TypeError)):
        s['abc'] = 123


def test_structure():
    s = ShareDict()
    s.update(a)
    s.update(b)
    s.update(c)

    assert all(any(d is x for d in s.dicts.values())
               for x in [a, b, c])


@pytest.mark.skip
def test_structure_2():
    s = ShareDict()
    s.update_with_key(a, key='a')
    s.update_with_key(b, key='b')
    s.update_with_key(c, key='c')

    assert s.order == ['a', 'b', 'c']

    s.update_with_key(b, key='b')

    assert s.order == ['a', 'c', 'b']


def test_keys_items():
    s = ShareDict()
    s.update_with_key(a, key='a')
    s.update_with_key(b, key='b')
    s.update_with_key(c, key='c')

    d = merge(a, b, c)

    for fn in [dict, set, len]:
        assert fn(s) == fn(d)

    for fn in [lambda x: x.values(), lambda x: x.keys(), lambda x: x.items()]:
        assert set(fn(s)) == set(fn(d))


def test_update_with_sharedict():
    s = ShareDict()
    s.update_with_key(a, key='a')
    s.update_with_key(b, key='b')
    s.update_with_key(c, key='c')

    d = {'z': 5}

    s2 = ShareDict()
    s2.update_with_key(a, key='a')
    s2.update_with_key(d, key='d')

    s.update(s2)

    assert s.dicts['a'] is s.dicts['a']
