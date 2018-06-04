from dask.cache import Cache
from dask.local import get_sync
from dask.threaded import get
from operator import add
from dask.context import _globals
from time import sleep
import pytest

cachey = pytest.importorskip('cachey')


flag = []


def inc(x):
    flag.append(x)
    return x + 1


def test_cache():
    c = cachey.Cache(10000)
    cc = Cache(c)

    with cc:
        assert get({'x': (inc, 1)}, 'x') == 2

    assert flag == [1]
    assert c.data['x'] == 2

    assert not cc.starttimes
    assert not cc.durations

    while flag:
        flag.pop()
    dsk = {'x': (inc, 1), 'y': (inc, 2), 'z': (add, 'x', 'y')}
    with cc:
        assert get(dsk, 'z') == 5

    assert flag == [2]  # no x present

    assert not _globals['callbacks']


def test_cache_with_number():
    c = Cache(10000, limit=1)
    assert isinstance(c.cache, cachey.Cache)
    assert c.cache.available_bytes == 10000
    assert c.cache.limit == 1


def f(duration, size, *args):
    sleep(duration)
    return [0] * size


def test_prefer_cheap_dependent():
    dsk = {'x': (f, 0.01, 10), 'y': (f, 0.000001, 1, 'x')}
    c = Cache(10000)
    with c:
        get_sync(dsk, 'y')

    assert c.cache.scorer.cost['x'] < c.cache.scorer.cost['y']
