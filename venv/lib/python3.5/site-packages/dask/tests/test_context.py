from dask.context import set_options, _globals, globalmethod
import dask.array as da
import dask


def test_with_get():
    var = [0]

    def myget(dsk, keys, **kwargs):
        var[0] = var[0] + 1
        return dask.get(dsk, keys, **kwargs)

    x = da.ones(10, chunks=(5,))

    assert x.sum().compute() == 10
    assert var[0] == 0

    with set_options(get=myget):
        assert x.sum().compute() == 10
    assert var[0] == 1

    # Make sure we've cleaned up
    assert x.sum().compute() == 10
    assert var[0] == 1


def test_set_options_context_manger():
    with set_options(foo='bar'):
        assert _globals['foo'] == 'bar'
    assert _globals['foo'] is None

    try:
        set_options(foo='baz')
        assert _globals['foo'] == 'baz'
    finally:
        del _globals['foo']


def foo():
    return 'foo'


def bar():
    return 'bar'


class Foo(object):
    @globalmethod(key='f')
    def f():
        return 1

    g = globalmethod(foo, key='g', falsey=bar)


def test_globalmethod():
    x = Foo()

    assert x.f() == 1

    with dask.set_options(f=lambda: 2):
        assert x.f() == 2

    with dask.set_options(f=foo):
        assert x.f is foo
        assert x.f() == 'foo'

    assert x.g is foo
    assert x.g() == 'foo'

    with dask.set_options(g=False):
        assert x.g is bar
        assert x.g() == 'bar'
