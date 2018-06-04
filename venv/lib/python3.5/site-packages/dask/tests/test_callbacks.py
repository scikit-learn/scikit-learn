from dask.local import get_sync
from dask.context import _globals
from dask.threaded import get as get_threaded
from dask.callbacks import Callback
from dask.utils_test import add


def test_start_callback():
    flag = [False]

    class MyCallback(Callback):
        def _start(self, dsk):
            flag[0] = True

    with MyCallback():
        get_sync({'x': 1}, 'x')

    assert flag[0] is True


def test_start_state_callback():
    flag = [False]

    class MyCallback(Callback):
        def _start_state(self, dsk, state):
            flag[0] = True
            assert dsk['x'] == 1
            assert len(state['cache']) == 1

    with MyCallback():
        get_sync({'x': 1}, 'x')

    assert flag[0] is True


def test_finish_always_called():
    flag = [False]

    class MyCallback(Callback):
        def _finish(self, dsk, state, errored):
            flag[0] = True
            assert errored

    dsk = {'x': (lambda: 1 / 0,)}

    # `raise_on_exception=True`
    try:
        with MyCallback():
            get_sync(dsk, 'x')
    except Exception as e:
        assert isinstance(e, ZeroDivisionError)
    assert flag[0]

    # `raise_on_exception=False`
    flag[0] = False
    try:
        with MyCallback():
            get_threaded(dsk, 'x')
    except Exception as e:
        assert isinstance(e, ZeroDivisionError)
    assert flag[0]

    # KeyboardInterrupt
    def raise_keyboard():
        raise KeyboardInterrupt()

    dsk = {'x': (raise_keyboard,)}
    flag[0] = False
    try:
        with MyCallback():
            get_sync(dsk, 'x')
    except BaseException as e:
        assert isinstance(e, KeyboardInterrupt)
    assert flag[0]


def test_nested_schedulers():

    class MyCallback(Callback):
        def _start(self, dsk):
            self.dsk = dsk

        def _pretask(self, key, dsk, state):
            assert key in self.dsk

    inner_callback = MyCallback()
    inner_dsk = {'x': (add, 1, 2),
                 'y': (add, 'x', 3)}

    def nested_call(x):
        assert not _globals['callbacks']
        with inner_callback:
            return get_threaded(inner_dsk, 'y') + x

    outer_callback = MyCallback()
    outer_dsk = {'a': (nested_call, 1),
                 'b': (add, 'a', 2)}

    with outer_callback:
        get_threaded(outer_dsk, 'b')

    assert not _globals['callbacks']
    assert outer_callback.dsk == outer_dsk
    assert inner_callback.dsk == inner_dsk
    assert not _globals['callbacks']


def test_add_remove_mutates_not_replaces():
    g = _globals.copy()

    assert not g['callbacks']

    with Callback():
        pass

    assert not g['callbacks']
