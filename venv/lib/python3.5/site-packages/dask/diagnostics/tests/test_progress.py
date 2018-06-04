from operator import add, mul

import pytest

from dask.local import get_sync
from dask.diagnostics import ProgressBar
from dask.diagnostics.progress import format_time
from dask.threaded import get as get_threaded
from dask.context import _globals


dsk = {'a': 1,
       'b': 2,
       'c': (add, 'a', 'b'),
       'd': (mul, 'a', 'b'),
       'e': (mul, 'c', 'd')}


def check_bar_completed(capsys, width=40):
    out, err = capsys.readouterr()
    bar, percent, time = [i.strip() for i in out.split('\r')[-1].split('|')]
    assert bar == '[' + '#' * width + ']'
    assert percent == '100% Completed'


def test_progressbar(capsys):
    with ProgressBar():
        out = get_threaded(dsk, 'e')
    assert out == 6
    check_bar_completed(capsys)
    with ProgressBar(width=20):
        out = get_threaded(dsk, 'e')
    check_bar_completed(capsys, 20)


def test_minimum_time(capsys):
    with ProgressBar(1.0):
        out = get_threaded(dsk, 'e')
    out, err = capsys.readouterr()
    assert out == '' and err == ''


@pytest.mark.parametrize('get', [get_threaded, get_sync])
def test_clean_exit(get):
    dsk = {'a': (lambda: 1 / 0, )}
    try:
        with ProgressBar() as pbar:
            get_threaded(dsk, 'a')
    except ZeroDivisionError:
        pass
    assert not pbar._running
    assert not pbar._timer.is_alive()


def test_format_time():
    assert format_time(1.4) == ' 1.4s'
    assert format_time(10.4) == '10.4s'
    assert format_time(100.4) == ' 1min 40.4s'
    assert format_time(1000.4) == '16min 40.4s'
    assert format_time(10000.4) == ' 2hr 46min 40.4s'


def test_register(capsys):
    try:
        p = ProgressBar()
        p.register()

        assert _globals['callbacks']

        get_threaded(dsk, 'e')
        check_bar_completed(capsys)

        p.unregister()

        assert not _globals['callbacks']
    finally:
        _globals['callbacks'].clear()


def test_no_tasks(capsys):
    with ProgressBar():
        get_threaded({'x': 1}, 'x')
    check_bar_completed(capsys)


def test_with_cache(capsys):
    cachey = pytest.importorskip('cachey')
    from dask.cache import Cache
    c = cachey.Cache(10000)
    cc = Cache(c)

    with cc:
        with ProgressBar():
            assert get_threaded({'x': (mul, 1, 2)}, 'x') == 2
    check_bar_completed(capsys)
    assert c.data['x'] == 2

    with cc:
        with ProgressBar():
            assert get_threaded({'x': (mul, 1, 2), 'y': (mul, 'x', 3)}, 'y') == 6
    check_bar_completed(capsys)


def test_with_alias(capsys):
    dsk = {'a': 1,
           'b': 2,
           'c': (add, 'a', 'b'),
           'd': (add, 1, 2),
           'e': 'd',
           'f': (mul, 'e', 'c')}
    with ProgressBar():
        get_threaded(dsk, 'f')
    check_bar_completed(capsys)


def test_store_time():
    p = ProgressBar()
    with p:
        get_threaded({'x': 1}, 'x')

    assert isinstance(p.last_duration, float)
