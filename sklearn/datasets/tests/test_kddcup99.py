"""Test  kddcup99 loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job).

Only 'percent10' mode is tested, as the full data
is too big to use in unit-testing.
"""

from functools import partial
from sklearn.datasets.tests.test_common import check_return_X_y

def test_percent10(fetch_kddcup99_fxt):
    data = fetch_kddcup99_fxt()

    assert data.data.shape == (494021, 41)
    assert data.target.shape == (494021,)

    data_shuffled = fetch_kddcup99_fxt(shuffle=True, random_state=0)
    assert data.data.shape == data_shuffled.data.shape
    assert data.target.shape == data_shuffled.target.shape

    data = fetch_kddcup99_fxt('SA')
    assert data.data.shape == (100655, 41)
    assert data.target.shape == (100655,)

    data = fetch_kddcup99_fxt('SF')
    assert data.data.shape == (73237, 4)
    assert data.target.shape == (73237,)

    data = fetch_kddcup99_fxt('http')
    assert data.data.shape == (58725, 3)
    assert data.target.shape == (58725,)

    data = fetch_kddcup99_fxt('smtp')
    assert data.data.shape == (9571, 3)
    assert data.target.shape == (9571,)

    fetch_func = partial(fetch_kddcup99_fxt, 'smtp')
    check_return_X_y(data, fetch_func)


def test_shuffle(fetch_kddcup99_fxt):
    dataset = fetch_kddcup99_fxt(random_state=0, subset='SA', shuffle=True,
                                 percent10=True)
    assert(any(dataset.target[-100:] == b'normal.'))

def test_fetch_kddcup99_check_as_frame_shape(fetch_kddcup99_fxt):

    data = fetch_kddcup99_fxt(as_frame=True)

    assert data.data.shape == (494021, 41)
    assert data.target.shape == (494021,)
    assert data.frame.shape == (494021, 41+1)

    data = fetch_kddcup99_fxt('SA', as_frame=True)
    assert data.data.shape == (100655, 41)
    assert data.target.shape == (100655,)
    assert data.frame.shape == (100655, 41+1)

    data = fetch_kddcup99_fxt('SF', as_frame=True)
    assert data.data.shape == (73237, 4)
    assert data.target.shape == (73237,)
    assert data.frame.shape == (73237, 4+1)

    data = fetch_kddcup99_fxt('http', as_frame=True)
    assert data.data.shape == (58725, 3)
    assert data.target.shape == (58725,)
    assert data.frame.shape == (58725, 3+1)

    data = fetch_kddcup99_fxt('smtp', as_frame=True)
    assert data.data.shape == (9571, 3)
    assert data.target.shape == (9571,)
    assert data.frame.shape == (9571, 3+1)
