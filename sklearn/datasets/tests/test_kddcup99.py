"""Test  kddcup99 loader. Only 'percent10' mode is tested, as the full data
is too big to use in unit-testing.

The test is skipped if the data wasn't previously fetched and saved to
scikit-learn data folder.
"""

from sklearn.datasets import fetch_kddcup99
from sklearn.utils.testing import assert_equal, SkipTest


def test_percent10():
    try:
        data = fetch_kddcup99(download_if_missing=False)
    except IOError:
        raise SkipTest("kddcup99 dataset can not be loaded.")

    assert_equal(data.data.shape, (494021, 41))
    assert_equal(data.target.shape, (494021,))

    data_shuffled = fetch_kddcup99(shuffle=True, random_state=0)
    assert_equal(data.data.shape, data_shuffled.data.shape)
    assert_equal(data.target.shape, data_shuffled.target.shape)

    data = fetch_kddcup99('SA')
    assert_equal(data.data.shape, (100655, 41))
    assert_equal(data.target.shape, (100655,))

    data = fetch_kddcup99('SF')
    assert_equal(data.data.shape, (73237, 4))
    assert_equal(data.target.shape, (73237,))

    data = fetch_kddcup99('http')
    assert_equal(data.data.shape, (58725, 3))
    assert_equal(data.target.shape, (58725,))

    data = fetch_kddcup99('smtp')
    assert_equal(data.data.shape, (9571, 3))
    assert_equal(data.target.shape, (9571,))
