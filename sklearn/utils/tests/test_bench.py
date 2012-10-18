
import datetime

from sklearn.utils.bench import total_seconds
from nose.tools import assert_equal


def test_total_seconds():
    delta = datetime.datetime(2012, 1, 1, 5, 5, 1) - datetime.datetime(2012, 1,
                                                                    1, 5, 5, 4)
    assert_equal(86397, total_seconds(delta))
