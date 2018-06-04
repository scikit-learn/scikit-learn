import pytest

from pandas.tests.series.common import TestData


@pytest.fixture(scope='module')
def test_data():
    return TestData()
