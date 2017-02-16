from sklearn.utils.testing import assert_true
from sklearn.utils.version import parse_version


def test_pre_release():
    assert_true(parse_version('1.12b1') < parse_version('1.12'))
    assert_true(parse_version('1.12b1') < (1, 12, 0))
    assert_true(parse_version('1.12.0b1') < (1, 12))


def test_cmp():
    assert_true(parse_version('1.12') < parse_version('1.13'))
    assert_true(parse_version('1.02') == parse_version('1.2'))
    assert_true(parse_version('1.12') < '1.13')
