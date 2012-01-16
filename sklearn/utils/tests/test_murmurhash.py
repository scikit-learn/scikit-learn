from sklearn.utils.murmurhash import murmurhash3
from nose.tools import assert_equal


def test_mmhash3():
    assert_equal(murmurhash3('foo', 0), 3)
    assert_equal(murmurhash3('foo', 42), 3)
    assert_equal(murmurhash3(u'foo', 0), 3)
    assert_equal(murmurhash3(u'foo', 42), 3)
    assert_equal(murmurhash3(3, 0), 3)
