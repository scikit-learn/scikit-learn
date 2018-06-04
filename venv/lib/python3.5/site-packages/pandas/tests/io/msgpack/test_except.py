# coding: utf-8

from datetime import datetime
from pandas.io.msgpack import packb, unpackb

import pytest
import pandas.util.testing as tm


class DummyException(Exception):
    pass


class TestExceptions(object):

    def test_raise_on_find_unsupported_value(self):
        msg = "can\'t serialize datetime"
        with tm.assert_raises_regex(TypeError, msg):
            packb(datetime.now())

    def test_raise_from_object_hook(self):
        def hook(_):
            raise DummyException()

        pytest.raises(DummyException, unpackb, packb({}), object_hook=hook)
        pytest.raises(DummyException, unpackb, packb({'fizz': 'buzz'}),
                      object_hook=hook)
        pytest.raises(DummyException, unpackb, packb({'fizz': 'buzz'}),
                      object_pairs_hook=hook)
        pytest.raises(DummyException, unpackb,
                      packb({'fizz': {'buzz': 'spam'}}), object_hook=hook)
        pytest.raises(DummyException, unpackb,
                      packb({'fizz': {'buzz': 'spam'}}),
                      object_pairs_hook=hook)

    def test_invalid_value(self):
        msg = "Unpack failed: error"
        with tm.assert_raises_regex(ValueError, msg):
            unpackb(b"\xd9\x97#DL_")
