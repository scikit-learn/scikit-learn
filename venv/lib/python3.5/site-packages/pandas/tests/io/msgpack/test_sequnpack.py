# coding: utf-8

from pandas import compat
from pandas.io.msgpack import Unpacker, BufferFull
from pandas.io.msgpack import OutOfData

import pytest
import pandas.util.testing as tm


class TestPack(object):

    def test_partial_data(self):
        unpacker = Unpacker()
        msg = "No more data to unpack"

        for data in [b"\xa5", b"h", b"a", b"l", b"l"]:
            unpacker.feed(data)
            with tm.assert_raises_regex(StopIteration, msg):
                next(iter(unpacker))

        unpacker.feed(b"o")
        assert next(iter(unpacker)) == b"hallo"

    def test_foobar(self):
        unpacker = Unpacker(read_size=3, use_list=1)
        unpacker.feed(b'foobar')
        assert unpacker.unpack() == ord(b'f')
        assert unpacker.unpack() == ord(b'o')
        assert unpacker.unpack() == ord(b'o')
        assert unpacker.unpack() == ord(b'b')
        assert unpacker.unpack() == ord(b'a')
        assert unpacker.unpack() == ord(b'r')
        pytest.raises(OutOfData, unpacker.unpack)

        unpacker.feed(b'foo')
        unpacker.feed(b'bar')

        k = 0
        for o, e in zip(unpacker, 'foobarbaz'):
            assert o == ord(e)
            k += 1
        assert k == len(b'foobar')

    def test_foobar_skip(self):
        unpacker = Unpacker(read_size=3, use_list=1)
        unpacker.feed(b'foobar')
        assert unpacker.unpack() == ord(b'f')
        unpacker.skip()
        assert unpacker.unpack() == ord(b'o')
        unpacker.skip()
        assert unpacker.unpack() == ord(b'a')
        unpacker.skip()
        pytest.raises(OutOfData, unpacker.unpack)

    def test_maxbuffersize(self):
        pytest.raises(ValueError, Unpacker, read_size=5, max_buffer_size=3)
        unpacker = Unpacker(read_size=3, max_buffer_size=3, use_list=1)
        unpacker.feed(b'fo')
        pytest.raises(BufferFull, unpacker.feed, b'ob')
        unpacker.feed(b'o')
        assert ord('f') == next(unpacker)
        unpacker.feed(b'b')
        assert ord('o') == next(unpacker)
        assert ord('o') == next(unpacker)
        assert ord('b') == next(unpacker)

    def test_readbytes(self):
        unpacker = Unpacker(read_size=3)
        unpacker.feed(b'foobar')
        assert unpacker.unpack() == ord(b'f')
        assert unpacker.read_bytes(3) == b'oob'
        assert unpacker.unpack() == ord(b'a')
        assert unpacker.unpack() == ord(b'r')

        # Test buffer refill
        unpacker = Unpacker(compat.BytesIO(b'foobar'), read_size=3)
        assert unpacker.unpack() == ord(b'f')
        assert unpacker.read_bytes(3) == b'oob'
        assert unpacker.unpack() == ord(b'a')
        assert unpacker.unpack() == ord(b'r')

    def test_issue124(self):
        unpacker = Unpacker()
        unpacker.feed(b'\xa1?\xa1!')
        assert tuple(unpacker) == (b'?', b'!')
        assert tuple(unpacker) == ()
        unpacker.feed(b"\xa1?\xa1")
        assert tuple(unpacker) == (b'?', )
        assert tuple(unpacker) == ()
        unpacker.feed(b"!")
        assert tuple(unpacker) == (b'!', )
        assert tuple(unpacker) == ()
