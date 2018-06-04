from __future__ import absolute_import, division, print_function

import pytest

from dask.hashing import hashers, hash_buffer, hash_buffer_hex


np = pytest.importorskip('numpy')

buffers = [
    b'abc',
    bytearray(b'123'),
    memoryview(b'456'),
    np.array(42),
    np.ones((100, 100)),
    np.zeros((100, 100), dtype=[('a', 'i4'), ('b', 'i2')]),
    np.ones(10000, dtype=np.int8)[1:],  # unaligned
]


@pytest.mark.parametrize('x', buffers)
def test_hash_buffer(x):
    for hasher in [None] + hashers:
        h = hash_buffer(x, hasher=hasher)
        assert isinstance(h, bytes)
        assert 8 <= len(h) < 32
        assert h == hash_buffer(x, hasher=hasher)


@pytest.mark.parametrize('x', buffers)
def test_hash_buffer_hex(x):
    for hasher in [None] + hashers:
        h = hash_buffer_hex(x, hasher=hasher)
        assert isinstance(h, str)
        assert 16 <= len(h) < 64
        assert h == hash_buffer_hex(x, hasher=hasher)


@pytest.mark.parametrize('hasher', hashers)
def test_hashers(hasher):
    # Sanity check
    x = b'x'
    h = hasher(x)
    assert isinstance(h, bytes)
    assert 8 <= len(h) < 32
