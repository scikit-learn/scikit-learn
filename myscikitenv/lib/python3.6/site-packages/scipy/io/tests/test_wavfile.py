from __future__ import division, print_function, absolute_import

import os
import sys
import tempfile
from io import BytesIO

import numpy as np
from numpy.testing import assert_equal, assert_, assert_array_equal
from pytest import raises as assert_raises
from scipy._lib._numpy_compat import suppress_warnings

from scipy.io import wavfile


def datafile(fn):
    return os.path.join(os.path.dirname(__file__), 'data', fn)


def test_read_1():
    for mmap in [False, True]:
        rate, data = wavfile.read(datafile('test-44100Hz-le-1ch-4bytes.wav'),
                                  mmap=mmap)

        assert_equal(rate, 44100)
        assert_(np.issubdtype(data.dtype, np.int32))
        assert_equal(data.shape, (4410,))

        del data


def test_read_2():
    for mmap in [False, True]:
        rate, data = wavfile.read(datafile('test-8000Hz-le-2ch-1byteu.wav'),
                                  mmap=mmap)
        assert_equal(rate, 8000)
        assert_(np.issubdtype(data.dtype, np.uint8))
        assert_equal(data.shape, (800, 2))

        del data

def test_read_3():
    for mmap in [False, True]:
        rate, data = wavfile.read(datafile('test-44100Hz-2ch-32bit-float-le.wav'),
                                  mmap=mmap)
        assert_equal(rate, 44100)
        assert_(np.issubdtype(data.dtype, np.float32))
        assert_equal(data.shape, (441, 2))

        del data

def test_read_4():
    for mmap in [False, True]:
        with suppress_warnings() as sup:
            sup.filter(wavfile.WavFileWarning,
                       "Chunk .non-data. not understood, skipping it")
            rate, data = wavfile.read(datafile('test-48000Hz-2ch-64bit-float-le-wavex.wav'),
                                      mmap=mmap)

        assert_equal(rate, 48000)
        assert_(np.issubdtype(data.dtype, np.float64))
        assert_equal(data.shape, (480, 2))

        del data


def test_read_5():
    for mmap in [False, True]:
        rate, data = wavfile.read(datafile('test-44100Hz-2ch-32bit-float-be.wav'),
                                  mmap=mmap)
        assert_equal(rate, 44100)
        assert_(np.issubdtype(data.dtype, np.float32))
        assert_(data.dtype.byteorder == '>' or (sys.byteorder == 'big' and
                                                data.dtype.byteorder == '='))
        assert_equal(data.shape, (441, 2))

        del data


def test_read_fail():
    for mmap in [False, True]:
        fp = open(datafile('example_1.nc'), 'rb')
        assert_raises(ValueError, wavfile.read, fp, mmap=mmap)
        fp.close()


def test_read_early_eof():
    for mmap in [False, True]:
        fp = open(datafile('test-44100Hz-le-1ch-4bytes-early-eof.wav'), 'rb')
        assert_raises(ValueError, wavfile.read, fp, mmap=mmap)
        fp.close()


def test_read_incomplete_chunk():
    for mmap in [False, True]:
        fp = open(datafile('test-44100Hz-le-1ch-4bytes-incomplete-chunk.wav'), 'rb')
        assert_raises(ValueError, wavfile.read, fp, mmap=mmap)
        fp.close()


def _check_roundtrip(realfile, rate, dtype, channels):
    if realfile:
        fd, tmpfile = tempfile.mkstemp(suffix='.wav')
        os.close(fd)
    else:
        tmpfile = BytesIO()
    try:
        data = np.random.rand(100, channels)
        if channels == 1:
            data = data[:,0]
        if dtype.kind == 'f':
            # The range of the float type should be in [-1, 1]
            data = data.astype(dtype)
        else:
            data = (data*128).astype(dtype)

        wavfile.write(tmpfile, rate, data)

        for mmap in [False, True]:
            rate2, data2 = wavfile.read(tmpfile, mmap=mmap)

            assert_equal(rate, rate2)
            assert_(data2.dtype.byteorder in ('<', '=', '|'), msg=data2.dtype)
            assert_array_equal(data, data2)

            del data2
    finally:
        if realfile:
            os.unlink(tmpfile)


def test_write_roundtrip():
    for realfile in (False, True):
        for dtypechar in ('i', 'u', 'f', 'g', 'q'):
            for size in (1, 2, 4, 8):
                if size == 1 and dtypechar == 'i':
                    # signed 8-bit integer PCM is not allowed
                    continue
                if size > 1 and dtypechar == 'u':
                    # unsigned > 8-bit integer PCM is not allowed
                    continue
                if (size == 1 or size == 2) and dtypechar == 'f':
                    # 8- or 16-bit float PCM is not expected
                    continue
                if dtypechar in 'gq':
                    # no size allowed for these types
                    if size == 1:
                        size = ''
                    else:
                        continue

                for endianness in ('>', '<'):
                    if size == 1 and endianness == '<':
                        continue
                    for rate in (8000, 32000):
                        for channels in (1, 2, 5):
                            dt = np.dtype('%s%s%s' % (endianness, dtypechar, size))
                            _check_roundtrip(realfile, rate, dt, channels)

