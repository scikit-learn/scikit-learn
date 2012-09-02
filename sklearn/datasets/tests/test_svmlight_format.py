from bz2 import BZ2File
import gzip
from io import BytesIO
import numpy as np
import os
import shutil
import tempfile

from numpy.testing import assert_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from nose.tools import assert_raises, raises

import sklearn
from sklearn.datasets import (load_svmlight_file, load_svmlight_files,
                              dump_svmlight_file)
from sklearn.utils.testing import assert_in

currdir = os.path.dirname(os.path.abspath(__file__))
datafile = os.path.join(currdir, "data", "svmlight_classification.txt")
multifile = os.path.join(currdir, "data", "svmlight_multilabel.txt")
invalidfile = os.path.join(currdir, "data", "svmlight_invalid.txt")
invalidfile2 = os.path.join(currdir, "data", "svmlight_invalid_order.txt")


def test_load_svmlight_file():
    X, y = load_svmlight_file(datafile)

    # test X's shape
    assert_equal(X.indptr.shape[0], 5)
    assert_equal(X.shape[0], 4)
    assert_equal(X.shape[1], 21)
    assert_equal(y.shape[0], 4)

    # test X's non-zero values
    for i, j, val in ((0, 2, 2.5), (0, 10, -5.2), (0, 15, 1.5),
                     (1, 5, 1.0), (1, 12, -3),
                     (2, 20, 27)):

        assert_equal(X[i, j], val)

    # tests X's zero values
    assert_equal(X[0, 3], 0)
    assert_equal(X[0, 5], 0)
    assert_equal(X[1, 8], 0)
    assert_equal(X[1, 16], 0)
    assert_equal(X[2, 18], 0)

    # test can change X's values
    X[0, 2] *= 2
    assert_equal(X[0, 2], 5)

    # test y
    assert_array_equal(y, [1, 2, 3, 4])


def test_load_svmlight_file_fd():
    # test loading from file descriptor
    X1, y1 = load_svmlight_file(datafile)

    fd = os.open(datafile, os.O_RDONLY)
    try:
        X2, y2 = load_svmlight_file(fd)
        assert_equal(X1.data, X2.data)
        assert_equal(y1, y2)
    finally:
        os.close(fd)


def test_load_svmlight_file_multilabel():
    X, y = load_svmlight_file(multifile, multilabel=True)
    assert_equal(y, [(0, 1), (2,), (1, 2)])


def test_load_svmlight_files():
    X_train, y_train, X_test, y_test = load_svmlight_files([datafile] * 2,
                                                           dtype=np.float32)
    assert_array_equal(X_train.toarray(), X_test.toarray())
    assert_array_equal(y_train, y_test)
    assert_equal(X_train.dtype, np.float32)
    assert_equal(X_test.dtype, np.float32)

    X1, y1, X2, y2, X3, y3 = load_svmlight_files([datafile] * 3,
                                                 dtype=np.float64)
    assert_equal(X1.dtype, X2.dtype)
    assert_equal(X2.dtype, X3.dtype)
    assert_equal(X3.dtype, np.float64)


def test_load_svmlight_file_n_features():
    X, y = load_svmlight_file(datafile, n_features=20)

    # test X'shape
    assert_equal(X.indptr.shape[0], 5)
    assert_equal(X.shape[0], 4)
    assert_equal(X.shape[1], 20)

    # test X's non-zero values
    for i, j, val in ((0, 2, 2.5), (0, 10, -5.2),
                     (1, 5, 1.0), (1, 12, -3)):

        assert_equal(X[i, j], val)


def test_load_compressed():
    X, y = load_svmlight_file(datafile)

    try:
        tempdir = tempfile.mkdtemp(prefix="sklearn-test")

        tmpgz = os.path.join(tempdir, "datafile.gz")
        shutil.copyfileobj(open(datafile, "rb"), gzip.open(tmpgz, "wb"))
        Xgz, ygz = load_svmlight_file(tmpgz)
        assert_array_equal(X.toarray(), Xgz.toarray())
        assert_array_equal(y, ygz)

        tmpbz = os.path.join(tempdir, "datafile.bz2")
        shutil.copyfileobj(open(datafile, "rb"), BZ2File(tmpbz, "wb"))
        Xbz, ybz = load_svmlight_file(tmpgz)
        assert_array_equal(X.toarray(), Xbz.toarray())
        assert_array_equal(y, ybz)
    except:
        shutil.rmtree(tempdir)
        raise


@raises(ValueError)
def test_load_invalid_file():
    load_svmlight_file(invalidfile)


@raises(ValueError)
def test_load_invalid_order_file():
    load_svmlight_file(invalidfile2)


@raises(ValueError)
def test_load_zero_based():
    f = BytesIO("-1 4:1.\n1 0:1\n")
    load_svmlight_file(f, zero_based=False)


def test_load_zero_based_auto():
    data1 = "-1 1:1 2:2 3:3\n"
    data2 = "-1 0:0 1:1\n"

    f1 = BytesIO(data1)
    X, y = load_svmlight_file(f1, zero_based="auto")
    assert_equal(X.shape, (1, 3))

    f1 = BytesIO(data1)
    f2 = BytesIO(data2)
    X1, y1, X2, y2 = load_svmlight_files([f1, f2], zero_based="auto")
    assert_equal(X1.shape, (1, 4))
    assert_equal(X2.shape, (1, 4))


@raises(ValueError)
def test_load_invalid_file2():
    load_svmlight_files([datafile, invalidfile, datafile])


@raises(TypeError)
def test_not_a_filename():
    # in python 3 integers are valid file opening arguments (taken as unix
    # file descriptors)
    load_svmlight_file(.42)


@raises(IOError)
def test_invalid_filename():
    load_svmlight_file("trou pic nic douille")


def test_dump():
    Xs, y = load_svmlight_file(datafile)
    Xd = Xs.toarray()

    for X in (Xs, Xd):
        for zero_based in (True, False):
            for dtype in [np.float32, np.float64]:
                f = BytesIO()
                dump_svmlight_file(X.astype(dtype), y, f,
                                   zero_based=zero_based)
                f.seek(0)

                comment = f.readline()
                assert_in("scikit-learn %s" % sklearn.__version__, comment)
                comment = f.readline()
                assert_in(["one", "zero"][zero_based] + "-based", comment)

                X2, y2 = load_svmlight_file(f, dtype=dtype,
                                            zero_based=zero_based)
                assert_equal(X2.dtype, dtype)
                if dtype == np.float32:
                    assert_array_almost_equal(
                        # allow a rounding error at the last decimal place
                        Xd.astype(dtype), X2.toarray(), 4)
                else:
                    assert_array_almost_equal(
                        # allow a rounding error at the last decimal place
                        Xd.astype(dtype), X2.toarray(), 15)
                assert_array_equal(y, y2)


def test_dump_comment():
    X, y = load_svmlight_file(datafile)
    X = X.toarray()

    f = BytesIO()
    ascii_comment = "This is a comment\nspanning multiple lines."
    dump_svmlight_file(X, y, f, comment=ascii_comment, zero_based=False)
    f.seek(0)

    X2, y2 = load_svmlight_file(f, zero_based=False)
    assert_array_almost_equal(X, X2.toarray())
    assert_array_equal(y, y2)

    # XXX we have to update this to support Python 3.x
    utf8_comment = "It is true that\n\xc2\xbd\xc2\xb2 = \xc2\xbc"
    f = BytesIO()
    assert_raises(UnicodeDecodeError,
                  dump_svmlight_file, X, y, f, comment=utf8_comment)

    unicode_comment = utf8_comment.decode("utf-8")
    f = BytesIO()
    dump_svmlight_file(X, y, f, comment=unicode_comment, zero_based=False)
    f.seek(0)

    X2, y2 = load_svmlight_file(f, zero_based=False)
    assert_array_almost_equal(X, X2.toarray())
    assert_array_equal(y, y2)

    f = BytesIO()
    assert_raises(ValueError,
                  dump_svmlight_file, X, y, f, comment="I've got a \0.")


def test_dump_invalid():
    X, y = load_svmlight_file(datafile)

    f = BytesIO()
    y2d = [y]
    assert_raises(ValueError, dump_svmlight_file, X, y2d, f)

    f = BytesIO()
    assert_raises(ValueError, dump_svmlight_file, X, y[:-1], f)
