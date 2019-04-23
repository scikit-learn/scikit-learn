import os
import pytest
from sklearn.utils._joblib import Parallel, delayed, Memory, parallel_backend


def test_old_pickle(tmpdir):
    import joblib

    # Check that a pickle that references sklearn.external.joblib can load
    f = tmpdir.join('foo.pkl')
    f.write(b'\x80\x02csklearn.externals.joblib.numpy_pickle\nNumpyArrayWrappe'
            b'r\nq\x00)\x81q\x01}q\x02(U\x05dtypeq\x03cnumpy\ndtype\nq\x04U'
            b'\x02i8q\x05K\x00K\x01\x87q\x06Rq\x07(K\x03U\x01<q\x08NNNJ\xff'
            b'\xff\xff\xffJ\xff\xff\xff\xffK\x00tq\tbU\x05shapeq\nK\x01\x85q'
            b'\x0bU\x05orderq\x0cU\x01Cq\rU\x08subclassq\x0ecnumpy\nndarray\nq'
            b'\x0fU\nallow_mmapq\x10\x88ub\x01\x00\x00\x00\x00\x00\x00\x00.',
            mode='wb')

    joblib.load(str(f))
