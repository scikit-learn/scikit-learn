
from sklearn.externals.joblib import load


def test_old_pickle(tmpdir):
    # Check that a pickle that references sklearn.external.joblib can load
    f = tmpdir.join('foo.pkl')
    f.write(b'\x80\x03csklearn.externals.joblib.numpy_pickle\nNumpyArrayWrappe'
            b'r\nq\x00)\x81q\x01}q\x02(X\n\x00\x00\x00allow_mmapq\x03\x88X\x05'
            b'\x00\x00\x00shapeq\x04K\x01\x85q\x05X\x05\x00\x00\x00dtypeq\x06c'
            b'numpy\ndtype\nq\x07X\x02\x00\x00\x00i8q\x08K\x00K\x01\x87q\tRq\n'
            b'(K\x03X\x01\x00\x00\x00<q\x0bNNNJ\xff\xff\xff\xffJ\xff\xff\xff'
            b'\xffK\x00tq\x0cbX\x08\x00\x00\x00subclassq\rcnumpy\nndarray\nq'
            b'\x0eX\x05\x00\x00\x00orderq\x0fX\x01\x00\x00\x00Cq\x10ub\x01\x00'
            b'\x00\x00\x00\x00\x00\x00.', mode='wb')

    load(str(f))
