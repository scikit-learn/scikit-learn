import os
import pytest
from sklearn import externals
from sklearn.externals import joblib as joblib_vendored
from sklearn.utils._joblib import Parallel, delayed, Memory, parallel_backend

if os.environ.get('SKLEARN_SITE_JOBLIB', False):
    import joblib as joblib_site
else:
    joblib_site = None


def test_old_pickle(tmpdir):
    vendored_joblib_home = os.path.dirname(joblib_vendored.__file__)
    sklearn_externals_home = os.path.dirname(externals.__file__)
    if not vendored_joblib_home.startswith(sklearn_externals_home):
        pytest.skip("joblib is physically unvendored (e.g. as in debian)")

    # Check that a pickle that references sklearn.external.joblib can load
    f = tmpdir.join('foo.pkl')
    f.write(b'\x80\x02csklearn.externals.joblib.numpy_pickle\nNumpyArrayWrappe'
            b'r\nq\x00)\x81q\x01}q\x02(U\x05dtypeq\x03cnumpy\ndtype\nq\x04U'
            b'\x02i8q\x05K\x00K\x01\x87q\x06Rq\x07(K\x03U\x01<q\x08NNNJ\xff'
            b'\xff\xff\xffJ\xff\xff\xff\xffK\x00tq\tbU\x05shapeq\nK\x01\x85q'
            b'\x0bU\x05orderq\x0cU\x01Cq\rU\x08subclassq\x0ecnumpy\nndarray\nq'
            b'\x0fU\nallow_mmapq\x10\x88ub\x01\x00\x00\x00\x00\x00\x00\x00.',
            mode='wb')

    joblib_vendored.load(str(f))


def test_site_joblib_dispatch():
    if os.environ.get('SKLEARN_SITE_JOBLIB', False):
        assert Parallel is joblib_site.Parallel
        assert delayed is joblib_site.delayed
        assert parallel_backend is joblib_site.parallel_backend
        assert Memory is joblib_site.Memory

        assert joblib_vendored.Parallel is not joblib_site.Parallel
        assert joblib_vendored.delayed is not joblib_site.delayed
        assert joblib_vendored.parallel_backend is not \
            joblib_site.parallel_backend
        assert joblib_vendored.Memory is not joblib_site.Memory

    else:
        assert Parallel is joblib_vendored.Parallel
        assert delayed is joblib_vendored.delayed
        assert parallel_backend is joblib_vendored.parallel_backend
        assert Memory is joblib_vendored.Memory
