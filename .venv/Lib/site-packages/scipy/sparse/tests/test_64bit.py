""" Test functions involving 64bit or 32bit indexing """
import pytest
import numpy as np
from scipy.sparse import (
    bsr_array, coo_array, csc_array, csr_array, dia_array,
    bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dia_matrix,
)

# rename to avoid pytest collecting them in this module
from .test_base import (
    TestBSR as _TestBSR,
    TestCOO as _TestCOO,
    TestCSC as _TestCSC,
    TestCSR as _TestCSR,
    TestDIA as _TestDIA,
    TestDOK as _TestDOK,
    TestLIL as _TestLIL,
    TestBSRMatrix as _TestBSRMatrix,
    TestCOOMatrix as _TestCOOMatrix,
    TestCSCMatrix as _TestCSCMatrix,
    TestCSRMatrix as _TestCSRMatrix,
    TestDIAMatrix as _TestDIAMatrix,
    TestDOKMatrix as _TestDOKMatrix,
    TestLILMatrix as _TestLILMatrix,
    with_64bit_maxval_limit,
)


# name : reason not tested here
SKIP_TESTS = {
    'test_expm': 'expm for 64-bit indices not available',
    'test_inv': 'linsolve for 64-bit indices not available',
    'test_solve': 'linsolve for 64-bit indices not available',
    'test_scalar_idx_dtype': 'test implemented in base class',
    'test_large_dimensions_reshape': 'test actually requires 64-bit to work',
    'test_constructor_smallcol': 'test verifies int32 indexes',
    'test_constructor_largecol': 'test verifies int64 indexes',
    'test_tocoo_tocsr_tocsc_gh19245': 'test verifies int32 indexes',
}


def cases_64bit(sp_api):
    """Yield all tests for all formats

    This is more than testing get_index_dtype. It allows checking whether upcasting
    or downcasting the index dtypes affects test results. The approach used here
    does not try to figure out which tests might fail due to 32/64-bit issues.
    We just run them all.
    So, each test method in that uses cases_64bit reruns most of the test suite!
    """
    if sp_api == "sparray":
        TEST_CLASSES = [_TestBSR, _TestCOO, _TestCSC, _TestCSR, _TestDIA]
    elif sp_api == "sparray-extra":
        # lil/dok->other conversion operations use get_index_dtype
        # so we include lil & dok test suite even though they do not
        # use get_index_dtype within the class. That means many of
        # these tests are superfluous, but it's hard to pick which
        TEST_CLASSES = [_TestDOK, _TestLIL]
    elif sp_api == "spmatrix":
        TEST_CLASSES = [_TestBSRMatrix, _TestCOOMatrix, _TestCSCMatrix,
                        _TestCSRMatrix, _TestDIAMatrix]
    elif sp_api == "spmatrix-extra":
        # lil/dok->other conversion operations use get_index_dtype
        TEST_CLASSES = [_TestDOKMatrix, _TestLILMatrix]
    else:
        raise ValueError(f"parameter {sp_api=} is not valid")

    for cls in TEST_CLASSES:
        for method_name in sorted(dir(cls)):
            method = getattr(cls, method_name)
            if (method_name.startswith('test_') and
                    not getattr(method, 'slow', False)):
                marks = []

                msg = SKIP_TESTS.get(method_name)
                if msg:
                    marks.append(pytest.mark.skip(reason=msg))

                markers = getattr(method, 'pytestmark', [])
                for mark in markers:
                    if mark.name in ('skipif', 'skip', 'xfail', 'xslow'):
                        marks.append(mark)

                yield pytest.param(cls, method_name, marks=marks)


@pytest.mark.thread_unsafe(reason="fails in parallel")
class RunAll64Bit:
    def _check_resiliency(self, cls, method_name, **kw):
        # Resiliency test, to check that sparse matrices deal reasonably
        # with varying index data types.

        @with_64bit_maxval_limit(**kw)
        def check(cls, method_name):
            instance = cls()
            if hasattr(instance, 'setup_method'):
                instance.setup_method()
            try:
                getattr(instance, method_name)()
            finally:
                if hasattr(instance, 'teardown_method'):
                    instance.teardown_method()

        check(cls, method_name)


class Test64BitArray(RunAll64Bit):
    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray"))
    def test_resiliency_random(self, cls, method_name):
        self._check_resiliency(cls, method_name)


class Test64BitMatrix(RunAll64Bit):
    # assert_32bit=True only for spmatrix cuz sparray does not check index content
    @pytest.mark.fail_slow(5)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_no_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, assert_32bit=True)


class Test64BitMatrixSameAsArray(RunAll64Bit):
    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix"))
    def test_resiliency_random(self, cls, method_name):
        # Resiliency check that sparse deals with varying index data types.
        self._check_resiliency(cls, method_name)

# Extra: LIL and DOK classes. no direct get_index_dtype, but convert to classes that do
@pytest.mark.xslow
class Test64BitArrayExtra(RunAll64Bit):
    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("sparray-extra"))
    def test_resiliency_random(self, cls, method_name):
        # Resiliency check that sparse deals with varying index data types.
        self._check_resiliency(cls, method_name)


# Extra: LIL and DOK classes. no direct get_index_dtype, but convert to classes that do
@pytest.mark.xslow
class Test64BitMatrixExtra(RunAll64Bit):
    # assert_32bit=True only for spmatrix cuz sparray does not check index content
    @pytest.mark.fail_slow(5)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_no_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, assert_32bit=True)

    # inheritance of pytest test classes does not separate marks for subclasses.
    # So we define these functions in both Array and Matrix versions.
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_limit_10(self, cls, method_name):
        self._check_resiliency(cls, method_name, maxval_limit=10)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_all_32(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int32)

    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_all_64(self, cls, method_name):
        self._check_resiliency(cls, method_name, fixed_dtype=np.int64)

    @pytest.mark.fail_slow(2)
    @pytest.mark.parametrize('cls,method_name', cases_64bit("spmatrix-extra"))
    def test_resiliency_random(self, cls, method_name):
        # Resiliency check that sparse deals with varying index data types.
        self._check_resiliency(cls, method_name)


@pytest.mark.thread_unsafe(reason="Fails in parallel for unknown reasons")
class Test64BitTools:
    # classes that use get_index_dtype
    MAT_CLASSES = [
        bsr_matrix, coo_matrix, csc_matrix, csr_matrix, dia_matrix,
        bsr_array, coo_array, csc_array, csr_array, dia_array,
    ]

    def _compare_index_dtype(self, m, dtype):
        dtype = np.dtype(dtype)
        if m.format in ['csc', 'csr', 'bsr']:
            return (m.indices.dtype == dtype) and (m.indptr.dtype == dtype)
        elif m.format == 'coo':
            return (m.row.dtype == dtype) and (m.col.dtype == dtype)
        elif m.format == 'dia':
            return (m.offsets.dtype == dtype)
        else:
            raise ValueError(f"matrix {m!r} has no integer indices")

    def test_decorator_maxval_limit(self):
        # Test that the with_64bit_maxval_limit decorator works

        @with_64bit_maxval_limit(maxval_limit=10)
        def check(mat_cls):
            m = mat_cls(np.random.rand(10, 1))
            assert self._compare_index_dtype(m, np.int32)
            m = mat_cls(np.random.rand(11, 1))
            assert self._compare_index_dtype(m, np.int64)

        for mat_cls in self.MAT_CLASSES:
            check(mat_cls)

    def test_decorator_maxval_random(self):
        # Test that the with_64bit_maxval_limit decorator works (2)

        @with_64bit_maxval_limit(random=True)
        def check(mat_cls):
            seen_32 = False
            seen_64 = False
            for k in range(100):
                m = mat_cls(np.random.rand(9, 9))
                seen_32 = seen_32 or self._compare_index_dtype(m, np.int32)
                seen_64 = seen_64 or self._compare_index_dtype(m, np.int64)
                if seen_32 and seen_64:
                    break
            else:
                raise AssertionError("both 32 and 64 bit indices not seen")

        for mat_cls in self.MAT_CLASSES:
            check(mat_cls)

    def test_downcast_intp(self):
        # Check that bincount and ufunc.reduceat intp downcasts are
        # dealt with. The point here is to trigger points in the code
        # that can fail on 32-bit systems when using 64-bit indices,
        # due to use of functions that only work with intp-size indices.

        @with_64bit_maxval_limit(fixed_dtype=np.int64, downcast_maxval=1)
        def check_limited(csc_container, csr_container, coo_container):
            # These involve indices larger than `downcast_maxval`
            a = csc_container([[1, 2], [3, 4], [5, 6]])
            pytest.raises(AssertionError, a.count_nonzero, axis=1)
            pytest.raises(AssertionError, a.sum, axis=0)

            a = csr_container([[1, 2, 3], [3, 4, 6]])
            pytest.raises(AssertionError, a.count_nonzero, axis=0)
            pytest.raises(AssertionError, a.sum, axis=1)

            a = coo_container([[1, 2, 3], [3, 4, 5]])
            pytest.raises(AssertionError, a.count_nonzero, axis=0)
            a.has_canonical_format = False
            pytest.raises(AssertionError, a.sum_duplicates)

        @with_64bit_maxval_limit(fixed_dtype=np.int64)
        def check_unlimited(csc_container, csr_container, coo_container):
            # These involve indices smaller than `downcast_maxval`
            a = csc_container([[1, 2], [3, 4], [5, 6]])
            a.count_nonzero(axis=1)
            a.sum(axis=0)

            a = csr_container([[1, 2, 3], [3, 4, 6]])
            a.count_nonzero(axis=0)
            a.sum(axis=1)

            a = coo_container([[1, 2, 3], [3, 4, 5]])
            a.count_nonzero(axis=0)
            a.has_canonical_format = False
            a.sum_duplicates()

        check_limited(csc_array, csr_array, coo_array)
        check_unlimited(csc_array, csr_array, coo_array)
        check_limited(csc_matrix, csr_matrix, coo_matrix)
        check_unlimited(csc_matrix, csr_matrix, coo_matrix)
