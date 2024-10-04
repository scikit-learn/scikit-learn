import numpy as np
from numpy.testing import assert_equal
import pytest
from scipy.sparse import coo_array


def test_shape_constructor():
    empty1d = coo_array((3,))
    assert empty1d.shape == (3,)
    assert_equal(empty1d.toarray(), np.zeros((3,)))

    empty2d = coo_array((3, 2))
    assert empty2d.shape == (3, 2)
    assert_equal(empty2d.toarray(), np.zeros((3, 2)))

    with pytest.raises(TypeError, match='invalid input format'):
        coo_array((3, 2, 2))


def test_dense_constructor():
    res1d = coo_array([1, 2, 3])
    assert res1d.shape == (3,)
    assert_equal(res1d.toarray(), np.array([1, 2, 3]))

    res2d = coo_array([[1, 2, 3], [4, 5, 6]])
    assert res2d.shape == (2, 3)
    assert_equal(res2d.toarray(), np.array([[1, 2, 3], [4, 5, 6]]))

    with pytest.raises(ValueError, match='shape must be a 1- or 2-tuple'):
        coo_array([[[3]], [[4]]])


def test_dense_constructor_with_shape():
    res1d = coo_array([1, 2, 3], shape=(3,))
    assert res1d.shape == (3,)
    assert_equal(res1d.toarray(), np.array([1, 2, 3]))

    res2d = coo_array([[1, 2, 3], [4, 5, 6]], shape=(2, 3))
    assert res2d.shape == (2, 3)
    assert_equal(res2d.toarray(), np.array([[1, 2, 3], [4, 5, 6]]))

    with pytest.raises(ValueError, match='shape must be a 1- or 2-tuple'):
        coo_array([[[3]], [[4]]], shape=(2, 1, 1))


def test_dense_constructor_with_inconsistent_shape():
    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([1, 2, 3], shape=(4,))

    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([1, 2, 3], shape=(3, 1))

    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([[1, 2, 3]], shape=(3,))

    with pytest.raises(ValueError,
                       match='axis 0 index 2 exceeds matrix dimension 2'):
        coo_array(([1], ([2],)), shape=(2,))

    with pytest.raises(ValueError, match='negative axis 0 index: -1'):
        coo_array(([1], ([-1],)))


def test_1d_sparse_constructor():
    empty1d = coo_array((3,))
    res = coo_array(empty1d)
    assert res.shape == (3,)
    assert_equal(res.toarray(), np.zeros((3,)))


def test_1d_tuple_constructor():
    res = coo_array(([9,8], ([1,2],)))
    assert res.shape == (3,)
    assert_equal(res.toarray(), np.array([0, 9, 8]))


def test_1d_tuple_constructor_with_shape():
    res = coo_array(([9,8], ([1,2],)), shape=(4,))
    assert res.shape == (4,)
    assert_equal(res.toarray(), np.array([0, 9, 8, 0]))

def test_non_subscriptability():
    coo_2d = coo_array((2, 2))

    with pytest.raises(TypeError,
                        match="'coo_array' object does not support item assignment"):
        coo_2d[0, 0] = 1

    with pytest.raises(TypeError,
                       match="'coo_array' object is not subscriptable"):
        coo_2d[0, :]

def test_reshape():
    arr1d = coo_array([1, 0, 3])
    assert arr1d.shape == (3,)

    col_vec = arr1d.reshape((3, 1))
    assert col_vec.shape == (3, 1)
    assert_equal(col_vec.toarray(), np.array([[1], [0], [3]]))

    row_vec = arr1d.reshape((1, 3))
    assert row_vec.shape == (1, 3)
    assert_equal(row_vec.toarray(), np.array([[1, 0, 3]]))

    arr2d = coo_array([[1, 2, 0], [0, 0, 3]])
    assert arr2d.shape == (2, 3)

    flat = arr2d.reshape((6,))
    assert flat.shape == (6,)
    assert_equal(flat.toarray(), np.array([1, 2, 0, 0, 0, 3]))


def test_nnz():
    arr1d = coo_array([1, 0, 3])
    assert arr1d.shape == (3,)
    assert arr1d.nnz == 2

    arr2d = coo_array([[1, 2, 0], [0, 0, 3]])
    assert arr2d.shape == (2, 3)
    assert arr2d.nnz == 3


def test_transpose():
    arr1d = coo_array([1, 0, 3]).T
    assert arr1d.shape == (3,)
    assert_equal(arr1d.toarray(), np.array([1, 0, 3]))

    arr2d = coo_array([[1, 2, 0], [0, 0, 3]]).T
    assert arr2d.shape == (3, 2)
    assert_equal(arr2d.toarray(), np.array([[1, 0], [2, 0], [0, 3]]))


def test_transpose_with_axis():
    arr1d = coo_array([1, 0, 3]).transpose(axes=(0,))
    assert arr1d.shape == (3,)
    assert_equal(arr1d.toarray(), np.array([1, 0, 3]))

    arr2d = coo_array([[1, 2, 0], [0, 0, 3]]).transpose(axes=(0, 1))
    assert arr2d.shape == (2, 3)
    assert_equal(arr2d.toarray(), np.array([[1, 2, 0], [0, 0, 3]]))

    with pytest.raises(ValueError, match="axes don't match matrix dimensions"):
        coo_array([1, 0, 3]).transpose(axes=(0, 1))

    with pytest.raises(ValueError, match="repeated axis in transpose"):
        coo_array([[1, 2, 0], [0, 0, 3]]).transpose(axes=(1, 1))


def test_1d_row_and_col():
    res = coo_array([1, -2, -3])
    assert_equal(res.col, np.array([0, 1, 2]))
    assert_equal(res.row, np.zeros_like(res.col))
    assert res.row.dtype == res.col.dtype
    assert res.row.flags.writeable is False

    res.col = [1, 2, 3]
    assert len(res.coords) == 1
    assert_equal(res.col, np.array([1, 2, 3]))
    assert res.row.dtype == res.col.dtype

    with pytest.raises(ValueError, match="cannot set row attribute"):
        res.row = [1, 2, 3]


def test_1d_toformats():
    res = coo_array([1, -2, -3])
    for f in [res.tobsr, res.tocsc, res.todia, res.tolil]:
        with pytest.raises(ValueError, match='Cannot convert'):
            f()
    for f in [res.tocoo, res.tocsr, res.todok]:
        assert_equal(f().toarray(), res.toarray())


@pytest.mark.parametrize('arg', [1, 2, 4, 5, 8])
def test_1d_resize(arg: int):
    den = np.array([1, -2, -3])
    res = coo_array(den)
    den.resize(arg, refcheck=False)
    res.resize(arg)
    assert res.shape == den.shape
    assert_equal(res.toarray(), den)


@pytest.mark.parametrize('arg', zip([1, 2, 3, 4], [1, 2, 3, 4]))
def test_1d_to_2d_resize(arg: tuple[int, int]):
    den = np.array([1, 0, 3])
    res = coo_array(den)

    den.resize(arg, refcheck=False)
    res.resize(arg)
    assert res.shape == den.shape
    assert_equal(res.toarray(), den)


@pytest.mark.parametrize('arg', [1, 4, 6, 8])
def test_2d_to_1d_resize(arg: int):
    den = np.array([[1, 0, 3], [4, 0, 0]])
    res = coo_array(den)
    den.resize(arg, refcheck=False)
    res.resize(arg)
    assert res.shape == den.shape
    assert_equal(res.toarray(), den)


def test_sum_duplicates():
    arr1d = coo_array(([2, 2, 2], ([1, 0, 1],)))
    assert arr1d.nnz == 3
    assert_equal(arr1d.toarray(), np.array([2, 4]))
    arr1d.sum_duplicates()
    assert arr1d.nnz == 2
    assert_equal(arr1d.toarray(), np.array([2, 4]))


def test_eliminate_zeros():
    arr1d = coo_array(([0, 0, 1], ([1, 0, 1],)))
    assert arr1d.nnz == 3
    assert arr1d.count_nonzero() == 1
    assert_equal(arr1d.toarray(), np.array([0, 1]))
    arr1d.eliminate_zeros()
    assert arr1d.nnz == 1
    assert arr1d.count_nonzero() == 1
    assert_equal(arr1d.toarray(), np.array([0, 1]))
    assert_equal(arr1d.col, np.array([1]))
    assert_equal(arr1d.row, np.array([0]))


def test_1d_add_dense():
    den_a = np.array([0, -2, -3, 0])
    den_b = np.array([0, 1, 2, 3])
    exp = den_a + den_b
    res = coo_array(den_a) + den_b
    assert type(res) == type(exp)
    assert_equal(res, exp)


def test_1d_add_sparse():
    den_a = np.array([0, -2, -3, 0])
    den_b = np.array([0, 1, 2, 3])
    dense_sum = den_a + den_b
    # this routes through CSR format
    sparse_sum = coo_array(den_a) + coo_array(den_b)
    assert_equal(dense_sum, sparse_sum.toarray())


def test_1d_matmul_vector():
    den_a = np.array([0, -2, -3, 0])
    den_b = np.array([0, 1, 2, 3])
    exp = den_a @ den_b
    res = coo_array(den_a) @ den_b
    assert np.ndim(res) == 0
    assert_equal(res, exp)


def test_1d_matmul_multivector():
    den = np.array([0, -2, -3, 0])
    other = np.array([[0, 1, 2, 3], [3, 2, 1, 0]]).T
    exp = den @ other
    res = coo_array(den) @ other
    assert type(res) == type(exp)
    assert_equal(res, exp)


def test_2d_matmul_multivector():
    den = np.array([[0, 1, 2, 3], [3, 2, 1, 0]])
    arr2d = coo_array(den)
    exp = den @ den.T
    res = arr2d @ arr2d.T
    assert_equal(res.toarray(), exp)


def test_1d_diagonal():
    den = np.array([0, -2, -3, 0])
    with pytest.raises(ValueError, match='diagonal requires two dimensions'):
        coo_array(den).diagonal()
