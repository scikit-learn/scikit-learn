import numpy as np
from numpy.testing import assert_equal
import pytest
from scipy.linalg import block_diag
from scipy.sparse import coo_array, random_array
from .._coo import _block_diag, _extract_block_diag


def test_shape_constructor():
    empty1d = coo_array((3,))
    assert empty1d.shape == (3,)
    assert_equal(empty1d.toarray(), np.zeros((3,)))

    empty2d = coo_array((3, 2))
    assert empty2d.shape == (3, 2)
    assert_equal(empty2d.toarray(), np.zeros((3, 2)))

    empty_nd = coo_array((2,3,4,6,7))
    assert empty_nd.shape == (2,3,4,6,7)
    assert_equal(empty_nd.toarray(), np.zeros((2,3,4,6,7)))


def test_dense_constructor():
    # 1d
    res1d = coo_array([1, 2, 3])
    assert res1d.shape == (3,)
    assert_equal(res1d.toarray(), np.array([1, 2, 3]))

    # 2d
    res2d = coo_array([[1, 2, 3], [4, 5, 6]])
    assert res2d.shape == (2, 3)
    assert_equal(res2d.toarray(), np.array([[1, 2, 3], [4, 5, 6]]))

    # 4d
    arr4d = np.array([[[[3, 7], [1, 0]], [[6, 5], [9, 2]]],
                      [[[4, 3], [2, 8]], [[7, 5], [1, 6]]],
                      [[[0, 9], [4, 3]], [[2, 1], [7, 8]]]])
    res4d = coo_array(arr4d)
    assert res4d.shape == (3, 2, 2, 2)
    assert_equal(res4d.toarray(), arr4d)

    # 9d
    np.random.seed(12)
    arr9d = np.random.randn(2,3,4,7,6,5,3,2,4)
    res9d = coo_array(arr9d)
    assert res9d.shape == (2,3,4,7,6,5,3,2,4)
    assert_equal(res9d.toarray(), arr9d)

    # storing nan as element of sparse array
    nan_3d = coo_array([[[1, np.nan]], [[3, 4]], [[5, 6]]])
    assert nan_3d.shape == (3, 1, 2)
    assert_equal(nan_3d.toarray(), np.array([[[1, np.nan]], [[3, 4]], [[5, 6]]]))


def test_dense_constructor_with_shape():
    res1d = coo_array([1, 2, 3], shape=(3,))
    assert res1d.shape == (3,)
    assert_equal(res1d.toarray(), np.array([1, 2, 3]))

    res2d = coo_array([[1, 2, 3], [4, 5, 6]], shape=(2, 3))
    assert res2d.shape == (2, 3)
    assert_equal(res2d.toarray(), np.array([[1, 2, 3], [4, 5, 6]]))

    res3d = coo_array([[[3]], [[4]]], shape=(2, 1, 1))
    assert res3d.shape == (2, 1, 1)
    assert_equal(res3d.toarray(), np.array([[[3]], [[4]]]))

    np.random.seed(12)
    arr7d = np.random.randn(2,4,1,6,5,3,2)
    res7d = coo_array((arr7d), shape=(2,4,1,6,5,3,2))
    assert res7d.shape == (2,4,1,6,5,3,2)
    assert_equal(res7d.toarray(), arr7d)


def test_dense_constructor_with_inconsistent_shape():
    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([1, 2, 3], shape=(4,))

    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([1, 2, 3], shape=(3, 1))

    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([[1, 2, 3]], shape=(3,))

    with pytest.raises(ValueError, match='inconsistent shapes'):
        coo_array([[[3]], [[4]]], shape=(1, 1, 1))

    with pytest.raises(ValueError,
                       match='axis 0 index 2 exceeds matrix dimension 2'):
        coo_array(([1], ([2],)), shape=(2,))

    with pytest.raises(ValueError,
                       match='axis 1 index 3 exceeds matrix dimension 3'):
        coo_array(([1,3], ([0, 1], [0, 3], [1, 1])), shape=(2, 3, 2))

    with pytest.raises(ValueError, match='negative axis 0 index: -1'):
        coo_array(([1], ([-1],)))

    with pytest.raises(ValueError, match='negative axis 2 index: -1'):
        coo_array(([1], ([0], [2], [-1])))


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

def test_reshape_overflow():
    # see gh-22353 : new idx_dtype can need to be int64 instead of int32
    M, N = (1045507, 523266)
    coords = (np.array([M - 1], dtype='int32'), np.array([N - 1], dtype='int32'))
    A = coo_array(([3.3], coords), shape=(M, N))

    # need new idx_dtype to not overflow
    B = A.reshape((M * N, 1))
    assert B.coords[0].dtype == np.dtype('int64')
    assert B.coords[0][0] == (M * N) - 1

    # need idx_dtype to stay int32 if before and after can be int32
    C = A.reshape(N, M)
    assert C.coords[0].dtype == np.dtype('int32')
    assert C.coords[0][0] == N - 1

def test_reshape():
    arr1d = coo_array([1, 0, 3])
    assert arr1d.shape == (3,)

    col_vec = arr1d.reshape((3, 1))
    assert col_vec.shape == (3, 1)
    assert_equal(col_vec.toarray(), np.array([[1], [0], [3]]))

    row_vec = arr1d.reshape((1, 3))
    assert row_vec.shape == (1, 3)
    assert_equal(row_vec.toarray(), np.array([[1, 0, 3]]))

    # attempting invalid reshape
    with pytest.raises(ValueError, match="cannot reshape array"):
        arr1d.reshape((3,3))

    # attempting reshape with a size 0 dimension
    with pytest.raises(ValueError, match="cannot reshape array"):
        arr1d.reshape((3,0))

    arr2d = coo_array([[1, 2, 0], [0, 0, 3]])
    assert arr2d.shape == (2, 3)

    flat = arr2d.reshape((6,))
    assert flat.shape == (6,)
    assert_equal(flat.toarray(), np.array([1, 2, 0, 0, 0, 3]))

    # 2d to 3d
    to_3d_arr = arr2d.reshape((2, 3, 1))
    assert to_3d_arr.shape == (2, 3, 1)
    assert_equal(to_3d_arr.toarray(), np.array([[[1], [2], [0]], [[0], [0], [3]]]))

    # attempting invalid reshape
    with pytest.raises(ValueError, match="cannot reshape array"):
        arr2d.reshape((1,3))


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
    # 1d case
    arr1d = coo_array(([2, 2, 2], ([1, 0, 1],)))
    assert arr1d.nnz == 3
    assert_equal(arr1d.toarray(), np.array([2, 4]))
    arr1d.sum_duplicates()
    assert arr1d.nnz == 2
    assert_equal(arr1d.toarray(), np.array([2, 4]))

    # 4d case
    arr4d = coo_array(([2, 3, 7], ([1, 0, 1], [0, 2, 0], [1, 2, 1], [1, 0, 1])))
    assert arr4d.nnz == 3
    expected = np.array(  # noqa: E501
        [[[[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [3, 0]]],
         [[[0, 0], [0, 9], [0, 0]], [[0, 0], [0, 0], [0, 0]], [[0, 0], [0, 0], [0, 0]]]]
    )
    assert_equal(arr4d.toarray(), expected)
    arr4d.sum_duplicates()
    assert arr4d.nnz == 2
    assert_equal(arr4d.toarray(), expected)

    # when there are no duplicates
    arr_nodups = coo_array(([1, 2, 3, 4], ([0, 1, 2, 3],)))
    assert arr_nodups.nnz == 4
    arr_nodups.sum_duplicates()
    assert arr_nodups.nnz == 4


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
    assert type(res) is type(exp)
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
    assert type(res) is type(exp)
    assert_equal(res, exp)


def test_2d_matmul_multivector():
    # sparse-sparse matmul
    den = np.array([[0, 1, 2, 3], [3, 2, 1, 0]])
    arr2d = coo_array(den)
    exp = den @ den.T
    res = arr2d @ arr2d.T
    assert_equal(res.toarray(), exp)

    # sparse-dense matmul for self.ndim = 2
    den = np.array([[0, 4, 3, 0, 5], [1, 0, 7, 3, 4]])
    arr2d = coo_array(den)
    exp = den @ den.T
    res = arr2d @ den.T
    assert_equal(res, exp)

    # sparse-dense matmul for self.ndim = 1
    den_a = np.array([[0, 4, 3, 0, 5], [1, 0, 7, 3, 4]])
    den_b = np.array([0, 1, 6, 0, 4])
    arr1d = coo_array(den_b)
    exp = den_b @ den_a.T
    res = arr1d @ den_a.T
    assert_equal(res, exp)

    # sparse-dense matmul for self.ndim = 1 and other.ndim = 2
    den_a = np.array([1, 0, 2])
    den_b = np.array([[3], [4], [0]])
    exp = den_a @ den_b
    res = coo_array(den_a) @ den_b
    assert_equal(res, exp)
    res = coo_array(den_a) @ list(den_b)
    assert_equal(res, exp)


def test_1d_diagonal():
    den = np.array([0, -2, -3, 0])
    with pytest.raises(ValueError, match='diagonal requires two dimensions'):
        coo_array(den).diagonal()


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_todense(shape):
    np.random.seed(12)
    arr = np.random.randint(low=0, high=5, size=shape)
    assert_equal(coo_array(arr).todense(), arr)


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_sparse_constructor(shape):
    empty_arr = coo_array(shape)
    res = coo_array(empty_arr)
    assert res.shape == (shape)
    assert_equal(res.toarray(), np.zeros(shape))


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_tuple_constructor(shape):
    np.random.seed(12)
    arr = np.random.randn(*shape)
    res = coo_array(arr)
    assert res.shape == shape
    assert_equal(res.toarray(), arr)


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_tuple_constructor_with_shape(shape):
    np.random.seed(12)
    arr = np.random.randn(*shape)
    res = coo_array(arr, shape=shape)
    assert res.shape == shape
    assert_equal(res.toarray(), arr)


def test_tuple_constructor_for_dim_size_zero():
    # arrays with a dimension of size 0
    with pytest.raises(ValueError, match='exceeds matrix dimension'):
        coo_array(([9, 8], ([1, 2], [1, 0], [2, 1])), shape=(3,4,0))

    empty_arr = coo_array(([], ([], [], [], [])), shape=(4,0,2,3))
    assert_equal(empty_arr.toarray(), np.empty((4,0,2,3)))


@pytest.mark.parametrize(('shape', 'new_shape'), [((4,9,6,5), (3,6,15,4)),
                                                  ((4,9,6,5), (36,30)),
                                                  ((4,9,6,5), (1080,)),
                                                  ((4,9,6,5), (2,3,2,2,3,5,3)),])
def test_nd_reshape(shape, new_shape):
    # reshaping a 4d sparse array
    rng = np.random.default_rng(23409823)

    arr4d = random_array(shape, density=0.6, rng=rng, dtype=int)
    assert arr4d.shape == shape
    den4d = arr4d.toarray()

    exp_arr = den4d.reshape(new_shape)
    res_arr = arr4d.reshape(new_shape)
    assert res_arr.shape == new_shape
    assert_equal(res_arr.toarray(), exp_arr)


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_nnz(shape):
    rng = np.random.default_rng(23409823)

    arr = random_array(shape, density=0.6, rng=rng, dtype=int)
    assert arr.nnz == np.count_nonzero(arr.toarray())


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_transpose(shape):
    rng = np.random.default_rng(23409823)

    arr = random_array(shape, density=0.6, rng=rng, dtype=int)
    exp_arr = arr.toarray().T
    trans_arr = arr.transpose()
    assert trans_arr.shape == shape[::-1]
    assert_equal(exp_arr, trans_arr.toarray())


@pytest.mark.parametrize(('shape', 'axis_perm'), [((3,), (0,)),
                                                  ((2,3), (0,1)),
                                                  ((2,4,3,6,5,3), (1,2,0,5,3,4)),])
def test_nd_transpose_with_axis(shape, axis_perm):
    rng = np.random.default_rng(23409823)

    arr = random_array(shape, density=0.6, rng=rng, dtype=int)
    trans_arr = arr.transpose(axes=axis_perm)
    assert_equal(trans_arr.toarray(), np.transpose(arr.toarray(), axes=axis_perm))


def test_transpose_with_inconsistent_axis():
    with pytest.raises(ValueError, match="axes don't match matrix dimensions"):
        coo_array([1, 0, 3]).transpose(axes=(0, 1))

    with pytest.raises(ValueError, match="repeated axis in transpose"):
        coo_array([[1, 2, 0], [0, 0, 3]]).transpose(axes=(1, 1))


def test_nd_eliminate_zeros():
    # for 3d sparse arrays
    arr3d = coo_array(([1, 0, 0, 4], ([0, 1, 1, 2], [0, 1, 0, 1], [1, 1, 2, 0])))
    assert arr3d.nnz == 4
    assert arr3d.count_nonzero() == 2
    assert_equal(arr3d.toarray(), np.array([[[0, 1, 0], [0, 0, 0]],
                                    [[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [4, 0, 0]]]))
    arr3d.eliminate_zeros()
    assert arr3d.nnz == 2
    assert arr3d.count_nonzero() == 2
    assert_equal(arr3d.toarray(), np.array([[[0, 1, 0], [0, 0, 0]],
                                    [[0, 0, 0], [0, 0, 0]], [[0, 0, 0], [4, 0, 0]]]))

    # for a 5d sparse array when all elements of data array are 0
    coords = ([0, 1, 1, 2], [0, 1, 0, 1], [1, 1, 2, 0], [0, 0, 2, 3], [1, 0, 0, 2])
    arr5d = coo_array(([0, 0, 0, 0], coords))
    assert arr5d.nnz == 4
    assert arr5d.count_nonzero() == 0
    arr5d.eliminate_zeros()
    assert arr5d.nnz == 0
    assert arr5d.count_nonzero() == 0
    assert_equal(arr5d.col, np.array([]))
    assert_equal(arr5d.row, np.array([]))
    assert_equal(arr5d.coords, ([], [], [], [], []))


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_add_dense(shape):
    rng = np.random.default_rng(23409823)
    sp_x = random_array(shape, density=0.6, rng=rng, dtype=int)
    sp_y = random_array(shape, density=0.6, rng=rng, dtype=int)
    den_x, den_y = sp_x.toarray(), sp_y.toarray()
    exp = den_x + den_y
    res = sp_x + den_y
    assert type(res) is type(exp)
    assert_equal(res, exp)


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_add_sparse(shape):
    rng = np.random.default_rng(23409823)
    sp_x = random_array((shape), density=0.6, rng=rng, dtype=int)
    sp_y = random_array((shape), density=0.6, rng=rng, dtype=int)
    den_x, den_y = sp_x.toarray(), sp_y.toarray()

    dense_sum = den_x + den_y
    sparse_sum = sp_x + sp_y
    assert_equal(dense_sum, sparse_sum.toarray())


def test_add_sparse_with_inf():
    # addition of sparse arrays with an inf element
    den_a = np.array([[[0], [np.inf]], [[-3], [0]]])
    den_b = np.array([[[0], [1]], [[2], [3]]])
    dense_sum = den_a + den_b
    sparse_sum = coo_array(den_a) + coo_array(den_b)
    assert_equal(dense_sum, sparse_sum.toarray())


@pytest.mark.parametrize(('a_shape', 'b_shape'), [((7,), (12,)),
                                                  ((6,4), (6,5)),
                                                  ((5,9,3,2), (9,5,2,3)),])
def test_nd_add_sparse_with_inconsistent_shapes(a_shape, b_shape):
    rng = np.random.default_rng(23409823)

    arr_a = random_array((a_shape), density=0.6, rng=rng, dtype=int)
    arr_b = random_array((b_shape), density=0.6, rng=rng, dtype=int)
    with pytest.raises(ValueError,
                       match="(Incompatible|inconsistent) shapes|cannot be broadcast"):
        arr_a + arr_b


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_sub_dense(shape):
    rng = np.random.default_rng(23409823)
    sp_x = random_array(shape, density=0.6, rng=rng, dtype=int)
    sp_y = random_array(shape, density=0.6, rng=rng, dtype=int)
    den_x, den_y = sp_x.toarray(), sp_y.toarray()
    exp = den_x - den_y
    res = sp_x - den_y
    assert type(res) is type(exp)
    assert_equal(res, exp)


@pytest.mark.parametrize('shape', [(0,), (7,), (4,7), (0,0,0), (3,6,2),
                                   (1,0,3), (7,9,3,2,4,5)])
def test_nd_sub_sparse(shape):
    rng = np.random.default_rng(23409823)

    sp_x = random_array(shape, density=0.6, rng=rng, dtype=int)
    sp_y = random_array(shape, density=0.6, rng=rng, dtype=int)
    den_x, den_y = sp_x.toarray(), sp_y.toarray()

    dense_sum = den_x - den_y
    sparse_sum = sp_x - sp_y
    assert_equal(dense_sum, sparse_sum.toarray())


def test_nd_sub_sparse_with_nan():
    # subtraction of sparse arrays with a nan element
    den_a = np.array([[[0], [np.nan]], [[-3], [0]]])
    den_b = np.array([[[0], [1]], [[2], [3]]])
    dense_sum = den_a - den_b
    sparse_sum = coo_array(den_a) - coo_array(den_b)
    assert_equal(dense_sum, sparse_sum.toarray())


@pytest.mark.parametrize(('a_shape', 'b_shape'), [((7,), (12,)),
                                                  ((6,4), (6,5)),
                                                  ((5,9,3,2), (9,5,2,3)),])
def test_nd_sub_sparse_with_inconsistent_shapes(a_shape, b_shape):
    rng = np.random.default_rng(23409823)

    arr_a = random_array((a_shape), density=0.6, rng=rng, dtype=int)
    arr_b = random_array((b_shape), density=0.6, rng=rng, dtype=int)
    with pytest.raises(ValueError, match="inconsistent shapes"):
        arr_a - arr_b


mat_vec_shapes = [
    ((2, 3, 4, 5), (5,)),
    ((0, 0), (0,)),
    ((2, 3, 4, 7, 8), (8,)),
    ((4, 4, 2, 0), (0,)),
    ((6, 5, 3, 2, 4), (4, 1)),
    ((2,5), (5,)),
    ((2, 5), (5, 1)),
    ((3,), (3, 1)),
    ((4,), (4,))
]
@pytest.mark.parametrize(('mat_shape', 'vec_shape'), mat_vec_shapes)
def test_nd_matmul_vector(mat_shape, vec_shape):
    rng = np.random.default_rng(23409823)

    sp_x = random_array(mat_shape, density=0.6, rng=rng, dtype=int)
    sp_y = random_array(vec_shape, density=0.6, rng=rng, dtype=int)
    den_x, den_y = sp_x.toarray(), sp_y.toarray()
    exp = den_x @ den_y
    res = sp_x @ den_y
    assert_equal(res,exp)
    res = sp_x @ list(den_y)
    assert_equal(res,exp)


mat_mat_shapes = [
    ((2, 3, 4, 5), (2, 3, 5, 7)),
    ((0, 0), (0,)),
    ((4, 4, 2, 0), (0,)),
    ((7, 8, 3), (3,)),
    ((7, 8, 3), (3, 1)),
    ((6, 5, 3, 2, 4), (4, 3)),
    ((1, 3, 2, 4), (6, 5, 1, 4, 3)),
    ((6, 1, 1, 2, 4), (1, 3, 4, 3)),
    ((4,), (2, 4, 3)),
    ((3,), (5, 6, 7, 3, 2)),
    ((4,), (4, 3)),
    ((2, 5), (5, 1)),
]
@pytest.mark.parametrize(('mat_shape1', 'mat_shape2'), mat_mat_shapes)
def test_nd_matmul(mat_shape1, mat_shape2):
    rng = np.random.default_rng(23409823)

    sp_x = random_array(mat_shape1, density=0.6, random_state=rng, dtype=int)
    sp_y = random_array(mat_shape2, density=0.6, random_state=rng, dtype=int)
    den_x, den_y = sp_x.toarray(), sp_y.toarray()
    exp = den_x @ den_y
    # sparse-sparse
    res = sp_x @ sp_y
    assert_equal(res.toarray(), exp)
    # sparse-dense
    res = sp_x @ den_y
    assert_equal(res, exp)
    res = sp_x @ list(den_y)
    assert_equal(res, exp)

    # dense-sparse
    res = den_x @ sp_y
    assert_equal(res, exp)


def test_nd_matmul_sparse_with_inconsistent_arrays():
    rng = np.random.default_rng(23409823)

    sp_x = random_array((4,5,7,6,3), density=0.6, random_state=rng, dtype=int)
    sp_y = random_array((1,5,3,2,5), density=0.6, random_state=rng, dtype=int)
    with pytest.raises(ValueError, match="matmul: dimension mismatch with signature"):
        sp_x @ sp_y
    with pytest.raises(ValueError, match="matmul: dimension mismatch with signature"):
        sp_x @ (sp_y.toarray())

    sp_z = random_array((1,5,3,2), density=0.6, random_state=rng, dtype=int)
    with pytest.raises(ValueError, match="Batch dimensions are not broadcastable"):
        sp_x @ sp_z
    with pytest.raises(ValueError, match="Batch dimensions are not broadcastable"):
        sp_x @ (sp_z.toarray())


def test_dot_1d_1d(): # 1-D inner product
    a = coo_array([1,2,3])
    b = coo_array([4,5,6])
    exp = np.dot(a.toarray(), b.toarray())
    res = a.dot(b)
    assert_equal(res, exp)
    res = a.dot(b.toarray())
    assert_equal(res, exp)


def test_dot_sparse_scalar():
    a = coo_array([[1, 2], [3, 4], [5, 6]])
    b = 3
    res = a.dot(b)
    exp = np.dot(a.toarray(), b)
    assert_equal(res.toarray(), exp)


def test_dot_with_inconsistent_shapes():
    arr_a = coo_array([[[1, 2]], [[3, 4]]])
    arr_b = coo_array([4, 5, 6])
    with pytest.raises(ValueError, match="not aligned for n-D dot"):
        arr_a.dot(arr_b)


def test_matmul_dot_not_implemented():
    arr_a = coo_array([[1, 2], [3, 4]])
    with pytest.raises(TypeError, match="argument not supported type"):
        arr_a.dot(None)
    with pytest.raises(TypeError, match="arg not supported type"):
        arr_a.tensordot(None)
    with pytest.raises(TypeError, match="unsupported operand type"):
        arr_a @ None
    with pytest.raises(TypeError, match="unsupported operand type"):
        None @ arr_a


dot_shapes = [
    ((3,3), (3,3)), ((4,6), (6,7)), ((1,4), (4,1)), # matrix multiplication 2-D
    ((3,2,4,7), (7,)), ((5,), (6,3,5,2)), # dot of n-D and 1-D arrays
    ((3,2,4,7), (7,1)), ((1,5,), (6,3,5,2)),
    ((4,6), (3,2,6,4)), ((2,8,7), (4,5,7,7,2)), # dot of n-D and m-D arrays
    ((4,5,7,6), (3,2,6,4)),
]
@pytest.mark.parametrize(('a_shape', 'b_shape'), dot_shapes)
def test_dot_nd(a_shape, b_shape):
    rng = np.random.default_rng(23409823)

    arr_a = random_array(a_shape, density=0.6, random_state=rng, dtype=int)
    arr_b = random_array(b_shape, density=0.6, random_state=rng, dtype=int)

    exp = np.dot(arr_a.toarray(), arr_b.toarray())
    # sparse-dense
    res = arr_a.dot(arr_b.toarray())
    assert_equal(res, exp)
    res = arr_a.dot(list(arr_b.toarray()))
    assert_equal(res, exp)
    # sparse-sparse
    res = arr_a.dot(arr_b)
    assert_equal(res.toarray(), exp)


tensordot_shapes_and_axes = [
    ((4,6), (6,7), ([1], [0])),
    ((3,2,4,7), (7,), ([3], [0])),
    ((5,), (6,3,5,2), ([0], [2])),
    ((4,5,7,6), (3,2,6,4), ([0, 3], [3, 2])),
    ((2,8,7), (4,5,7,8,2), ([0, 1, 2], [4, 3, 2])),
    ((4,5,3,2,6), (3,2,6,7,8), 3),
    ((4,5,7), (7,3,7), 1),
    ((2,3,4), (2,3,4), ([0, 1, 2], [0, 1, 2])),
]
@pytest.mark.parametrize(('a_shape', 'b_shape', 'axes'), tensordot_shapes_and_axes)
def test_tensordot(a_shape, b_shape, axes):
    rng = np.random.default_rng(23409823)

    arr_a = random_array(a_shape, density=0.6, random_state=rng, dtype=int)
    arr_b = random_array(b_shape, density=0.6, random_state=rng, dtype=int)

    exp = np.tensordot(arr_a.toarray(), arr_b.toarray(), axes=axes)

    # sparse-dense
    res = arr_a.tensordot(arr_b.toarray(), axes=axes)
    assert_equal(res, exp)
    res = arr_a.tensordot(list(arr_b.toarray()), axes=axes)
    assert_equal(res, exp)

    # sparse-sparse
    res = arr_a.tensordot(arr_b, axes=axes)
    if type(res) is coo_array:
        assert_equal(res.toarray(), exp)
    else:
        assert_equal(res, exp)


def test_tensordot_with_invalid_args():
    rng = np.random.default_rng(23409823)

    arr_a = random_array((3,4,5), density=0.6, random_state=rng, dtype=int)
    arr_b = random_array((3,4,6), density=0.6, random_state=rng, dtype=int)

    axes = ([2], [2]) # sizes of 2nd axes of both shapes do not match
    with pytest.raises(ValueError, match="sizes of the corresponding axes must match"):
        arr_a.tensordot(arr_b, axes=axes)

    arr_a = random_array((5,4,2,3,7), density=0.6, random_state=rng, dtype=int)
    arr_b = random_array((4,6,3,2), density=0.6, random_state=rng, dtype=int)

    axes = ([2,0,1], [1,3]) # lists have different lengths
    with pytest.raises(ValueError, match="axes lists/tuples must be of the"
                       " same length"):
        arr_a.tensordot(arr_b, axes=axes)


@pytest.mark.parametrize(('actual_shape', 'broadcast_shape'),
                         [((1,3,5,4), (2,3,5,4)), ((2,1,5,4), (6,2,3,5,4)),
                          ((1,1,7,8,9), (4,5,6,7,8,9)), ((1,3), (4,5,3)),
                          ((7,8,1), (7,8,5)), ((3,1), (3,4)), ((1,), (5,)),
                          ((1,1,1), (4,5,6)), ((1,3,1,5,4), (8,2,3,9,5,4)),])
def test_broadcast_to(actual_shape, broadcast_shape):
    rng = np.random.default_rng(23409823)

    arr = random_array(actual_shape, density=0.6, random_state=rng, dtype=int)
    res = arr._broadcast_to(broadcast_shape)
    exp = np.broadcast_to(arr.toarray(), broadcast_shape)
    assert_equal(res.toarray(), exp)


@pytest.mark.parametrize(('shape'), [(4,5,6,7,8), (6,4),
                                     (5,9,3,2), (9,5,2,3,4),])
def test_block_diag(shape):
    rng = np.random.default_rng(23409823)
    sp_x = random_array(shape, density=0.6, random_state=rng, dtype=int)
    den_x = sp_x.toarray()

    # converting n-d numpy array to an array of slices of 2-D matrices,
    # to pass as argument into scipy.linalg.block_diag
    num_slices = int(np.prod(den_x.shape[:-2]))
    reshaped_array = den_x.reshape((num_slices,) + den_x.shape[-2:])
    matrices = [reshaped_array[i, :, :] for i in range(num_slices)]
    exp = block_diag(*matrices)

    res = _block_diag(sp_x)

    assert_equal(res.toarray(), exp)


@pytest.mark.parametrize(('shape'), [(4,5,6,7,8), (6,4),
                                     (5,9,3,2), (9,5,2,3,4),])
def test_extract_block_diag(shape):
    rng = np.random.default_rng(23409823)
    sp_x = random_array(shape, density=0.6, random_state=rng, dtype=int)
    res = _extract_block_diag(_block_diag(sp_x), shape)

    assert_equal(res.toarray(), sp_x.toarray())
