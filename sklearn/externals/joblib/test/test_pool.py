import os
import shutil
import tempfile

from nose import SkipTest
from nose.tools import with_setup
from nose.tools import assert_equal
from nose.tools import assert_raises
from nose.tools import assert_false
from nose.tools import assert_true
from .common import with_numpy, np
from .common import setup_autokill
from .common import teardown_autokill

from .._multiprocessing_helpers import mp
if mp is not None:
    from ..pool import MemmapingPool
    from ..pool import has_shareable_memory
    from ..pool import ArrayMemmapReducer
    from ..pool import reduce_memmap


TEMP_FOLDER = None


def setup_module():
    setup_autokill(__name__, timeout=30)


def teardown_module():
    teardown_autokill(__name__)


def check_multiprocessing():
    if mp is None:
        raise SkipTest('Need multiprocessing to run')


with_multiprocessing = with_setup(check_multiprocessing)


def setup_temp_folder():
    global TEMP_FOLDER
    TEMP_FOLDER = tempfile.mkdtemp(prefix='joblib_test_pool_')


def teardown_temp_folder():
    global TEMP_FOLDER
    if TEMP_FOLDER is not None:
        shutil.rmtree(TEMP_FOLDER)
        TEMP_FOLDER = None


with_temp_folder = with_setup(setup_temp_folder, teardown_temp_folder)


def setup_if_has_dev_shm():
    if not os.path.exists('/dev/shm'):
        raise SkipTest("This test requires the /dev/shm shared memory fs.")


with_dev_shm = with_setup(setup_if_has_dev_shm)


def check_array(args):
    """Dummy helper function to be executed in subprocesses

    Check that the provided array has the expected values in the provided
    range.

    """
    assert_array_equal = np.testing.assert_array_equal
    data, position, expected = args
    assert_equal(data[position], expected)


def inplace_double(args):
    """Dummy helper function to be executed in subprocesses


    Check that the input array has the right values in the provided range
    and perform an inplace modification to double the values in the range by
    two.

    """
    assert_array_equal = np.testing.assert_array_equal
    data, position, expected = args
    assert_equal(data[position], expected)
    data[position] *= 2
    assert_equal(data[position], 2 * expected)


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_memmap_based_array_reducing():
    """Check that it is possible to reduce a memmap backed array"""
    assert_array_equal = np.testing.assert_array_equal
    filename = os.path.join(TEMP_FOLDER, 'test.mmap')

    # Create a file larger than what will be used by a
    buffer = np.memmap(filename, dtype=np.float64, shape=500, mode='w+')

    # Fill the original buffer with negative markers to detect over of
    # underflow in case of test failures
    buffer[:] = - 1.0 * np.arange(buffer.shape[0], dtype=buffer.dtype)
    buffer.flush()

    # Memmap a 2D fortran array on a offseted subsection of the previous
    # buffer
    a = np.memmap(filename, dtype=np.float64, shape=(3, 5, 4),
                  mode='r+', order='F', offset=4)
    a[:] = np.arange(60).reshape(a.shape)

    # Build various views that share the buffer with the original memmap

    # b is an memmap sliced view on an memmap instance
    b = a[1:-1, 2:-1, 2:4]

    # c and d are array views
    c = np.asarray(b)
    d = c.T

    # Array reducer with auto dumping disabled
    reducer = ArrayMemmapReducer(None, TEMP_FOLDER, 'c')

    def reconstruct_array(x):
        cons, args = reducer(x)
        return cons(*args)

    def reconstruct_memmap(x):
        cons, args = reduce_memmap(x)
        return cons(*args)

    # Reconstruct original memmap
    a_reconstructed = reconstruct_memmap(a)
    assert_true(has_shareable_memory(a_reconstructed))
    assert_true(isinstance(a_reconstructed, np.memmap))
    assert_array_equal(a_reconstructed, a)

    # Reconstruct strided memmap view
    b_reconstructed = reconstruct_memmap(b)
    assert_true(has_shareable_memory(b_reconstructed))
    assert_array_equal(b_reconstructed, b)

    # Reconstruct arrays views on memmap base
    c_reconstructed = reconstruct_array(c)
    assert_false(isinstance(c_reconstructed, np.memmap))
    assert_true(has_shareable_memory(c_reconstructed))
    assert_array_equal(c_reconstructed, c)

    d_reconstructed = reconstruct_array(d)
    assert_false(isinstance(d_reconstructed, np.memmap))
    assert_true(has_shareable_memory(d_reconstructed))
    assert_array_equal(d_reconstructed, d)

    # Test graceful degradation on fake memmap instances with in-memory
    # buffers
    a3 = a * 3
    assert_false(has_shareable_memory(a3))
    a3_reconstructed = reconstruct_memmap(a3)
    assert_false(has_shareable_memory(a3_reconstructed))
    assert_false(isinstance(a3_reconstructed, np.memmap))
    assert_array_equal(a3_reconstructed, a * 3)

    # Test graceful degradation on arrays derived from fake memmap instances
    b3 = np.asarray(a3)
    assert_false(has_shareable_memory(b3))

    b3_reconstructed = reconstruct_array(b3)
    assert_true(isinstance(b3_reconstructed, np.ndarray))
    assert_false(has_shareable_memory(b3_reconstructed))
    assert_array_equal(b3_reconstructed, b3)


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_high_dimension_memmap_array_reducing():
    assert_array_equal = np.testing.assert_array_equal

    filename = os.path.join(TEMP_FOLDER, 'test.mmap')

    # Create a high dimensional memmap
    a = np.memmap(filename, dtype=np.float64, shape=(100, 15, 15, 3), mode='w+')
    a[:] = np.arange(100 * 15 * 15 * 3).reshape(a.shape)

    # Create some slices/indices at various dimensions
    b = a[0:10]
    c = a[:, 5:10]
    d = a[:, :, :, 0]
    e = a[1:3:4]

    def reconstruct_memmap(x):
        cons, args = reduce_memmap(x)
        res = cons(*args)
        return res

    a_reconstructed = reconstruct_memmap(a)
    assert_true(has_shareable_memory(a_reconstructed))
    assert_true(isinstance(a_reconstructed, np.memmap))
    assert_array_equal(a_reconstructed, a)

    b_reconstructed = reconstruct_memmap(b)
    assert_true(has_shareable_memory(b_reconstructed))
    assert_array_equal(b_reconstructed, b)

    c_reconstructed = reconstruct_memmap(c)
    assert_true(has_shareable_memory(c_reconstructed))
    assert_array_equal(c_reconstructed, c)

    d_reconstructed = reconstruct_memmap(d)
    assert_true(has_shareable_memory(d_reconstructed))
    assert_array_equal(d_reconstructed, d)

    e_reconstructed = reconstruct_memmap(e)
    assert_true(has_shareable_memory(e_reconstructed))
    assert_array_equal(e_reconstructed, e)


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_pool_with_memmap():
    """Check that subprocess can access and update shared memory memmap"""
    assert_array_equal = np.testing.assert_array_equal

    # Fork the subprocess before allocating the objects to be passed
    pool_temp_folder = os.path.join(TEMP_FOLDER, 'pool')
    os.makedirs(pool_temp_folder)
    p = MemmapingPool(10, max_nbytes=2, temp_folder=pool_temp_folder)
    try:
        filename = os.path.join(TEMP_FOLDER, 'test.mmap')
        a = np.memmap(filename, dtype=np.float32, shape=(3, 5), mode='w+')
        a.fill(1.0)

        p.map(inplace_double, [(a, (i, j), 1.0)
                               for i in range(a.shape[0])
                               for j in range(a.shape[1])])

        assert_array_equal(a, 2 * np.ones(a.shape))

        # Open a copy-on-write view on the previous data
        b = np.memmap(filename, dtype=np.float32, shape=(5, 3), mode='c')

        p.map(inplace_double, [(b, (i, j), 2.0)
                               for i in range(b.shape[0])
                               for j in range(b.shape[1])])

        # Passing memmap instances to the pool should not trigger the creation
        # of new files on the FS
        assert_equal(os.listdir(pool_temp_folder), [])

        # the original data is untouched
        assert_array_equal(a, 2 * np.ones(a.shape))
        assert_array_equal(b, 2 * np.ones(b.shape))

        # readonly maps can be read but not updated
        c = np.memmap(filename, dtype=np.float32, shape=(10,), mode='r',
                      offset=5 * 4)

        assert_raises(AssertionError, p.map, check_array,
                      [(c, i, 3.0) for i in range(c.shape[0])])

        # depending on the version of numpy one can either get a RuntimeError
        # or a ValueError
        assert_raises((RuntimeError, ValueError), p.map, inplace_double,
                      [(c, i, 2.0) for i in range(c.shape[0])])
    finally:
        # Clean all filehandlers held by the pool
        p.terminate()
        del p


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_pool_with_memmap_array_view():
    """Check that subprocess can access and update shared memory array"""
    assert_array_equal = np.testing.assert_array_equal

    # Fork the subprocess before allocating the objects to be passed
    pool_temp_folder = os.path.join(TEMP_FOLDER, 'pool')
    os.makedirs(pool_temp_folder)
    p = MemmapingPool(10, max_nbytes=2, temp_folder=pool_temp_folder)
    try:

        filename = os.path.join(TEMP_FOLDER, 'test.mmap')
        a = np.memmap(filename, dtype=np.float32, shape=(3, 5), mode='w+')
        a.fill(1.0)

        # Create an ndarray view on the memmap instance
        a_view = np.asarray(a)
        assert_false(isinstance(a_view, np.memmap))
        assert_true(has_shareable_memory(a_view))

        p.map(inplace_double, [(a_view, (i, j), 1.0)
                               for i in range(a.shape[0])
                               for j in range(a.shape[1])])

        # Both a and the a_view have been updated
        assert_array_equal(a, 2 * np.ones(a.shape))
        assert_array_equal(a_view, 2 * np.ones(a.shape))

        # Passing memmap array view to the pool should not trigger the
        # creation of new files on the FS
        assert_equal(os.listdir(pool_temp_folder), [])

    finally:
        p.terminate()
        del p


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_memmaping_pool_for_large_arrays():
    """Check that large arrays are not copied in memory"""
    assert_array_equal = np.testing.assert_array_equal

    # Check that the tempfolder is empty
    assert_equal(os.listdir(TEMP_FOLDER), [])

    # Build an array reducers that automaticaly dump large array content
    # to filesystem backed memmap instances to avoid memory explosion
    p = MemmapingPool(3, max_nbytes=40, temp_folder=TEMP_FOLDER)
    try:
        # The tempory folder for the pool is not provisioned in advance
        assert_equal(os.listdir(TEMP_FOLDER), [])
        assert_false(os.path.exists(p._temp_folder))

        small = np.ones(5, dtype=np.float32)
        assert_equal(small.nbytes, 20)
        p.map(check_array, [(small, i, 1.0) for i in range(small.shape[0])])

        # Memory has been copied, the pool filesystem folder is unused
        assert_equal(os.listdir(TEMP_FOLDER), [])

        # Try with a file larger than the memmap threshold of 40 bytes
        large = np.ones(100, dtype=np.float64)
        assert_equal(large.nbytes, 800)
        p.map(check_array, [(large, i, 1.0) for i in range(large.shape[0])])

        # The data has been dumped in a temp folder for subprocess to share it
        # without per-child memory copies
        assert_true(os.path.isdir(p._temp_folder))
        dumped_filenames = os.listdir(p._temp_folder)
        assert_equal(len(dumped_filenames), 2)

    finally:
        # check FS garbage upon pool termination
        p.terminate()
        assert_false(os.path.exists(p._temp_folder))
        del p


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_memmaping_pool_for_large_arrays_disabled():
    """Check that large arrays memmaping can be disabled"""
    # Set max_nbytes to None to disable the auto memmaping feature
    p = MemmapingPool(3, max_nbytes=None, temp_folder=TEMP_FOLDER)
    try:

        # Check that the tempfolder is empty
        assert_equal(os.listdir(TEMP_FOLDER), [])

        # Try with a file largish than the memmap threshold of 40 bytes
        large = np.ones(100, dtype=np.float64)
        assert_equal(large.nbytes, 800)
        p.map(check_array, [(large, i, 1.0) for i in range(large.shape[0])])

        # Check that the tempfolder is still empty
        assert_equal(os.listdir(TEMP_FOLDER), [])

    finally:
        # Cleanup open file descriptors
        p.terminate()
        del p


@with_numpy
@with_multiprocessing
@with_dev_shm
def test_memmaping_on_dev_shm():
    """Check that large arrays memmaping can be disabled"""
    p = MemmapingPool(3, max_nbytes=10)
    try:
        # Check that the pool has correctly detected the presence of the
        # shared memory filesystem.
        pool_temp_folder = p._temp_folder
        folder_prefix = '/dev/shm/joblib_memmaping_pool_'
        assert_true(pool_temp_folder.startswith(folder_prefix))
        assert_true(os.path.exists(pool_temp_folder))

        # Try with a file larger than the memmap threshold of 10 bytes
        a = np.ones(100, dtype=np.float64)
        assert_equal(a.nbytes, 800)
        p.map(id, [a] * 10)
        # a should have been memmaped to the pool temp folder: the joblib
        # pickling procedure generate a .pkl and a .npy file:
        assert_equal(len(os.listdir(pool_temp_folder)), 2)

        b = np.ones(100, dtype=np.float64)
        assert_equal(b.nbytes, 800)
        p.map(id, [b] * 10)
        # A copy of both a and b are not stored in the shared memory folder
        assert_equal(len(os.listdir(pool_temp_folder)), 4)

    finally:
        # Cleanup open file descriptors
        p.terminate()
        del p

    # The temp folder is cleaned up upon pool termination
    assert_false(os.path.exists(pool_temp_folder))


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_memmaping_pool_for_large_arrays_in_return():
    """Check that large arrays are not copied in memory in return"""
    assert_array_equal = np.testing.assert_array_equal

    # Build an array reducers that automaticaly dump large array content
    # but check that the returned datastructure are regular arrays to avoid
    # passing a memmap array pointing to a pool controlled temp folder that
    # might be confusing to the user

    # The MemmapingPool user can always return numpy.memmap object explicitly
    # to avoid memory copy
    p = MemmapingPool(3, max_nbytes=10, temp_folder=TEMP_FOLDER)
    try:
        res = p.apply_async(np.ones, args=(1000,))
        large = res.get()
        assert_false(has_shareable_memory(large))
        assert_array_equal(large, np.ones(1000))
    finally:
        p.terminate()
        del p


def _worker_multiply(a, n_times):
    """Multiplication function to be executed by subprocess"""
    assert_true(has_shareable_memory(a))
    return a * n_times


@with_numpy
@with_multiprocessing
@with_temp_folder
def test_workaround_against_bad_memmap_with_copied_buffers():
    """Check that memmaps with a bad buffer are returned as regular arrays

    Unary operations and ufuncs on memmap instances return a new memmap
    instance with an in-memory buffer (probably a numpy bug).
    """
    assert_array_equal = np.testing.assert_array_equal

    p = MemmapingPool(3, max_nbytes=10, temp_folder=TEMP_FOLDER)
    try:
        # Send a complex, large-ish view on a array that will be converted to
        # a memmap in the worker process
        a = np.asarray(np.arange(6000).reshape((1000, 2, 3)),
                       order='F')[:, :1, :]

        # Call a non-inplace multiply operation on the worker and memmap and
        # send it back to the parent.
        b = p.apply_async(_worker_multiply, args=(a, 3)).get()
        assert_false(has_shareable_memory(b))
        assert_array_equal(b, 3 * a)
    finally:
        p.terminate()
        del p
