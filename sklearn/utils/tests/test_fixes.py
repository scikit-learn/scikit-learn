# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import io
import tarfile

import numpy as np
import pytest
import scipy as sp

from sklearn.utils._testing import assert_array_equal
from sklearn.utils.fixes import (
    _ensure_sparse_index_int32,
    _object_dtype_isnan,
    _smallest_admissible_index_dtype,
    _validate_tar_members,
)


def _make_tar(members):
    """Create an in-memory tar archive with the given members.

    Each member is a (name, content) or (name, content, kind) tuple.
    kind is one of "file" (default), "symlink", or "hardlink".
    """
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        for item in members:
            name, content = item[0], item[1]
            kind = item[2] if len(item) > 2 else "file"
            info = tarfile.TarInfo(name=name)
            if kind == "symlink":
                info.type = tarfile.SYMTYPE
                info.linkname = content
                tf.addfile(info)
            elif kind == "hardlink":
                info.type = tarfile.LNKTYPE
                info.linkname = content
                tf.addfile(info)
            else:
                data = content if isinstance(content, bytes) else content.encode()
                info.size = len(data)
                tf.addfile(info, io.BytesIO(data))
    buf.seek(0)
    return buf


def test_validate_tar_members_safe(tmp_path):
    """A well-formed archive passes validation without error."""
    buf = _make_tar([("subdir/file.txt", b"hello")])
    with tarfile.open(fileobj=buf) as tf:
        _validate_tar_members(tf, str(tmp_path))  # should not raise


@pytest.mark.parametrize(
    "members, match",
    [
        ([("../../evil.txt", b"pwned")], "outside destination"),
        ([("link", "../../outside", "symlink")], "links outside destination"),
        ([("link", "/etc/passwd", "symlink")], "absolute link target"),
        ([("link", "../../outside", "hardlink")], "links outside destination"),
    ],
)
def test_validate_tar_members_rejects_unsafe(members, match, tmp_path):
    """Path traversal, absolute paths, and unsafe links are rejected."""
    buf = _make_tar(members)
    with tarfile.open(fileobj=buf) as tf:
        with pytest.raises(ValueError, match=match):
            _validate_tar_members(tf, str(tmp_path))


def test_validate_tar_members_special_file(tmp_path):
    """Special file types (e.g. FIFO) are rejected."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tf:
        info = tarfile.TarInfo(name="fifo")
        info.type = tarfile.FIFOTYPE
        tf.addfile(info)
    buf.seek(0)
    with tarfile.open(fileobj=buf) as tf:
        with pytest.raises(ValueError, match="special file type"):
            _validate_tar_members(tf, str(tmp_path))


@pytest.mark.parametrize("dtype, val", ([object, 1], [object, "a"], [float, 1]))
def test_object_dtype_isnan(dtype, val):
    X = np.array([[val, np.nan], [np.nan, val]], dtype=dtype)

    expected_mask = np.array([[False, True], [True, False]])

    mask = _object_dtype_isnan(X)

    assert_array_equal(mask, expected_mask)


@pytest.mark.parametrize(
    "params, expected_dtype",
    [
        ({}, np.int32),  # default behaviour
        ({"maxval": np.iinfo(np.int32).max}, np.int32),
        ({"maxval": np.iinfo(np.int32).max + 1}, np.int64),
    ],
)
def test_smallest_admissible_index_dtype_max_val(params, expected_dtype):
    """Check the behaviour of `smallest_admissible_index_dtype` depending only on the
    `max_val` parameter.
    """
    assert _smallest_admissible_index_dtype(**params) == expected_dtype


@pytest.mark.parametrize(
    "params, expected_dtype",
    [
        # Arrays dtype is int64 and thus should not be downcasted to int32 without
        # checking the content of providing maxval.
        ({"arrays": np.array([1, 2], dtype=np.int64)}, np.int64),
        # One of the array is int64 and should not be downcasted to int32
        # for the same reasons.
        (
            {
                "arrays": (
                    np.array([1, 2], dtype=np.int32),
                    np.array([1, 2], dtype=np.int64),
                )
            },
            np.int64,
        ),
        # Both arrays are already int32: we can just keep this dtype.
        (
            {
                "arrays": (
                    np.array([1, 2], dtype=np.int32),
                    np.array([1, 2], dtype=np.int32),
                )
            },
            np.int32,
        ),
        # Arrays should be upcasted to at least int32 precision.
        ({"arrays": np.array([1, 2], dtype=np.int8)}, np.int32),
        # Check that `maxval` takes precedence over the arrays and thus upcast to
        # int64.
        (
            {
                "arrays": np.array([1, 2], dtype=np.int32),
                "maxval": np.iinfo(np.int32).max + 1,
            },
            np.int64,
        ),
    ],
)
def test_smallest_admissible_index_dtype_without_checking_contents(
    params, expected_dtype
):
    """Check the behaviour of `smallest_admissible_index_dtype` using the passed
    arrays but without checking the contents of the arrays.
    """
    assert _smallest_admissible_index_dtype(**params) == expected_dtype


@pytest.mark.parametrize(
    "params, expected_dtype",
    [
        # empty arrays should always be converted to int32 indices
        (
            {
                "arrays": (np.array([], dtype=np.int64), np.array([], dtype=np.int64)),
                "check_contents": True,
            },
            np.int32,
        ),
        # arrays respecting np.iinfo(np.int32).min < x < np.iinfo(np.int32).max should
        # be converted to int32,
        (
            {"arrays": np.array([1], dtype=np.int64), "check_contents": True},
            np.int32,
        ),
        # otherwise, it should be converted to int64. We need to create a uint32
        # arrays to accommodate a value > np.iinfo(np.int32).max
        (
            {
                "arrays": np.array([np.iinfo(np.int32).max + 1], dtype=np.uint32),
                "check_contents": True,
            },
            np.int64,
        ),
        # maxval should take precedence over the arrays contents and thus upcast to
        # int64.
        (
            {
                "arrays": np.array([1], dtype=np.int32),
                "check_contents": True,
                "maxval": np.iinfo(np.int32).max + 1,
            },
            np.int64,
        ),
        # when maxval is small, but check_contents is True and the contents
        # require np.int64, we still require np.int64 indexing in the end.
        (
            {
                "arrays": np.array([np.iinfo(np.int32).max + 1], dtype=np.uint32),
                "check_contents": True,
                "maxval": 1,
            },
            np.int64,
        ),
    ],
)
def test_smallest_admissible_index_dtype_by_checking_contents(params, expected_dtype):
    """Check the behaviour of `smallest_admissible_index_dtype` using the dtype of the
    arrays but as well the contents.
    """
    assert _smallest_admissible_index_dtype(**params) == expected_dtype


@pytest.mark.parametrize(
    "params, err_type, err_msg",
    [
        (
            {"maxval": np.iinfo(np.int64).max + 1},
            ValueError,
            "is to large to be represented as np.int64",
        ),
        (
            {"arrays": np.array([1, 2], dtype=np.float64)},
            ValueError,
            "Array dtype float64 is not supported",
        ),
        ({"arrays": [1, 2]}, TypeError, "Arrays should be of type np.ndarray"),
    ],
)
def test_smallest_admissible_index_dtype_error(params, err_type, err_msg):
    """Check that we raise the proper error message."""
    with pytest.raises(err_type, match=err_msg):
        _smallest_admissible_index_dtype(**params)


INDEX_CONSTRUCTORS = [
    sp.sparse.csc_array,
    sp.sparse.csr_array,
    sp.sparse.coo_array,
    sp.sparse.csc_matrix,
    sp.sparse.csr_matrix,
    sp.sparse.coo_matrix,
]
NO_INDEX_TEST_CONSTRUCTORS = [
    sp.sparse.bsr_array,
    sp.sparse.bsr_matrix,
    sp.sparse.dia_array,
    sp.sparse.dok_array,
    sp.sparse.lil_array,
    sp.sparse.dia_matrix,
    sp.sparse.dok_matrix,
    sp.sparse.lil_matrix,
]
SPARSE_CONSTRUCTORS = INDEX_CONSTRUCTORS + NO_INDEX_TEST_CONSTRUCTORS


@pytest.mark.parametrize("constructor", SPARSE_CONSTRUCTORS)
def test_ensure_sparse_index_int32(constructor):
    A = constructor(np.array([[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]))
    _ensure_sparse_index_int32(A)


@pytest.mark.parametrize("constructor", INDEX_CONSTRUCTORS)
def test_ensure_int32_raises(constructor):
    with pytest.raises(ValueError, match="too large"):
        rows, cols = [2, 0], [1, np.iinfo(np.int32).max + 1]
        if "csc" in constructor.__name__:
            rows, cols = cols, rows
        A = sp.sparse.coo_array(([1.0, 2.0], (rows, cols)))
        _ensure_sparse_index_int32(constructor(A))
