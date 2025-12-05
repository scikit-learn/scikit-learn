# Optimized inner loop of load_svmlight_file.
#
# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

cimport cython
from libc.string cimport strchr

import numpy as np


cdef bytes COMMA = u','.encode('ascii')
cdef bytes COLON = u':'.encode('ascii')


# Small dynamic array helpers backed by NumPy arrays. These provide
# an append() method with amortized O(1) behaviour and expose the
# underlying ndarray slice via get_array(). They use typed memoryviews
# internally for fast access.


cdef class _DynDouble:
    cdef double[:] view
    cdef object arr
    cdef Py_ssize_t size, capacity

    def __cinit__(self, Py_ssize_t init=16):
        if init <= 0:
            init = 16
        self.capacity = init
        self.arr = np.empty(self.capacity, dtype=np.float64)
        self.view = self.arr
        self.size = 0

    cpdef append(self, double val):
        if self.size >= self.capacity:
            newcap = self.capacity * 2
            newarr = np.empty(newcap, dtype=np.float64)
            newarr[:self.size] = self.arr[:self.size]
            self.arr = newarr
            self.view = self.arr
            self.capacity = newcap
        self.view[self.size] = val
        self.size += 1

    cpdef object get_array(self):
        return self.arr[:self.size]


cdef class _DynFloat:
    cdef float[:] view
    cdef object arr
    cdef Py_ssize_t size, capacity

    def __cinit__(self, Py_ssize_t init=16):
        if init <= 0:
            init = 16
        self.capacity = init
        self.arr = np.empty(self.capacity, dtype=np.float32)
        self.view = self.arr
        self.size = 0

    cpdef append(self, float val):
        if self.size >= self.capacity:
            newcap = self.capacity * 2
            newarr = np.empty(newcap, dtype=np.float32)
            newarr[:self.size] = self.arr[:self.size]
            self.arr = newarr
            self.view = self.arr
            self.capacity = newcap
        self.view[self.size] = val
        self.size += 1

    cpdef object get_array(self):
        return self.arr[:self.size]


cdef class _DynLongLong:
    cdef long long[:] view
    cdef object arr
    cdef Py_ssize_t size, capacity

    def __cinit__(self, Py_ssize_t init=16):
        if init <= 0:
            init = 16
        self.capacity = init
        self.arr = np.empty(self.capacity, dtype=np.int64)
        self.view = self.arr
        self.size = 0

    cpdef append(self, long long val):
        if self.size >= self.capacity:
            newcap = self.capacity * 2
            newarr = np.empty(newcap, dtype=np.int64)
            newarr[:self.size] = self.arr[:self.size]
            self.arr = newarr
            self.view = self.arr
            self.capacity = newcap
        self.view[self.size] = val
        self.size += 1

    cpdef object get_array(self):
        return self.arr[:self.size]


def _load_svmlight_file(f, dtype, bint multilabel, bint zero_based,
                        bint query_id, long long offset, long long length):
    cdef object data, indices, indptr
    cdef bytes line
    cdef char *hash_ptr
    cdef char *line_cstr
    cdef int idx, prev_idx
    cdef Py_ssize_t i
    cdef bytes qid_prefix = b'qid'
    cdef Py_ssize_t n_features
    cdef long long offset_max = offset + length if length > 0 else -1

    # Special-case float32 but use float64 for everything else; the Python
    # code will do further conversions. Use small dynamic arrays backed by
    # NumPy arrays for performance and limited-API compatibility.
    if dtype == np.float32:
        data_dyn = _DynFloat()
        dtype = np.float32
    else:
        data_dyn = _DynDouble()
        dtype = np.float64

    indices_dyn = _DynLongLong()
    indptr_dyn = _DynLongLong()
    # start indptr with 0
    indptr_dyn.append(0)
    query = np.arange(0, dtype=np.int64)

    if multilabel:
        labels = []
    else:
        labels = _DynDouble()

    if offset > 0:
        f.seek(offset)
        # drop the current line that might be truncated and is to be
        # fetched by another call
        f.readline()

    for line in f:
        # skip comments
        line_cstr = line
        hash_ptr = strchr(line_cstr, 35)  # ASCII value of '#' is 35
        if hash_ptr != NULL:
            line = line[:hash_ptr - line_cstr]

        line_parts = line.split()
        if len(line_parts) == 0:
            continue

        target, features = line_parts[0], line_parts[1:]
        if multilabel:
            if COLON in target:
                target, features = [], line_parts[0:]
            else:
                target = [float(y) for y in target.split(COMMA)]
            target.sort()
            labels.append(tuple(target))
        else:
            # Use a dynamic array for labels for performance.
            labels.append(float(target))

        prev_idx = -1
        n_features = len(features)
        if n_features and features[0].startswith(qid_prefix):
            _, value = features[0].split(COLON, 1)
            if query_id:
                query = np.append(query, np.int64(value))
            features.pop(0)
            n_features -= 1

        for i in range(0, n_features):
            idx_s, value = features[i].split(COLON, 1)
            idx = int(idx_s)
            if idx < 0 or not zero_based and idx == 0:
                raise ValueError(
                    "Invalid index %d in SVMlight/LibSVM data file." % idx)
            if idx <= prev_idx:
                raise ValueError("Feature indices in SVMlight/LibSVM data "
                                 "file should be sorted and unique.")

            # append into the dynamic int array
            indices_dyn.append(idx)

            # append into the dynamic float array
            data_dyn.append(float(value))

            prev_idx = idx

        # increment index pointer array size
        # data_dyn.size is a cdef attribute not visible to Python-level
        # code paths; obtain the current length via the NumPy view.
        indptr_dyn.append(data_dyn.get_array().shape[0])

        if offset_max != -1 and f.tell() > offset_max:
            # Stop here and let another call deal with the following.
            break

    # Convert dynamic arrays to NumPy arrays for the return value. The
    # calling Python code expects array-like sequences; NumPy arrays are
    # efficient and consistent with the rest of scikit-learn's code.
    data = data_dyn.get_array()
    indices = indices_dyn.get_array()
    indptr = indptr_dyn.get_array()
    if not multilabel:
        labels = labels.get_array()

    return (dtype, data, indices, indptr, labels, query)


# Two fused types are defined to be able to
# use all possible combinations of parameters.
ctypedef fused int_or_float:
    cython.integral
    cython.floating
    signed long long

ctypedef fused double_or_longlong:
    double
    signed long long

ctypedef fused int_or_longlong:
    cython.integral
    signed long long


def get_dense_row_string(
    const int_or_float[:, :] X,
    Py_ssize_t[:] x_inds,
    double_or_longlong[:] x_vals,
    Py_ssize_t row,
    str value_pattern,
    bint one_based,
):
    cdef:
        Py_ssize_t row_length = X.shape[1]
        Py_ssize_t x_nz_used = 0
        Py_ssize_t k
        int_or_float val

    for k in range(row_length):
        val = X[row, k]
        if val == 0:
            continue
        x_inds[x_nz_used] = k
        x_vals[x_nz_used] = <double_or_longlong> val
        x_nz_used += 1

    reprs = [
        value_pattern % (x_inds[i] + one_based, x_vals[i])
        for i in range(x_nz_used)
    ]

    return " ".join(reprs)


def get_sparse_row_string(
    int_or_float[:] X_data,
    int[:] X_indptr,
    int[:] X_indices,
    Py_ssize_t row,
    str value_pattern,
    bint one_based,
):
    cdef:
        Py_ssize_t row_start = X_indptr[row]
        Py_ssize_t row_end = X_indptr[row+1]

    reprs = [
        value_pattern % (X_indices[i] + one_based, X_data[i])
        for i in range(row_start, row_end)
    ]

    return " ".join(reprs)


def _dump_svmlight_file(
    X,
    y,
    f,
    bint multilabel,
    bint one_based,
    int_or_longlong[:] query_id,
    bint X_is_sp,
    bint y_is_sp,
):
    cdef bint X_is_integral
    cdef bint query_id_is_not_empty = query_id.size > 0
    X_is_integral = X.dtype.kind == "i"
    if X_is_integral:
        value_pattern = "%d:%d"
    else:
        value_pattern = "%d:%.16g"
    if y.dtype.kind == "i":
        label_pattern = "%d"
    else:
        label_pattern = "%.16g"

    line_pattern = "%s"
    if query_id_is_not_empty:
        line_pattern += " qid:%d"
    line_pattern += " %s\n"

    cdef:
        Py_ssize_t num_labels = y.shape[1]
        Py_ssize_t x_len = X.shape[0]
        Py_ssize_t row_length = X.shape[1]
        Py_ssize_t i
        Py_ssize_t j
        Py_ssize_t col_start
        Py_ssize_t col_end
        Py_ssize_t[:] x_inds = np.empty(row_length, dtype=np.intp)
        signed long long[:] x_vals_int
        double[:] x_vals_float

    if not X_is_sp:
        if X_is_integral:
            x_vals_int = np.zeros(row_length, dtype=np.longlong)
        else:
            x_vals_float = np.zeros(row_length, dtype=np.float64)

    for i in range(x_len):
        if not X_is_sp:
            if X_is_integral:
                s = get_dense_row_string(X, x_inds, x_vals_int, i, value_pattern, one_based)
            else:
                s = get_dense_row_string(X, x_inds, x_vals_float, i, value_pattern, one_based)
        else:
            s = get_sparse_row_string(X.data, X.indptr, X.indices, i, value_pattern, one_based)
        if multilabel:
            if y_is_sp:
                col_start = y.indptr[i]
                col_end = y.indptr[i+1]
                labels_str = ','.join(tuple(label_pattern % y.indices[j] for j in range(col_start, col_end) if y.data[j] != 0))
            else:
                labels_str = ','.join(label_pattern % j for j in range(num_labels) if y[i, j] != 0)
        else:
            if y_is_sp:
                labels_str = label_pattern % y.data[i]
            else:
                labels_str = label_pattern % y[i, 0]

        if query_id_is_not_empty:
            feat = (labels_str, query_id[i], s)
        else:
            feat = (labels_str, s)

        f.write((line_pattern % feat).encode("utf-8"))
