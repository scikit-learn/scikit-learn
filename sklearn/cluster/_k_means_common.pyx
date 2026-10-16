# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
from cython cimport floating
from cython.parallel cimport parallel, prange
from libc.math cimport sqrt
from libc.stdlib cimport calloc, free, malloc

from sklearn.utils._cython_blas cimport _gemm
from sklearn.utils._cython_blas cimport RowMajor, Trans, NoTrans
from sklearn.utils._typedefs cimport intp_t
from sklearn.utils.extmath import row_norms


# Number of samples per data chunk defined as a global constant.
CHUNK_SIZE = 256


cdef floating _euclidean_dense_dense(
        const floating* a,  # IN
        const floating* b,  # IN
        int n_features,
        bint squared
) noexcept nogil:
    """Euclidean distance between a dense and b dense"""
    cdef:
        int i
        int n = n_features // 4
        int rem = n_features % 4
        floating result = 0

    # We manually unroll the loop for better cache optimization.
    for _ in range(n):
        result += (
            (a[0] - b[0]) * (a[0] - b[0]) +
            (a[1] - b[1]) * (a[1] - b[1]) +
            (a[2] - b[2]) * (a[2] - b[2]) +
            (a[3] - b[3]) * (a[3] - b[3])
        )
        a += 4
        b += 4

    for i in range(rem):
        result += (a[i] - b[i]) * (a[i] - b[i])

    return result if squared else sqrt(result)


def _euclidean_dense_dense_wrapper(
    const floating[::1] a,
    const floating[::1] b,
    bint squared
):
    """Wrapper of _euclidean_dense_dense for testing purpose"""
    return _euclidean_dense_dense(&a[0], &b[0], a.shape[0], squared)


cdef floating _euclidean_sparse_dense(
        const floating[::1] a_data,  # IN
        const int[::1] a_indices,    # IN
        const floating[::1] b,       # IN
        floating b_squared_norm,
        bint squared
) noexcept nogil:
    """Euclidean distance between a sparse and b dense"""
    cdef:
        int nnz = a_indices.shape[0]
        int i
        floating tmp, bi
        floating result = 0.0

    for i in range(nnz):
        bi = b[a_indices[i]]
        tmp = a_data[i] - bi
        result += tmp * tmp - bi * bi

    result += b_squared_norm

    if result < 0:
        result = 0.0

    return result if squared else sqrt(result)


def _euclidean_sparse_dense_wrapper(
        const floating[::1] a_data,
        const int[::1] a_indices,
        const floating[::1] b,
        floating b_squared_norm,
        bint squared
):
    """Wrapper of _euclidean_sparse_dense for testing purpose"""
    return _euclidean_sparse_dense(
        a_data, a_indices, b, b_squared_norm, squared)


cpdef floating _inertia_dense(
        const floating[:, ::1] X,           # IN
        const floating[::1] sample_weight,  # IN
        const floating[:, ::1] centers,     # IN
        const int[::1] labels,              # IN
        int n_threads,
        int single_label=-1,
):
    """Compute inertia for dense input data

    Sum of squared distance between each sample and its assigned center.

    If single_label is >= 0, the inertia is computed only for that label.
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int i, j

        floating sq_dist = 0.0
        floating inertia = 0.0

    for i in prange(n_samples, nogil=True, num_threads=n_threads,
                    schedule='static'):
        j = labels[i]
        if single_label < 0 or single_label == j:
            sq_dist = _euclidean_dense_dense(&X[i, 0], &centers[j, 0],
                                             n_features, True)
            inertia += sq_dist * sample_weight[i]

    return inertia


cpdef floating _inertia_sparse(
        X,                                  # IN
        const floating[::1] sample_weight,  # IN
        const floating[:, ::1] centers,     # IN
        const int[::1] labels,              # IN
        int n_threads,
        int single_label=-1,
):
    """Compute inertia for sparse input data

    Sum of squared distance between each sample and its assigned center.

    If single_label is >= 0, the inertia is computed only for that label.
    """
    cdef:
        floating[::1] X_data = X.data
        int[::1] X_indices = X.indices
        int[::1] X_indptr = X.indptr

        int n_samples = X.shape[0]
        int i, j

        floating sq_dist = 0.0
        floating inertia = 0.0

        floating[::1] centers_squared_norms = row_norms(centers, squared=True)

    for i in prange(n_samples, nogil=True, num_threads=n_threads,
                    schedule='static'):
        j = labels[i]
        if single_label < 0 or single_label == j:
            sq_dist = _euclidean_sparse_dense(
                X_data[X_indptr[i]: X_indptr[i + 1]],
                X_indices[X_indptr[i]: X_indptr[i + 1]],
                centers[j], centers_squared_norms[j], True)
            inertia += sq_dist * sample_weight[i]

    return inertia


cdef inline floating _recompute_threshold(floating dummy) noexcept nogil:
    """Relative threshold below which a squared distance is recomputed directly.

    ||x||² - 2 x.c + ||c||² has an absolute error of a few eps * (||x||² + ||c||²),
    which is large relative to the distance when x is close to c, in particular in
    float32. Distances below `threshold * (||x||² + ||c||²)` are recomputed, which
    bounds the relative error by about eps / threshold. In particular, a sample
    identical to an existing center gets a zero distance and is never sampled
    again. The threshold is lower in float64 because recomputing is costly when the
    data is far from the origin, where most distances are small relative to the
    norms.
    """
    if floating is float:
        return 1e-2
    else:
        return 1e-8


def _kmeans_plusplus_trials_dense(
        const floating[:, ::1] X,                       # IN
        const floating[::1] x_squared_norms,            # IN
        const floating[::1] sample_weight,              # IN
        const floating[:, ::1] candidates,              # IN
        const floating[::1] candidates_squared_norms,   # IN
        const floating[::1] closest_dist_sq,            # IN
        floating[:, ::1] trials_dist_sq,                # OUT
        double[:, ::1] trials_chunk_pot,                # OUT
        int n_threads):
    """Evaluate k-means++ candidate centers for dense input data.

    For each candidate `j` and sample `i`, compute the squared distance from
    sample `i` to its closest center if candidate `j` was added, i.e.
    `trials_dist_sq[j, i] = min(closest_dist_sq[i], ||X[i] - candidates[j]||²)`.

    `trials_chunk_pot[j, chunk_idx]` is the sum of these distances weighted by
    `sample_weight` over the chunk of `CHUNK_SIZE` samples `chunk_idx`. This is
    the parallel part of the sum over all the samples: each chunk is summed by a
    single thread, without synchronization, and the caller finishes the sum over
    the chunks (only `n_samples / CHUNK_SIZE` values per candidate). The chunks
    have a fixed size, so the result doesn't depend on the number of threads.
    The caller uses these sums to:

    - pick the candidate with the smallest potential `trials_chunk_pot[j].sum()`,
    - sample the next candidates without a cumulative sum over all the samples
      (see `_kmeans_plusplus_sample_candidates`).
    """
    cdef:
        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int n_trials = candidates.shape[0]
        int n_chunks = trials_chunk_pot.shape[1]
        int chunk_size = CHUNK_SIZE
        int ld_trials = trials_dist_sq.strides[0] // sizeof(floating)
        # Each thread processes a contiguous block of chunks with a single BLAS
        # call: many small concurrent BLAS calls scale badly with the number of
        # threads (at least with OpenBLAS).
        int n_blocks = min(n_threads, n_chunks)
        int block_idx, first_chunk, last_chunk, chunk_idx, start, end, i, j
        floating dist, norms
        floating threshold = _recompute_threshold(<floating> 0)
        double pot

    for block_idx in prange(n_blocks, nogil=True, num_threads=n_blocks,
                            schedule='static'):
        # Samples X[start:end] of the chunks first_chunk to last_chunk - 1
        first_chunk = block_idx * n_chunks // n_blocks
        last_chunk = (block_idx + 1) * n_chunks // n_blocks
        start = first_chunk * chunk_size
        end = min(last_chunk * chunk_size, n_samples)

        # trials_dist_sq[:, start:end] = -2 candidates @ X[start:end].T
        _gemm(RowMajor, NoTrans, Trans, n_trials, end - start, n_features,
              -2.0, &candidates[0, 0], n_features, &X[start, 0], n_features,
              0.0, &trials_dist_sq[0, start], ld_trials)

        for chunk_idx in range(first_chunk, last_chunk):
            start = chunk_idx * chunk_size
            end = min(start + chunk_size, n_samples)
            for j in range(n_trials):
                # pot = sum over the chunk of sample_weight[i] * trials_dist_sq[j, i]
                pot = 0.0
                for i in range(start, end):
                    # dist = ||x||² - 2 x.c + ||c||²
                    norms = x_squared_norms[i] + candidates_squared_norms[j]
                    dist = trials_dist_sq[j, i] + norms
                    if dist < threshold * norms:
                        # The sum of squared differences has no cancellation: its
                        # relative error is a few eps, whatever the norms.
                        dist = _euclidean_dense_dense(
                            &X[i, 0], &candidates[j, 0], n_features, True)
                    # trials_dist_sq[j, i] = min(closest_dist_sq[i], dist)
                    if closest_dist_sq[i] < dist:
                        dist = closest_dist_sq[i]
                    trials_dist_sq[j, i] = dist
                    pot = pot + sample_weight[i] * dist
                trials_chunk_pot[j, chunk_idx] = pot


def _kmeans_plusplus_trials_sparse(
        X,                                              # IN
        const floating[::1] x_squared_norms,            # IN
        const floating[::1] sample_weight,              # IN
        const floating[:, ::1] candidates,              # IN
        const floating[::1] candidates_squared_norms,   # IN
        const floating[::1] closest_dist_sq,            # IN
        floating[:, ::1] trials_dist_sq,                # OUT
        double[:, ::1] trials_chunk_pot,                # OUT
        int n_threads):
    """Evaluate k-means++ candidate centers for sparse input data.

    See `_kmeans_plusplus_trials_dense` for details.
    """
    cdef:
        const floating[::1] X_data = X.data
        const int[::1] X_indices = X.indices
        const int[::1] X_indptr = X.indptr
        # candidates_t[k, j] = candidates[j, k], so that the inner loop over the
        # candidates reads contiguous memory
        const floating[:, ::1] candidates_t = np.ascontiguousarray(
            np.asarray(candidates).T)

        int n_samples = X.shape[0]
        int n_features = X.shape[1]
        int n_trials = candidates.shape[0]
        int n_chunks = trials_chunk_pot.shape[1]
        int chunk_size = CHUNK_SIZE
        int chunk_idx, start, end, i, j, k
        floating dist, norms, x_k
        floating threshold = _recompute_threshold(<floating> 0)
        # Dense copy of a sample, for the direct computation of small distances
        floating *x_dense
        # Dot products of a sample with the candidates
        floating *dot

    with nogil, parallel(num_threads=min(n_threads, n_chunks)):
        x_dense = <floating*> calloc(n_features, sizeof(floating))
        dot = <floating*> malloc(n_trials * sizeof(floating))

        for chunk_idx in prange(n_chunks, schedule='static'):
            start = chunk_idx * chunk_size
            end = min(start + chunk_size, n_samples)
            for j in range(n_trials):
                trials_chunk_pot[j, chunk_idx] = 0.0
            for i in range(start, end):
                # dot[j] = x.c_j, summed over the non-zeros of x. The loop over the
                # candidates is the inner one: the n_trials sums are independent,
                # while a single sum per candidate would be one long chain of
                # dependent additions (much slower, especially on one thread).
                for j in range(n_trials):
                    dot[j] = 0.0
                for k in range(X_indptr[i], X_indptr[i + 1]):
                    x_k = X_data[k]
                    for j in range(n_trials):
                        dot[j] += x_k * candidates_t[X_indices[k], j]
                for j in range(n_trials):
                    # dist = ||x||² - 2 x.c + ||c||²
                    norms = x_squared_norms[i] + candidates_squared_norms[j]
                    dist = norms - 2 * dot[j]
                    if dist < threshold * norms:
                        # Same direct computation as for dense input, with
                        # x_dense = X[i].toarray()
                        for k in range(X_indptr[i], X_indptr[i + 1]):
                            x_dense[X_indices[k]] = X_data[k]
                        dist = _euclidean_dense_dense(
                            x_dense, &candidates[j, 0], n_features, True)
                        # Reset x_dense to zeros for the next sample
                        for k in range(X_indptr[i], X_indptr[i + 1]):
                            x_dense[X_indices[k]] = 0.0
                    # trials_dist_sq[j, i] = min(closest_dist_sq[i], dist)
                    if closest_dist_sq[i] < dist:
                        dist = closest_dist_sq[i]
                    trials_dist_sq[j, i] = dist
                    trials_chunk_pot[j, chunk_idx] = (
                        trials_chunk_pot[j, chunk_idx] + sample_weight[i] * dist)

        free(x_dense)
        free(dot)


def _kmeans_plusplus_sample_candidates(
        const double[::1] rand_vals,            # IN
        const double[::1] chunk_pot,            # IN
        const floating[::1] closest_dist_sq,    # IN
        const floating[::1] sample_weight,      # IN
        intp_t[::1] candidate_ids):             # OUT
    """Sample k-means++ candidate centers.

    The probability of a sample is proportional to its weighted squared
    distance to the closest center. For each value of `rand_vals`, drawn
    uniformly in [0, 1), a chunk of `CHUNK_SIZE` samples is first sampled from
    the potentials per chunk, then a sample within that chunk.
    """
    cdef:
        const double[::1] chunk_pot_cumsum = np.cumsum(chunk_pot)
        int n_samples = closest_dist_sq.shape[0]
        int n_chunks = chunk_pot.shape[0]
        int n_candidates = rand_vals.shape[0]
        int chunk_size = CHUNK_SIZE
        int j, chunk_idx, lo, hi, start, end, i
        intp_t candidate
        double rand_val, cumsum, pot

    with nogil:
        for j in range(n_candidates):
            rand_val = rand_vals[j] * chunk_pot_cumsum[n_chunks - 1]

            # First chunk whose cumulative potential reaches rand_val
            lo, hi = 0, n_chunks - 1
            while lo < hi:
                chunk_idx = (lo + hi) // 2
                if chunk_pot_cumsum[chunk_idx] < rand_val:
                    lo = chunk_idx + 1
                else:
                    hi = chunk_idx
            if lo > 0:
                rand_val -= chunk_pot_cumsum[lo - 1]

            # First sample of the chunk with a positive potential whose cumulative
            # potential reaches rand_val. Rounding errors can make rand_val larger
            # than the potential of the chunk: then take its last sample with a
            # positive potential.
            start = lo * chunk_size
            end = min(start + chunk_size, n_samples)
            candidate = start
            cumsum = 0.0
            for i in range(start, end):
                pot = sample_weight[i] * closest_dist_sq[i]
                if pot > 0:
                    candidate = i
                    cumsum += pot
                    if cumsum >= rand_val:
                        break
            candidate_ids[j] = candidate


cpdef void _relocate_empty_clusters_dense(
        const floating[:, ::1] X,            # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # INOUT
        floating[::1] weight_in_clusters,    # INOUT
        const int[::1] labels                # IN
):
    """Relocate centers which have no sample assigned to them."""
    cdef:
        int[::1] empty_clusters = np.where(np.equal(weight_in_clusters, 0))[0].astype(np.int32)
        int n_empty = empty_clusters.shape[0]

    if n_empty == 0:
        return

    cdef:
        int n_features = X.shape[1]

        floating[::1] distances = ((np.asarray(X) - np.asarray(centers_old)[labels])**2).sum(axis=1)
        int[::1] far_from_centers = np.argpartition(distances, -n_empty)[:-n_empty-1:-1].astype(np.int32)

        int new_cluster_id, old_cluster_id, far_idx, idx, k
        floating weight

    if np.max(distances) == 0:
        # Happens when there are more clusters than non-duplicate samples. Relocating
        # is pointless in this case.
        return

    for idx in range(n_empty):

        new_cluster_id = empty_clusters[idx]

        far_idx = far_from_centers[idx]
        weight = sample_weight[far_idx]

        old_cluster_id = labels[far_idx]

        for k in range(n_features):
            centers_new[old_cluster_id, k] -= X[far_idx, k] * weight
            centers_new[new_cluster_id, k] = X[far_idx, k] * weight

        weight_in_clusters[new_cluster_id] = weight
        weight_in_clusters[old_cluster_id] -= weight


cpdef void _relocate_empty_clusters_sparse(
        const floating[::1] X_data,          # IN
        const int[::1] X_indices,            # IN
        const int[::1] X_indptr,             # IN
        const floating[::1] sample_weight,   # IN
        const floating[:, ::1] centers_old,  # IN
        floating[:, ::1] centers_new,        # INOUT
        floating[::1] weight_in_clusters,    # INOUT
        const int[::1] labels                # IN
):
    """Relocate centers which have no sample assigned to them."""
    cdef:
        int[::1] empty_clusters = np.where(np.equal(weight_in_clusters, 0))[0].astype(np.int32)
        int n_empty = empty_clusters.shape[0]

    if n_empty == 0:
        return

    cdef:
        int n_samples = X_indptr.shape[0] - 1
        int i, j, k

        floating[::1] distances = np.zeros(n_samples, dtype=X_data.base.dtype)
        floating[::1] centers_squared_norms = row_norms(centers_old, squared=True)

    for i in range(n_samples):
        j = labels[i]
        distances[i] = _euclidean_sparse_dense(
            X_data[X_indptr[i]: X_indptr[i + 1]],
            X_indices[X_indptr[i]: X_indptr[i + 1]],
            centers_old[j], centers_squared_norms[j], True)

    if np.max(distances) == 0:
        # Happens when there are more clusters than non-duplicate samples. Relocating
        # is pointless in this case.
        return

    cdef:
        int[::1] far_from_centers = np.argpartition(distances, -n_empty)[:-n_empty-1:-1].astype(np.int32)

        int new_cluster_id, old_cluster_id, far_idx, idx
        floating weight

    for idx in range(n_empty):

        new_cluster_id = empty_clusters[idx]

        far_idx = far_from_centers[idx]
        weight = sample_weight[far_idx]

        old_cluster_id = labels[far_idx]

        for k in range(X_indptr[far_idx], X_indptr[far_idx + 1]):
            centers_new[old_cluster_id, X_indices[k]] -= X_data[k] * weight
            centers_new[new_cluster_id, X_indices[k]] = X_data[k] * weight

        weight_in_clusters[new_cluster_id] = weight
        weight_in_clusters[old_cluster_id] -= weight


cdef void _average_centers(
        floating[:, ::1] centers,               # INOUT
        const floating[::1] weight_in_clusters  # IN
):
    """Average new centers wrt weights."""
    cdef:
        int n_clusters = centers.shape[0]
        int n_features = centers.shape[1]
        int j, k
        floating alpha
        int argmax_weight = np.argmax(weight_in_clusters)

    for j in range(n_clusters):
        if weight_in_clusters[j] > 0:
            alpha = 1.0 / weight_in_clusters[j]
            for k in range(n_features):
                centers[j, k] *= alpha
        else:
            # For convenience, we avoid setting empty clusters at the origin but place
            # them at the location of the biggest cluster.
            for k in range(n_features):
                centers[j, k] = centers[argmax_weight, k]


cdef void _center_shift(
        const floating[:, ::1] centers_old,  # IN
        const floating[:, ::1] centers_new,  # IN
        floating[::1] center_shift           # OUT
):
    """Compute shift between old and new centers."""
    cdef:
        int n_clusters = centers_old.shape[0]
        int n_features = centers_old.shape[1]
        int j

    for j in range(n_clusters):
        center_shift[j] = _euclidean_dense_dense(
            &centers_new[j, 0], &centers_old[j, 0], n_features, False)


def _is_same_clustering(
    const int[::1] labels1,
    const int[::1] labels2,
    n_clusters
):
    """Check if two arrays of labels are the same up to a permutation of the labels"""
    cdef int[::1] mapping = np.full(fill_value=-1, shape=(n_clusters,), dtype=np.int32)
    cdef int i

    for i in range(labels1.shape[0]):
        if mapping[labels1[i]] == -1:
            mapping[labels1[i]] = labels2[i]
        elif mapping[labels1[i]] != labels2[i]:
            return False
    return True
