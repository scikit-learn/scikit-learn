from cython cimport floating, integral
from cython.parallel cimport parallel, prange
from libcpp.map cimport map as cmap
from libc.stdlib cimport free

cimport numpy as cnp

cnp.import_array()

from ._argkmin cimport PairwiseDistancesArgKmin64
from ._datasets_pair cimport DatasetsPair, DenseDenseDatasetsPair
from ...utils._typedefs cimport ITYPE_t, DTYPE_t
from ...utils._typedefs import ITYPE, DTYPE
from ...utils._sorting cimport simultaneous_sort
import numpy as np

from sklearn.utils.fixes import threadpool_limits


cpdef enum Weight:
    uniform = 0
    distance = 1
    other = 2

cdef class PairwiseDistancesArgKminLabels64(PairwiseDistancesArgKmin64):
    """
    64bit implementation of PairwiseDistancesArgKmin that also keeps of labels.
    """
    cdef ITYPE_t[:] argkmin_labels, labels, unique_labels
    cdef DTYPE_t[:,:] label_counts
    cdef cmap[ITYPE_t, ITYPE_t] labels_to_index
    cdef Weight weight_type

    @classmethod
    def compute(
        cls,
        X,
        Y,
        ITYPE_t k,
        weights,
        labels,
        str metric="euclidean",
        chunk_size=None,
        dict metric_kwargs=None,
        str strategy=None,
        bint return_distance=False,
    ):
        """Compute the argkmin reduction.

        This classmethod is responsible for introspecting the arguments
        values to dispatch to the most appropriate implementation of
        :class:`PairwiseDistancesArgKmin64`.

        This allows decoupling the API entirely from the implementation details
        whilst maintaining RAII: all temporarily allocated datastructures necessary
        for the concrete implementation are therefore freed when this classmethod
        returns.

        No instance should directly be created outside of this class method.
        """
        pda = PairwiseDistancesArgKminLabels64(
            datasets_pair=DatasetsPair.get_for(X, Y, metric, metric_kwargs),
            k=k,
            chunk_size=chunk_size,
            strategy=strategy,
            weights=weights,
            labels=labels,
        )
        # Limit the number of threads in second level of nested parallelism for BLAS
        # to avoid threads over-subscription (in GEMM for instance).
        with threadpool_limits(limits=1, user_api="blas"):
            if pda.execute_in_parallel_on_Y:
                pda._parallel_on_Y()
            else:
                pda._parallel_on_X()

        return pda._finalize_results(return_distance)

    def __init__(
        self,
        DatasetsPair datasets_pair,
        chunk_size=None,
        strategy=None,
        ITYPE_t k=1,
        weights=None,
        cnp.ndarray[ndim=1, dtype=ITYPE_t] labels=None,
    ):
        super().__init__(
            datasets_pair=datasets_pair,
            chunk_size=chunk_size,
            strategy=strategy,
            k=k,
        )

        #TODO: Remove after adding implementation
        self.execute_in_parallel_on_Y = False

        if weights == "uniform":
            self.weight_type = Weight.uniform
        elif weights == "distance":
            self.weight_type = Weight.distance
        else:
            self.weight_type = Weight.other
        self.labels = labels
        self.unique_labels = np.unique(labels)

        # Heap to be returned
        self.argkmin_labels = np.full(self.n_samples_X, -1, dtype=ITYPE)

        # Map from set of unique labels to their indices in `label_counts`
        self.labels_to_index = {label:idx for idx, label in enumerate(self.unique_labels)}
        # Buffer used in building a histogram for one-pass weighted mode
        self.label_counts = np.zeros((self.n_samples_X,  len(self.unique_labels)), dtype=DTYPE)

    def _finalize_results(self, bint return_distance=False):
        if return_distance:
            # We need to recompute distances because we relied on
            # surrogate distances for the reduction.
            self.compute_exact_distances()

            # Values are returned identically to the way `KNeighborsMixin.kneighbors`
            # returns values. This is counter-intuitive but this allows not using
            # complex adaptations where `PairwiseDistancesArgKmin64.compute` is called.
            return np.asarray(self.argkmin_distances), np.asarray(self.argkmin_indices), np.asarray(self.argkmin_labels)

        return np.asarray(self.argkmin_indices), np.asarray(self.argkmin_labels)

    cdef void _parallel_on_X_prange_iter_finalize(
        self,
        ITYPE_t thread_num,
        ITYPE_t X_start,
        ITYPE_t X_end,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, y_idx, label_index, sample_index
            DTYPE_t max_weight
            ITYPE_t label, max_label
            DTYPE_t val, weight
        # Sorting the main heaps portion associated to `X[X_start:X_end]`
        # in ascending order w.r.t the distances.
        for idx in range(X_end - X_start):
            simultaneous_sort(
                self.heaps_r_distances_chunks[thread_num] + idx * self.k,
                self.heaps_indices_chunks[thread_num] + idx * self.k,
                self.k
            )
            # One-pass top-one weighted mode
            sample_index = thread_num * self.X_n_samples_chunk + idx
            max_weight = -1
            for jdx in range(self.k): # Iterate through k-nearest neighbors
                if self.weight_type == Weight.uniform:
                    weight = 1.
                elif self.weight_type == Weight.distance:
                    weight = 1. / self.heaps_r_distances_chunks[thread_num][idx*self.k+jdx]

                y_idx = self.heaps_indices_chunks[thread_num][idx*self.k+jdx]
                label = self.labels[y_idx]
                label_index = self.labels_to_index[label]
                self.label_counts[sample_index][label_index] += weight
                val = self.label_counts[sample_index][label_index]
                if max_weight < val or (max_weight == val and label < max_label):
                    max_label = label
                    max_weight = val
            self.argkmin_labels[thread_num*self.X_n_samples_chunk + idx] = max_label

    cdef void _parallel_on_Y_finalize(
        self,
    ) nogil:
        cdef:
            ITYPE_t idx, jdx, y_idx, label_index, sample_index, thread_num
            DTYPE_t max_weight
            ITYPE_t label, max_label
            DTYPE_t val, weight

        with nogil, parallel(num_threads=self.chunks_n_threads):
            # Deallocating temporary datastructures
            for thread_num in prange(self.chunks_n_threads, schedule='static'):
                free(self.heaps_r_distances_chunks[thread_num])
                free(self.heaps_indices_chunks[thread_num])

            # Sorting the main in ascending order w.r.t the distances.
            # This is done in parallel sample-wise (no need for locks).
            for idx in prange(self.n_samples_X, schedule='static'):
                simultaneous_sort(
                    &self.argkmin_distances[idx, 0],
                    &self.argkmin_indices[idx, 0],
                    self.k,
                )
                # One-pass top-one weighted mode
                max_weight = -1
                for jdx in range(self.k):
                    if self.weight_type == Weight.uniform:
                        weight = 1.
                    elif self.weight_type == Weight.distance:
                        weight = 1. / self.argkmin_distances[idx][jdx]

                    y_idx = self.argkmin_indices[idx][jdx]
                    label = self.labels[y_idx]
                    label_index = self.labels_to_index[label]
                    self.label_counts[idx][label_index] += weight
                    val = self.label_counts[idx][label_index]
                    if max_weight < val or (max_weight == val and label < max_label):
                        max_label = label
                        max_weight = val
                self.argkmin_labels[idx] = max_label
        return
