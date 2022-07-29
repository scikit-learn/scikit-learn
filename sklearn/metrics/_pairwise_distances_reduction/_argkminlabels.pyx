from cython cimport floating, integral
from libcpp.map cimport map as cmap

cimport numpy as cnp

cnp.import_array()

from ._argkmin cimport PairwiseDistancesArgKmin64
from ._datasets_pair cimport DatasetsPair
from ...utils._typedefs cimport ITYPE_t, DTYPE_t
from ...utils._typedefs import ITYPE, DTYPE
import numpy as np

from sklearn.utils.fixes import threadpool_limits


cdef class PairwiseDistancesArgKminLabels64(PairwiseDistancesArgKmin64):
    """
    64bit implementation of PairwiseDistancesArgKmin that also keeps of labels.
    """
    cdef ITYPE_t[:] argkmin_labels
    cdef cmap[int, int] labels_to_index

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
        cnp.ndarray[ndim=1, dtype=DTYPE_t] weights=None,
        cnp.ndarray[ndim=1, dtype=ITYPE_t] labels=None,
    ):
        super().__init__(
            datasets_pair=datasets_pair,
            chunk_size=chunk_size,
            strategy=strategy,
            k=k,
        )
        self.weights = weights
        self.labels = labels
        self.argkmin_labels = cnp.empty(self.n_samples_X, dtype=ITYPE_t)
        self.unique_labels = labels.unique()
        self.labels_to_index = {label:idx for idx, label in enumerate(self.unique_labels)}
        self.max_label = labels.max()

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
            ITYPE_t idx, jdx, kdx
            DTYPE_t[:] label_counts = cnp.zeros((len(self.unique_labels),), dtype=DTYPE_t)
            DTYPE_t max_weight = 0
            ITYPE_t max_label
            DTYPE_t val
        # Sorting the main heaps portion associated to `X[X_start:X_end]`
        # in ascending order w.r.t the distances.
        for idx in range(X_end - X_start):
            simultaneous_sort(
                self.heaps_r_distances_chunks[thread_num] + idx * self.k,
                self.heaps_indices_chunks[thread_num] + idx * self.k,
                self.k
            )
            # One-pass top-one weighted mode
            for jdx in range(k):
                kdx = self.heaps_indices_chunks[thread_num][idx*self.k+jdx]
                label_counts[kdx] += self.weights[kdx]
                val = label_counts[kdx]
                if max_weight < val:
                    max_label = kdx
                    max_weight = val
            self.argkmin_labels[thread_num*self.X_n_samples_chunk + idk] = max_label
