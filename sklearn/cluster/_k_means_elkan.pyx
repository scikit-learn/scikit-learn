# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Andreas Mueller
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport sqrt

from ..metrics import euclidean_distances
from ._k_means import _centers_dense


cdef double *getfloatpointer(np.ndarray[np.float64_t, ndim=2, mode='c'] data):
    return <double *>(data.data) 


cdef double d(double* a, double* b, int n_features) nogil:
    cdef double result, tmp
    result = 0
    cdef int i
    for i in range(n_features):
        tmp = (a[0] - b[0])
        a += 1
        b += 1
        result += tmp * tmp
    return sqrt(result)


cdef assign_labels(double* X, double* centers, double[:, :] center_distances,
                   int[:] labels, double[:, :] lower_bounds, double[:]
                   distances, int n_samples, int n_features, int n_clusters):
    # assigns closest center to X
    # uses triangle inequality
    cdef double* x
    cdef double* c
    cdef double d_c, dist
    cdef int c_x, j, sample
    for sample in range(n_samples):
        # assign first cluster center
        c_x = 0
        x = X + sample * n_features
        d_c = d(x, centers, n_features)
        lower_bounds[sample, 0] = d_c
        for j in range(n_clusters):
            if d_c > center_distances[c_x, j]:
                c = centers + j * n_features
                dist = d(x, c, n_features)
                lower_bounds[sample, j] = dist
                if dist < d_c:
                    d_c = dist
                    c_x = j
        labels[sample] = c_x
        distances[sample] = d_c


def k_means_elkan(X_, int n_clusters, init, float tol=1e-4, int max_iter=30, verbose=False):
    #initialize
    centers_ = init
    cdef double* centers_p = getfloatpointer(centers_)
    cdef double* X_p = getfloatpointer(X_)
    cdef double* x_p
    cdef int n_samples = X_.shape[0]
    cdef int n_features = X_.shape[1]
    cdef int point_index, center_index, label
    cdef float upper_bound, distance
    cdef double[:, :] center_distances = euclidean_distances(centers_) / 2.
    cdef double[:, :] lower_bounds = np.zeros((n_samples, n_clusters))
    cdef double[:] distance_next_center
    labels_ = np.empty(n_samples, dtype=np.int32)
    cdef int[:] labels = labels_
    upper_bounds_ = np.empty(n_samples, dtype=np.float)
    cdef double[:] upper_bounds = upper_bounds_
    assign_labels(X_p, centers_p, center_distances, labels, lower_bounds,
                  upper_bounds, n_samples, n_features, n_clusters)
    cdef np.uint8_t[:] bounds_tight = np.ones(n_samples, dtype=np.uint8)
    cdef np.uint8_t[:] points_to_update = np.zeros(n_samples, dtype=np.uint8)
    for iteration in range(max_iter):
        if verbose:
            print("start iteration")
        # we could find the closest center in O(n) but
        # this does not seem to be the bottleneck
        distance_next_center = np.sort(center_distances, axis=0)[1]
        if verbose:
            print("done sorting")
        for point_index in range(n_samples):
            upper_bound = upper_bounds[point_index]
            label = labels[point_index]
            if distance_next_center[label] >= upper_bound:
                continue
            x_p = X_p + point_index * n_features
            # check other update conditions
            # TODO: get pointer to lower_bounds[point_index, center_index]
            for center_index in range(n_clusters):
                if (center_index != label
                        and (upper_bound > lower_bounds[point_index, center_index])
                        and (upper_bound > center_distances[center_index, label])):
                    # update distance to center
                    if not bounds_tight[point_index]:
                        upper_bound = d(x_p, centers_p + label * n_features, n_features)
                        lower_bounds[point_index, label] = upper_bound
                        bounds_tight[point_index] = 1
                    # check for relabels
                    if (upper_bound > lower_bounds[point_index, center_index]
                            or (upper_bound > center_distances[label, center_index])):
                        distance = d(x_p, centers_p + center_index * n_features, n_features)
                        lower_bounds[point_index, center_index] = distance
                        if distance < upper_bound:
                            label = center_index
                            upper_bound = distance
            labels[point_index] = label
            upper_bounds[point_index] = upper_bound
        if verbose:
            print("end inner loop")
        # compute new centers
        new_centers = _centers_dense(X_, labels_, n_clusters, upper_bounds_)
        bounds_tight = np.zeros(n_samples, dtype=np.uint8)

        # compute distance each center moved
        center_shift = np.sqrt(np.sum((centers_ - new_centers) ** 2, axis=1))
        # update bounds accordingly
        lower_bounds = np.maximum(lower_bounds - center_shift, 0)
        upper_bounds = upper_bounds + center_shift[labels_]
        # reassign centers
        centers_ = new_centers
        centers_p = getfloatpointer(new_centers)
        # update between-center distances
        center_distances = euclidean_distances(centers_) / 2.
        if verbose:
            print('Iteration %i, inertia %s'
                  % (iteration, np.sum((X_ - centers_[labels]) ** 2)))
        center_shift_total = np.sum(center_shift)
        if center_shift_total ** 2 < tol:
            if verbose:
                print("center shift %e within tolerance %e"
                      % (center_shift_total, tol))
            break
    return centers_, labels_, iteration
