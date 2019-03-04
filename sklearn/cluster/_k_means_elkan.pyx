# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3
#
# Author: Andreas Mueller
#
# Licence: BSD 3 clause

import numpy as np
cimport numpy as np
cimport cython
from cython cimport floating

from libc.math cimport sqrt

from ..metrics import euclidean_distances
from ._k_means import _centers_dense


cdef floating euclidean_dist(floating* a, floating* b, int n_features) nogil:
    cdef floating result, tmp
    result = 0
    cdef int i
    for i in range(n_features):
        tmp = (a[i] - b[i])
        result += tmp * tmp
    return sqrt(result)


cdef update_labels_distances_inplace(
        floating* X, floating* centers, floating[:, :] center_half_distances,
        int[:] labels, floating[:, :] lower_bounds, floating[:] upper_bounds,
        int n_samples, int n_features, int n_clusters):
    """
    Calculate upper and lower bounds for each sample.

    Given X, centers and the pairwise distances divided by 2.0 between the
    centers this calculates the upper bounds and lower bounds for each sample.
    The upper bound for each sample is set to the distance between the sample
    and the closest center.

    The lower bound for each sample is a one-dimensional array of n_clusters.
    For each sample i assume that the previously assigned cluster is c1 and the
    previous closest distance is dist, for a new cluster c2, the
    lower_bound[i][c2] is set to distance between the sample and this new
    cluster, if and only if dist > center_half_distances[c1][c2]. This prevents
    computation of unnecessary distances for each sample to the clusters that
    it is unlikely to be assigned to.

    Parameters
    ----------
    X : nd-array, shape (n_samples, n_features)
        The input data.

    centers : nd-array, shape (n_clusters, n_features)
        The cluster centers.

    center_half_distances : nd-array, shape (n_clusters, n_clusters)
        The half of the distance between any 2 clusters centers.

    labels : nd-array, shape(n_samples)
        The label for each sample. This array is modified in place.

    lower_bounds : nd-array, shape(n_samples, n_clusters)
        The lower bound on the distance between a sample and each cluster
        center. It is modified in place.

    upper_bounds : nd-array, shape(n_samples,)
        The distance of each sample from its closest cluster center.  This is
        modified in place by the function.

    n_samples : int
        The number of samples.

    n_features : int
        The number of features.

    n_clusters : int
        The number of clusters.
    """
    # assigns closest center to X
    # uses triangle inequality
    cdef floating* x
    cdef floating* c
    cdef floating d_c, dist
    cdef int c_x, j, sample
    for sample in range(n_samples):
        # assign first cluster center
        c_x = 0
        x = X + sample * n_features
        d_c = euclidean_dist(x, centers, n_features)
        lower_bounds[sample, 0] = d_c
        for j in range(1, n_clusters):
            if d_c > center_half_distances[c_x, j]:
                c = centers + j * n_features
                dist = euclidean_dist(x, c, n_features)
                lower_bounds[sample, j] = dist
                if dist < d_c:
                    d_c = dist
                    c_x = j
        labels[sample] = c_x
        upper_bounds[sample] = d_c


def k_means_elkan(np.ndarray[floating, ndim=2, mode='c'] X_,
                  np.ndarray[floating, ndim=1, mode='c'] sample_weight,
                  int n_clusters,
                  np.ndarray[floating, ndim=2, mode='c'] init,
                  float tol=1e-4, int max_iter=30, verbose=False):
    """Run Elkan's k-means.

    Parameters
    ----------
    X_ : nd-array, shape (n_samples, n_features)

    sample_weight : nd-array, shape (n_samples,)
        The weights for each observation in X.

    n_clusters : int
        Number of clusters to find.

    init : nd-array, shape (n_clusters, n_features)
        Initial position of centers.

    tol : float, default=1e-4
        The relative increment in cluster means before declaring convergence.

    max_iter : int, default=30
    Maximum number of iterations of the k-means algorithm.

    verbose : bool, default=False
        Whether to be verbose.

    """
    if floating is float:
        dtype = np.float32
    else:
        dtype = np.float64

    # initialize
    cdef np.ndarray[floating, ndim=2, mode='c'] centers_ = init
    cdef floating* centers_p = <floating*>centers_.data
    cdef floating* X_p = <floating*>X_.data
    cdef floating* x_p
    cdef Py_ssize_t n_samples = X_.shape[0]
    cdef Py_ssize_t n_features = X_.shape[1]
    cdef int point_index, center_index, label
    cdef floating upper_bound, distance
    cdef floating[:, :] center_half_distances = euclidean_distances(centers_) / 2.
    cdef floating[:, :] lower_bounds = np.zeros((n_samples, n_clusters), dtype=dtype)
    cdef floating[:] distance_next_center
    labels_ = np.empty(n_samples, dtype=np.int32)
    cdef int[:] labels = labels_
    upper_bounds_ = np.empty(n_samples, dtype=dtype)
    cdef floating[:] upper_bounds = upper_bounds_

    # Get the initial set of upper bounds and lower bounds for each sample.
    update_labels_distances_inplace(X_p, centers_p, center_half_distances,
                                    labels, lower_bounds, upper_bounds,
                                    n_samples, n_features, n_clusters)
    cdef np.uint8_t[:] bounds_tight = np.ones(n_samples, dtype=np.uint8)
    cdef np.uint8_t[:] points_to_update = np.zeros(n_samples, dtype=np.uint8)
    cdef np.ndarray[floating, ndim=2, mode='c'] new_centers

    if max_iter <= 0:
        raise ValueError('Number of iterations should be a positive number'
        ', got %d instead' % max_iter)

    col_indices = np.arange(center_half_distances.shape[0], dtype=np.int)
    for iteration in range(max_iter):
        if verbose:
            print("start iteration")

        cd =  np.asarray(center_half_distances)
        distance_next_center = np.partition(cd, kth=1, axis=0)[1]

        if verbose:
            print("done sorting")

        for point_index in range(n_samples):
            upper_bound = upper_bounds[point_index]
            label = labels[point_index]

            # This means that the next likely center is far away from the
            # currently assigned center and the sample is unlikely to be
            # reassigned.
            if distance_next_center[label] >= upper_bound:
                continue
            x_p = X_p + point_index * n_features

            # TODO: get pointer to lower_bounds[point_index, center_index]
            for center_index in range(n_clusters):

                # If this holds, then center_index is a good candidate for the
                # sample to be relabelled, and we need to confirm this by
                # recomputing the upper and lower bounds.
                if (center_index != label
                        and (upper_bound > lower_bounds[point_index, center_index])
                        and (upper_bound > center_half_distances[center_index, label])):

                    # Recompute the upper bound by calculating the actual distance
                    # between the sample and label.
                    if not bounds_tight[point_index]:
                        upper_bound = euclidean_dist(x_p, centers_p + label * n_features, n_features)
                        lower_bounds[point_index, label] = upper_bound
                        bounds_tight[point_index] = 1

                    # If the condition still holds, then compute the actual distance between
                    # the sample and center_index. If this is still lesser than the previous
                    # distance, reassign labels.
                    if (upper_bound > lower_bounds[point_index, center_index]
                            or (upper_bound > center_half_distances[label, center_index])):
                        distance = euclidean_dist(x_p, centers_p + center_index * n_features, n_features)
                        lower_bounds[point_index, center_index] = distance
                        if distance < upper_bound:
                            label = center_index
                            upper_bound = distance

            labels[point_index] = label
            upper_bounds[point_index] = upper_bound

        if verbose:
            print("end inner loop")

        # compute new centers
        new_centers = _centers_dense(X_, sample_weight, labels_,
                                     n_clusters, upper_bounds_)
        bounds_tight[:] = 0

        # compute distance each center moved
        center_shift = np.sqrt(np.sum((centers_ - new_centers) ** 2, axis=1))

        # update bounds accordingly
        lower_bounds = np.maximum(lower_bounds - center_shift, 0)
        upper_bounds = upper_bounds + center_shift[labels_]

        # reassign centers
        centers_ = new_centers
        centers_p = <floating*>new_centers.data

        # update between-center distances
        center_half_distances = euclidean_distances(centers_) / 2.
        if verbose:
            print('Iteration %i, inertia %s'
                    % (iteration, np.sum((X_ - centers_[labels]) ** 2 *
                                         sample_weight[:,np.newaxis])))
        center_shift_total = np.sum(center_shift)
        if center_shift_total ** 2 < tol:
            if verbose:
                print("center shift %e within tolerance %e"
                      % (center_shift_total, tol))
            break

    # We need this to make sure that the labels give the same output as
    # predict(X)
    if center_shift_total > 0:
        update_labels_distances_inplace(X_p, centers_p, center_half_distances,
                                        labels, lower_bounds, upper_bounds,
                                        n_samples, n_features, n_clusters)
    return centers_, labels_, iteration + 1
