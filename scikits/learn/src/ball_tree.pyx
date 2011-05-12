""" Cython bindings for the C++ BallTree code.

A Ball Tree is a data structure which can be used
to perform fast neighbor searches in data sets of
low to medium dimensionality.
"""
# Author: Thouis Jones
# License: BSD

from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

# C++ functions
cdef extern from "BallTreePoint.h":
   cdef cppclass Point:
       Point(size_t)
       size_t size()
   void SET(Point *, size_t, double)

ctypedef Point *Point_p


cdef extern from "BallTree.h":
    cdef cppclass cBallTree "BallTree<Point>":
        cBallTree(vector[Point_p] *, size_t)
        double query(Point *, vector[size_t] &) except +
    double Euclidean_Dist(Point *, Point *) except +


cdef Point *make_point(vals):
    pt = new Point(vals.size)
    for idx, v in enumerate(vals.flat):
        SET(pt, idx, v)
    return pt


################################################################################
# Cython wrapper
cdef class BallTree:
    """
    Ball Tree for fast nearest-neighbor searches :

    BallTree(M, leafsize=20)

    Parameters
    ----------
    M : array-like, shape = [N,D]
            N is the number of points in the data set, and
            D is the dimension of the parameter space.
            Note: if M is an aligned array of doubles (not
            necessarily contiguous) then data will not be
            copied. Otherwise, an internal copy will be made.

    leafsize : positive integer (default = 20)
        number of points at which to switch to brute-force. Currently not
        implemented.

    Notes
    -----
    brute-force search was removed. docs should be accordingly.
    """
    cdef cBallTree *bt_ptr
    cdef vector[Point_p] *ptdata
    cdef size_t num_points
    cdef size_t num_dims
    cdef public object data

    def __cinit__(self, arr, leafsize=20):
        # copy points into ptdata
        arr = np.atleast_2d(arr).astype(np.double)
        assert arr.ndim == 2, "input points must be 2 dimensional (points x dimensions)"
        num_points, num_dims = self.num_points, self.num_dims = arr.shape
        self.ptdata = new vector[Point_p]()
        for i in range(num_points):
            self.ptdata.push_back(make_point(arr[i, :]))
        self.bt_ptr = new cBallTree(self.ptdata, leafsize)
        self.data = arr.copy()

    def __dealloc__(self):
        cdef Point *temp
        # __dealloc__ is called if __cinit__ fails at any point
        if self.ptdata:
            for idx in range(self.ptdata.size()):
                # Cython won't allow the more direct form
                temp = self.ptdata.at(idx)
                del temp
            del self.ptdata
        if self.bt_ptr:
            del self.bt_ptr

    def query(self, x, k=1, return_distance=True):
        """
        query(x, k=1, return_distance=True)

        query the Ball Tree for the k nearest neighbors

        Parameters
        ----------
        x : array-like, last dimension self.dim
              An array of points to query
        k : integer  (default = 1)
              The number of nearest neighbors to return
        return_distance : boolean (default = True)
              if True, return a tuple (d,i)
              if False, return array i

        Returns
        -------
        i    : if return_distance == False
        (d,i) : if return_distance == True

        d : array of doubles - shape: x.shape[:-1] + (k,)
            each entry gives the list of distances to the
            neighbors of the corresponding point
            (note that distances are not sorted)

        i : array of integers - shape: x.shape[:-1] + (k,)
            each entry gives the list of indices of
            neighbors of the corresponding point
            (note that neighbors are not sorted)
        """
        x = np.atleast_2d(x)
        assert x.shape[-1] == self.num_dims
        assert k <= self.num_points

        cdef Point *temp
        cdef vector[size_t] results = vector[size_t](<size_t>k)

        # almost-flatten x for iteration
        orig_shape = x.shape
        x.reshape((-1, self.num_dims))

        # allocate output
        out_indices = np.zeros((x.shape[0], k), np.int64)
        if return_distance:
            out_distances = np.zeros((x.shape[0], k), np.float64)
        for pt_idx, pt in enumerate(x):
            temp = make_point(pt)
            self.bt_ptr.query(temp, results)
            for neighbor_idx in range(k):
                out_indices[pt_idx, neighbor_idx] = results[neighbor_idx]
                if return_distance:
                    # It would be better to use self.bt_ptr.Dist(), but it's private.
                    out_distances[pt_idx, neighbor_idx] = Euclidean_Dist(<Point *> self.ptdata.at(results[neighbor_idx]),
                                                                         temp)
            del temp

        # deflatten results
        if return_distance:
            return (out_distances.reshape((orig_shape[:-1]) + (k,)),
                    out_indices.reshape((orig_shape[:-1]) + (k,)))
        else:
            return out_indices.reshape((orig_shape[:-1]) + (k,))
