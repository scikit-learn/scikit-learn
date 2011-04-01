""" Cython bindings for the C++ BallTree code.
"""
# Author: Thouis Jones
# License: BSD

from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

# C++ functions
cdef extern from "BallTreePoint.h":
   cdef cppclass Point:
       Point(int)
       int size()
   void SET(Point *, int, double)

ctypedef Point *Point_p


cdef extern from "BallTree.h":
    cdef cppclass cBallTree "BallTree<Point>":
        cBallTree(vector[Point_p] *, int)
        double query(Point *, vector[long int] &) except +
    double Euclidean_Dist(Point *, Point *) except +


cdef Point *make_point(vals):
    pt = new Point(vals.size)
    for idx, v in enumerate(vals.flat):
        SET(pt, idx, v)
    return pt
    

################################################################################
# Cython wrapper
cdef class BallTree:
    cdef cBallTree *bt_ptr
    cdef vector[Point_p] *ptdata
    cdef int num_points
    cdef int num_dims
    cdef public object data

    def __cinit__(self, arr, leafsize=20):
        # copy points into ptdata
        num_points, num_dims = self.num_points, self.num_dims = arr.shape
        self.ptdata = new vector[Point_p]()
        for i in range(num_points):
            self.ptdata.push_back(make_point(arr[i, :]))
        self.bt_ptr = new cBallTree(self.ptdata, leafsize)
        self.data = arr.copy()

    def __dealloc__(self):
        cdef Point *temp
        for idx in range(self.ptdata.size()):
            # Cython won't allow the more direct form
            temp = self.ptdata.at(idx)
            del temp
        del self.ptdata
        del self.bt_ptr

    def query(self, x, k=1, return_distance=True):
        x = np.atleast_2d(x)
        assert x.shape[-1] == self.num_dims
        assert k <= self.num_points

        cdef Point *temp
        cdef vector[long] results = vector[long](<int>k)

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
