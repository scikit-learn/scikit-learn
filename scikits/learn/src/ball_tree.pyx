""" Cython bindings for the C++ BallTree code.

A Ball Tree is a data structure which can be used
to perform fast neighbor searches in data sets of
low to medium dimensionality.
"""
# Authors: Thouis Jones
#          Jake Vanderplas
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

ctypedef double Point_dtype
ctypedef Point *Point_p


cdef extern from "BallTree.h":
    cdef cppclass cBallTree "BallTree<Point>":
        cBallTree(vector[Point_p] *, size_t)
        void query(Point *, vector[size_t] &) except +
        void query(Point *, vector[size_t] &, vector[Point_dtype] &) except +
        double query_radius(Point *, Point_dtype) except +
        double query_radius(Point *, Point_dtype, vector[size_t] &) except +
        double query_radius(Point *, Point_dtype,
                            vector[size_t] &, vector[Point_dtype] &) except +
    double Euclidean_Dist(Point *, Point *) except +


cdef Point *make_point(vals):
    pt = new Point(vals.size)
    for idx, v in enumerate(vals.flat):
        SET(pt, idx, v)
    return pt


##############################################################################
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

    def __cinit__(self, arr, size_t leafsize=20):
        # copy points into ptdata
        arr = np.atleast_2d(arr).astype(np.double)
        assert arr.ndim == 2, ("input points must be 2 dimensional "
                               "(points x dimensions)")
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
        cdef vector[size_t] ind = vector[size_t](<size_t>k)
        cdef vector[Point_dtype] dist = vector[Point_dtype](<size_t>k)

        # almost-flatten x for iteration
        orig_shape = x.shape
        x = x.reshape((-1, self.num_dims))

        # allocate output
        out_indices = np.zeros((x.shape[0], k), np.int64)
        if return_distance:
            out_distances = np.zeros((x.shape[0], k), np.float64)
        for pt_idx, pt in enumerate(x):
            temp = make_point(pt)
            self.bt_ptr.query(temp, ind, dist)
            for neighbor_idx in range(k):
                out_indices[pt_idx, neighbor_idx] = ind[neighbor_idx]
                if return_distance:
                    out_distances[pt_idx, neighbor_idx] = dist[neighbor_idx]
            del temp

        # deflatten results
        if return_distance:
            return (out_distances.reshape((orig_shape[:-1]) + (k,)),
                    out_indices.reshape((orig_shape[:-1]) + (k,)))
        else:
            return out_indices.reshape((orig_shape[:-1]) + (k,))

    def query_radius(self, x, r, return_distance=False, count_only=False):
        """
        query_radius(self, x, r, count_only = False):

        query the Ball Tree for neighbors within a ball of size r

        Parameters
        ----------
        x : array-like, last dimension self.dim
              An array of points to query
        r : distance within which neighbors are returned
              r can be a single value, or an array of values
              of shape x.shape[:-1] if different radii are
              desired for each point.
        return_distance : boolean (default = False)
              if True,  return distances to neighbors of each point
              if False, return only neighbors
              Note that unlike query() above, setting return_distance=True
              adds to the computation time.  Not all distances must be
              calculated for return_distance=False.
        count_only : boolean (default = False)
              if True,  return only the count of points
                         within distance r
              if False, return the indices of all points
                         within distance r
              If return_distance==True, setting count_only=True will
              raise an error.

        Returns
        -------
        n     : if count_only == True
        i     : if count_only == False and return_distance == False
        (i,d) : if count_only == False and return_distance == True

        n : array of integers - shape: x.shape[:-1]
            each entry gives the number of neighbors within
            a distance r of the corresponding point.

        i : array of objects  - shape: x.shape[:-1]
            each element is a numpy integer array
            listing the indices of neighbors
            of the corresponding point
            (note that neighbors are not sorted by distance)

        d : array of objects  - shape: x.shape[:-1]
            each element is a numpy double array
            listing the distances corresponding to indices in i.
        """
        if count_only and return_distance:
            raise ValueError("count_only and return_distance "
                             "cannot both be true")
        x = np.atleast_2d(x)
        assert x.shape[-1] == self.num_dims

        r = np.atleast_1d(r)
        if r.shape==(1,):
            r = r[0]*np.ones(x.shape[:-1], dtype=np.double)
        else:
            assert r.shape == x.shape[:-1]

        cdef Point *temp
        cdef vector[size_t] ind_vec = vector[size_t](0)
        cdef vector[Point_dtype] dist_vec = vector[Point_dtype](0)

        # almost-flatten x for iteration
        orig_shape = x.shape
        x = x.reshape((-1, self.num_dims))
        r = r.reshape(-1)

        # allocate output
        if count_only:
            count = np.zeros(x.shape[0], np.int64)
        elif return_distance:
            indices = np.empty(x.shape[0], dtype='object')
            distances = np.empty(x.shape[0], dtype='object')
        else:
            indices = np.empty(x.shape[0], dtype='object')

        for pt_idx, pt in enumerate(x):
            temp = make_point(pt)
            if count_only:
                count[pt_idx] = self.bt_ptr.query_radius(temp, r[pt_idx])
            elif return_distance:
                self.bt_ptr.query_radius(temp, r[pt_idx], ind_vec, dist_vec)
                indices[pt_idx] = np.zeros(ind_vec.size(), dtype=np.int)
                distances[pt_idx] = np.zeros(dist_vec.size(), dtype=np.double)
                for neighbor_idx in range(ind_vec.size()):
                    indices[pt_idx][neighbor_idx] = ind_vec[neighbor_idx]
                    distances[pt_idx][neighbor_idx] = dist_vec[neighbor_idx]
                ind_vec.resize(0)
                dist_vec.resize(0)
            else:
                self.bt_ptr.query_radius(temp, r[pt_idx], ind_vec)
                indices[pt_idx] = np.zeros(ind_vec.size(), dtype=np.int)
                for neighbor_idx in range(ind_vec.size()):
                    indices[pt_idx][neighbor_idx] = ind_vec[neighbor_idx]
                ind_vec.resize(0)
            del temp

        # deflatten results
        if count_only:
            return count.reshape(orig_shape[:-1])
        elif return_distance:
            return (indices.reshape(orig_shape[:-1]),
                    distances.reshape(orig_shape[:-1]))
        else:
            return indices.reshape(orig_shape[:-1])
