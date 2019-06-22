# Copyright Anne M. Archibald 2008
# Released under the scipy license

from __future__ import division, print_function, absolute_import

from numpy.testing import (assert_equal, assert_array_equal, assert_,
                           assert_almost_equal, assert_array_almost_equal)
from pytest import raises as assert_raises
import pytest
from platform import python_implementation
import numpy as np
from scipy.spatial import KDTree, Rectangle, distance_matrix, cKDTree
from scipy.spatial.ckdtree import cKDTreeNode
from scipy.spatial import minkowski_distance

import itertools

def distance_box(a, b, p, boxsize):
    diff = a - b
    diff[diff > 0.5 * boxsize] -= boxsize
    diff[diff < -0.5 * boxsize] += boxsize
    d = minkowski_distance(diff, 0, p)
    return d

class ConsistencyTests:
    def distance(self, a, b, p):
        return minkowski_distance(a, b, p)

    def test_nearest(self):
        x = self.x
        d, i = self.kdtree.query(x, 1)
        assert_almost_equal(d**2,np.sum((x-self.data[i])**2))
        eps = 1e-8
        assert_(np.all(np.sum((self.data-x[np.newaxis,:])**2,axis=1) > d**2-eps))

    def test_m_nearest(self):
        x = self.x
        m = self.m
        dd, ii = self.kdtree.query(x, m)
        d = np.amax(dd)
        i = ii[np.argmax(dd)]
        assert_almost_equal(d**2,np.sum((x-self.data[i])**2))
        eps = 1e-8
        assert_equal(np.sum(np.sum((self.data-x[np.newaxis,:])**2,axis=1) < d**2+eps),m)

    def test_points_near(self):
        x = self.x
        d = self.d
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd,ii):
            if near_d == np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d**2,np.sum((x-self.data[near_i])**2))
            assert_(near_d < d+eps, "near_d=%g should be less than %g" % (near_d,d))
        assert_equal(np.sum(self.distance(self.data,x,2) < d**2+eps),hits)

    def test_points_near_l1(self):
        x = self.x
        d = self.d
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, p=1, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd,ii):
            if near_d == np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d,self.distance(x,self.data[near_i],1))
            assert_(near_d < d+eps, "near_d=%g should be less than %g" % (near_d,d))
        assert_equal(np.sum(self.distance(self.data,x,1) < d+eps),hits)

    def test_points_near_linf(self):
        x = self.x
        d = self.d
        dd, ii = self.kdtree.query(x, k=self.kdtree.n, p=np.inf, distance_upper_bound=d)
        eps = 1e-8
        hits = 0
        for near_d, near_i in zip(dd,ii):
            if near_d == np.inf:
                continue
            hits += 1
            assert_almost_equal(near_d,self.distance(x,self.data[near_i],np.inf))
            assert_(near_d < d+eps, "near_d=%g should be less than %g" % (near_d,d))
        assert_equal(np.sum(self.distance(self.data,x,np.inf) < d+eps),hits)

    def test_approx(self):
        x = self.x
        k = self.k
        eps = 0.1
        d_real, i_real = self.kdtree.query(x, k)
        d, i = self.kdtree.query(x, k, eps=eps)
        assert_(np.all(d <= d_real*(1+eps)))


class Test_random(ConsistencyTests):
    def setup_method(self):
        self.n = 100
        self.m = 4
        np.random.seed(1234)
        self.data = np.random.randn(self.n, self.m)
        self.kdtree = KDTree(self.data,leafsize=2)
        self.x = np.random.randn(self.m)
        self.d = 0.2
        self.k = 10

class Test_random_far(Test_random):
    def setup_method(self):
        Test_random.setup_method(self)
        self.x = np.random.randn(self.m)+10


class Test_small(ConsistencyTests):
    def setup_method(self):
        self.data = np.array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1,0,0],
                              [1,0,1],
                              [1,1,0],
                              [1,1,1]])
        self.kdtree = KDTree(self.data)
        self.n = self.kdtree.n
        self.m = self.kdtree.m
        np.random.seed(1234)
        self.x = np.random.randn(3)
        self.d = 0.5
        self.k = 4

    def test_nearest(self):
        assert_array_equal(
                self.kdtree.query((0,0,0.1), 1),
                (0.1,0))

    def test_nearest_two(self):
        assert_array_equal(
                self.kdtree.query((0,0,0.1), 2),
                ([0.1,0.9],[0,1]))


class Test_small_nonleaf(Test_small):
    def setup_method(self):
        Test_small.setup_method(self)
        self.kdtree = KDTree(self.data,leafsize=1)


class Test_small_compiled(Test_small):
    def setup_method(self):
        Test_small.setup_method(self)
        self.kdtree = cKDTree(self.data)


class Test_small_nonleaf_compiled(Test_small):
    def setup_method(self):
        Test_small.setup_method(self)
        self.kdtree = cKDTree(self.data,leafsize=1)


class Test_random_compiled(Test_random):
    def setup_method(self):
        Test_random.setup_method(self)
        self.kdtree = cKDTree(self.data)


class Test_random_far_compiled(Test_random_far):
    def setup_method(self):
        Test_random_far.setup_method(self)
        self.kdtree = cKDTree(self.data)


class Test_vectorization:
    def setup_method(self):
        self.data = np.array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1,0,0],
                              [1,0,1],
                              [1,1,0],
                              [1,1,1]])
        self.kdtree = KDTree(self.data)

    def test_single_query(self):
        d, i = self.kdtree.query(np.array([0,0,0]))
        assert_(isinstance(d,float))
        assert_(np.issubdtype(i, np.signedinteger))

    def test_vectorized_query(self):
        d, i = self.kdtree.query(np.zeros((2,4,3)))
        assert_equal(np.shape(d),(2,4))
        assert_equal(np.shape(i),(2,4))

    def test_single_query_multiple_neighbors(self):
        s = 23
        kk = self.kdtree.n+s
        d, i = self.kdtree.query(np.array([0,0,0]),k=kk)
        assert_equal(np.shape(d),(kk,))
        assert_equal(np.shape(i),(kk,))
        assert_(np.all(~np.isfinite(d[-s:])))
        assert_(np.all(i[-s:] == self.kdtree.n))

    def test_vectorized_query_multiple_neighbors(self):
        s = 23
        kk = self.kdtree.n+s
        d, i = self.kdtree.query(np.zeros((2,4,3)),k=kk)
        assert_equal(np.shape(d),(2,4,kk))
        assert_equal(np.shape(i),(2,4,kk))
        assert_(np.all(~np.isfinite(d[:,:,-s:])))
        assert_(np.all(i[:,:,-s:] == self.kdtree.n))

    def test_single_query_all_neighbors(self):
        d, i = self.kdtree.query([0,0,0],k=None,distance_upper_bound=1.1)
        assert_(isinstance(d,list))
        assert_(isinstance(i,list))

    def test_vectorized_query_all_neighbors(self):
        d, i = self.kdtree.query(np.zeros((2,4,3)),k=None,distance_upper_bound=1.1)
        assert_equal(np.shape(d),(2,4))
        assert_equal(np.shape(i),(2,4))

        assert_(isinstance(d[0,0],list))
        assert_(isinstance(i[0,0],list))


class Test_vectorization_compiled:
    def setup_method(self):
        self.data = np.array([[0,0,0],
                              [0,0,1],
                              [0,1,0],
                              [0,1,1],
                              [1,0,0],
                              [1,0,1],
                              [1,1,0],
                              [1,1,1]])
        self.kdtree = cKDTree(self.data)

    def test_single_query(self):
        d, i = self.kdtree.query([0,0,0])
        assert_(isinstance(d,float))
        assert_(isinstance(i,int))

    def test_vectorized_query(self):
        d, i = self.kdtree.query(np.zeros((2,4,3)))
        assert_equal(np.shape(d),(2,4))
        assert_equal(np.shape(i),(2,4))

    def test_vectorized_query_noncontiguous_values(self):
        np.random.seed(1234)
        qs = np.random.randn(3,1000).T
        ds, i_s = self.kdtree.query(qs)
        for q, d, i in zip(qs,ds,i_s):
            assert_equal(self.kdtree.query(q),(d,i))

    def test_single_query_multiple_neighbors(self):
        s = 23
        kk = self.kdtree.n+s
        d, i = self.kdtree.query([0,0,0],k=kk)
        assert_equal(np.shape(d),(kk,))
        assert_equal(np.shape(i),(kk,))
        assert_(np.all(~np.isfinite(d[-s:])))
        assert_(np.all(i[-s:] == self.kdtree.n))

    def test_vectorized_query_multiple_neighbors(self):
        s = 23
        kk = self.kdtree.n+s
        d, i = self.kdtree.query(np.zeros((2,4,3)),k=kk)
        assert_equal(np.shape(d),(2,4,kk))
        assert_equal(np.shape(i),(2,4,kk))
        assert_(np.all(~np.isfinite(d[:,:,-s:])))
        assert_(np.all(i[:,:,-s:] == self.kdtree.n))

class ball_consistency:
    tol = 0.0

    def distance(self, a, b, p):
        return minkowski_distance(a * 1.0, b * 1.0, p)

    def test_in_ball(self):
        x = np.atleast_2d(self.x)
        d = np.broadcast_to(self.d, x.shape[:-1])
        l = self.T.query_ball_point(x, self.d, p=self.p, eps=self.eps)
        for i, ind in enumerate(l):
            dist = self.distance(self.data[ind], x[i],self.p) - d[i]*(1.+self.eps)
            norm = self.distance(self.data[ind], x[i],self.p) + d[i]*(1.+self.eps)
            assert_array_equal(dist < self.tol * norm, True)

    def test_found_all(self):
        x = np.atleast_2d(self.x)
        d = np.broadcast_to(self.d, x.shape[:-1])
        l = self.T.query_ball_point(x, self.d, p=self.p, eps=self.eps)
        for i, ind in enumerate(l):
            c = np.ones(self.T.n, dtype=bool)
            c[ind] = False
            dist = self.distance(self.data[c], x[i],self.p) - d[i]/(1.+self.eps)
            norm = self.distance(self.data[c], x[i],self.p) + d[i]/(1.+self.eps)
            assert_array_equal(dist > -self.tol * norm, True)

class Test_random_ball(ball_consistency):

    def setup_method(self):
        n = 100
        m = 4
        np.random.seed(1234)
        self.data = np.random.randn(n,m)
        self.T = KDTree(self.data,leafsize=2)
        self.x = np.random.randn(m)
        self.p = 2.
        self.eps = 0
        self.d = 0.2


class Test_random_ball_compiled(ball_consistency):

    def setup_method(self):
        n = 100
        m = 4
        np.random.seed(1234)
        self.data = np.random.randn(n,m)
        self.T = cKDTree(self.data,leafsize=2)
        self.x = np.random.randn(m)
        self.p = 2.
        self.eps = 0
        self.d = 0.2

class Test_random_ball_compiled_periodic(ball_consistency):
    def distance(self, a, b, p):
        return distance_box(a, b, p, 1.0)

    def setup_method(self):
        n = 10000
        m = 4
        np.random.seed(1234)
        self.data = np.random.uniform(size=(n,m))
        self.T = cKDTree(self.data,leafsize=2, boxsize=1)
        self.x = np.ones(m) * 0.1
        self.p = 2.
        self.eps = 0
        self.d = 0.2

    def test_in_ball_outside(self):
        l = self.T.query_ball_point(self.x + 1.0, self.d, p=self.p, eps=self.eps)
        for i in l:
            assert_(self.distance(self.data[i],self.x,self.p) <= self.d*(1.+self.eps))
        l = self.T.query_ball_point(self.x - 1.0, self.d, p=self.p, eps=self.eps)
        for i in l:
            assert_(self.distance(self.data[i],self.x,self.p) <= self.d*(1.+self.eps))

    def test_found_all_outside(self):
        c = np.ones(self.T.n,dtype=bool)
        l = self.T.query_ball_point(self.x + 1.0, self.d, p=self.p, eps=self.eps)
        c[l] = False
        assert_(np.all(self.distance(self.data[c],self.x,self.p) >= self.d/(1.+self.eps)))

        l = self.T.query_ball_point(self.x - 1.0, self.d, p=self.p, eps=self.eps)
        c[l] = False
        assert_(np.all(self.distance(self.data[c],self.x,self.p) >= self.d/(1.+self.eps)))

class Test_random_ball_compiled_largep_issue9890(ball_consistency):

    # allow some roundoff errors due to numerical issues
    tol = 1e-13

    def setup_method(self):
        n = 1000
        m = 2
        np.random.seed(123)
        self.data = np.random.randint(100, 1000, size=(n, m))
        self.T = cKDTree(self.data)
        self.x = self.data
        self.p = 100
        self.eps = 0
        self.d = 10

class Test_random_ball_approx(Test_random_ball):

    def setup_method(self):
        Test_random_ball.setup_method(self)
        self.eps = 0.1


class Test_random_ball_approx_compiled(Test_random_ball_compiled):

    def setup_method(self):
        Test_random_ball_compiled.setup_method(self)
        self.eps = 0.1

class Test_random_ball_approx_compiled_periodic(Test_random_ball_compiled_periodic):

    def setup_method(self):
        Test_random_ball_compiled_periodic.setup_method(self)
        self.eps = 0.1


class Test_random_ball_far(Test_random_ball):

    def setup_method(self):
        Test_random_ball.setup_method(self)
        self.d = 2.


class Test_random_ball_far_compiled(Test_random_ball_compiled):

    def setup_method(self):
        Test_random_ball_compiled.setup_method(self)
        self.d = 2.

class Test_random_ball_far_compiled_periodic(Test_random_ball_compiled_periodic):

    def setup_method(self):
        Test_random_ball_compiled_periodic.setup_method(self)
        self.d = 2.


class Test_random_ball_l1(Test_random_ball):

    def setup_method(self):
        Test_random_ball.setup_method(self)
        self.p = 1


class Test_random_ball_l1_compiled(Test_random_ball_compiled):

    def setup_method(self):
        Test_random_ball_compiled.setup_method(self)
        self.p = 1

class Test_random_ball_l1_compiled_periodic(Test_random_ball_compiled_periodic):

    def setup_method(self):
        Test_random_ball_compiled_periodic.setup_method(self)
        self.p = 1


class Test_random_ball_linf(Test_random_ball):

    def setup_method(self):
        Test_random_ball.setup_method(self)
        self.p = np.inf

class Test_random_ball_linf_compiled_periodic(Test_random_ball_compiled_periodic):

    def setup_method(self):
        Test_random_ball_compiled_periodic.setup_method(self)
        self.p = np.inf


def test_random_ball_vectorized():

    n = 20
    m = 5
    T = KDTree(np.random.randn(n,m))

    r = T.query_ball_point(np.random.randn(2,3,m),1)
    assert_equal(r.shape,(2,3))
    assert_(isinstance(r[0,0],list))


def test_random_ball_vectorized_compiled():

    n = 20
    m = 5
    np.random.seed(1234)
    T = cKDTree(np.random.randn(n,m))

    r = T.query_ball_point(np.random.randn(2,3,m),1)
    assert_equal(r.shape,(2,3))
    assert_(isinstance(r[0,0],list))


def test_query_ball_point_multithreading():
    np.random.seed(0)
    n = 5000
    k = 2
    points = np.random.randn(n,k)
    T = cKDTree(points)
    l1 = T.query_ball_point(points,0.003,n_jobs=1)
    l2 = T.query_ball_point(points,0.003,n_jobs=64)
    l3 = T.query_ball_point(points,0.003,n_jobs=-1)

    for i in range(n):
        if l1[i] or l2[i]:
            assert_array_equal(l1[i],l2[i])

    for i in range(n):
        if l1[i] or l3[i]:
            assert_array_equal(l1[i],l3[i])


class two_trees_consistency:

    def distance(self, a, b, p):
        return minkowski_distance(a, b, p)

    def test_all_in_ball(self):
        r = self.T1.query_ball_tree(self.T2, self.d, p=self.p, eps=self.eps)
        for i, l in enumerate(r):
            for j in l:
                assert_(self.distance(self.data1[i],self.data2[j],self.p) <= self.d*(1.+self.eps))

    def test_found_all(self):
        r = self.T1.query_ball_tree(self.T2, self.d, p=self.p, eps=self.eps)
        for i, l in enumerate(r):
            c = np.ones(self.T2.n,dtype=bool)
            c[l] = False
            assert_(np.all(self.distance(self.data2[c],self.data1[i],self.p) >= self.d/(1.+self.eps)))


class Test_two_random_trees(two_trees_consistency):

    def setup_method(self):
        n = 50
        m = 4
        np.random.seed(1234)
        self.data1 = np.random.randn(n,m)
        self.T1 = KDTree(self.data1,leafsize=2)
        self.data2 = np.random.randn(n,m)
        self.T2 = KDTree(self.data2,leafsize=2)
        self.p = 2.
        self.eps = 0
        self.d = 0.2


class Test_two_random_trees_compiled(two_trees_consistency):

    def setup_method(self):
        n = 50
        m = 4
        np.random.seed(1234)
        self.data1 = np.random.randn(n,m)
        self.T1 = cKDTree(self.data1,leafsize=2)
        self.data2 = np.random.randn(n,m)
        self.T2 = cKDTree(self.data2,leafsize=2)
        self.p = 2.
        self.eps = 0
        self.d = 0.2

class Test_two_random_trees_compiled_periodic(two_trees_consistency):
    def distance(self, a, b, p):
        return distance_box(a, b, p, 1.0)

    def setup_method(self):
        n = 50
        m = 4
        np.random.seed(1234)
        self.data1 = np.random.uniform(size=(n,m))
        self.T1 = cKDTree(self.data1,leafsize=2, boxsize=1.0)
        self.data2 = np.random.uniform(size=(n,m))
        self.T2 = cKDTree(self.data2,leafsize=2, boxsize=1.0)
        self.p = 2.
        self.eps = 0
        self.d = 0.2

class Test_two_random_trees_far(Test_two_random_trees):

    def setup_method(self):
        Test_two_random_trees.setup_method(self)
        self.d = 2


class Test_two_random_trees_far_compiled(Test_two_random_trees_compiled):

    def setup_method(self):
        Test_two_random_trees_compiled.setup_method(self)
        self.d = 2

class Test_two_random_trees_far_compiled_periodic(Test_two_random_trees_compiled_periodic):

    def setup_method(self):
        Test_two_random_trees_compiled_periodic.setup_method(self)
        self.d = 2


class Test_two_random_trees_linf(Test_two_random_trees):

    def setup_method(self):
        Test_two_random_trees.setup_method(self)
        self.p = np.inf


class Test_two_random_trees_linf_compiled(Test_two_random_trees_compiled):

    def setup_method(self):
        Test_two_random_trees_compiled.setup_method(self)
        self.p = np.inf

class Test_two_random_trees_linf_compiled_periodic(Test_two_random_trees_compiled_periodic):

    def setup_method(self):
        Test_two_random_trees_compiled_periodic.setup_method(self)
        self.p = np.inf


class Test_rectangle:

    def setup_method(self):
        self.rect = Rectangle([0,0],[1,1])

    def test_min_inside(self):
        assert_almost_equal(self.rect.min_distance_point([0.5,0.5]),0)

    def test_min_one_side(self):
        assert_almost_equal(self.rect.min_distance_point([0.5,1.5]),0.5)

    def test_min_two_sides(self):
        assert_almost_equal(self.rect.min_distance_point([2,2]),np.sqrt(2))

    def test_max_inside(self):
        assert_almost_equal(self.rect.max_distance_point([0.5,0.5]),1/np.sqrt(2))

    def test_max_one_side(self):
        assert_almost_equal(self.rect.max_distance_point([0.5,1.5]),np.hypot(0.5,1.5))

    def test_max_two_sides(self):
        assert_almost_equal(self.rect.max_distance_point([2,2]),2*np.sqrt(2))

    def test_split(self):
        less, greater = self.rect.split(0,0.1)
        assert_array_equal(less.maxes,[0.1,1])
        assert_array_equal(less.mins,[0,0])
        assert_array_equal(greater.maxes,[1,1])
        assert_array_equal(greater.mins,[0.1,0])


def test_distance_l2():
    assert_almost_equal(minkowski_distance([0,0],[1,1],2),np.sqrt(2))


def test_distance_l1():
    assert_almost_equal(minkowski_distance([0,0],[1,1],1),2)


def test_distance_linf():
    assert_almost_equal(minkowski_distance([0,0],[1,1],np.inf),1)


def test_distance_vectorization():
    np.random.seed(1234)
    x = np.random.randn(10,1,3)
    y = np.random.randn(1,7,3)
    assert_equal(minkowski_distance(x,y).shape,(10,7))


class count_neighbors_consistency:
    def test_one_radius(self):
        r = 0.2
        assert_equal(self.T1.count_neighbors(self.T2, r),
                np.sum([len(l) for l in self.T1.query_ball_tree(self.T2,r)]))

    def test_large_radius(self):
        r = 1000
        assert_equal(self.T1.count_neighbors(self.T2, r),
                np.sum([len(l) for l in self.T1.query_ball_tree(self.T2,r)]))

    def test_multiple_radius(self):
        rs = np.exp(np.linspace(np.log(0.01),np.log(10),3))
        results = self.T1.count_neighbors(self.T2, rs)
        assert_(np.all(np.diff(results) >= 0))
        for r,result in zip(rs, results):
            assert_equal(self.T1.count_neighbors(self.T2, r), result)

class Test_count_neighbors(count_neighbors_consistency):

    def setup_method(self):
        n = 50
        m = 2
        np.random.seed(1234)
        self.T1 = KDTree(np.random.randn(n,m),leafsize=2)
        self.T2 = KDTree(np.random.randn(n,m),leafsize=2)


class Test_count_neighbors_compiled(count_neighbors_consistency):

    def setup_method(self):
        n = 50
        m = 2
        np.random.seed(1234)
        self.T1 = cKDTree(np.random.randn(n,m),leafsize=2)
        self.T2 = cKDTree(np.random.randn(n,m),leafsize=2)


class sparse_distance_matrix_consistency:

    def distance(self, a, b, p):
        return minkowski_distance(a, b, p)

    def test_consistency_with_neighbors(self):
        M = self.T1.sparse_distance_matrix(self.T2, self.r)
        r = self.T1.query_ball_tree(self.T2, self.r)
        for i,l in enumerate(r):
            for j in l:
                assert_almost_equal(M[i,j],
                                    self.distance(self.T1.data[i], self.T2.data[j], self.p),
                                    decimal=14)
        for ((i,j),d) in M.items():
            assert_(j in r[i])

    def test_zero_distance(self):
        # raises an exception for bug 870 (FIXME: Does it?)
        self.T1.sparse_distance_matrix(self.T1, self.r)

class Test_sparse_distance_matrix(sparse_distance_matrix_consistency):

    def setup_method(self):
        n = 50
        m = 4
        np.random.seed(1234)
        data1 = np.random.randn(n,m)
        data2 = np.random.randn(n,m)
        self.T1 = cKDTree(data1,leafsize=2)
        self.T2 = cKDTree(data2,leafsize=2)
        self.r = 0.5
        self.p = 2
        self.data1 = data1
        self.data2 = data2
        self.n = n
        self.m = m

class Test_sparse_distance_matrix_compiled(sparse_distance_matrix_consistency):

    def setup_method(self):
        n = 50
        m = 4
        np.random.seed(0)
        data1 = np.random.randn(n,m)
        data2 = np.random.randn(n,m)
        self.T1 = cKDTree(data1,leafsize=2)
        self.T2 = cKDTree(data2,leafsize=2)
        self.ref_T1 = KDTree(data1, leafsize=2)
        self.ref_T2 = KDTree(data2, leafsize=2)
        self.r = 0.5
        self.n = n
        self.m = m
        self.data1 = data1
        self.data2 = data2
        self.p = 2

    def test_consistency_with_python(self):
        M1 = self.T1.sparse_distance_matrix(self.T2, self.r)
        M2 = self.ref_T1.sparse_distance_matrix(self.ref_T2, self.r)
        assert_array_almost_equal(M1.todense(), M2.todense(), decimal=14)

    def test_against_logic_error_regression(self):
        # regression test for gh-5077 logic error
        np.random.seed(0)
        too_many = np.array(np.random.randn(18, 2), dtype=int)
        tree = cKDTree(too_many, balanced_tree=False, compact_nodes=False)
        d = tree.sparse_distance_matrix(tree, 3).todense()
        assert_array_almost_equal(d, d.T, decimal=14)

    def test_ckdtree_return_types(self):
        # brute-force reference
        ref = np.zeros((self.n,self.n))
        for i in range(self.n):
            for j in range(self.n):
                v = self.data1[i,:] - self.data2[j,:]
                ref[i,j] = np.dot(v,v)
        ref = np.sqrt(ref)
        ref[ref > self.r] = 0.
        # test return type 'dict'
        dist = np.zeros((self.n,self.n))
        r = self.T1.sparse_distance_matrix(self.T2, self.r, output_type='dict')
        for i,j in r.keys():
            dist[i,j] = r[(i,j)]
        assert_array_almost_equal(ref, dist, decimal=14)
        # test return type 'ndarray'
        dist = np.zeros((self.n,self.n))
        r = self.T1.sparse_distance_matrix(self.T2, self.r,
            output_type='ndarray')
        for k in range(r.shape[0]):
            i = r['i'][k]
            j = r['j'][k]
            v = r['v'][k]
            dist[i,j] = v
        assert_array_almost_equal(ref, dist, decimal=14)
        # test return type 'dok_matrix'
        r = self.T1.sparse_distance_matrix(self.T2, self.r,
            output_type='dok_matrix')
        assert_array_almost_equal(ref, r.todense(), decimal=14)
        # test return type 'coo_matrix'
        r = self.T1.sparse_distance_matrix(self.T2, self.r,
            output_type='coo_matrix')
        assert_array_almost_equal(ref, r.todense(), decimal=14)


def test_distance_matrix():
    m = 10
    n = 11
    k = 4
    np.random.seed(1234)
    xs = np.random.randn(m,k)
    ys = np.random.randn(n,k)
    ds = distance_matrix(xs,ys)
    assert_equal(ds.shape, (m,n))
    for i in range(m):
        for j in range(n):
            assert_almost_equal(minkowski_distance(xs[i],ys[j]),ds[i,j])


def test_distance_matrix_looping():
    m = 10
    n = 11
    k = 4
    np.random.seed(1234)
    xs = np.random.randn(m,k)
    ys = np.random.randn(n,k)
    ds = distance_matrix(xs,ys)
    dsl = distance_matrix(xs,ys,threshold=1)
    assert_equal(ds,dsl)


def check_onetree_query(T,d):
    r = T.query_ball_tree(T, d)
    s = set()
    for i, l in enumerate(r):
        for j in l:
            if i < j:
                s.add((i,j))

    assert_(s == T.query_pairs(d))

def test_onetree_query():
    np.random.seed(0)
    n = 50
    k = 4
    points = np.random.randn(n,k)
    T = KDTree(points)
    check_onetree_query(T, 0.1)

    points = np.random.randn(3*n,k)
    points[:n] *= 0.001
    points[n:2*n] += 2
    T = KDTree(points)
    check_onetree_query(T, 0.1)
    check_onetree_query(T, 0.001)
    check_onetree_query(T, 0.00001)
    check_onetree_query(T, 1e-6)


def test_onetree_query_compiled():
    np.random.seed(0)
    n = 100
    k = 4
    points = np.random.randn(n,k)
    T = cKDTree(points)
    check_onetree_query(T, 0.1)

    points = np.random.randn(3*n,k)
    points[:n] *= 0.001
    points[n:2*n] += 2
    T = cKDTree(points)
    check_onetree_query(T, 0.1)
    check_onetree_query(T, 0.001)
    check_onetree_query(T, 0.00001)
    check_onetree_query(T, 1e-6)


def test_query_pairs_single_node():
    tree = KDTree([[0, 1]])
    assert_equal(tree.query_pairs(0.5), set())


def test_query_pairs_single_node_compiled():
    tree = cKDTree([[0, 1]])
    assert_equal(tree.query_pairs(0.5), set())


def test_ckdtree_query_pairs():
    np.random.seed(0)
    n = 50
    k = 2
    r = 0.1
    r2 = r**2
    points = np.random.randn(n,k)
    T = cKDTree(points)
    # brute force reference
    brute = set()
    for i in range(n):
        for j in range(i+1,n):
            v = points[i,:] - points[j,:]
            if np.dot(v,v) <= r2:
                brute.add((i,j))
    l0 = sorted(brute)
    # test default return type
    s = T.query_pairs(r)
    l1 = sorted(s)
    assert_array_equal(l0,l1)
    # test return type 'set'
    s = T.query_pairs(r, output_type='set')
    l1 = sorted(s)
    assert_array_equal(l0,l1)
    # test return type 'ndarray'
    s = set()
    arr = T.query_pairs(r, output_type='ndarray')
    for i in range(arr.shape[0]):
        s.add((int(arr[i,0]),int(arr[i,1])))
    l2 = sorted(s)
    assert_array_equal(l0,l2)


def test_ball_point_ints():
    # Regression test for #1373.
    x, y = np.mgrid[0:4, 0:4]
    points = list(zip(x.ravel(), y.ravel()))
    tree = KDTree(points)
    assert_equal(sorted([4, 8, 9, 12]),
                 sorted(tree.query_ball_point((2, 0), 1)))
    points = np.asarray(points, dtype=float)
    tree = KDTree(points)
    assert_equal(sorted([4, 8, 9, 12]),
                 sorted(tree.query_ball_point((2, 0), 1)))


def test_kdtree_comparisons():
    # Regression test: node comparisons were done wrong in 0.12 w/Py3.
    nodes = [KDTree.node() for _ in range(3)]
    assert_equal(sorted(nodes), sorted(nodes[::-1]))


def test_ckdtree_build_modes():
    # check if different build modes for cKDTree give
    # similar query results
    np.random.seed(0)
    n = 5000
    k = 4
    points = np.random.randn(n, k)
    T1 = cKDTree(points).query(points, k=5)[-1]
    T2 = cKDTree(points, compact_nodes=False).query(points, k=5)[-1]
    T3 = cKDTree(points, balanced_tree=False).query(points, k=5)[-1]
    T4 = cKDTree(points, compact_nodes=False, balanced_tree=False).query(points, k=5)[-1]
    assert_array_equal(T1, T2)
    assert_array_equal(T1, T3)
    assert_array_equal(T1, T4)

def test_ckdtree_pickle():
    # test if it is possible to pickle
    # a cKDTree
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    np.random.seed(0)
    n = 50
    k = 4
    points = np.random.randn(n, k)
    T1 = cKDTree(points)
    tmp = pickle.dumps(T1)
    T2 = pickle.loads(tmp)
    T1 = T1.query(points, k=5)[-1]
    T2 = T2.query(points, k=5)[-1]
    assert_array_equal(T1, T2)

def test_ckdtree_pickle_boxsize():
    # test if it is possible to pickle a periodic
    # cKDTree
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    np.random.seed(0)
    n = 50
    k = 4
    points = np.random.uniform(size=(n, k))
    T1 = cKDTree(points, boxsize=1.0)
    tmp = pickle.dumps(T1)
    T2 = pickle.loads(tmp)
    T1 = T1.query(points, k=5)[-1]
    T2 = T2.query(points, k=5)[-1]
    assert_array_equal(T1, T2)

def test_ckdtree_copy_data():
    # check if copy_data=True makes the kd-tree
    # impervious to data corruption by modification of
    # the data arrray
    np.random.seed(0)
    n = 5000
    k = 4
    points = np.random.randn(n, k)
    T = cKDTree(points, copy_data=True)
    q = points.copy()
    T1 = T.query(q, k=5)[-1]
    points[...] = np.random.randn(n, k)
    T2 = T.query(q, k=5)[-1]
    assert_array_equal(T1, T2)

def test_ckdtree_parallel():
    # check if parallel=True also generates correct
    # query results
    np.random.seed(0)
    n = 5000
    k = 4
    points = np.random.randn(n, k)
    T = cKDTree(points)
    T1 = T.query(points, k=5, n_jobs=64)[-1]
    T2 = T.query(points, k=5, n_jobs=-1)[-1]
    T3 = T.query(points, k=5)[-1]
    assert_array_equal(T1, T2)
    assert_array_equal(T1, T3)

def test_ckdtree_view():
    # Check that the nodes can be correctly viewed from Python.
    # This test also sanity checks each node in the cKDTree, and
    # thus verifies the internal structure of the kd-tree.
    np.random.seed(0)
    n = 100
    k = 4
    points = np.random.randn(n, k)
    kdtree = cKDTree(points)

    # walk the whole kd-tree and sanity check each node
    def recurse_tree(n):
        assert_(isinstance(n, cKDTreeNode))
        if n.split_dim == -1:
            assert_(n.lesser is None)
            assert_(n.greater is None)
            assert_(n.indices.shape[0] <= kdtree.leafsize)
        else:
            recurse_tree(n.lesser)
            recurse_tree(n.greater)
            x = n.lesser.data_points[:, n.split_dim]
            y = n.greater.data_points[:, n.split_dim]
            assert_(x.max() < y.min())

    recurse_tree(kdtree.tree)
    # check that indices are correctly retrieved
    n = kdtree.tree
    assert_array_equal(np.sort(n.indices), range(100))
    # check that data_points are correctly retrieved
    assert_array_equal(kdtree.data[n.indices, :], n.data_points)

# cKDTree is specialized to type double points, so no need to make
# a unit test corresponding to test_ball_point_ints()

def test_ckdtree_list_k():
    # check ckdtree periodic boundary
    n = 200
    m = 2
    klist = [1, 2, 3]
    kint = 3

    np.random.seed(1234)
    data = np.random.uniform(size=(n, m))
    kdtree = cKDTree(data, leafsize=1)

    # check agreement between arange(1,k+1) and k
    dd, ii = kdtree.query(data, klist)
    dd1, ii1 = kdtree.query(data, kint)
    assert_equal(dd, dd1)
    assert_equal(ii, ii1)

    # now check skipping one element
    klist = np.array([1, 3])
    kint = 3
    dd, ii = kdtree.query(data, kint)
    dd1, ii1 = kdtree.query(data, klist)
    assert_equal(dd1, dd[..., klist - 1])
    assert_equal(ii1, ii[..., klist - 1])

    # check k == 1 special case
    # and k == [1] non-special case
    dd, ii = kdtree.query(data, 1)
    dd1, ii1 = kdtree.query(data, [1])
    assert_equal(len(dd.shape), 1)
    assert_equal(len(dd1.shape), 2)
    assert_equal(dd, np.ravel(dd1))
    assert_equal(ii, np.ravel(ii1))

def test_ckdtree_box():
    # check ckdtree periodic boundary
    n = 2000
    m = 3
    k = 3
    np.random.seed(1234)
    data = np.random.uniform(size=(n, m))
    kdtree = cKDTree(data, leafsize=1, boxsize=1.0)

    # use the standard python KDTree for the simulated periodic box
    kdtree2 = cKDTree(data, leafsize=1)

    for p in [1, 2, 3.0, np.inf]:
        dd, ii = kdtree.query(data, k, p=p)

        dd1, ii1 = kdtree.query(data + 1.0, k, p=p)
        assert_almost_equal(dd, dd1)
        assert_equal(ii, ii1)

        dd1, ii1 = kdtree.query(data - 1.0, k, p=p)
        assert_almost_equal(dd, dd1)
        assert_equal(ii, ii1)

        dd2, ii2 = simulate_periodic_box(kdtree2, data, k, boxsize=1.0, p=p)
        assert_almost_equal(dd, dd2)
        assert_equal(ii, ii2)

def test_ckdtree_box_0boxsize():
    # check ckdtree periodic boundary that mimics non-periodic
    n = 2000
    m = 2
    k = 3
    np.random.seed(1234)
    data = np.random.uniform(size=(n, m))
    kdtree = cKDTree(data, leafsize=1, boxsize=0.0)

    # use the standard python KDTree for the simulated periodic box
    kdtree2 = cKDTree(data, leafsize=1)

    for p in [1, 2, np.inf]:
        dd, ii = kdtree.query(data, k, p=p)

        dd1, ii1 = kdtree2.query(data, k, p=p)
        assert_almost_equal(dd, dd1)
        assert_equal(ii, ii1)

def test_ckdtree_box_upper_bounds():
    data = np.linspace(0, 2, 10).reshape(-1, 2)
    data[:, 1] += 10
    assert_raises(ValueError, cKDTree, data, leafsize=1, boxsize=1.0)
    assert_raises(ValueError, cKDTree, data, leafsize=1, boxsize=(0.0, 2.0))
    # skip a dimension.
    cKDTree(data, leafsize=1, boxsize=(2.0, 0.0))

def test_ckdtree_box_lower_bounds():
    data = np.linspace(-1, 1, 10)
    assert_raises(ValueError, cKDTree, data, leafsize=1, boxsize=1.0)

def simulate_periodic_box(kdtree, data, k, boxsize, p):
    dd = []
    ii = []
    x = np.arange(3 ** data.shape[1])
    nn = np.array(np.unravel_index(x, [3] * data.shape[1])).T
    nn = nn - 1.0
    for n in nn:
        image = data + n * 1.0 * boxsize
        dd2, ii2 = kdtree.query(image, k, p=p)
        dd2 = dd2.reshape(-1, k)
        ii2 = ii2.reshape(-1, k)
        dd.append(dd2)
        ii.append(ii2)
    dd = np.concatenate(dd, axis=-1)
    ii = np.concatenate(ii, axis=-1)

    result = np.empty([len(data), len(nn) * k], dtype=[
            ('ii', 'i8'),
            ('dd', 'f8')])
    result['ii'][:] = ii
    result['dd'][:] = dd
    result.sort(order='dd')
    return result['dd'][:, :k], result['ii'][:,:k]


@pytest.mark.skipif(python_implementation() == 'PyPy',
                    reason="Fails on PyPy CI runs. See #9507")
def test_ckdtree_memuse():
    # unit test adaptation of gh-5630

    # NOTE: this will fail when run via valgrind,
    # because rss is no longer a reliable memory usage indicator.

    try:
        import resource
    except ImportError:
        # resource is not available on Windows with Python 2.6
        return
    # Make some data
    dx, dy = 0.05, 0.05
    y, x = np.mgrid[slice(1, 5 + dy, dy),
                    slice(1, 5 + dx, dx)]
    z = np.sin(x)**10 + np.cos(10 + y*x) * np.cos(x)
    z_copy = np.empty_like(z)
    z_copy[:] = z
    # Place FILLVAL in z_copy at random number of random locations
    FILLVAL = 99.
    mask = np.random.randint(0, z.size, np.random.randint(50) + 5)
    z_copy.flat[mask] = FILLVAL
    igood = np.vstack(np.nonzero(x != FILLVAL)).T
    ibad = np.vstack(np.nonzero(x == FILLVAL)).T
    mem_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # burn-in
    for i in range(10):
        tree = cKDTree(igood)
    # count memleaks while constructing and querying cKDTree
    num_leaks = 0
    for i in range(100):
        mem_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        tree = cKDTree(igood)
        dist, iquery = tree.query(ibad, k=4, p=2)
        new_mem_use = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if new_mem_use > mem_use:
            num_leaks += 1
    # ideally zero leaks, but errors might accidentally happen
    # outside cKDTree
    assert_(num_leaks < 10)

def test_ckdtree_weights():

    data = np.linspace(0, 1, 4).reshape(-1, 1)
    tree1 = cKDTree(data, leafsize=1)
    weights = np.ones(len(data), dtype='f4')

    nw = tree1._build_weights(weights)
    assert_array_equal(nw, [4, 2, 1, 1, 2, 1, 1])

    assert_raises(ValueError, tree1._build_weights, weights[:-1])

    for i in range(10):
        # since weights are uniform, these shall agree:
        c1 = tree1.count_neighbors(tree1, np.linspace(0, 10, i))
        c2 = tree1.count_neighbors(tree1, np.linspace(0, 10, i),
                weights=(weights, weights))
        c3 = tree1.count_neighbors(tree1, np.linspace(0, 10, i),
                weights=(weights, None))
        c4 = tree1.count_neighbors(tree1, np.linspace(0, 10, i),
                weights=(None, weights))
        c5 = tree1.count_neighbors(tree1, np.linspace(0, 10, i),
                weights=weights)

        assert_array_equal(c1, c2)
        assert_array_equal(c1, c3)
        assert_array_equal(c1, c4)

    for i in range(len(data)):
        # this tests removal of one data point by setting weight to 0
        w1 = weights.copy()
        w1[i] = 0
        data2 = data[w1 != 0]
        w2 = weights[w1 != 0]
        tree2 = cKDTree(data2)

        c1 = tree1.count_neighbors(tree1, np.linspace(0, 10, 100),
                weights=(w1, w1))
        # "c2 is correct"
        c2 = tree2.count_neighbors(tree2, np.linspace(0, 10, 100))

        assert_array_equal(c1, c2)

        #this asserts for two different trees, singular weights
        # crashes
        assert_raises(ValueError, tree1.count_neighbors,
            tree2, np.linspace(0, 10, 100), weights=w1)

def test_ckdtree_count_neighbous_multiple_r():
    n = 2000
    m = 2
    np.random.seed(1234)
    data = np.random.normal(size=(n, m))
    kdtree = cKDTree(data, leafsize=1)
    r0 = [0, 0.01, 0.01, 0.02, 0.05]
    i0 = np.arange(len(r0))
    n0 = kdtree.count_neighbors(kdtree, r0)
    nnc = kdtree.count_neighbors(kdtree, r0, cumulative=False)
    assert_equal(n0, nnc.cumsum())

    for i, r in zip(itertools.permutations(i0),
                    itertools.permutations(r0)):
        # permute n0 by i and it shall agree
        n = kdtree.count_neighbors(kdtree, r)
        assert_array_equal(n, n0[list(i)])

def test_len0_arrays():
    # make sure len-0 arrays are handled correctly
    # in range queries (gh-5639)
    np.random.seed(1234)
    X = np.random.rand(10,2)
    Y = np.random.rand(10,2)
    tree = cKDTree(X)
    # query_ball_point (single)
    d,i = tree.query([.5, .5], k=1)
    z = tree.query_ball_point([.5, .5], 0.1*d)
    assert_array_equal(z, [])
    # query_ball_point (multiple)
    d,i = tree.query(Y, k=1)
    mind = d.min()
    z = tree.query_ball_point(Y, 0.1*mind)
    y = np.empty(shape=(10,), dtype=object)
    y.fill([])
    assert_array_equal(y, z)
    # query_ball_tree
    other = cKDTree(Y)
    y = tree.query_ball_tree(other, 0.1*mind)
    assert_array_equal(10*[[]], y)
    # count_neighbors
    y = tree.count_neighbors(other, 0.1*mind)
    assert_(y == 0)
    # sparse_distance_matrix
    y = tree.sparse_distance_matrix(other, 0.1*mind, output_type='dok_matrix')
    assert_array_equal(y == np.zeros((10,10)), True)
    y = tree.sparse_distance_matrix(other, 0.1*mind, output_type='coo_matrix')
    assert_array_equal(y == np.zeros((10,10)), True)
    y = tree.sparse_distance_matrix(other, 0.1*mind, output_type='dict')
    assert_equal(y, {})
    y = tree.sparse_distance_matrix(other,0.1*mind, output_type='ndarray')
    _dtype = [('i',np.intp), ('j',np.intp), ('v',np.float64)]
    res_dtype = np.dtype(_dtype, align=True)
    z = np.empty(shape=(0,), dtype=res_dtype)
    assert_array_equal(y, z)
    # query_pairs
    d,i = tree.query(X, k=2)
    mind = d[:,-1].min()
    y = tree.query_pairs(0.1*mind, output_type='set')
    assert_equal(y, set())
    y = tree.query_pairs(0.1*mind, output_type='ndarray')
    z = np.empty(shape=(0,2), dtype=np.intp)
    assert_array_equal(y, z)

def test_ckdtree_duplicated_inputs():
    # check ckdtree with duplicated inputs
    n = 1024
    for m in range(1, 8):
        data = np.concatenate([
            np.ones((n // 2, m)) * 1,
            np.ones((n // 2, m)) * 2], axis=0)

        # it shall not divide more than 3 nodes.
        # root left (1), and right (2)
        kdtree = cKDTree(data, leafsize=1)
        assert_equal(kdtree.size, 3)

        kdtree = cKDTree(data)
        assert_equal(kdtree.size, 3)

        # if compact_nodes are disabled, the number
        # of nodes is n (per leaf) + (m - 1)* 2 (splits per dimension) + 1
        # and the root
        kdtree = cKDTree(data, compact_nodes=False, leafsize=1)
        assert_equal(kdtree.size, n + m * 2 - 1)

def test_ckdtree_noncumulative_nondecreasing():
    # check ckdtree with duplicated inputs

    # it shall not divide more than 3 nodes.
    # root left (1), and right (2)
    kdtree = cKDTree([[0]], leafsize=1)

    assert_raises(ValueError, kdtree.count_neighbors,
        kdtree, [0.1, 0], cumulative=False)

def test_short_knn():

    # The test case is based on github: #6425 by @SteveDoyle2

    xyz = np.array([
        [0., 0., 0.],
        [1.01, 0., 0.],
        [0., 1., 0.],
        [0., 1.01, 0.],
        [1., 0., 0.],
        [1., 1., 0.],],
    dtype='float64')

    ckdt = cKDTree(xyz)

    deq, ieq = ckdt.query(xyz, k=4, distance_upper_bound=0.2)

    assert_array_almost_equal(deq,
            [[0., np.inf, np.inf, np.inf],
            [0., 0.01, np.inf, np.inf],
            [0., 0.01, np.inf, np.inf],
            [0., 0.01, np.inf, np.inf],
            [0., 0.01, np.inf, np.inf],
            [0., np.inf, np.inf, np.inf]])

def test_query_ball_point_vector_r():

    np.random.seed(1234)
    data = np.random.normal(size=(100, 3))
    query = np.random.normal(size=(100, 3))
    tree = cKDTree(data)
    d = np.random.uniform(0, 0.3, size=len(query))

    rvector = tree.query_ball_point(query, d)
    rscalar = [tree.query_ball_point(qi, di) for qi, di in zip(query, d)]
    for a, b in zip(rvector, rscalar):
        assert_array_equal(sorted(a), sorted(b))

def test_query_ball_point_length():

    np.random.seed(1234)
    data = np.random.normal(size=(100, 3))
    query = np.random.normal(size=(100, 3))
    tree = cKDTree(data)
    d = 0.3

    length = tree.query_ball_point(query, d, return_length=True)
    length2 = [len(ind) for ind in tree.query_ball_point(query, d, return_length=False)]
    length3 = [len(tree.query_ball_point(qi, d)) for qi in query]
    length4 = [tree.query_ball_point(qi, d, return_length=True) for qi in query]
    assert_array_equal(length, length2)
    assert_array_equal(length, length3)
    assert_array_equal(length, length4)

class Test_sorted_query_ball_point(object):

    def setup_method(self):
        np.random.seed(1234)
        self.x = np.random.randn(100, 1)
        self.ckdt = cKDTree(self.x)

    def test_return_sorted_True(self):
        idxs_list = self.ckdt.query_ball_point(self.x, 1., return_sorted=True)
        for idxs in idxs_list:
            assert_array_equal(idxs, sorted(idxs))

        for xi in self.x:
            idxs = self.ckdt.query_ball_point(xi, 1., return_sorted=True)
            assert_array_equal(idxs, sorted(idxs))

    def test_return_sorted_None(self):
        """Previous behavior was to sort the returned indices if there were
        multiple points per query but not sort them if there was a single point
        per query."""
        idxs_list = self.ckdt.query_ball_point(self.x, 1.)
        for idxs in idxs_list:
            assert_array_equal(idxs, sorted(idxs))

        idxs_list_single = [self.ckdt.query_ball_point(xi, 1.) for xi in self.x]
        idxs_list_False = self.ckdt.query_ball_point(self.x, 1., return_sorted=False)
        for idxs0, idxs1 in zip(idxs_list_False, idxs_list_single):
            assert_array_equal(idxs0, idxs1)
