# 3dproj.py
#
"""
Various transforms used for by the 3D code
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import zip

import numpy as np
import numpy.linalg as linalg



def line2d(p0, p1):
    """
    Return 2D equation of line in the form ax+by+c = 0
    """
    # x + x1  = 0
    x0, y0 = p0[:2]
    x1, y1 = p1[:2]
    #
    if x0 == x1:
        a = -1
        b = 0
        c = x1
    elif y0 == y1:
        a = 0
        b = 1
        c = -y1
    else:
        a = (y0-y1)
        b = (x0-x1)
        c = (x0*y1 - x1*y0)
    return a, b, c

def line2d_dist(l, p):
    """
    Distance from line to point
    line is a tuple of coefficients a,b,c
    """
    a, b, c = l
    x0, y0 = p
    return abs((a*x0 + b*y0 + c)/np.sqrt(a**2+b**2))


def line2d_seg_dist(p1, p2, p0):
    """distance(s) from line defined by p1 - p2 to point(s) p0

    p0[0] = x(s)
    p0[1] = y(s)

    intersection point p = p1 + u*(p2-p1)
    and intersection point lies within segment if u is between 0 and 1
    """

    x21 = p2[0] - p1[0]
    y21 = p2[1] - p1[1]
    x01 = np.asarray(p0[0]) - p1[0]
    y01 = np.asarray(p0[1]) - p1[1]

    u = (x01*x21 + y01*y21) / (x21**2 + y21**2)
    u = np.clip(u, 0, 1)
    d = np.sqrt((x01 - u*x21)**2 + (y01 - u*y21)**2)

    return d


def mod(v):
    """3d vector length"""
    return np.sqrt(v[0]**2+v[1]**2+v[2]**2)

def world_transformation(xmin, xmax,
                         ymin, ymax,
                         zmin, zmax):
    dx, dy, dz = (xmax-xmin), (ymax-ymin), (zmax-zmin)
    return np.array([
        [1.0/dx,0,0,-xmin/dx],
        [0,1.0/dy,0,-ymin/dy],
        [0,0,1.0/dz,-zmin/dz],
        [0,0,0,1.0]])


def view_transformation(E, R, V):
    n = (E - R)
    ## new
#    n /= mod(n)
#    u = np.cross(V,n)
#    u /= mod(u)
#    v = np.cross(n,u)
#    Mr = np.diag([1.]*4)
#    Mt = np.diag([1.]*4)
#    Mr[:3,:3] = u,v,n
#    Mt[:3,-1] = -E
    ## end new

    ## old
    n = n / mod(n)
    u = np.cross(V, n)
    u = u / mod(u)
    v = np.cross(n, u)
    Mr = [[u[0],u[1],u[2],0],
          [v[0],v[1],v[2],0],
          [n[0],n[1],n[2],0],
          [0,   0,   0,   1],
          ]
    #
    Mt = [[1, 0, 0, -E[0]],
          [0, 1, 0, -E[1]],
          [0, 0, 1, -E[2]],
          [0, 0, 0, 1]]
    ## end old

    return np.dot(Mr, Mt)

def persp_transformation(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)
    return np.array([[1,0,0,0],
                     [0,1,0,0],
                     [0,0,a,b],
                     [0,0,-1,0]
                     ])

def ortho_transformation(zfront, zback):
    # note: w component in the resulting vector will be (zback-zfront), not 1
    a = -(zfront + zback)
    b = -(zfront - zback)
    return np.array([[2,0,0,0],
                     [0,2,0,0],
                     [0,0,-2,0],
                     [0,0,a,b]
                     ])

def proj_transform_vec(vec, M):
    vecw = np.dot(M, vec)
    w = vecw[3]
    # clip here..
    txs, tys, tzs = vecw[0]/w, vecw[1]/w, vecw[2]/w
    return txs, tys, tzs

def proj_transform_vec_clip(vec, M):
    vecw = np.dot(M, vec)
    w = vecw[3]
    # clip here.
    txs, tys, tzs = vecw[0] / w, vecw[1] / w, vecw[2] / w
    tis = (0 <= vecw[0]) & (vecw[0] <= 1) & (0 <= vecw[1]) & (vecw[1] <= 1)
    if np.any(tis):
        tis = vecw[1] < 1
    return txs, tys, tzs, tis

def inv_transform(xs, ys, zs, M):
    iM = linalg.inv(M)
    vec = vec_pad_ones(xs, ys, zs)
    vecr = np.dot(iM, vec)
    try:
        vecr = vecr/vecr[3]
    except OverflowError:
        pass
    return vecr[0], vecr[1], vecr[2]

def vec_pad_ones(xs, ys, zs):
    return np.array([xs, ys, zs, np.ones_like(xs)])

def proj_transform(xs, ys, zs, M):
    """
    Transform the points by the projection matrix
    """
    vec = vec_pad_ones(xs, ys, zs)
    return proj_transform_vec(vec, M)

def proj_transform_clip(xs, ys, zs, M):
    """
    Transform the points by the projection matrix
    and return the clipping result
    returns txs,tys,tzs,tis
    """
    vec = vec_pad_ones(xs, ys, zs)
    return proj_transform_vec_clip(vec, M)
transform = proj_transform

def proj_points(points, M):
    return np.column_stack(proj_trans_points(points, M))

def proj_trans_points(points, M):
    xs, ys, zs = zip(*points)
    return proj_transform(xs, ys, zs, M)

def proj_trans_clip_points(points, M):
    xs, ys, zs = zip(*points)
    return proj_transform_clip(xs, ys, zs, M)


def rot_x(V, alpha):
    cosa, sina = np.cos(alpha), np.sin(alpha)
    M1 = np.array([[1,0,0,0],
                   [0,cosa,-sina,0],
                   [0,sina,cosa,0],
                   [0,0,0,1]])

    return np.dot(M1, V)
