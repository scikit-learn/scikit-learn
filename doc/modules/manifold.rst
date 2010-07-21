=================
Manifold Learning
=================

Introduction
============

Manifold learning is a set of dimensionality reduction techniques. They try to
find an embedded space where the original data can be efficiently described.
Then, other tools can create a mapping between the original space and the
embedded one so as to allow classification, for instance, in the embedded
space.

Dimensionality Reduction
========================

Dimensionality reduction is the set of tools that will actually search the
embedded space. Each technique has its advantages and drawbacks, so you will
have to select the correct one for your purpose.

Multidimensional Scaling (MDS)
------------------------------

MSD is a technique that tries to find a lower dimension space where input
distances are equals to the distances in the embedded space. Input distances
are usually Euclidian distances and in this case, MDS is tantamount to PCA.

Isomap
------

Isomap[1] embeds data in a lower dimension space by computing geodesic
distances on the manifold and then applying MDS on it. Geodesic distances are
approximated by distances on a neighbors graph.

Locally Linear Embedding (LLE)
------------------------------

Laplacian Eigenmaps
-------------------

Diffusion Maps
--------------

NonLinear Mapping (NLM)
-----------------------


Embedding New Data
==================

Barycenter
----------

Kernel Projection
-----------------

Linear Mapping
--------------

Piecewise Linear Mapping
------------------------

Examples
========

Notes
=====
    .. [1] Tenenbaum, J. B., de Silva, V. and Langford, J. C.,
           "A Global Geometric Framework for Nonlinear Dimensionality 
           Reduction,"
           Science, 290(5500), pp. 2319-2323, 2000
