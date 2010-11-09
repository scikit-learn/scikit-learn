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

Not yet implemented.

MSD is a technique that tries to find a lower dimension space where input
distances are equals to the distances in the embedded space. Input distances
are usually Euclidian distances and in this case, MDS is tantamount to PCA.

Isomap
------

Not yet implemented.

Isomap[1] embeds data in a lower dimension space by computing geodesic
distances on the manifold and then applying MDS on it. Geodesic distances are
approximated by distances on a neighbors graph.

Locally Linear Embedding (LLE)
------------------------------

LLE[2] creates an embedded space by computing barycenter weights in the
original space and by trying to reproduce them in the embedded space.

Laplacian Eigenmaps
-------------------

Laplacian Eigenmaps[3] are an algorithm that will reproduce diffusion distances
inside the reduced space. Only the closest distances will be taken into account.

Diffusion Maps
--------------

Diffusion Maps[4] are based on the same principle than Laplacian Eigenmaps,
safe that all distances are used and are whittened before being processed.

Hessian Eigenmaps
-----------------

Hessian Maps[5] (also known as Hessian LLE) solve a quadratic form that averages
the Froebenius norm of the Hessian of the manifold.

NonLinear Mapping (NLM)
-----------------------

Not yet implemented.

NLM[6] is one of the first embedding algorithms. Its goal is to minimize a
quadratic cost function weighted with the inverse of the original distances.
This allows stable estimation of an embedding space.

Curvilinear component analysis (CCA)
------------------------------------

Not yet implemented.

Originally implemented as a Self Organizing Neural Network, CCA[7] tries to
find an embedding space where small distances are preserved. It is highly
instable as badly estimated distances (distances in the embedding space far
greater than in the original space) lead to a null cost.

Robust Embedding
----------------

Not yet implemented.

Embedding New Data
==================

Embedding new data is done through mapping the original space into the embedding
space. Several algorithms are available, depending on whether the original data
can be kept and on whether the associated kernel[9] is known.

Barycenter
----------

Kernel Projection
-----------------

Not yet implemented.

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
           doi: 10.1126/science.290.5500.2319

    .. [2] Sam T. Roweis, and Lawrence K. Saul,
           "Nonlinear Dimensionality Reduction by Locally Linear Embedding,"
           Science, 290(5500), pp. 2323 - 2326
           doi: 10.1126/science.290.5500.2323

    .. [3] Partha Niyogi, andMikhail Belkin,
           "Laplacian Eigenmaps for Dimensionality Reduction and Data
           Representation",
           Neural Computation, Vol. 15, No. 6, Pages 1373-1396
           doi: 10.1162/089976603321780317

    .. [4] R.R. Coifman, and S. Lafon, "Diffusion maps",
           Applied and Computational Harmonic Analysis,
           Vol 21, July 2006, pp 5-30

    .. [5] David L. Donoho, and Carrie Grimes,
           "Hessian eigenmaps: Locally linear embedding techniques for
           high-dimensional data",
           PNAS, Vol 100, May 2003, n°10, p5591-5596,
           doi: 10.1073/pnas.1031596100

    .. [6] JR. J. Sammon.,
           "Data reduction with NonLinear Mapping algorithm"
           IEEE Transactions on Computers, C-18(No. 5):401--409, May 1969

    .. [7] doi: 10.1109/72.554199

    .. [9] Schölkopf, B. and Smola, A.J.,
           "Learning with kernels"
