# Pairwise Distances Reductions
# =============================
#
#    Author: Julien Jerphanion <git@jjerphan.xyz>
#
# Overview
# --------
#
#    This module provides routines to compute pairwise distances between a set
#    of row vectors of X and another set of row vectors of Y and apply a
#    reduction on top. The canonical example is the brute-force computation
#    of the top k nearest neighbors by leveraging the arg-k-min reduction.
#
#    The reduction takes a matrix of pairwise distances between rows of X and Y
#    as input and outputs an aggregate data-structure for each row of X. The
#    aggregate values are typically smaller than the number of rows in Y, hence
#    the term reduction.
#
#    For computational reasons, the reduction are performed on the fly on chunks
#    of rows of X and Y so as to keep intermediate data-structures in CPU cache
#    and avoid unnecessary round trips of large distance arrays with the RAM
#    that would otherwise severely degrade the speed by making the overall
#    processing memory-bound.
#
#    Finally, the routines follow a generic parallelization template to process
#    chunks of data with OpenMP loops (via Cython prange), either on rows of X
#    or rows of Y depending on their respective sizes.
#
#
# Dispatching to specialized implementations
# ------------------------------------------
#
#    Dispatchers are meant to be used in the Python code. Under the hood, a
#    dispatcher must only define the logic to choose at runtime to the correct
#    dtype-specialized :class:`PairwiseDistancesReduction` implementation based
#    on the dtype of X and of Y.
#
#
# High-level diagram
# ------------------
#
#    Legend:
#
#      A ---⊳ B: A inherits from B
#      A ---x B: A dispatches to B
#
#
#                               (base dispatcher)
#                           PairwiseDistancesReduction
#                                       ∆
#                                       |
#                                       |
#                     +-----------------+-----------------+
#                     |                                   |
#               (dispatcher)                        (dispatcher)
#         PairwiseDistancesArgKmin           PairwiseDistancesRadiusNeighbors
#               |                                              |
#               |                                              |
#               |                                              |
#               |                  (64bit implem.)             |
#               |          PairwiseDistancesReduction64        |
#               |                       ∆                      |
#               |                       |                      |
#               |                       |                      |
#               |     +-----------------+-----------------+    |
#               |     |                                   |    |
#               |     |                                   |    |
#               x     |                                   |    x
#        PairwiseDistancesArgKmin64       PairwiseDistancesRadiusNeighbors64
#               |     ∆                                   ∆    |
#               |     |                                   |    |
#               x     |                                   |    |
#     FastEuclideanPairwiseDistancesArgKmin64             |    |
#                                                         |    |
#                                                         |    x
#                                  FastEuclideanPairwiseDistancesRadiusNeighbors64
#
#    For instance :class:`PairwiseDistancesArgKmin`, dispatches to
#    :class:`PairwiseDistancesArgKmin64` if X and Y are both dense NumPy arrays
#    with a float64 dtype.
#
#    In addition, if the metric parameter is set to "euclidean" or "sqeuclidean",
#    :class:`PairwiseDistancesArgKmin64` further dispatches to
#    :class:`FastEuclideanPairwiseDistancesArgKmin64` a specialized subclass
#    to optimally handle the Euclidean distance case using the Generalized Matrix
#    Multiplication (see the docstring of :class:`GEMMTermComputer64` for details).


from ._dispatcher import (
    PairwiseDistancesReduction,
    PairwiseDistancesArgKmin,
    PairwiseDistancesRadiusNeighborhood,
    sqeuclidean_row_norms,
)

__all__ = [
    "PairwiseDistancesReduction",
    "PairwiseDistancesArgKmin",
    "PairwiseDistancesRadiusNeighborhood",
    "sqeuclidean_row_norms",
]
