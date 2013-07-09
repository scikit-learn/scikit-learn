.. _biclustering:

============
Biclustering
============

`Biclustering <http://en.wikipedia.org/wiki/Biclustering>`__ of
unlabeled data can be performed with the module
:mod:`sklearn.bicluster`.

After fitting a model, row and column cluster membership can be found
in the `rows_` and `columns_` attributes. `rows_[i]` is a binary vector
with nonzero entries corresponding to rows that belong to bicluster
`i`.

Some models also have `row_labels_` and `column_labels_` attributes.
These models partition the rows and columns, such as in the block
diagonal and checkerboard bicluster structures.

.. currentmodule:: sklearn.bicluster


.. _spectral_coclustering:

Spectral Co-Clustering
======================

The :class:`SpectralCoclustering` algorithm treats the input data
matrix as a bipartite graph: the rows and columns of the matrix
correspond to the two sets of vertices, and each entry corresponds to
an edge between a row and a column. The algorithm approximates the
normalized cut of this graph to find heavy subgraphs.

The input matrix :math:`A` is preprocessed as follows:

.. math::
    A_n = R^{−1/2} A C^{−1/2}

Where :math:`R` is the diagonal matrix with entry :math:`i` equal to
:math:`\sum_{j} A_{ij}` and :math:`C` is the diagonal matrix with
entry :math:`j` equal to :math:`\sum_{i} A_{ij}`.

The singular value decomposition, :math:`A_n = U \Sigma V^\top`,
provides the partitions of the rows and columns of :math:`A`. A subset
of the left singular vectors gives the word partitions, and a subset
of the right singular vectors gives the song partitions.

The :math:`\ell = \lceil \log_2 k \rceil` singular vectors, starting
from the second, provide the desired partitioning information. They
are used to form the matrix :math:`Z`:

.. math::
    Z = \begin{bmatrix} R^{-1/2} U \\\\
                        C^{-1/2} V
          \end{bmatrix}

where the the columns of :math:`U` are :math:`u_2, \dots, u_{\ell +
1}`, and similarly for :math:`V`.

Then the rows of :math:`Z` are clustered using k-means. The first
n_rows labels provide the row partitioning, and the remaining n_column
labels provide the column partitioning.

.. topic:: Examples:

 * :ref:`example_bicluster_spectral_coclustering.py`: A simple example
   showing how to generate a data matrix with biclusters and apply
   this method to it.


.. topic:: References:

 * Dhillon, Inderjit S, 2001. `Co-clustering documents and words using
   bipartite spectral graph partitioning
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.140.3011>`__.


.. _spectral_biclustering:

Spectral Biclustering
=====================

The :class:`SpectralBiclustering` algorithm assumes that the input
data matrix has a hidden checkerboard structure. For instance, if
there are two row partitions and three column partitions, each row
will belong to three biclusters, and each column will belong to two
biclusters. The outer product of the corresponding row and column
label vectors gives this checkerboard structure. The algorithm tries
to partition the rows and columns to uncover this structure.

The input matrix :math:`A` is first preprocessed to make the
checkerboard pattern more obvious. There are three possible
preprocessing methods:

1. Independent row and column rescaling, as in Spectral Co-Clustering:

    .. math::
        A_n = R^{−1/2} A C^{−1/2}

    This method makes the rows sum to a constant and the columns sum
    to a different constant.

2. Bistochastization

    Repeated row and column rescaling until convergence. This method
    makes both rows and columns sum to the same constant.

3. Log scaling

    The log of the data matrix is computed: :math:`L = \log A`. Then
    the column mean :math:`\overline{L_{i *}}`, row mean
    :math:`\overline{L_{* j}}`, and overall mean
    :math:`\overline{L_{* *}}` of :math:`L` are computed. The
    final matrix is computed according to the formula

    .. math::
        L_{ij} = L_{ij} - \overline{L_{i *}} - \overline{L_{* j}} + \overline{L_{* *}}

The best of the first few right singular vectors of the resulting
matrix are used to project the rows to a lower dimensional space,
where k-means finds the row partitions. Similarly, the best of the
first few left singular vectors are used to project and cluster the
columns.

.. topic:: Examples:

 * :ref:`example_bicluster_spectral_biclustering.py`: a simple example
   showing how to generate a checkerboard matrix and bicluster it.

.. topic:: References:

 * Kluger, Yuval, et al, 2003. `Spectral biclustering of microarray
   data: coclustering genes and conditions
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.135.1608>`__.
