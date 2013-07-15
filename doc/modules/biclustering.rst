.. _biclustering:

============
Biclustering
============

Biclustering can be performed with the module
:mod:`sklearn.cluster.bicluster`. Biclustering algorithms simultaneously
cluster rows and columns of a data matrix. These clusters of rows and
columns are known as biclusters. Each determines a submatrix of the
original data matrix with some desired properties.

For instance, given a matrix of shape `(10, 10)`, one possible bicluster
with three rows and two columns induces a submatrix of shape `(3, 2)`::

    >>> import numpy as np
    >>> data = np.arange(100).reshape(10, 10)
    >>> rows = np.array([0, 2, 3])
    >>> columns = np.array([1, 2])
    >>> data[rows][:, columns]
    array([[ 1,  2],
           [21, 22],
           [31, 32]])

For visualization purposes, given a bicluster, the rows and columns of
the data matrix may be rearranged to make the bicluster contiguous.

.. figure:: ./images/bicluster_submatrix.png
   :align: center
   :scale: 75%

Algorithms differ in how they define biclusters. Some of the
common types include:

* constant values, constant rows, or constant columns
* unusually high or low values
* submatrices with low variance
* correlated rows or columns

Algorithms also differ in how rows and columns may be assigned to
biclusters, which leads to different bicluster structures.

Block diagonal or checkerboard structures occur when rows and columns
are divided into partitions. If each row and each column belongs to
exactly one bicluster, then rearranging the data matrix reveals the
biclusters on the diagonal:

.. figure:: ./images/partitioned_biclusters.png
   :align: center
   :scale: 75%

In the checkerboard case, each row belongs to all column clusters, and
each column belongs to all row clusters:

.. figure:: ./images/checkerboard_biclusters.png
   :align: center
   :scale: 75%

Many other structures have been created. In the general case,
each row and column may belong to any number of biclusters:

.. figure:: ./images/unrestricted_biclusters.png
   :align: center
   :scale: 75%

After fitting a model, row and column cluster membership can be found
in the `rows_` and `columns_` attributes. `rows_[i]` is a binary vector
with nonzero entries corresponding to rows that belong to bicluster
`i`. Similarly, `columns_[i]` indicates which columns belong to
bicluster `i`.

Some models also have `row_labels_` and `column_labels_` attributes.
These models partition the rows and columns, such as in the block
diagonal and checkerboard bicluster structures.

.. currentmodule:: sklearn.cluster.bicluster


.. _spectral_coclustering:

Spectral Co-Clustering
======================

The :class:`SpectralCoclustering` algorithm treats the input data
matrix as a bipartite graph: the rows and columns of the matrix
correspond to the two sets of vertices, and each entry corresponds to
an edge between a row and a column. The algorithm approximates the
normalized cut of this graph to find heavy subgraphs.

For the following data matrix:

.. figure:: ./images/bipartite_matrix.png
   :align: center
   :scale: 45%

the corresponding bipartite graph is:

.. figure:: ./images/bipartite_graph.png
   :align: center
   :scale: 45%

The normalized cut of this bipartite graph is the partitioning of
nodes into equal sizes so that edges within partitions have large
weighs and edges between partitions have light weights. After
clustering, the nodes of the graph may be rearranged to show these
partitions:

.. figure:: ./images/bipartite_graph_clustered.png
   :align: center
   :scale: 45%

The corresponding rearrangement of the data matrix shows the
characteristic block diagonal structure:

.. figure:: ./images/bipartite_matrix_clustered.png
   :align: center
   :scale: 45%


Mathematical formulation
------------------------

An approximate solution to the optimal normalized cut may be found via
the generalized eigenvalue decomposition of the Laplacian of the
graph. Usually this would mean working directly with the Laplacian
matrix. If the original data matrix :math:`A` has shape :math:`m
\times n`, the Laplacian matrix for the corresponding bipartite graph
has shape :math:`(m + n) \times (m + n)`. However, in this case it is
possible to work directly with :math:`A`, which is smaller and more
efficient.

The input matrix :math:`A` is preprocessed as follows:

.. math::
    A_n = R^{−1/2} A C^{−1/2}

Where :math:`R` is the diagonal matrix with entry :math:`i` equal to
:math:`\sum_{j} A_{ij}` and :math:`C` is the diagonal matrix with
entry :math:`j` equal to :math:`\sum_{i} A_{ij}`.

The singular value decomposition, :math:`A_n = U \Sigma V^\top`,
provides the partitions of the rows and columns of :math:`A`. A subset
of the left singular vectors gives the row partitions, and a subset
of the right singular vectors gives the column partitions.

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
biclusters. This partitioning reveals the checkerboard structure:

.. figure:: ./images/checkerboard_example.png
   :align: center
   :scale: 30%

The algorithm tries to partition the rows and columns to uncover this
structure.


Mathematical formulation
------------------------

The input matrix :math:`A` is first normalized to make the
checkerboard pattern more obvious. There are three possible methods:

1. *Independent row and column normalization*, as in Spectral
   Co-Clustering. This method makes the rows sum to a constant and the
   columns sum to a different constant.

2. **Bistochastization**: repeated row and column normalization until
   convergence. This method makes both rows and columns sum to the
   same constant.

3. **Log normalization**: the log of the data matrix is computed: :math:`L =
   \log A`. Then the column mean :math:`\overline{L_{i \cdot}}`, row mean
   :math:`\overline{L_{\cdot j}}`, and overall mean :math:`\overline{L_{\cdot
   \cdot}}` of :math:`L` are computed. The final matrix is computed
   according to the formula

.. math::
    K_{ij} = L_{ij} - \overline{L_{i \cdot}} - \overline{L_{\cdot
    j}} + \overline{L_{\cdot \cdot}}

After normalizing, the first few singular vectors are computed, just
as in the Spectral Co-Clustering algorithm.

If log normalization was used, all the singular vectors are
meaningful. However, if independent normalization or bistochastization
were used, the first singular vectors, :math:`u_1` and :math:`v_1`.
are discarded. From now on, the "first" singular vectors refers to
:math:`u_2 \dots u_{p+1}` and :math:`v_2 \dots v_{p+1}` except in the
case of log normalization.

Given these singular vectors, they are ranked according to which can
be best approximated by a piecewise-constant vector. The
approximations for each vector are found using one-dimensional k-means
and scored using the Euclidean distance. Some subset of the best left
and right singular vector are selected. Next, the data is projected to
this best subset of singular vectors and clustered.

For instance, if :math:`p` singular vectors were calculated, the
:math:`q` best are found as described, where :math:`q<p`. Let
:math:`U_b` be the matrix with columns the :math:`q` best left
singular vectors, and similarly :math:`V_b` for the right. To
partition the rows, the rows of :math:`A` are projected to a :math:`q`
dimensional space: :math:`A * V_{b}`. Treating the :math:`m` rows of
this :math:`m \times q` matrix as samples and clustering using k-means
yields the row labels. Similarly, projecting the columns to
:math:`A^{\top} * U_{b}` and clustering this :math:`n \times q` matrix
yields the column labels.

.. topic:: Examples:

 * :ref:`example_bicluster_spectral_biclustering.py`: a simple example
   showing how to generate a checkerboard matrix and bicluster it.

.. topic:: References:

 * Kluger, Yuval, et al, 2003. `Spectral biclustering of microarray
   data: coclustering genes and conditions
   <http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.135.1608>`__.
