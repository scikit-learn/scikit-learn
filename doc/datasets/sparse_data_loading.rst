.. _sparse_data_loading

Loading sparse data
===================

When dealing with large problems it is often interesting to load data into
sparse matrices to save memory and speed up computations (see the 
`performance documentation <../modules/computational_performance.html>`_ 
on this asect).

A common and efficient format is `Compressed Sparse Row
<http://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html#scipy.sparse.csr_matrix>`_.

The following snippet shows how to read a dataset in CSV format and build a `CSR` 
matrix directly::

  import array
  import csv
  import numpy as np
  from scipy.sparse import csr_matrix

  def csv_to_csr(f):
      """Read content of CSV file f, return as CSR matrix."""
      data = array.array("f")
      indices = array.array("i")
      indptr = array.array("i", [0])

      for i, row in enumerate(csv.reader(f), 1):
          row = np.array(map(float, row))
          shape1 = len(row)
          nonzero = np.where(row)[0]
          data.extend(row[nonzero])
          indices.extend(nonzero)
          indptr.append(indptr[-1]+len(nonzero))

      return csr_matrix((data, indices, indptr),
                        dtype=float, shape=(i, shape1))

A similar code may be easily developed to build such a matrix from other text
formats or from data read from the network.
