
.. _pandas:

=======================
Pandas Interoperability
=======================

This is a first section test.

Sparse DataFrames Handling
=============================

**Issue:**
Sparse DataFrames are not automatically converted to scipy.sparse matrices.

This is an issue which has vastly improved from pandas version 0.21.1 onwards. The conversation from dataframes has been largely optimized and are much faster to convert.

In general, Sparse datastructures (i.e. DataFrames, Series, Arrays) are memory optimised structures of their standard counterparts. They work on the principle that they contain a lot of NaN, 0, or another repeating value (this can be specified), and as such a lot of memory can be saved, which means one can potentially work with datasets that would otherwise be too large to fit into available memory. However one has to be careful they don't get converted into the dense format by mistake.

In Pandas, the sparse datastructrures are: :class:`~pandas.SparseDataFrame`, :class:`~pandas.SparseSeries` and :class:`~pandas.SparseArray`.
The methods: :meth:`.to_sparse(fill_value=0)` and :meth:`.to_dense()` can be used to convert between normal and sparse data structures.
The `.density` property can be called on the sparse structures to report sparseness.

In scipy.sparse we have a number of various sparse matrix classes:

==========  =====================================
Class
==========  =====================================
bsr_matrix  Block Sparse Row matrix
coo_matrix  Sparse matrix in COOrdinate format
csc_matrix  Compresed Sparse Column matrix
csr_matrix  Compresed Row matrix
dia_matrix  Sparse matrix with diagonal storage
dok_matrix  Dictionary of Keys based sparse matrix
lil_matrix  Row-based linked list sparse matrix
==========  =====================================

Example Usage
-------------

  >>> import numpy as np
  >>> import pandas as pd
  >>> from scipy.sparse import coo_matrix, csr_matrix, csc_matrix, issparse
  >>>
  >>> arr = np.random.random(size=(1000, 1000))
  >>> arr[arr < .9] = 0
  >>>
  >>> sparse_df = pd.SparseDataFrame(arr, default_fill_value=0)
  >>> print('Density: {:.2%}'.format(sparse_df.density))
  >>> # Output: Density: 10.00%
  >>>
  >>> coo = sparse_df.to_coo()
  >>> #or
  >>> coo = coo_matrix(sparse_df)
  >>>
  >>> csr = coo.tocsr()
  >>> csc = coo.tocsc()
  >>>
  >>> print('Confirm both are sparse:', issparse(coo) == issparse(csr) == issparse(csc) == True)
  >>> # Output: Confirm both are sparse: True
  >>> print('Confirm same amount of non-empty values:', coo.nnz == csr.nnz == csc.nnz)
  >>> # Output: Confirm same amount of non-empty values: True


The code above highlights the following three elements:

1) If your sparse value is not NaN then it is important to specify *default_fill_value* property when creating your pandas DataFrame, otherwise no space saving will occur. Check this using the :meth:`.density` property, which should be less than 100% if successful. When creating the scipy sparse matrix, this *default_fill_value* will be used for use as the sparse value (nnz).

2) Either the :meth:`.to_coo()` method on the pandas dataframe, or :meth:`coo_matrix()` constructor are alternative ways you can convert to a scipy sparse datastructure.

3) It is generally better to convert from your pandas Dataframe first to a :class:`coo_matrix`, as this is far quicker to construct, and from this to then convert to a Compressed Row :class:`csr_matrix`, or Compressed Column :class:`csc_matrix` sparse matrix using the :meth:`.tocsr()` or :meth:`.tocsc()` methods respectively.
