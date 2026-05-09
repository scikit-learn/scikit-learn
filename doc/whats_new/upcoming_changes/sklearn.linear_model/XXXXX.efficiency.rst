- :func:`linear_model._base._rescale_data` avoids constructing an intermediate
  ``n_samples x n_samples`` diagonal sparse matrix when scaling sparse ``X`` or
  ``y`` by sample weights. Row-wise scaling is now performed directly via
  :meth:`scipy.sparse.csr_array.multiply`, reducing memory usage and improving
  performance for large sparse datasets.
  By :user:`Vishal Chaurasiya <vishalchaurasiya>`.
