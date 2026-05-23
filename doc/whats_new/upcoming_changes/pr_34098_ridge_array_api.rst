:mod:`sklearn.linear_model`
...........................

- Fixes a bug where :class:`linear_model.Ridge` and :func:`linear_model.ridge_regression` 
  failed parameter validation or scalar detection when using alternative Array API backends 
  (such as PyTorch or CuPy) with array-valued ``alpha`` hyperparameters.
  by :user:`prem-479`
