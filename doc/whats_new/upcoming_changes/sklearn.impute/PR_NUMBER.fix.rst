- :class:`impute.IterativeImputer` with `keep_empty_features=True` no longer
  modifies the input array in place and no longer fails on pandas DataFrames
  with `ValueError: assignment destination is read-only`.
  By :user:`digit50 <digit50>`.
