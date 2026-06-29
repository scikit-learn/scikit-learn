- :class:`neighbors.KDTree` and :class:`neighbors.BallTree` can now be used
  with joblib parallelism (``n_jobs > 1``) when the loky backend memory-maps
  internal arrays as read-only buffers. :issue:`34344` :pr:`XXXXX`
  by :user:`imann_128`.
