- Optimize parallel prediction in :class:`ensemble.RandomForestClassifier`
  and :class:`ensemble.RandomForestRegressor` by removing threading lock
  contention. Predictions are now collected without locks and accumulated
  in the main thread, providing 10-15x faster prediction accumulation on
  multi-core systems when using ``n_jobs > 1``.
  By :user:`Karamvir Singh <karamvirsingh1998>`.

