- :class:`ensemble.RandomForestClassifier`, :class:`ensemble.RandomForestRegressor`,
  :class:`ensemble.ExtraTreesClassifier`, :class:`ensemble.ExtraTreesRegressor`,
  :class:`ensemble.GradientBoostingClassifier` and
  :class:`ensemble.GradientBoostingRegressor` no longer check `y` and
  `sample_weight` again in each tree. With many parallel jobs, this makes fitting
  forests up to ~1.3x faster.
  By :user:`Arthur Lacote <cakedev0>`.
