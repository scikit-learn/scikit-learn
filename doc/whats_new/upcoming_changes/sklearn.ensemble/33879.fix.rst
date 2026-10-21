- :class:`ensemble.RandomForestClassifier`,
  :class:`ensemble.RandomForestRegressor`,
  :class:`ensemble.ExtraTreesClassifier`,
  :class:`ensemble.ExtraTreesRegressor`,
  :class:`ensemble.BaggingClassifier`,
  :class:`ensemble.BaggingRegressor`: Fix OOB predictions and scores to
  use ``NaN`` for samples that were never in the out-of-bag set (which
  can happen when ``n_estimators`` is small), and exclude those samples
  from ``oob_score_`` computation.
  :pr:`XXXXX` by :user:`Dhruv Sharma <dhruv7477>`.
