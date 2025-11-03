- :class:`ensemble.RandomForestClassifier` and
  :class:`ensemble.RandomForestRegressor` now support early stopping via the
  new ``early_stopping_rounds`` parameter. When enabled with ``oob_score=True``,
  training automatically stops when the out-of-bag score stops improving for
  the specified number of consecutive rounds, enabling automatic determination
  of the optimal forest size and reducing training time.
  By :user:`Karamvir Singh <karamvirsingh1998>`.

