- :class:`ensemble.GradientBoostingClassifier` and
  :class:`ensemble.GradientBoostingRegressor` now support missing values in the
  data matrix `X`, for dense inputs. During training, the trees learn at each split
  whether samples with missing values should go to the left or right child.
  By :user:`Arthur Lacote <cakedev0>`.
