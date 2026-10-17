- :class:`ensemble.GradientBoostingClassifier` and
  :class:`ensemble.GradientBoostingRegressor` now accept the `categorical_features`
  parameter, including `"from_dtype"`. The default remains `None`. `"from_dtype"`
  will become the default in version 1.13. Up to 255 categories per feature are
  supported, for all losses, including multi-class classification. Unknown
  categories at prediction time are treated as missing values.
  ``warm_start`` is not supported together with categorical features.
  By :user:`Arthur Lacote <cakedev0>`.
