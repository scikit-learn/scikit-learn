- Negative values for numeric-encoded categorical features in
  :class:`ensemble.HistGradientBoostingClassifier` and
  :class:`ensemble.HistGradientBoostingRegressor` are no longer treated as
  missing values. This convention (borrowed from LightGBM) had already become
  unreachable through the public API since version 1.4, when categorical
  preprocessing started going through an internal :class:`~sklearn.preprocessing.OrdinalEncoder`
  that always remaps categories to non-negative codes; only categories unseen
  at fit time (regardless of sign) are still treated as missing, as documented.
  This change only affects the internal, private ``_BinMapper`` and
  ``TreePredictor`` classes when used directly, and removes now-dead code.
  By :user:`Arthur Lacote <cakedev0>`.
