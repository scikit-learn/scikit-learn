- The `criterion` parameter is now deprecated for classes
  :class:`ensemble.GradientBoostingRegressor`
  and :class:`ensemble.GradientBoostingClassifier`, as both options
  (`"friedman_mse"` and `"squared_error"`) were producing the same results,
  up to floating-point rounding discrepancies and a bug in `"friedman_mse"`.
  By :user:`Arthur Lacote <cakedev0>`
