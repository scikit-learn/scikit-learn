- `criterion="friedman_mse"` is now deprecated. This criterion was intended for
  gradient boosting but was incorrectly implemented in scikit-learn's trees and
  was actually behaving identically to `criterion="squared_error"`. Use
  `criterion="squared_error"` instead. This affects:
  - :class:`tree.DecisionTreeRegressor`
  - :class:`tree.ExtraTreeRegressor`
  - :class:`ensemble.RandomForestRegressor`
  - :class:`ensemble.ExtraTreesRegressor`
  By :user:`Arthur Lacote <cakedev0>`
