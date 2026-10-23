- :class:`ensemble.HistGradientBoostingClassifier` and
  :class:`ensemble.HistGradientBoostingRegressor` use much less memory when early
  stopping holds out a validation set from the training data: the data is now split
  after binning instead of before. For instance, fitting on 1 million samples with
  100 features uses about 750 MB less memory. The fitted models are unchanged.
  By :user:`Arthur Lacote <cakedev0>`.
