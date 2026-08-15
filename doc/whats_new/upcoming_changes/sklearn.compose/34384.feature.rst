- Added a ``bias_correction`` parameter to
  :class:`compose.TransformedTargetRegressor` to correct the systematic
  prediction bias caused by Jensen's inequality when using nonlinear target
  transformations. Supported modes are ``"additive"``, ``"multiplicative"``,
  and ``"taylor"`` (second-order Taylor expansion).
  By :user:`Imran Khan <iamimrankhan1997>`.
