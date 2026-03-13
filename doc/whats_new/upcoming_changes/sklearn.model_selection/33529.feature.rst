- :func:`model_selection.cross_validate` now supports a ``return_predictions``
  parameter (``'predict'``, ``'predict_proba'``, ``'predict_log_proba'``,
  ``'decision_function'``). When set, out-of-fold predictions are returned in
  original sample order under the key ``'predictions'``, analogous to
  :func:`model_selection.cross_val_predict` but combined with scoring in a
  single pass.
  By :user:`aakgna <aakgna>`.
