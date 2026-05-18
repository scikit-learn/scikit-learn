- The `probability` parameter of :class:`sklearn.svm.SVC` and :class:`sklearn.svm.NuSVC`
  is deprecated due to not being thread-safe and will be removed in 1.11. Use
  :class:`sklearn.calibration.CalibratedClassifierCV` with the respective estimator and
  `ensemble=False` instead.
  By :user:`Shruti Nath <snath-xoc>`.
