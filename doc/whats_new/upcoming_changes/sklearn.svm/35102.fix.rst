- Fixed the fitted attributes of :class:`svm.SVC`, :class:`svm.NuSVC`,
  :class:`svm.SVR`, :class:`svm.NuSVR` and :class:`svm.OneClassSVM` when some
  training samples have a null or negative `sample_weight`: `support_` now
  contains indices into the training data passed to `fit` instead of the data
  without these samples. In :class:`svm.SVC`, when all the samples of a class are
  ignored, `n_support_`, `dual_coef_`, `intercept_` and the output of
  `predict_proba` and `decision_function` now describe all the classes in
  `classes_`, instead of being inconsistent with it and producing wrong
  predictions. Such a class has no support vector and is never predicted.
  By :user:`bodapatisaikrishna`.
