- Fixed truncated :class:`~sklearn.exceptions.UndefinedMetricWarning` message in
  :func:`metrics.class_likelihood_ratios` when ``fp == 0``, caused by implicit
  string statement bug (missing parentheses around multi-line string).
  Also fixed ``__ALL__`` typo to ``__all__`` in :mod:`sklearn.metrics._regression`.
  By :user:`naveen-111q <naveen-111q>`.
