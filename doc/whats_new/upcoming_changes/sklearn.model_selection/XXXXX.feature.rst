- :func:`model_selection.permutation_test_score` now has a ``permutation``
  parameter. ``permutation="blocks"`` permutes the values of ``y`` between the
  groups given in ``groups``, for data where ``y`` is constant within each group.
  By :user:`Dongjin Choi <jinchoi7651507>`.
