- :func:`model_selection.permutation_test_score` now raises a ``UserWarning``
  when each group in ``groups`` has a single value of ``y``. Permuting within
  groups then leaves ``y`` unchanged, so the permutation test is not
  meaningful.
  By :user:`Dongjin Choi <jinchoi7651507>`.
