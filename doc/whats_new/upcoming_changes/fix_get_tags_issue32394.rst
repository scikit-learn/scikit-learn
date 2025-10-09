FIX get_tags handles missing __sklearn_tags__ (issue #32394)
------------------------------------------------------------

`get_tags` now safely handles third-party estimators that do not
implement the ``__sklearn_tags__`` method. Previously, such estimators
could raise an ``AttributeError`` when passed to scikit-learn utilities
that query tags (e.g., ``is_regressor``).

This change ensures consistent behavior across both built-in and
external estimators.

:pr:`32447` by :user:`Bruno César Maymone Galvão <bcmaymonegalvao>`.
