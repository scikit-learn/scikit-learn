- Fix division by zero bug in :func:`utils.compute_class_weight` when using
  ``class_weight="balanced"`` with sample weights where some classes have zero
  total weight. Classes with zero weight now correctly receive a weight of 0.0
  instead of producing ``inf`` values.
  By :user:`Karamvir Singh <karamvirsingh1998>`.

