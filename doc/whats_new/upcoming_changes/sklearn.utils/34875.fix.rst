- :func:`utils.sparsefuncs.count_nonzero` now raises a :class:`TypeError` for
  non-CSR sparse input when a negative ``axis`` (``-1`` or ``-2``) is passed,
  consistent with the positive-axis behavior, instead of silently returning
  incorrect counts.
  By :user:`Viraj Mishra <VirajMishra1>`.
