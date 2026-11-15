- :func:`utils.check_array` now avoids an unnecessary memory copy for
  :class:`narwhals.DataFrame` backed by pandas when `copy=False` and
  `force_writeable=True` under pandas Copy-on-Write.
  By :user:`Direk Kakkar <direkkakkar319-ops>`.
