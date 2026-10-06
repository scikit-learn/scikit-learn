- :class:`utils.parallel.Parallel` with a thread-based backend no longer resets
  the warning filters in each task. With many threads, this could corrupt the
  process-wide warning filters (losing the user's filters and raising spurious
  warnings) on Python without context-aware warnings, and was costly on
  free-threaded Python.
  By :user:`Arthur Lacote <cakedev0>`.
