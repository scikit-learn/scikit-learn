- :class:`manifold.MDS` no longer accepts the `dissimilarity` parameter.
  Use `metric` instead. Boolean values of `metric` (`True` or `False`) are no
  longer accepted; use `metric_mds` to select metric or non-metric MDS.
  The default value of `init` is now ``"classical_mds"``. With that default,
  `n_init` is ignored. Pass ``init="random"`` to keep the previous random
  initialization and to run several SMACOF starts via `n_init`.
  By :user:`François Paugam <FrancoisPgm>`.
