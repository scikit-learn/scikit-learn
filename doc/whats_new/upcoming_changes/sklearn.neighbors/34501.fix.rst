- :class:`neighbors.LocalOutlierFactor` now detects duplicated training samples
  directly through zero k-distances and warns also when ``novelty=True``. The
  previous score-based heuristic could miss moderate score explosions and warn
  spuriously on genuine extreme outliers.
  By :user:`Emanuele Andaloro <EmaAnd8>`
