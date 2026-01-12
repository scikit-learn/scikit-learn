RANSACRegressor: skip iterations with non-finite score
------------------------------------------------------

Fixed a bug in ``RANSACRegressor`` where iterations were not skipped when
``estimator.score`` returned a non-finite value (e.g. ``NaN`` or ``inf``),
which could lead to warnings and invalid score comparisons.
