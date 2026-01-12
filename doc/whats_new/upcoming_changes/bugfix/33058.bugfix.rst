RANSACRegressor now skips iterations where ``estimator.score`` returns a
non-finite value (e.g. ``NaN``) without emitting warnings.
This fixes :pr:`33058`.