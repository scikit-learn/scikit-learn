RANSACRegressor
---------------

- Fixed a bug where ``RANSACRegressor`` did not properly skip iterations
  when ``estimator.score`` returned a non-finite value (e.g. NaN), which
  could previously emit warnings.
