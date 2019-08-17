import pytest

from sklearn.feature_selection import CovarianceThreshold


# Clustered Features by covariance are removed
# Uncorrelated features are kept
# All collreated features -> One feature
# ordering by changing threshold
