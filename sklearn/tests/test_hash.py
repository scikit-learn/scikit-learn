import pytest
from sklearn.feature_extraction import FeatureHasher

def test_feature_hasher_requires_fit_tag():
    hasher = FeatureHasher()
    assert hasher._get_tags()["requires_fit"] is False
