import pytest
from sklearn.feature_extraction import FeatureHasher
from sklearn.feature_extraction.text import HashingVectorizer

def test_feature_hasher_requires_fit_tag():
    """Test that FeatureHasher has requires_fit=False tag."""
    hasher = FeatureHasher()
    assert hasher._more_tags()['requires_fit'] is False

def test_hashing_vectorizer_requires_fit_tag():
    """Test that HashingVectorizer has requires_fit=False tag."""
    vectorizer = HashingVectorizer()
    assert vectorizer._more_tags()['requires_fit'] is False
