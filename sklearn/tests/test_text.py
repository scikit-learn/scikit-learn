import pytest
from sklearn.feature_extraction.text import HashingVectorizer

def test_hashing_vectorizer_requires_fit_tag():
    vectorizer = HashingVectorizer()
    assert vectorizer._get_tags()["requires_fit"] is False
