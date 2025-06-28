from sklearn.feature_extraction import FeatureHasher, HashingVectorizer


def test_feature_hasher_requires_fit():
    assert not FeatureHasher()._get_tags()["requires_fit"]


def test_hashing_vectorizer_requires_fit():
    assert not HashingVectorizer()._get_tags()["requires_fit"]
