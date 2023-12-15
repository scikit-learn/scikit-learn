from pathlib import Path

import numpy as np
import scipy.sparse as sp
from joblib import Memory

from sklearn.datasets import (
    fetch_20newsgroups,
    fetch_olivetti_faces,
    fetch_openml,
    load_digits,
    make_blobs,
    make_classification,
    make_regression,
)
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MaxAbsScaler, StandardScaler

# memory location for caching datasets
M = Memory(location=str(Path(__file__).resolve().parent / "cache"))


@M.cache
def _blobs_dataset(n_samples=500000, n_features=3, n_clusters=100, dtype=np.float32):
    X, _ = make_blobs(
        n_samples=n_samples, n_features=n_features, centers=n_clusters, random_state=0
    )
    X = X.astype(dtype, copy=False)

    X, X_val = train_test_split(X, test_size=0.1, random_state=0)
    return X, X_val, None, None


@M.cache
def _20newsgroups_highdim_dataset(n_samples=None, ngrams=(1, 1), dtype=np.float32):
    newsgroups = fetch_20newsgroups(random_state=0)
    vectorizer = TfidfVectorizer(ngram_range=ngrams, dtype=dtype)
    X = vectorizer.fit_transform(newsgroups.data[:n_samples])
    y = newsgroups.target[:n_samples]

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _20newsgroups_lowdim_dataset(n_components=100, ngrams=(1, 1), dtype=np.float32):
    newsgroups = fetch_20newsgroups()
    vectorizer = TfidfVectorizer(ngram_range=ngrams)
    X = vectorizer.fit_transform(newsgroups.data)
    X = X.astype(dtype, copy=False)
    svd = TruncatedSVD(n_components=n_components)
    X = svd.fit_transform(X)
    y = newsgroups.target

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _mnist_dataset(dtype=np.float32):
    X, y = fetch_openml("mnist_784", version=1, return_X_y=True, as_frame=False)
    X = X.astype(dtype, copy=False)
    X = MaxAbsScaler().fit_transform(X)

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _digits_dataset(n_samples=None, dtype=np.float32):
    X, y = load_digits(return_X_y=True)
    X = X.astype(dtype, copy=False)
    X = MaxAbsScaler().fit_transform(X)
    X = X[:n_samples]
    y = y[:n_samples]

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _synth_regression_dataset(n_samples=100000, n_features=100, dtype=np.float32):
    X, y = make_regression(
        n_samples=n_samples,
        n_features=n_features,
        n_informative=n_features // 10,
        noise=50,
        random_state=0,
    )
    X = X.astype(dtype, copy=False)
    X = StandardScaler().fit_transform(X)

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _synth_regression_sparse_dataset(
    n_samples=10000, n_features=10000, density=0.01, dtype=np.float32
):
    X = sp.random(
        m=n_samples, n=n_features, density=density, format="csr", random_state=0
    )
    X.data = np.random.RandomState(0).randn(X.getnnz())
    X = X.astype(dtype, copy=False)
    coefs = sp.random(m=n_features, n=1, density=0.5, random_state=0)
    coefs.data = np.random.RandomState(0).randn(coefs.getnnz())
    y = X.dot(coefs.toarray()).reshape(-1)
    y += 0.2 * y.std() * np.random.randn(n_samples)

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _synth_classification_dataset(
    n_samples=1000, n_features=10000, n_classes=2, dtype=np.float32
):
    X, y = make_classification(
        n_samples=n_samples,
        n_features=n_features,
        n_classes=n_classes,
        random_state=0,
        n_informative=n_features,
        n_redundant=0,
    )
    X = X.astype(dtype, copy=False)
    X = StandardScaler().fit_transform(X)

    X, X_val, y, y_val = train_test_split(X, y, test_size=0.1, random_state=0)
    return X, X_val, y, y_val


@M.cache
def _olivetti_faces_dataset():
    dataset = fetch_olivetti_faces(shuffle=True, random_state=42)
    faces = dataset.data
    n_samples, n_features = faces.shape
    faces_centered = faces - faces.mean(axis=0)
    # local centering
    faces_centered -= faces_centered.mean(axis=1).reshape(n_samples, -1)
    X = faces_centered

    X, X_val = train_test_split(X, test_size=0.1, random_state=0)
    return X, X_val, None, None


@M.cache
def _random_dataset(
    n_samples=1000, n_features=1000, representation="dense", dtype=np.float32
):
    if representation == "dense":
        X = np.random.RandomState(0).random_sample((n_samples, n_features))
        X = X.astype(dtype, copy=False)
    else:
        X = sp.random(
            n_samples,
            n_features,
            density=0.05,
            format="csr",
            dtype=dtype,
            random_state=0,
        )

    X, X_val = train_test_split(X, test_size=0.1, random_state=0)
    return X, X_val, None, None
