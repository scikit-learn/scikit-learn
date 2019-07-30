# Authors: Joel Nothman
#          Oliver Rausch
import numpy as np
from sklearn.base import BaseEstimator
from sklearn.datasets import make_classification
from sklearn.svm import SVC, OneClassSVM
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.compose import ResampledTrainer
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.validation import _num_samples, check_X_y_kwargs


class HalfSampler(BaseEstimator):
    "Train with every second sample"

    def fit_resample(self, X, y, **kws):
        X, y, kws = check_X_y_kwargs(X, y, kws, accept_sparse="csr")
        if _num_samples(X) > 1:
            return X[::2], y[::2], {kw: kws[kw][::2] for kw in kws}

        return X, y, kws


class DataSaver(BaseEstimator):
    "remembers the data that it was fitted with"

    def fit(self, X, y, **kws):
        self.X = X
        self.y = y
        self.kws = kws
        return self

    def predict(self, X):
        return np.zeros((X.shape[0],))

    def transform(self, X):
        return np.zeros((X.shape[0],))


def test_estimator_checks():
    check_estimator(ResampledTrainer(HalfSampler(), SVC()))


def test_correct_halfsampler():
    # check that the estimator is fitted with the correct data
    X = np.zeros((10, 2))
    y = np.arange(10)

    rt = ResampledTrainer(HalfSampler(), DataSaver())
    for method in [rt.fit, rt.fit_transform, rt.fit_predict]:
        method(X, y)

        np.testing.assert_array_equal(
            rt.estimator_.y, np.array([0, 2, 4, 6, 8])
        )
        assert rt.estimator_.kws == {}
        method(X, y, sample_weight=np.arange(10, 20),
               sample_prop=np.arange(20, 30))

        np.testing.assert_array_equal(
            rt.estimator_.y, np.array([0, 2, 4, 6, 8])
        )
        np.testing.assert_array_equal(
            rt.estimator_.kws['sample_weight'], np.array([10, 12, 14, 16, 18])
        )
        np.testing.assert_array_equal(
            rt.estimator_.kws['sample_prop'], np.array([20, 22, 24, 26, 28])
        )


def test_pca_outlier_svm():
    # Test the various methods of the pipeline (pca + svm).
    X, y = make_classification(
        n_classes=2,
        class_sep=2,
        weights=[0.1, 0.9],
        n_informative=3,
        n_redundant=1,
        flip_y=0,
        n_features=20,
        n_clusters_per_class=1,
        n_samples=500,
        random_state=0,
    )

    # Test with PCA + SVC
    clf = SVC(gamma="scale", probability=True, random_state=0)
    pca = PCA()
    outlier = OneClassSVM(gamma="scale")
    pipe = Pipeline([("pca", pca), ("clf", ResampledTrainer(outlier, clf))])
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)


def test_outlier_pca_svm():
    # Test the various methods of the pipeline (pca + svm).
    X, y = make_classification(
        n_classes=2,
        class_sep=2,
        weights=[0.1, 0.9],
        n_informative=3,
        n_redundant=1,
        flip_y=0,
        n_features=20,
        n_clusters_per_class=1,
        n_samples=500,
        random_state=0,
    )

    # Test with PCA + SVC
    clf = SVC(gamma="scale", probability=True, random_state=0)
    pca = PCA()
    outlier = OneClassSVM(gamma="scale")
    pipe = ResampledTrainer(outlier, Pipeline([("pca", pca), ("svc", clf)]))
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)
