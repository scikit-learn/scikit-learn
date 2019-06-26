from sklearn.base import BaseEstimator
from sklearn.datasets import make_classification
from sklearn.svm import SVC, OneClassSVM
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.compose import ResampledTrainer
from sklearn.utils.estimator_checks import check_estimator
from sklearn.utils.validation import _num_samples


class HalfSampler(BaseEstimator):
    "Train with every odd sample"
    def fit_resample(self, X, y, **kw):
        if _num_samples(X) > 1 and getattr(X, 'format', None) not in ['dia',
                                                                      'coo',
                                                                      'bsr']:
            return X[::2], y[::2]
        return X, y


def test_estimator_checks():
    check_estimator(ResampledTrainer(HalfSampler(), SVC()))


def test_pca_outlier_svm():
    # Test the various methods of the pipeline (pca + svm).
    X, y = make_classification(n_classes=2, class_sep=2, weights=[0.1, 0.9],
                               n_informative=3, n_redundant=1, flip_y=0,
                               n_features=20, n_clusters_per_class=1,
                               n_samples=500, random_state=0)

    # Test with PCA + SVC
    clf = SVC(gamma='scale', probability=True, random_state=0)
    pca = PCA()
    outlier = OneClassSVM(gamma='scale')
    pipe = Pipeline([('pca', pca), ('clf', ResampledTrainer(outlier, clf))])
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)


def test_outlier_pca_svm():
    # Test the various methods of the pipeline (pca + svm).
    X, y = make_classification(n_classes=2, class_sep=2, weights=[0.1, 0.9],
                               n_informative=3, n_redundant=1, flip_y=0,
                               n_features=20, n_clusters_per_class=1,
                               n_samples=500, random_state=0)

    # Test with PCA + SVC
    clf = SVC(gamma='scale', probability=True, random_state=0)
    pca = PCA()
    outlier = OneClassSVM(gamma='scale')
    pipe = ResampledTrainer(outlier, Pipeline([('pca', pca), ('svc', clf)]))
    pipe.fit(X, y)
    pipe.predict(X)
    pipe.predict_proba(X)
    pipe.predict_log_proba(X)
    pipe.score(X, y)
