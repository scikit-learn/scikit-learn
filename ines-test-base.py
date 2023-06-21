import pytest
from sklearn.base import BaseEstimator
from sklearn.naive_bayes import GaussianNB
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.cluster import KMeans
from packaging.version import parse as parse_version
import sklearn

class CustomEstimator(BaseEstimator):
    _doc_link = "https://custom_link/{major}.{minor}/docs/{estimator_name}.html"
    _doc_link_module = "custom_module"

    def _get_url_link(self):
        major = 2
        minor = 19
        estimator_name = self.__class__.__name__
        full_url = self._doc_link.format(**locals())
        return full_url
    def __init__(self, parameter1=None, parameter2=None):
        # Initialize any required parameters
        self.parameter1 = parameter1
        self.parameter2 = parameter2

    def fit(self, X, y=None):
        # Implement the fitting logic of your estimator
        # ...
        return self

    def predict(self, X):
        # Implement the prediction logic of your estimator
        # ...
        return X


class CustomEstimator_wrong(BaseEstimator):
    def __init__(self, parameter1=None, parameter2=None):
        # Initialize any required parameters
        self.parameter1 = parameter1
        self.parameter2 = parameter2

    def fit(self, X, y=None):
        # Implement the fitting logic of your estimator
        # ...
        return self

    def predict(self, X):
        # Implement the prediction logic of your estimator
        # ...
        return X

def _local_get_estimator_doc_url(estimator):
    """Generating a link to the API documentation for a given estimator and scikit-learn version."""
    version = parse_version(sklearn.__version__)
    major = version.major
    minor = version.minor
    estimator_name = estimator.__class__.__name__
    estimator_module = ".".join(
        [_ for _ in estimator.__class__.__module__.split(".") if not _.startswith("_")]
    )
    base_url = f"https://scikit-learn.org/{major}.{minor}/modules/generated/"
    full_url = f"{base_url}{estimator_module}.{estimator_name}.html"
    return full_url

@pytest.fixture(params=[BaseEstimator, GaussianNB, StandardScaler, LogisticRegression, SVC, DecisionTreeClassifier, KMeans, CustomEstimator, CustomEstimator_wrong])
def estimator(request):
    return request.param()

def test_get_url_link(estimator):
    # Call the _get_url_link function
    url_link = estimator._get_url_link()
    # Perform assertions to check the expected behavior
    assert isinstance(url_link, str)
    if isinstance(estimator, CustomEstimator):
        assert url_link.startswith("https://custom_link")
    elif isinstance(estimator, CustomEstimator_wrong):
        assert url_link == ""
    else:
        assert url_link.startswith("https://scikit-learn.org/")
        local_url = _local_get_estimator_doc_url(estimator)
        assert url_link == local_url
