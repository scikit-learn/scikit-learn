import pytest
from sklearn.base import BaseEstimator
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.compose import TransformedTargetRegressor
from sklearn.compose import ColumnTransformer

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
    
    
# creating a pipeline
# Sub-pipeline 1
preprocessing1 = FeatureUnion([
    ('numeric_pipeline', Pipeline([
        ('scaler', StandardScaler()),
        ('pca', PCA(n_components=2))
    ])),
    ('categorical_pipeline', Pipeline([
        ('encoder', OneHotEncoder())
    ]))
])

# Sub-pipeline 2
preprocessing2 = FeatureUnion([
    ('numeric_pipeline', Pipeline([
        ('scaler', StandardScaler()),
        ('pca', PCA(n_components=3))
    ])),
    ('categorical_pipeline', Pipeline([
        ('encoder', OneHotEncoder())
    ]))
])

# Main pipeline that will contain _VisualBlocks
pipeline = lambda: Pipeline([
    ('preprocess1', preprocessing1),
    ('preprocess2', preprocessing2),
    ('feature_selection', SelectKBest(k=5)),
    ('meta_estimator', TransformedTargetRegressor(
        regressor=GradientBoostingClassifier(),
        transformer=StandardScaler()
    )),
    ('column_transformer', ColumnTransformer([
        ('numeric_pipeline', Pipeline([
            ('scaler', StandardScaler()),
            ('pca', PCA(n_components=2))
        ]), ['numeric_feature1', 'numeric_feature2']),
        ('categorical_pipeline', Pipeline([
            ('encoder', OneHotEncoder())
        ]), ['categorical_feature1', 'categorical_feature2'])
    ])),
    ('estimator', RandomForestClassifier())
])

@pytest.fixture(params=[pipeline, GaussianNB, StandardScaler, LogisticRegression, SVC, DecisionTreeClassifier, KMeans, CustomEstimator, CustomEstimator_wrong])
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

def test_repr_html(estimator):
    html_repr = estimator._repr_html_()
    # Perform assertions to check the expected behavior
    if isinstance(estimator, CustomEstimator): 
        assert "https://custom_link" in html_repr
    elif isinstance(estimator, CustomEstimator_wrong):
        assert "<a" not in html_repr
    else:
        assert "https://scikit-learn.org/" in html_repr