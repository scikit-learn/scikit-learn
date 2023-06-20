from sklearn.preprocessing import StandardScaler, OneHotEncoder, PolynomialFeatures
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.feature_selection import SelectKBest
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.ensemble import StackingClassifier
from sklearn.preprocessing import FunctionTransformer
from sklearn.pipeline import Pipeline, FeatureUnion
from sklearn.base import BaseEstimator
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

pipeline0 = Pipeline([
    ('custom', CustomEstimator_wrong()),
    ('custom_right', CustomEstimator())
])

pipeline1 = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', LogisticRegression())
])

preprocessing = FeatureUnion([
    ('numeric_pipeline', Pipeline([
        ('scaler', StandardScaler()),
        ('pca', PCA(n_components=2))
    ])),
    ('categorical_pipeline', Pipeline([
        ('encoder', OneHotEncoder())
    ]))
])

pipeline2 = Pipeline([
    ('preprocess', preprocessing),
    ('classifier', LogisticRegression())
])

preprocessing = FeatureUnion([
    ('numeric_pipeline', Pipeline([
        ('scaler', StandardScaler()),
        ('polynomial', PolynomialFeatures(degree=2))
    ])),
    ('categorical_pipeline', Pipeline([
        ('encoder', OneHotEncoder())
    ]))
])

pipeline3 = Pipeline([
    ('preprocess', preprocessing),
    ('regressor', LinearRegression())
])

pipeline4 = Pipeline([
    ('feature_selection', SelectKBest(k=5)),
    ('classifier', LogisticRegression())
])

estimators = [
    ('rf', RandomForestClassifier()),
    ('svc', SVC())
]

stacking_classifier = StackingClassifier(
    estimators=estimators,
    final_estimator=LogisticRegression()
)

pipeline5 = Pipeline([
    ('classifier', stacking_classifier)
])

def custom_transformer(X):
    # Perform custom transformations on X
    # ...
    return X

pipeline6 = Pipeline([
    ('custom_transform', FunctionTransformer(custom_transformer)),
    ('pca', PCA(n_components=5)),
    ('classifier', LogisticRegression())
])

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

# Main pipeline
pipeline7 = Pipeline([
    ('preprocess1', preprocessing1),
    ('preprocess2', preprocessing2),
    ('feature_selection', SelectKBest(k=5)),
    ('classifier', StackingClassifier(
        estimators=[
            ('rf', RandomForestClassifier()),
            ('svc', SVC())
        ],
        final_estimator=LogisticRegression()
    ))
])

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

# Main pipeline
pipeline8 = Pipeline([
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

pipelines = [pipeline0, pipeline1, pipeline2, pipeline3, pipeline4, pipeline5, pipeline6, pipeline7, pipeline8]
for i, pipeline in enumerate(pipelines):
    with open(f'pipeline_repr{i}.html', 'w') as repr_file:
        repr_file.write(pipeline._repr_html_())