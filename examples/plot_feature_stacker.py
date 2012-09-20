from sklearn.pipeline import Pipeline, FeatureStacker
from sklearn.grid_search import GridSearchCV
from sklearn.svm import SVC
from sklearn.datasets import load_iris
from sklearn.decomposition import PCA
from sklearn.feature_selection import SelectKBest

iris = load_iris()

X, y = iris.data, iris.target

# this dataset is way to high-dimensional. Better do PCA:
pca = PCA(n_components=2)

# but maybe some original features where good?
selection = SelectKBest(k=1)

# build estimator from PCA and Univariate selection:

combined_features = FeatureStacker([("pca", pca), ("univ_select", selection)])

# use combined features to transform dataset:
X_features = combined_features.fit(X, y).transform(X)

# classify:
svm = SVC(kernel="linear")
svm.fit(X_features, y)

# Do grid search over k, n_components and C:

pipeline = Pipeline([("features", combined_features), ("svm", svm)])

param_grid = dict(features__pca__n_components=[1, 2, 3],
                  features__univ_select__k=[1, 2],
                  svm__C=[0.1, 1, 10])

grid_search = GridSearchCV(pipeline, param_grid=param_grid, verbose=10)
grid_search.fit(X, y)
print(grid_search.best_estimator_)
