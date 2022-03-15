import numpy as np

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import ObliqueDecisionTreeClassifier
from sklearn.datasets import load_iris
from sklearn.model_selection import cross_val_score

# set random seed
random_state = 123456

iris = load_iris()

X, y = iris.data, iris.target

# either axis-aligned
clf = DecisionTreeClassifier(
    random_state=random_state,
    # max_leaf_nodes=5,
)

cv_scores = cross_val_score(clf, X, y, scoring="accuracy", cv=10)

print(len(cv_scores))
print(cv_scores)
print(
    "10-Fold CV score for axis-aligned decision tree is: "
    f"{np.mean(cv_scores)} +/- {np.std(cv_scores)}"
)
# assert False
# or oblique
n_features = X.shape[1]
clf = ObliqueDecisionTreeClassifier(
    max_features=n_features,
    random_state=random_state,
    # max_leaf_nodes=5,
)

print("About to fit...")
clf = clf.fit(X, y)
print("now done...?")

cv_scores = cross_val_score(clf, X, y, scoring="accuracy", cv=10, error_score="raise")

# from scipy.sparse import issparse
# print(X.shape)
# print(X)
# print(y)
# print(X[29, :], y[29])
# print(X[57, :], y[57])
# # print(clf.predict(X))
# print(issparse(X))
# print(type(X))
print(len(cv_scores))
print(cv_scores)
print(
    "10-Fold CV score for oblique decision tree is: "
    f"{np.nanmean(cv_scores)} +/- {np.nanstd(cv_scores)}"
)


from sklearn.ensemble import ObliqueRandomForestClassifier
from sklearn.datasets import make_classification
X, y = make_classification(n_samples=1000, n_features=4,
    n_informative=2, n_redundant=0,
    random_state=0, shuffle=False)
clf = ObliqueRandomForestClassifier(max_depth=2, random_state=0)
clf.fit(X, y)

print(clf.predict([[0, 0, 0, 0]]))