"""
=====================================================================
Parameter estimation using grid search with a nested cross-validation
=====================================================================

The classifier is optimized by "nested" cross-validation using the
GridSearchCV object.

The performance of the selected parameters is evaluated using
cross-validation (different than the nested cross-validation that is used
to select the best classifier).

"""
print __doc__

from pprint import pprint

from sklearn import datasets
from sklearn.cross_validation import StratifiedKFold
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.svm import SVC

###############################################################################
# Loading the Digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

# split the dataset in two equal part respecting label proportions
train, test = iter(StratifiedKFold(y, 2)).next()

###############################################################################
# Set the parameters by cross-validation
tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                     'C': [1, 10, 100, 1000]},
                    {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

scores = [
    ('precision', precision_score),
    ('recall', recall_score),
]

for score_name, score_func in scores:
    clf = GridSearchCV(SVC(C=1), tuned_parameters, score_func=score_func)
    clf.fit(X[train], y[train], cv=StratifiedKFold(y[train], 5))
    y_true, y_pred = y[test], clf.predict(X[test])

    print "Classification report for the best estimator: "
    print clf.best_estimator_
    print "Tuned for '%s' with optimal value: %0.3f" % (
        score_name, score_func(y_true, y_pred))
    print classification_report(y_true, y_pred)
    print "Grid scores:"
    pprint(clf.grid_scores_)
    print

# Note the problem is too easy: the hyperparameter plateau is too flat and the
# output model is the same for precision and recall with ties in quality
