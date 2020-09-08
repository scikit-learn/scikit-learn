"""
============================================================
Parameter estimation using grid search with cross-validation
============================================================

This examples shows how a classifier is optimized by cross-validation,
which is done using the :class:`~sklearn.model_selection.GridSearchCV` object
on a development set that comprises only half of the available labeled data.

The performance of the selected hyper-parameters and trained model is
then measured on a dedicated evaluation set that was not used during
the model selection step.

More details on tools available for model selection can be found in the
sections on :ref:`cross_validation` and :ref:`grid_search`.

"""
from sklearn import datasets
from sklearn.base import clone
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report
from sklearn.svm import SVC

print(__doc__)

# Loading the Digits dataset
digits = datasets.load_digits()

# To apply an classifier on this data, we need to flatten the image, to
# turn the data in a (samples, feature) matrix:
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target

# Split the dataset in two equal parts
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.5, random_state=0)

# Set the parameters by cross-validation
tuned_parameters = [{'kernel': ['rbf'], 'gamma': [1e-3, 1e-4],
                     'C': [1, 10, 100, 1000]},
                    {'kernel': ['linear'], 'C': [1, 10, 100, 1000]}]

scores = ['precision_macro', 'recall_macro']

print("# Tuning hyper-parameters")
print()

# Returns the Best Index
def refit_best_index(results):

    # Set the Selected Index as -1
    # As such, if we have no Scoring, then an error will be thrown
    best_index = -1
    # Set the Best Score as 0
    best_score = 0
    for score in scores:
        print("Best parameters set found on development set %s:" % score)
        print()
        best_score_index = results["rank_test_%s" % score].argmin()
        best_score_params = results["params"][best_score_index]
        print()
        print("Grid scores on development set:")
        print()
        means = results['mean_test_%s' % score]
        stds = results['std_test_%s' % score]

        # Check if this is the Biggest Index
        if(best_score < means[best_score_index]):
            best_score = means[best_score_index]
            best_index = best_score_index

        for mean, std, params in zip(means, stds, results['params']):
            print("%0.3f (+/-%0.03f) for %r"
                  % (mean, std * 2, params))
        print()

        print("Detailed classification report:")
        print()
        print("The model is trained on the full development set.")
        print("The scores are computed on the full evaluation set.")
        print()

        best_estimator = clone(clf.estimator)
        best_estimator.set_params(**best_score_params)
        best_estimator.fit(X_train, y_train)
        print("Best parameters set found on development set")
        print(best_score_params)
        y_true, y_pred = y_test, best_estimator.predict(X_test)
        print(classification_report(y_true, y_pred))
        print()
    return best_index

clf = GridSearchCV(
    SVC(), tuned_parameters, scoring=scores, refit=refit_best_index
)
clf.fit(X_train, y_train)

# Note the problem is too easy: the hyperparameter plateau is too flat and the
# output model is the same for precision and recall with ties in quality.
