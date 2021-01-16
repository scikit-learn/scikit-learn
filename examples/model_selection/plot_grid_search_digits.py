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

# Macro scores are needed because target is multiclass
macro_scores = ['precision_macro', 'recall_macro']

# Instanciate a C-Support Vector Classification as an estimator for the
# GridSearch
estimator = SVC()

# For multi-metric scoring, the parameter refit must be set to a scorer key
# or a callable to refit an estimator with the best parameter setting on
# the whole data
clf = GridSearchCV(
    estimator, tuned_parameters, scoring=macro_scores, refit=False
)
clf.fit(X_train, y_train)

for score in macro_scores:
    print("# Tuning hyper-parameters for %s \n" % score)
    print("Best parameters set found on development set: \n")
    best_rank = clf.cv_results_["rank_test_%s" % score].argmin()
    best_params = clf.cv_results_["params"][best_rank]
    print(best_params)
    print()
    print("Grid scores on development set:\n")
    means = clf.cv_results_['mean_test_%s' % score]
    stds = clf.cv_results_['std_test_%s' % score]
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))
    print()
    print("Detailed classification report: \n")
    print("The model is trained on the full development set.")
    print("The scores are computed on the full evaluation set. \n")
    # Train the estimator with the best set of parameters found
    # by the GridSearch
    estimator.set_params(**best_params)
    estimator.fit(X_train, y_train)
    y_pred = estimator.predict(X_test)
    print(classification_report(y_test, y_pred))
    print()
# Note the problem is too easy: the hyperparameter plateau is too flat and the
# output model is the same for precision and recall with ties in quality.
