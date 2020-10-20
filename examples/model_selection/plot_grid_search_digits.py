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

print(__doc__)

# %%
# The dataset
# -----------
#
# We will work with the `digits` dataset. The goal is to classify handwritten
# digits images.
from sklearn import datasets

digits = datasets.load_digits()

# %%
# To apply a classifier on this data, we need to flatten the images.
# Each image of 8 by 8 pixels needs to be transformed to a vector of 64 pixels.
# Thus, we will get a final data array of shape (n_images, n_pixels)
n_samples = len(digits.images)
X = digits.images.reshape((n_samples, -1))
y = digits.target
print(
    f"The number of images is {X.shape[0]} and each image contains "
    f"{X.shape[1]} pixels"
)

# %%
# As presented in the introduction, the data will be split into a development
# and evaluation set of equal size.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.5, random_state=0
)

# %%
# Define our grid-search framework
# --------------------------------
#
# We will fine-tune a classifier by searching the best hyper-parameters on the
# development dataset. We need to define the values of the hyper-parameters
# and the scores to compute to later on select the best candidate.
tuned_parameters = [
    {"kernel": ["rbf"], "gamma": [1e-3, 1e-4], "C": [1, 10, 100, 1000]},
    {"kernel": ["linear"], "C": [1, 10, 100, 1000]},
]
scores = ['precision_macro', 'recall_macro']

# %%
# We can also define a function that will implement the strategy to use to
# select the best candidate from the `cv_results_` returns from the
# :class:`~sklearn.model_selection.GridSearchCV`. Once the candidate selected,
# it will be automatically refitted by the
# :class:`~sklearn.model_selection.GridSearchCV` instance. Here, the strategy
# is to short list the model which are the best in terms of precision and
# recall. From the selected models, we finally select the fastest model at
# predicting.
#
# Note that these choices are completely arbitrary.
import pandas as pd


def refit_strategy(cv_results):
    """Define the strategy to select the best estimator.

    The strategy defines here is to select the best ranking model in terms of
    precision and recall. Once these models are selected, we can select the
    fastest model to predict.

    Parameters
    ----------
    cv_results : dict of numpy (masked) ndarrays
        CV results as returned by the `GridSearchCV`.

    Returns
    -------
    best_index : int
        The index of the best estimator as it appears in `cv_results`.
    """
    # print the info about the grid-search for the different scores
    for score in scores:
        print("\nGrid scores on development set using {score}:\n")
        for mean, std, params in zip(
            cv_results[f"mean_test_{score}"],
            cv_results[f"std_test_{score}"],
            cv_results["params"],
        ):
            print(f"{mean:0.3f} (+/-{std * 2:0.03f}) for {params}")

    cv_results_ = pd.DataFrame(cv_results)
    cv_results_ = cv_results_[
        [
            "mean_score_time",
            "mean_test_recall_macro",
            "mean_test_precision_macro",
            "rank_test_recall_macro",
            "rank_test_precision_macro",
        ]
    ]

    # Select the most performant models in terms of precision/recall
    mask_best_recall = cv_results_["rank_test_recall_macro"] == 1
    mask_best_precision = cv_results_["rank_test_precision_macro"] == 1
    mask_combined = mask_best_recall & mask_best_precision
    print(
        f"The following models are the most performance in terms of "
        f"recall/precision:\n\n {cv_results_[mask_combined]}"
    )

    # From the best candidates, select the fastest model to predict
    best_index = cv_results_.loc[mask_combined, "mean_score_time"].idxmin()

    print(
        f"\nWe selected the model with the following parameters which maximize"
        f" the performance metric and minimize the scoring time:\n\n"
        f"{cv_results_.loc[best_index]}"
    )

    return cv_results_.index[best_index]


# %%
# Once we defined our strategy to select the best model, we can create the
# grid-search instance. Subsequently, we can check the best parameters found.
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC

clf = GridSearchCV(
    SVC(), tuned_parameters, scoring=scores, refit=refit_strategy
)

print("# Tuning hyper-parameters \n")
clf.fit(X_train, y_train)
print(f"\nThe best set of parameters found are:\n{clf.best_params_}")

# %%
# Finally, we evaluate the fine-tuned model on the left-out evaluation set.
from sklearn.metrics import classification_report

y_pred = clf.predict(X_test)
print(
    f"\nOur selected model will have the following performance on the "
    f"evaluation set:\n\n {classification_report(y_test, y_pred)}"
)

# %%
# .. note::
#    The problem is too easy: the hyperparameter plateau is too flat and the
#    output model is the same for precision and recall with ties in quality.
